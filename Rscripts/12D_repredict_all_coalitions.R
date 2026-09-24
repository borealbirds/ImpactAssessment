# ---
# title: Impact Assessment: re-predict bird densities for ALL 255 coalitions, one BCR per task
# author: Mannfred Boehm
# ---
# One SLURM array task = ONE species x BCR. coalition_task_table() (12F) lists every
# species' 06_bootstraps files in list.files() order and task i is row i, so the task list
# is the same in every task and in 12H. Each task builds the backfilled field ONCE over its
# BCR's superset and reduces all 255 coalitions internally (see "12F restructure
# invariants" in CLAUDE.md), then writes a single file:
#   density_tables/by_bcr/{species}_{year}_{bcr}.rds
# 12H_merge_bcr_tables.R binds those, in the BCR order the one-job-per-species 12D summed
# them in, into density_tables/{species}_{year}_coalition_{cid}.rds (cid 2..256) and
# density_tables/arrays/, which is what 14B and 15A read.
#
# Why per BCR (2026-09-24): the BCRs are independent, but they used to run in series inside
# one 16-core / 384 GB / multi-day job per species. Split, each task is sized to one BCR,
# queues as a small job, and a failure costs one BCR instead of the whole species.

suppressPackageStartupMessages({
  library(BAMexploreR); library(gbm); library(terra); library(tidyverse)
})

# set paths ------------------------------------------------------
nm_root <- "/home/mannfred/projects/def-ecknight/NationalModels"   # must not change

cc <- TRUE; local <- FALSE
if (cc)            { ia_dir <- "/home/mannfred/scratch/impact_assessment" }
if (!cc && local)  { ia_dir <- getwd() }
if (!cc && !local) { ia_dir <- file.path("G:/Shared drives/BAM_NationalModels5", "data", "Extras", "sandbox_data", "impactassessment_sandbox") }

.terra_tmp <- file.path(
  Sys.getenv("SLURM_TMPDIR", unset = tempdir()),
  paste0("terra_", Sys.getenv("SLURM_JOB_ID", unset = "local"),
         "_", Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "0")))
dir.create(.terra_tmp, recursive = TRUE, showWarnings = FALSE)
terra::terraOptions(tempdir = .terra_tmp)

# import data ------------------------------------------------------
load(file.path(ia_dir, "data", "raw_data", "SpeciesPredictionTruncationValues.Rdata"))
bam_boundary         <- terra::vect(file.path(ia_dir, "data", "raw_data", "Regions", "BAM_BCR_NationalModel_Unbuffered.shp"))
all_subbasins_subset <- terra::vect(file.path(ia_dir, "data", "raw_data", "hydrobasins_masked_merged_subset.gpkg"))

bcr_subbasins_ref <- {
  hits <- terra::relate(all_subbasins_subset, bam_boundary, relation = "intersects")
  ij   <- which(hits, arr.ind = TRUE)
  tibble(sub_index = ij[, 1],
         HYBAS_ID  = all_subbasins_subset$first_HYBAS_ID[ij[, 1]],
         bcr_label = paste(bam_boundary$country[ij[, 2]], bam_boundary$subUnit[ij[, 2]], sep = "_"),
         bcr_code  = gsub("_", "", bcr_label))
}

# define covariate types -------------------------------------------------------
categorical_responses <- c("ABoVE_1km", "NLCD_1km", "MODISLCC_1km", "MODISLCC_5x5", "SCANFI_1km", "VLCE_1km")
predictor_metadata <-
  dplyr::tibble(BAMexploreR::predictor_metadata) |>
  dplyr::filter(version == "v5") |>
  dplyr::select(predictor, definition, predictor_class) |>
  dplyr::mutate(dplyr::across('predictor', stringr::str_replace, 'Year', 'year')) |>
  dplyr::mutate(dplyr::across('predictor', stringr::str_replace, 'Method', 'method'))

actually_biotic_what <- c("Peatland_5x5", "Peatland_1km")
actually_biotic_df   <- tibble::tibble(predictor = actually_biotic_what, predictor_class = c("Wetland", "Wetland"))
abiotic_vars <- predictor_metadata |>
  dplyr::filter(predictor_class %in% c("Annual Climate", "Climate Normals", "Topography", "Wetland", "Disturbance", "Time", "Method")) |>
  dplyr::filter(!(predictor %in% actually_biotic_what))
biotic_continuous_vars <- predictor_metadata |>
  dplyr::filter(!(predictor_class %in% abiotic_vars$predictor_class)) |>
  dplyr::bind_rows(actually_biotic_df) |>
  dplyr::filter(!predictor %in% categorical_responses) |>
  dplyr::pull(predictor)

disturbance_vars <- dplyr::tibble(BAMexploreR::predictor_metadata) |>
  dplyr::filter(version == "v5") |>
  dplyr::select(predictor, definition, predictor_class) |>
  dplyr::filter(predictor_class == "Disturbance")

# Shapley utils + the all-coalitions predictor + the fast gbm tree walk ----------
source(file.path(ia_dir, "Rscripts", "12E_shapley_utils.R"))
source(file.path(ia_dir, "Rscripts", "12F_predict_species_all_coalitions.R"))
source(file.path(ia_dir, "Rscripts", "12G_gbm_tree_walk.R"))
Rcpp::sourceCpp(file.path(ia_dir, "Rscripts", "12G_gbm_tree_walk.cpp"))

# pick this task's species x BCR ---------------------------------------------------
species_vec <- c("CAWA", "OVEN")   # 12H_merge_bcr_tables.R must list the same species and year
year <- 2020

tasks   <- coalition_task_table(species_vec, nm_root)
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
message(Sys.time(), " | task table: ", nrow(tasks), " species x BCR task(s); this is task ", task_id)
if (is.na(task_id) || task_id < 1L || task_id > nrow(tasks)) {
  message(Sys.time(), " | no task ", task_id, " (table has ", nrow(tasks), ") - nothing to do")
  quit(save = "no", status = 0)
}
task    <- tasks[task_id, ]
species <- task$species
bcr     <- task$bcr
message(Sys.time(), " | running ALL coalitions for species=", species, " BCR=", bcr)

hirsh_dir  <- file.path(ia_dir, "data", "raw_data", "hirshpearson")

# preflight: this species x BCR must have weight.tif -------------------------------
# 12F itself stop()s on a missing, stale or all-zero/NA weight.tif; this is the cheap early
# copy of that check. Without the weight, totals are quietly wrong rather than obviously
# broken: since 08A median-imputes partial-NA BART predictors, water, out-of-range and
# out-of-extent pixels no longer drop out on their own, so weight.tif is the only masking.
if (!file.exists(file.path(ia_dir, "data", "derived_data", "predictions", species, bcr,
                           year, "weight.tif")))
  stop(species, " ", bcr, " | weight.tif missing - run 12C_build_prediction_weights.sh before 12D")

per_bcr <- predict_species_all_coalitions(species, year = year,
                                          all_subbasins_subset = all_subbasins_subset,
                                          hirsh_dir = hirsh_dir,
                                          save_arrays_ids = coalition_array_ids(),
                                          rdata_files = task$rdata_path,
                                          return_per_bcr = TRUE)

# write this BCR's result, atomically ----------------------------------------------
# A NULL result is a BCR 12F skipped on purpose (no subbasins, empty superset, no mosaic);
# it is still written so 12H can tell "skipped" from "never ran". The temp-then-rename means
# a task killed mid-write leaves no file, never a truncated one 12H would accept. code_md5
# lets 12H refuse to merge BCRs produced by different versions of the prediction code.
out_dir <- file.path(ia_dir, "data", "derived_data", "density_tables", "by_bcr")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_f <- file.path(out_dir, paste0(species, "_", year, "_", bcr, ".rds"))
code_files <- file.path(ia_dir, "Rscripts", c("12E_shapley_utils.R", "12F_predict_species_all_coalitions.R",
                                              "12G_gbm_tree_walk.R", "12G_gbm_tree_walk.cpp"))
saveRDS(list(species   = species, year = year, bcr = bcr, bcr_order = task$bcr_order,
             code_md5  = paste(unname(tools::md5sum(code_files)), collapse = ":"),
             job       = paste0(Sys.getenv("SLURM_ARRAY_JOB_ID", "local"), "_", task_id),
             created   = Sys.time(),
             test_n_boot = as.integer(Sys.getenv("TEST_N_BOOT", "0")),
             result    = if (length(per_bcr) > 0L) per_bcr[[1]] else NULL),
        paste0(out_f, ".tmp"))
if (!file.rename(paste0(out_f, ".tmp"), out_f)) stop("could not move ", out_f, ".tmp into place")
message(Sys.time(), " | wrote ", out_f, if (length(per_bcr) == 0L) " (BCR skipped by 12F)" else "")

message(Sys.time(), " nice.")
