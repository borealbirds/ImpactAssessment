# ---
# title: Impact Assessment: re-predict bird densities for ALL 255 coalitions in one job
# author: Mannfred Boehm
# ---
# Restructured entry point. One SLURM array job = ONE species (array task 1=CAWA,
# 2=OVEN). Unlike 12B_repredict_birds.R (RETIRED — a pre-2026-09-15 name, not
# today's 12B_v5_truncate.R; it ran one job per species x coalition), this computes
# the backfilled field ONCE per species x BCR over the superset and reduces all 255
# coalitions internally — see "12F restructure invariants" in CLAUDE.md.
#
# Output: density_tables/{species}_{year}_coalition_{cid}.rds  (cid 2..256).
# Observed bootstraps are always present (12A outputs staged on the cluster), so
# every table is complete on write — no _bf_only suffix, no separate combine step.

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

# Shapley utils + the all-coalitions predictor ---------------------------------
source(file.path(ia_dir, "Rscripts", "12E_shapley_utils.R"))
source(file.path(ia_dir, "Rscripts", "12F_predict_species_all_coalitions.R"))

# run one species ------------------------------------------------------
species_vec <- c("CAWA", "OVEN")
year <- 2020

task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
species <- species_vec[task_id]
message("running ALL coalitions for species: ", species)

hirsh_dir  <- file.path(ia_dir, "data", "raw_data", "hirshpearson")
sectors    <- canonical_sectors()
target_ids <- c(sectors_to_coalition_id(sectors, sectors),
                vapply(sectors, function(s) sectors_to_coalition_id(s, sectors), numeric(1L)))

# preflight: every species x BCR must have weight.tif --------------------------
# 12F itself stop()s on a missing, stale or all-zero/NA weight.tif (it used to
# fall back to an UNMASKED run, w == 1, and only warn). This preflight is the
# cheap early copy of that check: fail here rather than 20 h into a job.
# Without the weight, totals are quietly wrong rather than obviously broken — it
# silently counts birds on water, outside the species' range, and outside the V5
# data-limitation extent. That got more dangerous once 08A started median-imputing
# partial-NA BART predictors: water pixels used to fail complete.cases in 12F and
# drop out of BOTH obs and bf on their own, so weight.tif is now the only masking
# left.
# BCR codes are taken from the bootstrap filenames (12F reads them from the
# b.list "bcr" attribute, which would mean loading every .Rdata just to check).
weight_bcrs <- basename(list.files(file.path(nm_root, "output/06_bootstraps", species),
                                   pattern = "can.*\\.Rdata$"))
weight_bcrs <- unique(regmatches(weight_bcrs, regexpr("can[0-9]+", weight_bcrs)))
test_bcr_env <- Sys.getenv("TEST_BCR", "")   # honour the same smoke-test filter as 12C
if (nchar(test_bcr_env) > 0)
  weight_bcrs <- intersect(weight_bcrs, strsplit(test_bcr_env, ",")[[1]])

if (length(weight_bcrs) == 0) {
  warning(species, " | could not derive any BCR codes from 06_bootstraps filenames",
          " — skipping the weight.tif preflight")
} else {
  missing_w <- weight_bcrs[!file.exists(
    file.path(ia_dir, "data", "derived_data", "predictions", species, weight_bcrs,
              year, "weight.tif"))]
  if (length(missing_w) > 0)
    stop(species, " | weight.tif missing for ", length(missing_w), "/",
         length(weight_bcrs), " BCR(s): ", paste(missing_w, collapse = ", "),
         " — run 12C_build_prediction_weights.sh before 12D")
  message(Sys.time(), " | weight.tif present for all ", length(weight_bcrs), " BCR(s)")
}

res <- predict_species_all_coalitions(species, year = year,
                                      all_subbasins_subset = all_subbasins_subset,
                                      hirsh_dir = hirsh_dir, save_arrays_ids = target_ids)

dt_dir <- file.path(ia_dir, "data", "derived_data", "density_tables")
dir.create(dt_dir, showWarnings = FALSE)

for (cid in names(res$tables_by_cid)) {
  tbl <- res$tables_by_cid[[cid]]
  if (is.null(tbl) || nrow(tbl) == 0) {
    message(Sys.time(), " | WARNING: empty table for coalition ", cid, " — not written")
    next
  }
  saveRDS(tbl, file = file.path(dt_dir, paste0(species, "_", year, "_coalition_", cid, ".rds")))
}
message(Sys.time(), " | wrote ", length(res$tables_by_cid), " coalition tables")

# national bootstrap x scenario arrays for 15A (full coalition + 8 singletons)
arr_dir <- file.path(dt_dir, "arrays")
dir.create(arr_dir, showWarnings = FALSE)
for (cid in names(res$arrays_by_cid)) {
  arr <- res$arrays_by_cid[[cid]]
  if (is.null(arr)) {
    message(Sys.time(), " | WARNING: no arrays for coalition ", cid, " — not written")
    next
  }
  saveRDS(arr, file = file.path(arr_dir, paste0(species, "_", year, "_coalition_", cid, "_arrays.rds")))
}
message(Sys.time(), " | wrote ", length(Filter(Negate(is.null), res$arrays_by_cid)), " array files")

message(Sys.time(), " nice.")
