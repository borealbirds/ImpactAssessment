# ---
# title: Impact Assessment: re-predict bird densities using coalition-based backfilling
# author: Mannfred Boehm
# ---
# Entry point for SLURM array jobs. Each job processes one species for one
# coalition of sectors (specified by COALITION_ID env var). The coalition ID
# maps to a specific subset of sectors via 12E_shapley_utils.R.
#
# With 8 sectors there are 2^8 = 256 coalitions (each sector is either in or
# out: a binary choice). Running all 256 x N_species jobs produces one density
# table per species per coalition; 15B then reads all 256 tables per species to
# compute exact Shapley values that partition the total HF impact.
#
# For a given species x coalition job:
#   1. Build the coalition mask: pixels where any sector in the coalition has
#      footprint AND CanHF >= 1.
#   2. For each BCR the species occupies: read the canonical observed bootstraps
#      from 12A, run joint BRT x BART sampling via predict_species_bcr(), and
#      compute the counterfactual population under the coalition scenario.
#   3. Save the resulting density table to
#      density_tables/{species}_{year}_coalition_{id}.rds
#
# Phase 1: Run 12A_observed.R to produce canonical observed bootstraps.
# Phase 2: Run this script for each coalition x species combination.
# ---

suppressPackageStartupMessages({
  library(BAMexploreR)
  library(gbm)
  library(terra)
  library(tidyverse)
})


# set paths ------------------------------------------------------

# colleague's NationalModels directory — paths here must not change
nm_root <- "/home/mannfred/projects/def-ecknight/NationalModels"

cc    <- TRUE    # TRUE = Compute Canada cluster
local <- FALSE   # TRUE = local RProject machine

if (cc)            { ia_dir <- "/home/mannfred/scratch/impact_assessment" }
if (!cc && local)  { ia_dir <- getwd() }
if (!cc && !local) { ia_dir <- file.path("G:/Shared drives/BAM_NationalModels5", "data", "Extras", "sandbox_data", "impactassessment_sandbox") }

.terra_tmp <- file.path(
  Sys.getenv("SLURM_TMPDIR", unset = tempdir()),
  paste0("terra_", Sys.getenv("SLURM_JOB_ID", unset = "local"),
         "_", Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "0"))
)
dir.create(.terra_tmp, recursive = TRUE, showWarnings = FALSE)
terra::terraOptions(tempdir = .terra_tmp)


# import data ------------------------------------------------------

# prediction thresholds (q.out = high/qsp, l.out = low/denshthresh)
load(file.path(ia_dir, "data", "raw_data", "SpeciesPredictionTruncationValues.Rdata"))

# import BCR boundaries
bam_boundary <- terra::vect(file.path(ia_dir, "data", "raw_data", "Regions", "BAM_BCR_NationalModel_Unbuffered.shp"))

# import merged + subsetted subbasins
all_subbasins_subset <- terra::vect(file.path(ia_dir, "data", "raw_data", "hydrobasins_masked_merged_subset.gpkg"))

# ------------------------------------------------------
# create a reference table for which subbasins are in which BCRs
# some subbasins will be in multiple BCRs, and that's OK because
# we will ultimately crop to the BCR boundary to run the `gbm` model

bcr_subbasins_ref <- {

    # logical matrix: rows=subbasins, cols=BCRs
    hits <- terra::relate(all_subbasins_subset, bam_boundary, relation = "intersects")

    # row/col indexes of TRUE
    ij <- which(hits, arr.ind = TRUE)

    tibble(
       sub_index = ij[, 1],
       HYBAS_ID  = all_subbasins_subset$first_HYBAS_ID[ij[, 1]],
       bcr_label = paste(bam_boundary$country[ij[, 2]],
                    bam_boundary$subUnit[ij[, 2]], sep = "_"),
       bcr_code  = gsub("_", "", bcr_label))
}

# define covariate types -------------------------------------------------------------

# define categorical covariates
categorical_responses <- c("ABoVE_1km", "NLCD_1km","MODISLCC_1km", "MODISLCC_5x5","SCANFI_1km","VLCE_1km")

predictor_metadata <-
  dplyr::tibble(BAMexploreR::predictor_metadata) |>
  dplyr::filter(version == "v5") |>
  dplyr::select(predictor, definition, predictor_class) |>
  dplyr::mutate(dplyr::across('predictor', stringr::str_replace, 'Year', 'year')) |>
  dplyr::mutate(dplyr::across('predictor', stringr::str_replace, 'Method','method'))


# convert some abiotic variables to biotic variables
actually_biotic_what <- c("Peatland_5x5", "Peatland_1km")
actually_biotic_df <- tibble::tibble(predictor = actually_biotic_what, predictor_class = c("Wetland", "Wetland"))

# define abiotic variables
abiotic_vars <-
  predictor_metadata |>
  dplyr::filter(predictor_class %in% c("Annual Climate", "Climate Normals", "Topography", "Wetland", "Disturbance", "Time", "Method")) |>
  dplyr::filter(!(predictor %in% actually_biotic_what))

# define biotic variables
biotic_continuous_vars <-
  predictor_metadata |>
  dplyr::filter(!(predictor_class %in% abiotic_vars$predictor_class)) |>
  dplyr::bind_rows(actually_biotic_df) |>
  dplyr::filter(!predictor %in% categorical_responses) |>
  dplyr::pull(predictor)

# import human footprint rasters ------------------------------------------------------

industry_rasters <- list.files(file.path(ia_dir, "data", "raw_data", "hirshpearson"), pattern = "\\.tif$", full.names = TRUE)
industry_names <- sub(".tif", "", basename(industry_rasters))
industry_stack <- setNames(lapply(industry_rasters, terra::rast), industry_names)

# define disturbance variables, which we'll set to zero when re-predicting
disturbance_vars <-
  dplyr::tibble(BAMexploreR::predictor_metadata) |>
  dplyr::filter(version == "v5") |>
  dplyr::select(predictor, definition, predictor_class) |>
  dplyr::filter(predictor_class == "Disturbance")


# import Shapley utilities and density prediction script -------------------------------------------------

source(file.path(ia_dir, "Rscripts", "12E_shapley_utils.R"))
source(file.path(ia_dir, "Rscripts", "12C_predict_species_bcr.R"))



# run one species on one core ------------------------------------------------------

species_vec <- c("CAWA", "OVEN")
# species_vec <- sort(c("BANS", "BARS", "BOBO", "CAWA", "EAWP", "EVGR", "GCTH", "GRSP", "GWWA", "LEYE", "OSFL"))
# species_vec <- sort(list.dirs(file.path(nm_root, "output/06_bootstraps"), full.names = FALSE, recursive = FALSE))
year <- 2020

# get species index from SLURM
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
species <- species_vec[task_id]
message("running species: ", species)

# get coalition ID from environment — required
coalition_id <- Sys.getenv("COALITION_ID")
if (nchar(coalition_id) == 0) stop("COALITION_ID env var must be set (e.g. export COALITION_ID=42)")
coalition_id <- as.integer(coalition_id)

# map coalition ID to sector names
sectors  <- canonical_sectors()
coalition <- coalition_id_to_sectors(coalition_id, sectors)
message("Coalition ID: ", coalition_id, " = {",
        if (length(coalition) == 0) "empty" else paste(coalition, collapse = ", "), "}")

# skip empty coalition (v = 0 by definition)
if (length(coalition) == 0) {
  message("Empty coalition — nothing to compute. Exiting.")
  quit(save = "no", status = 0)
}

hirsh_dir <- file.path(ia_dir, "data", "raw_data", "hirshpearson")

# coalitions for which we save full bootstrap x scenario arrays (for distribution plots):
# full coalition (all sectors removed) + each of the 8 single-sector coalitions
target_ids  <- c(
  sectors_to_coalition_id(sectors, sectors),
  vapply(sectors, function(s) sectors_to_coalition_id(s, sectors), numeric(1L))
)
save_arrays <- coalition_id %in% target_ids

res <- predict_species_bcr(species, year = year, all_subbasins_subset = all_subbasins_subset,
                           coalition = coalition, coalition_id = coalition_id,
                           hirsh_dir = hirsh_dir, save_arrays = save_arrays)

dt_dir <- file.path(ia_dir, "data", "derived_data", "density_tables")
dir.create(dt_dir, showWarnings = FALSE)

is_bf_only <- isTRUE(res$is_bf_only)
suffix     <- if (is_bf_only) "_bf_only" else ""
saveRDS(res$table, file = file.path(dt_dir,
                                    paste0(species, "_", year, "_coalition_", coalition_id, suffix, ".rds")))
if (is_bf_only) message(Sys.time(), " | bf-only table saved — run 12D_combine.R after 12A completes")

if (save_arrays && !is_bf_only && !is.null(res$national_arrays)) {
  arr_dir <- file.path(dt_dir, "arrays")
  dir.create(arr_dir, showWarnings = FALSE)
  saveRDS(res$national_arrays,
          file = file.path(arr_dir,
                           paste0(species, "_", year, "_coalition_", coalition_id, "_arrays.rds")))
  message(Sys.time(), " | saved national bootstrap arrays for coalition ", coalition_id)
}
message(Sys.time(), " nice.")
