# ---
# title: Impact Assessment: pre-mosaic backfilled subbasin stacks per BCR
# author: Mannfred Boehm
# created: 2026-03-23
#
# Purpose: mosaic per-subbasin BART backfill rasters into one BCR-wide stack,
# resample to the BCR covariate grid, mask to the BCR polygon, and write to
#   data/derived_data/bart_models_mosaics/{year}/{bcr_code}_backfilled.tif
#
# Run as a SLURM array (one task per BCR) before 12A_repredict_birds.R.
# 12B_predict_species_bcr.R reads these files instead of re-mosaicing per
# species x sector, which eliminates the dominant ~11h bottleneck per BCR.
# ---

suppressPackageStartupMessages({
  library(terra)
  library(tidyverse)
})


# set paths ------------------------------------------------------

# NationalModels directory 
nm_root <- "/home/mannfred/projects/def-ecknight/NationalModels"

cc    <- TRUE    # TRUE = Compute Canada cluster
local <- FALSE   # TRUE = local RProject machine

if (cc)            { ia_dir <- "/home/mannfred/scratch/impact_assessment" }
if (!cc && local)  { ia_dir <- getwd() }
if (!cc && !local) { ia_dir <- file.path("G:/Shared drives/BAM_NationalModels5", "data", "Extras", "sandbox_data", "impactassessment_sandbox") }


# import data ------------------------------------------------------

bam_boundary       <- terra::vect(file.path(ia_dir, "data", "raw_data", "Regions", "BAM_BCR_NationalModel_Unbuffered.shp"))
all_subbasins_subset <- terra::vect(file.path(ia_dir, "data", "raw_data", "hydrobasins_masked_merged_subset.gpkg"))


# build BCR-subbasin reference table ------------------------------------------------------

bcr_subbasins_ref <- {
  hits <- terra::relate(all_subbasins_subset, bam_boundary, relation = "intersects")
  ij   <- which(hits, arr.ind = TRUE)
  tibble(
    sub_index = ij[, 1],
    HYBAS_ID  = all_subbasins_subset$first_HYBAS_ID[ij[, 1]],
    bcr_label = paste(bam_boundary$country[ij[, 2]],
                      bam_boundary$subUnit[ij[, 2]], sep = "_"),
    bcr_code  = gsub("_", "", bcr_label)
  )
}

bcr_vec <- sort(unique(bcr_subbasins_ref$bcr_code))
bcr_vec <- bcr_vec[startsWith(bcr_vec, "can")]
message("Total Canadian BCRs: ", length(bcr_vec))


# mosaic helper (same logic as in 12A) ------------------------------------------------------

mosaic_backfilled_stacks <- function(sub_ids, year) {

  paths <- file.path(
    file.path(ia_dir, "data", "derived_data", "bart_models", year),
    paste0("subbasin_", sub_ids),
    paste0("subbasin_", sub_ids, "_backfill.tif")
  )

  paths <- paths[file.exists(paths)]
  if (length(paths) == 0) return(NULL)

  stacks   <- lapply(paths, terra::rast)
  full_ext <- Reduce(terra::union, lapply(stacks, terra::ext))
  vars     <- sort(unique(unlist(lapply(stacks, names))))

  var_list <- lapply(vars, function(v) {
    available <- Filter(Negate(is.null),
                        lapply(stacks, function(r) if (v %in% names(r)) r[[v]] else NULL))
    if (length(available) == 0) return(NULL)
    terra::extend(Reduce(terra::merge, available), full_ext)
  })

  keep <- !sapply(var_list, is.null)
  out  <- terra::rast(var_list[keep])
  names(out) <- vars[keep]
  out
}


# get task from SLURM ------------------------------------------------------

year    <- 2020
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
bcr_code <- bcr_vec[task_id]
message(Sys.time(), " | task=", task_id, " BCR=", bcr_code)

out_dir  <- file.path(ia_dir, "data", "derived_data", "bart_models_mosaics", year)
out_path <- file.path(out_dir, paste0(bcr_code, "_backfilled.tif"))

if (file.exists(out_path)) {
  message(Sys.time(), " | ", bcr_code, " | output already exists — skipping")
  quit(save = "no", status = 0)
}

# subbasins for this BCR
sub_ids <- bcr_subbasins_ref |>
  dplyr::filter(bcr_code == !!bcr_code) |>
  dplyr::pull(sub_index) |>
  unique()

message(Sys.time(), " | ", bcr_code, " | subbasins=", length(sub_ids))

if (length(sub_ids) == 0) {
  message(Sys.time(), " | ", bcr_code, " | no subbasins — skipping")
  quit(save = "no", status = 0)
}


# mosaic ------------------------------------------------------

message(Sys.time(), " | ", bcr_code, " | starting mosaic")
stack_bf <- mosaic_backfilled_stacks(sub_ids, year)

if (is.null(stack_bf)) {
  message(Sys.time(), " | ", bcr_code, " | no backfilled rasters found — skipping")
  quit(save = "no", status = 0)
}
message(Sys.time(), " | ", bcr_code, " | mosaic done")


# resample to BCR covariate grid and mask to BCR polygon ------------------------------------------------------

stack_obs <- terra::rast(file.path(nm_root, "gis/stacks", paste0(bcr_code, "_", year, ".tif")))

stack_bf <- terra::resample(stack_bf, stack_obs, method = "near")
message(Sys.time(), " | ", bcr_code, " | resampled to BCR grid")

bam_bcr_codes <- gsub("_", "", paste(bam_boundary$country, bam_boundary$subUnit, sep = "_"))
bcr_poly      <- bam_boundary[bam_bcr_codes == bcr_code, ]
stack_bf      <- terra::mask(stack_bf, bcr_poly)
message(Sys.time(), " | ", bcr_code, " | masked to BCR polygon")


# write ------------------------------------------------------

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
terra::writeRaster(stack_bf, out_path, overwrite = TRUE)
message(Sys.time(), " | ", bcr_code, " | written to ", out_path)
message(Sys.time(), " | done.")
