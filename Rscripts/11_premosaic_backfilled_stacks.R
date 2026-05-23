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

mosaic_backfilled_stacks <- function(sub_ids, year, ref) {

  paths <- file.path(
    file.path(ia_dir, "data", "derived_data", "bart_models", year),
    paste0("subbasin_", sub_ids),
    paste0("subbasin_", sub_ids, "_backfill.tif")
  )

  paths <- paths[file.exists(paths)]
  if (length(paths) == 0) return(NULL)

  # Memory-frugal, IO-frugal mosaic (v3, 2026-05-23).
  # v2 cut peak memory during the variable loop but OOM-killed afterwards:
  # 2546 BCR-grid SpatRasters were kept alive in `var_list`, then stitched +
  # masked + written all at once, which materialized everything in RAM.
  # v2 was also slow (~14h for can10) because terra::resample was called
  # once per (variable, subbasin) pair (~200k times for can10).
  #
  # v3 changes:
  #   (1) cheap metadata sweep (unchanged from v2),
  #   (2) PRE-RESAMPLE each subbasin's full stack to the BCR grid ONCE,
  #       writing to its own tempfile. ~80 resamples instead of ~200k --
  #       this collapses the dominant walltime cost.
  #   (3) pre-open lightweight SpatRaster handles to the resampled tifs so
  #       the inner loop is a cheap name lookup + layer reference,
  #   (4) per variable, build a cover() chain across subbasin handles and
  #       FLUSH the accumulator to its own tempfile. The in-memory raster
  #       is then dropped -- only file paths persist across iterations.
  #   (5) the returned SpatRaster is file-backed (2546 single-layer tifs);
  #       downstream mask + writeRaster can stream block-by-block.
  ref1    <- ref[[1]]
  ref_ext <- terra::ext(ref1)

  # (1) cheap metadata sweep -- no pixels read
  meta <- lapply(paths, function(p) {
    r <- terra::rast(p)
    list(path = p, names = names(r), ext = terra::ext(r))
  })

  # overlap filter
  meta <- Filter(function(m) {
    e1 <- m$ext
    e1$xmin < ref_ext$xmax && e1$xmax > ref_ext$xmin &&
      e1$ymin < ref_ext$ymax && e1$ymax > ref_ext$ymin
  }, meta)
  if (length(meta) == 0) return(NULL)

  # (2) pre-resample each subbasin's full stack to BCR grid (on disk)
  rs_dir <- tempfile("premosaic_rs_")
  dir.create(rs_dir, recursive = TRUE, showWarnings = FALSE)
  for (mi in seq_along(meta)) {
    m     <- meta[[mi]]
    e_int <- tryCatch(terra::intersect(m$ext, ref_ext), error = function(e) NULL)
    if (is.null(e_int)) { meta[[mi]]$rs_path <- NA_character_; next }
    r_full <- terra::rast(m$path)
    r_c    <- terra::crop(r_full, e_int)
    pth    <- file.path(rs_dir, sprintf("sub_%04d.tif", mi))
    terra::resample(r_c, ref1, method = "near", filename = pth, overwrite = TRUE)
    meta[[mi]]$rs_path <- pth
    rm(r_full, r_c)
    gc(verbose = FALSE)
    message(Sys.time(), " | pre-resampled subbasin ", mi, "/", length(meta))
  }
  meta <- Filter(function(m) !is.na(m$rs_path), meta)
  if (length(meta) == 0) return(NULL)

  # (3) pre-open SpatRaster handles (metadata only -- no pixels loaded)
  sub_handles <- lapply(meta, function(m) terra::rast(m$rs_path))
  names_per   <- lapply(sub_handles, names)

  vars <- sort(unique(unlist(names_per)))

  # (4) per-variable accumulator, flushed to its own tempfile
  layer_dir <- tempfile("premosaic_layers_")
  dir.create(layer_dir, recursive = TRUE, showWarnings = FALSE)
  layer_paths <- character(length(vars))

  for (k in seq_along(vars)) {
    v   <- vars[k]
    acc <- NULL
    for (i in seq_along(sub_handles)) {
      if (!(v %in% names_per[[i]])) next
      r_layer <- sub_handles[[i]][[v]]
      acc <- if (is.null(acc)) r_layer else terra::cover(acc, r_layer)
    }
    if (!is.null(acc)) {
      lp <- file.path(layer_dir, sprintf("layer_%05d.tif", k))
      terra::writeRaster(acc, lp, overwrite = TRUE)
      layer_paths[k] <- lp
    }
    rm(acc)
    gc(verbose = FALSE)
    if (k %% 50 == 0 || k == length(vars)) {
      message(Sys.time(), " | var ", k, "/", length(vars), " (", v, ") flushed")
    }
  }

  # (5) assemble final stack from on-disk per-layer tifs
  keep <- nzchar(layer_paths)
  if (!any(keep)) return(NULL)
  out <- terra::rast(layer_paths[keep])
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

# load BCR observed stack first; it is the canonical resample target for the mosaic
stack_obs <- terra::rast(file.path(nm_root, "gis/stacks", paste0(bcr_code, "_", year, ".tif")))

message(Sys.time(), " | ", bcr_code, " | starting mosaic")
stack_bf <- mosaic_backfilled_stacks(sub_ids, year, stack_obs)

if (is.null(stack_bf)) {
  message(Sys.time(), " | ", bcr_code, " | no backfilled rasters found — skipping")
  quit(save = "no", status = 0)
}
message(Sys.time(), " | ", bcr_code, " | mosaic done (already on BCR grid)")


# mask + write ------------------------------------------------------
# Combined into one streaming pass: terra::mask(filename=...) block-processes
# the file-backed stack so peak memory stays bounded regardless of layer count.
# Earlier v2 split this into mask -> writeRaster, which materialized all 2546
# layers in RAM at once and OOM-killed the job at 128G.

bam_bcr_codes <- gsub("_", "", paste(bam_boundary$country, bam_boundary$subUnit, sep = "_"))
bcr_poly      <- bam_boundary[bam_bcr_codes == bcr_code, ]

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
terra::mask(stack_bf, bcr_poly, filename = out_path, overwrite = TRUE)
message(Sys.time(), " | ", bcr_code, " | masked and written to ", out_path)
message(Sys.time(), " | done.")
