# ---
# title: Impact Assessment: build per-species x BCR prediction weight rasters
# author: Mannfred Boehm
# ---
# Run ONCE per species (SLURM array) on the cluster, BEFORE 12B/12D.
#
# Builds a multiplicative prediction weight that replicates V5 10.Truncate's
# range / water / data-extent masking:
#
#     weight = range_membership  x  (not water)  x  (inside data-limit)
#
#   - range_membership : continuous V5 range raster (ranges/{spp}.tif), reprojected
#                        to the BCR prediction grid; outside-range (NA) -> 0.
#                        (V5 applies the range as a multiplicative weight, not a
#                        hard cutoff: mask.i <- truncate2.i * range.i.)
#   - not water        : 1 on land, 0 on WaterMask_Canada polygons
#   - inside data-limit: 1 inside DataLimitationsMask, 0 outside
#
# weight is coalition-independent, so it is built ONCE per species x BCR here and
# reused by every coalition job. 12C and 12D multiply BOTH the observed and the
# backfilled density by this weight, so the obs/bf symmetry is preserved
# exactly:  w*bf - w*obs = w*(bf - obs).
#
# Source masks live in data/raw_data/v5_gis (transferred from G: once); the
# cluster never reads G:. The prediction grid template is the same BCR stack
# (nm_root/gis/stacks/{bcr}_{year}.tif) that 12C/12D use as stack_obs, so the
# weight aligns with both stack_obs and observed_bootstraps.tif with no resample.
#
# Output: data/derived_data/predictions/{species}/{bcr_code}/{year}/weight.tif
# (single band, values in [0, 1], defined everywhere in the BCR extent).
#
# US water (WaterMask_US) is intentionally not used: only Canadian BCRs are
# processed and aggregation is over Canadian subbasins, so US-side water pixels
# never enter any subbasin zone.
# ---

suppressPackageStartupMessages({
  library(terra)
  library(tidyverse)
})

# ---- Paths -------------------------------------------------------------------

nm_root <- "/home/mannfred/projects/def-ecknight/NationalModels"

cc    <- TRUE
local <- FALSE

if (cc)            { ia_dir <- "/home/mannfred/scratch/impact_assessment" }
if (!cc && local)  { ia_dir <- getwd() }
if (!cc && !local) { ia_dir <- file.path("G:/Shared drives/BAM_NationalModels5", "data", "Extras",
                                          "sandbox_data", "impactassessment_sandbox") }

if (!cc) { nm_root <- "G:/Shared drives/BAM_NationalModels5" }

gis_dir <- file.path(ia_dir, "data", "raw_data", "v5_gis")
year    <- 2020

# ---- Species from SLURM ------------------------------------------------------

species_vec <- c("CAWA", "OVEN")
# species_vec <- sort(c("BANS", "BARS", "BOBO", "CAWA", "EAWP", "EVGR", "GCTH", "GRSP", "GWWA", "LEYE", "OSFL"))
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
species <- species_vec[task_id]
message(Sys.time(), " | building prediction weights for species=", species)

# ---- Static (species-independent) mask vectors -------------------------------

water <- terra::vect(file.path(gis_dir, "WaterMask_Canada.shp"))
limit <- terra::vect(file.path(gis_dir, "DataLimitationsMask.shp"))

range_path <- file.path(gis_dir, "ranges", paste0(species, ".tif"))
if (!file.exists(range_path)) stop("no range raster for ", species, " at ", range_path)
range_r <- terra::rast(range_path)

# ---- Find BCR models (same discovery as 12A/12C/12D) -------------------------

rdata_files <- list.files(file.path(nm_root, "output/06_bootstraps", species),
                          pattern = "can.*\\.Rdata$", full.names = TRUE)
message(Sys.time(), " | found ", length(rdata_files), " BCR models")

# crop a vector to the prediction grid extent, working in the vector's native CRS
# first so only the relevant subset is reprojected to the grid CRS.
crop_to_grid <- function(v, tmpl) {
  tmpl_poly <- terra::as.polygons(terra::ext(tmpl), crs = terra::crs(tmpl))
  v_sub     <- terra::crop(v, terra::project(tmpl_poly, terra::crs(v)))
  terra::project(v_sub, terra::crs(tmpl))
}

for (rdata_path in rdata_files) {

  e <- new.env(parent = emptyenv())
  load(rdata_path, envir = e)
  if (!exists("b.list", envir = e)) {
    message("  b.list not found in ", basename(rdata_path), " — skipping"); next
  }
  bcr_code <- attr(e$b.list[[1]], "bcr"); rm(e)

  out_dir  <- file.path(ia_dir, "data", "derived_data", "predictions", species, bcr_code, year)
  out_path <- file.path(out_dir, "weight.tif")
  if (file.exists(out_path)) {
    message(Sys.time(), " | ", bcr_code, " | weight.tif exists — skipping"); next
  }

  stack_path <- file.path(nm_root, "gis/stacks", paste0(bcr_code, "_", year, ".tif"))
  if (!file.exists(stack_path)) {
    message(Sys.time(), " | ", bcr_code, " | no prediction stack — skipping"); next
  }
  tmpl <- terra::rast(stack_path)[[1]]

  # range membership (continuous) on the grid; outside range (NA) -> 0
  w_range <- terra::project(range_r, tmpl, method = "bilinear")
  w_range <- terra::classify(w_range, cbind(NA, 0))

  # not-water (1 = land) and inside-data-limit (1 = inside)
  notwater <- 1 - terra::rasterize(crop_to_grid(water, tmpl), tmpl, field = 1, background = 0)
  inlim    <-     terra::rasterize(crop_to_grid(limit, tmpl), tmpl, field = 1, background = 0)

  weight <- w_range * notwater * inlim   # in [0, 1], defined everywhere in extent

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  terra::writeRaster(weight, out_path, overwrite = TRUE,
                     wopt = list(gdal = c("COMPRESS=DEFLATE")))

  message(Sys.time(), " | ", bcr_code, " | wrote weight.tif | mean=",
          round(terra::global(weight, "mean", na.rm = TRUE)[[1]], 3),
          " frac_excluded=",
          round(terra::global(weight == 0, "mean", na.rm = TRUE)[[1]], 3))
  rm(w_range, notwater, inlim, weight); gc()
}

message(Sys.time(), " | species=", species, " | prediction weights complete.")
