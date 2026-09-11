# ---
# title: Impact Assessment: build per-species x BCR prediction weight rasters
# author: Mannfred Boehm
# ---
# Run ONCE per species (SLURM array) on the cluster, BEFORE 12B/12D.
#
# Builds a multiplicative prediction weight that replicates V5 10.Truncate's
# range / water / data-extent masking:
#
#     weight = range_membership x (not water) x (inside data-limit) x (inside BCR)
#
#   - range_membership : continuous V5 range raster (ranges/{spp}.tif), reprojected
#                        to the BCR prediction grid; outside-range (NA) -> 0.
#                        (V5 applies the range as a multiplicative weight, not a
#                        hard cutoff: mask.i <- truncate2.i * range.i.)
#   - not water        : 1 on land, 0 on WaterMask_Canada polygons
#   - inside data-limit: 1 inside DataLimitationsMask, 0 outside
#   - inside BCR       : 1 inside this subunit's own polygon, 0 outside.
#
# All three vector terms are rasterized with touches = TRUE, which is what V5's
# mask()/crop() do and is NOT terra::rasterize()'s default -- see the note at the
# rasterize calls below.
#
# ADDED 2026-09-11 -- the BCR term was missing and it is the LARGEST of the four.
# V5 predicts each subunit on a BUFFERED grid and cuts it back at 10.Truncate.R:146
# (`crop(vect(sf.i), mask = TRUE)`, sf.i from Subregions_Mosaics_EPSG3978.shp)
# before mosaicking. Measured on our own staged stacks, 59-71% of non-NA pixels in
# {bcr}_2020.tif lie OUTSIDE the subunit's unbuffered polygon, carrying 41-84% of
# the raw density sum. Without this term those buffer pixels enter the density
# tables, and because 12B assigns a subbasin to EVERY BCR it intersects (329 of 674
# subbasins intersect >1 Canadian BCR; 59% of total area), a straddling subbasin was
# summed in full under each of them -- a straight double count at every BCR seam.
# Cropping here reproduces V5's mosaic cut: each BCR contributes only its own share
# and the per-BCR rows in density_tables sum to the subbasin once.
#
# Subregions_Mosaics_EPSG3978.shp and our staged Regions/BAM_BCR_NationalModel_
# Unbuffered.shp are the same geometry (verified: per-BCR areas agree to <1 km^2,
# IoU = 1.0000 on can10/can11/can60), so this reads the file already on the cluster
# and needs no extra Globus staging.
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

# Version stamp written as the band name of weight.tif and checked on re-runs, so a
# weight built by an older version of this script is rebuilt instead of skipped.
# Bump this whenever the weight definition changes.
WEIGHT_VERSION <- "weight_v3_touches"

# ---- Species from SLURM ------------------------------------------------------

species_vec <- c("CAWA", "OVEN")
# species_vec <- sort(c("BANS", "BARS", "BOBO", "CAWA", "EAWP", "EVGR", "GCTH", "GRSP", "GWWA", "LEYE", "OSFL"))
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
species <- species_vec[task_id]
message(Sys.time(), " | building prediction weights for species=", species)

# ---- Static (species-independent) mask vectors -------------------------------

water <- terra::vect(file.path(gis_dir, "WaterMask_Canada.shp"))
limit <- terra::vect(file.path(gis_dir, "DataLimitationsMask.shp"))

# V5's mosaic-cut polygons. Same geometry as V5 Subregions_Mosaics_EPSG3978.shp;
# `code` reconstructs the "can10"-style key the rest of the pipeline uses.
bcr_polys <- terra::vect(file.path(ia_dir, "data", "raw_data", "Regions",
                                   "BAM_BCR_NationalModel_Unbuffered.shp"))
bcr_polys$code <- gsub("_", "", paste(bcr_polys$country, bcr_polys$subUnit, sep = "_"))

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
  # Skip only a weight built by THIS version. The band name is the version stamp:
  # weights written before the BCR-cut fix (2026-09-11) are missing the dominant
  # mask term, and silently skipping them would leave the double count in place.
  if (file.exists(out_path)) {
    existing <- tryCatch(names(terra::rast(out_path))[1], error = function(e) NA_character_)
    if (isTRUE(existing == WEIGHT_VERSION)) {
      message(Sys.time(), " | ", bcr_code, " | weight.tif is current (", WEIGHT_VERSION,
              ") — skipping"); next
    }
    message(Sys.time(), " | ", bcr_code, " | weight.tif is stale (band name '",
            existing, "', want '", WEIGHT_VERSION, "') — REBUILDING")
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
  #
  # touches = TRUE is REQUIRED for V5 parity and is not terra's default. V5 applies
  # these as terra::mask(vect, inverse = TRUE) and terra::crop(vect, mask = TRUE),
  # both of which retain/remove any cell the polygon TOUCHES; terra::rasterize()
  # instead defaults to touches = FALSE (cell-centre rule). The gap is large wherever
  # the layer is fragmented: WaterMask_Canada has 29,545 polygons, most of them thin
  # rivers and small lakes that touch a cell without covering its centre, so the
  # centre rule under-masked water by 4.9-7.7% of total abundance. On the single-blob
  # BCR polygon the same discrepancy is only 0.7-2.2% (perimeter cells). With
  # touches = TRUE all three terms reproduce V5's vector masking exactly -- verified
  # to the digit on OVEN can10 and CAWA can71.
  notwater <- 1 - terra::rasterize(crop_to_grid(water, tmpl), tmpl, field = 1,
                                   background = 0, touches = TRUE)
  inlim    <-     terra::rasterize(crop_to_grid(limit, tmpl), tmpl, field = 1,
                                   background = 0, touches = TRUE)

  # inside this BCR's own polygon (V5 10.Truncate.R:146). The prediction grid is
  # buffered well past the subunit, so this is the dominant exclusion -- see header.
  poly_i <- bcr_polys[bcr_polys$code == bcr_code, ]
  if (nrow(poly_i) == 0)
    stop("no unbuffered polygon for bcr_code=", bcr_code,
         " in BAM_BCR_NationalModel_Unbuffered.shp")
  inbcr <- terra::rasterize(terra::project(poly_i, terra::crs(tmpl)), tmpl,
                            field = 1, background = 0, touches = TRUE)

  weight <- w_range * notwater * inlim * inbcr   # in [0, 1], defined everywhere

  # Every term is built with background = 0 / NA -> 0, so weight must be non-NA
  # everywhere in the extent. 12C multiplies density by it; an NA here would drop
  # the pixel silently instead of zeroing it (V5 10.Truncate.R:141 zeroes).
  if (terra::global(is.na(weight), "sum", na.rm = TRUE)[[1]] > 0)
    stop(bcr_code, ": weight.tif has NA cells - one of the mask terms returned NA ",
         "(an empty rasterize input will do this). Investigate before using it.")
  if (terra::global(weight, "sum", na.rm = TRUE)[[1]] == 0)
    stop(bcr_code, ": weight.tif is all zero - the BCR polygon may not overlap the ",
         "prediction grid. Check the bcr_code -> polygon match.")

  names(weight) <- WEIGHT_VERSION
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  terra::writeRaster(weight, out_path, overwrite = TRUE,
                     wopt = list(gdal = c("COMPRESS=DEFLATE")))

  message(Sys.time(), " | ", bcr_code, " | wrote weight.tif | mean=",
          round(terra::global(weight, "mean", na.rm = TRUE)[[1]], 3),
          " frac_excluded=",
          round(terra::global(weight == 0, "mean", na.rm = TRUE)[[1]], 3),
          " frac_outside_bcr=",
          round(terra::global(inbcr == 0, "mean", na.rm = TRUE)[[1]], 3))
  rm(w_range, notwater, inlim, inbcr, weight); gc()
}

message(Sys.time(), " | species=", species, " | prediction weights complete.")
