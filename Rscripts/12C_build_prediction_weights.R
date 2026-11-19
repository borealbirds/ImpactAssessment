# ---
# title: Impact Assessment: build per-species x BCR prediction weight rasters
# author: Mannfred Boehm
# ---

# V5's range, water and data-extent masking (10.Truncate.R) as one multiplicative weight, on
# the BCR's prediction grid:
#   weight = range (continuous, outside = 0) x not water x inside data limit x inside the BCR's polygon
#   -> data/derived_data/predictions/{species}/{bcr}/{year}/weight.tif (band name = version stamp)
# 12F multiplies the observed and the counterfactual densities by it, so w*bf - w*obs = w*(bf - obs).
# - the BCR term: V5 predicts each BCR on a grid buffered 100 km past it, and cuts back to the
#   BCR's own polygon only when packaging (10.Truncate.R:146). we read those buffered grids, so
#   the cut happens here; without it a subbasin straddling two BCRs was counted in full under
#   both (59-71% of each grid lies outside its BCR; OL #7a). our unbuffered BCR shapefile has
#   the same geometry as V5's Subregions_Mosaics_EPSG3978.shp
# - touches = TRUE throughout, as V5's mask() and crop() behave; rasterize()'s default (cell
#   centres) under-masked thin rivers and lakes by 4.9-7.7% of abundance (OL #7b)
# - only Canadian BCRs are run, so the US water mask is not needed
# once per species, before 12D: bash 12C_build_prediction_weights.sh OSFL GRSP

library(terra)

nm_root <- "/home/mannfred/projects/def-ecknight/NationalModels"
ia_dir  <- "/home/mannfred/scratch/impact_assessment"
gis_dir <- file.path(ia_dir, "data", "raw_data", "v5_gis")

# bump when the weight's definition changes: a weight with another stamp is rebuilt, not skipped
WEIGHT_VERSION <- "weight_v3_touches"

species_vec <- strsplit(Sys.getenv("SPECIES"), "[, ]+")[[1]]
species     <- species_vec[as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))]
years       <- as.integer(strsplit(Sys.getenv("YEARS", "2020"), ",")[[1]])
message(Sys.time(), " | prediction weights for ", species)

water <- terra::vect(file.path(gis_dir, "WaterMask_Canada.shp"))
limit <- terra::vect(file.path(gis_dir, "DataLimitationsMask.shp"))

bcr_polys <- terra::vect(file.path(ia_dir, "data", "raw_data", "Regions", "BAM_BCR_NationalModel_Unbuffered.shp"))
bcr_polys$code <- paste0(bcr_polys$country, bcr_polys$subUnit)

range_path <- file.path(gis_dir, "ranges", paste0(species, ".tif"))
if (!file.exists(range_path)) stop("no range raster for ", species, " at ", range_path)
range_r <- terra::rast(range_path)

# the species' BCRs, from its model files' names ({species}_{bcr}.Rdata, as 12A reads them)
rdata_files <- list.files(file.path(nm_root, "output/06_bootstraps", species), pattern = "can.*\\.Rdata$")
bcr_codes   <- sub(paste0("^", species, "_(.*)\\.Rdata$"), "\\1", rdata_files)
message(Sys.time(), " | ", length(bcr_codes), " BCR models")

# a vector cropped to the grid's extent in its own CRS first, so only that part is reprojected
crop_to_grid <- function(v, tmpl) {
  tmpl_poly <- terra::as.polygons(terra::ext(tmpl), crs = terra::crs(tmpl))
  terra::project(terra::crop(v, terra::project(tmpl_poly, terra::crs(v))), terra::crs(tmpl))
}

is_current <- function(f) isTRUE(tryCatch(names(terra::rast(f))[1], error = function(e) NA_character_) == WEIGHT_VERSION)

for (bcr_code in bcr_codes) {

  # a current weight for any year of this BCR is copied to the years on the same grid
  built  <- Filter(is_current, Sys.glob(file.path(ia_dir, "data", "derived_data", "predictions", species, bcr_code, "*", "weight.tif")))
  weight <- if (length(built) > 0L) terra::rast(built[1]) else NULL

  for (year in years) {

    out_dir  <- file.path(ia_dir, "data", "derived_data", "predictions", species, bcr_code, year)
    out_path <- file.path(out_dir, "weight.tif")
    if (file.exists(out_path) && is_current(out_path)) {
      message(Sys.time(), " | ", bcr_code, " ", year, " | weight.tif is current - skipping")
      next
    }

    stack_path <- file.path(nm_root, "gis/stacks", paste0(bcr_code, "_", year, ".tif"))
    if (!file.exists(stack_path)) {
      message(Sys.time(), " | ", bcr_code, " ", year, " | no prediction stack - skipping")
      next
    }
    tmpl <- terra::rast(stack_path)[[1]]
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    if (!is.null(weight) && isTRUE(terra::compareGeom(weight, tmpl, stopOnError = FALSE))) {
      terra::writeRaster(weight, out_path, overwrite = TRUE, wopt = list(gdal = c("COMPRESS=DEFLATE")))
      message(Sys.time(), " | ", bcr_code, " ", year, " | weight.tif copied from a current weight on the same grid")
      next
    }

    # every term is 0, never NA, outside its shape, so the weight is defined everywhere (V5
    # zeroes every NA, 10.Truncate.R:141); an NA would drop a pixel instead of zeroing it
    w_range  <- terra::classify(terra::project(range_r, tmpl, method = "bilinear"), cbind(NA, 0))
    notwater <- 1 - terra::rasterize(crop_to_grid(water, tmpl), tmpl, field = 1, background = 0, touches = TRUE)
    inlim    <- terra::rasterize(crop_to_grid(limit, tmpl), tmpl, field = 1, background = 0, touches = TRUE)
    poly_i   <- bcr_polys[bcr_polys$code == bcr_code, ]
    if (nrow(poly_i) == 0) stop("no polygon for ", bcr_code, " in BAM_BCR_NationalModel_Unbuffered.shp")
    inbcr    <- terra::rasterize(terra::project(poly_i, terra::crs(tmpl)), tmpl, field = 1, background = 0, touches = TRUE)

    weight <- w_range * notwater * inlim * inbcr
    if (terra::global(is.na(weight), "sum", na.rm = TRUE)[[1]] > 0)
      stop(bcr_code, ": weight.tif has NA cells (an empty rasterize input does this)")
    if (terra::global(weight, "sum", na.rm = TRUE)[[1]] == 0)
      stop(bcr_code, ": weight.tif is all 0 - does the BCR polygon overlap the prediction grid?")

    names(weight) <- WEIGHT_VERSION
    terra::writeRaster(weight, out_path, overwrite = TRUE, wopt = list(gdal = c("COMPRESS=DEFLATE")))
    message(Sys.time(), " | ", bcr_code, " ", year, " | wrote weight.tif | mean ",
            round(terra::global(weight, "mean", na.rm = TRUE)[[1]], 3), ", share outside the BCR ",
            round(terra::global(inbcr == 0, "mean", na.rm = TRUE)[[1]], 3))
    rm(w_range, notwater, inlim, inbcr)
    gc()
  }
  rm(weight)
  gc()
}
