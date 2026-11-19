# ---
# title: Impact Assessment: crop soil properties to study area (Canada)
# author: Mannfred Boehm
# created: May 7, 2025
# ---

# ISRIC soil organic carbon (soc) and cation exchange capacity (cec), six depths each, from
# https://files.isric.org/soilgrids/latest/data_aggregated/
# no aggregate() as in 02: these layers are already close to 1 km

library(terra)

root    <- "G:/Shared drives/BAM_NationalModels5"
raw_dir <- file.path(root, "data", "Extras", "sandbox_data", "impactassessment_sandbox")   # the downloaded tiles

soil_covariates <- terra::rast(list.files(raw_dir, pattern = "*cm_mean_*", full.names = TRUE))

# a rough crop to the BAM boundary first makes the reprojection much faster
bam_boundary <- terra::vect(file.path(root, "Regions", "BAM_BCR_NationalModel_Buffered.shp"))
bam_boundary_reproj <- terra::project(bam_boundary, terra::crs(soil_covariates))

soil_covariates_roughcrop <-
  terra::crop(soil_covariates, bam_boundary_reproj) |>
  terra::mask(bam_boundary_reproj)

# 1 km template in EPSG:5072 (projecting onto a template avoids GDAL errors)
bam_template <- terra::rast(terra::ext(bam_boundary), resolution = 1000, crs = terra::crs(bam_boundary))

soil_covariates_mask <-
  terra::project(soil_covariates_roughcrop, bam_template, method = "bilinear") |>
  terra::crop(bam_boundary) |>
  terra::mask(bam_boundary)

terra::writeRaster(soil_covariates_mask,
                   file.path("data", "raw_data", "isric_soil", "isric_soil_covariates_masked.tif"),
                   overwrite = TRUE)
