# ---
# title: Impact Assessment: crop soil properties to study area (Canada)
# author: Mannfred Boehm
# created: May 7, 2025
# ---

library(terra)
library(tidyverse)


# set root path
root <- "G:/Shared drives/BAM_NationalModels5"
ia_dir <- file.path(root, "data", "Extras", "sandbox_data", "impactassessment_sandbox")



# import soil carbon and pH data from ISRIC (International Soil Reference and Information Centre)
# https://files.isric.org/soilgrids/latest/data_aggregated/
# note: we do not use `aggregate` here as in "CAfire" because the resolution is already close to 1x1km
soil_path <- file.path(ia_dir)
soil_covariates <- terra::rast(list.files(soil_path, pattern = "*cm_mean_*", full.names=TRUE))


# reproject BAM boundary to roughly crop/mask the global soil layers (speeds up reprojecting the soil layers)
bam_boundary <- terra::vect(file.path(root, "Regions", "BAM_BCR_NationalModel_Buffered.shp"))
bam_boundary_reproj <- terra::project(x=bam_boundary, y=crs(soil_covariates))

# crop/mask the global soil layer to roughly BAM boundary
soil_covariates_roughcrop <-
   terra::crop(x=soil_covariates, y=bam_boundary_reproj) |> 
   terra::mask(x=_, mask=bam_boundary_reproj)

# now that we're approximately cropped to BAM boundary, we can reproject faster
# first, create a template raster from the BAM boundary in EPSG:5072 (avoids GDAL errors)
bam_template <- terra::rast(ext(bam_boundary), resolution=1000, crs=crs(bam_boundary)) 

soil_covariates_mask <- 
  terra::project(x=soil_covariates_roughcrop, y=bam_template, method = "bilinear") |> 
  terra::crop(x=_, y=bam_boundary) |> # crop to rectangle
  terra::mask(x=_, mask=bam_boundary) # set pixels in rectangle but outside of bam_boundary to NA


# save cropped/masked soil carbon and pH data 
terra::writeRaster(soil_covariates_mask, file.path(ia_dir, "isric_soil_covariates_masked.tif"), overwrite=TRUE)


