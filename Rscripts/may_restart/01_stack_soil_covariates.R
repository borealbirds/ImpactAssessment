# ---
# title: Impact Assessment: crop soil properties to study area (Canada)
# author: Mannfred Boehm
# created: May 7, 2025
# ---

library(terra)
library(tidyverse)


# set root path
root <- "G:/Shared drives/BAM_NationalModels5"

# define study area 
bam_boundary <- terra::vect(file.path(root, "Regions", "BAM_BCR_NationalModel_Buffered.shp"))

# import soil carbon and pH data from ISRIC (International Soil Reference and Information Centre)
# https://files.isric.org/soilgrids/latest/data_aggregated/
soil_path <- file.path(root, "gis", "other_landscape_covariates")

soil_covariates <- 
  terra::rast(list.files(soil_path, pattern = "*cm_mean_*", full.names=TRUE)) |> 
  terra::project(x=_, y=bam_boundary) 

# crop to rectangle
soil_covariates_crop <- terra::crop(x=soil_covariates, y=bam_boundary)

# set pixels in rectangle but outside of bam_boundary to NA
soil_covariates_mask <- terra::mask(x=soil_covariates_crop, mask=bam_boundary)

# save cropped/masked soil carbon and pH data 
terra::writeRaster(soil_covariates_mask, filename="soil_covariates_masked.tif")


