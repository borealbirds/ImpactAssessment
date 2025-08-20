# ---
# title: Impact Assessment: crop Level 6 HydroBasins to study area (Canada)
# author: Mannfred Boehm
# created: August 19, 2025
# ---

library(terra)
library(tidyverse)

root <- "G:/Shared drives/BAM_NationalModels5"

bam_boundary <- terra::vect(file.path(root, "Regions", "BAM_BCR_NationalModel_Buffered.shp"))

basins_na <- terra::vect("data/raw_data/hydrobasins/hybas_na_lev06_v1c.shp")
basins_ar <- terra::vect("data/raw_data/hydrobasins/hybas_ar_lev06_v1c.shp")

basins <- 
  rbind(basins_na, basins_ar) |> 
  terra::project(x=_, y=bam_boundary) |> 
  terra::crop(x=_, y=ext(bam_boundary)) |> 
  terra::mask(x=_, mask=ext(bam_boundary))

# save cropped/masked soil carbon and pH data 
terra::writeVector(basins, file.path(root, "gis", "other_landscape_covariates", "hydrobasins_masked.gpkg"), overwrite=TRUE)

