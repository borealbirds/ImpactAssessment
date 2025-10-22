# ---
# title: Impact Assessment: crop Level 5 HydroBasins to study area (Canada)
# author: Mannfred Boehm
# created: August 19, 2025
# ---

library(terra)
library(tidyverse)

root <- "G:/Shared drives/BAM_NationalModels5"
ia_dir <- file.path(root, "data", "Extras", "sandbox_data", "impactassessment_sandbox")

# import an arbitrary bam template for re-projecting HF layer
bam_template <- terra::rast(file.path(root, "PredictionRasters", "Biomass", "SCANFI", "1km", "SCANFIBalsamFir_1km_2020.tif"))

bam_boundary <- 
  terra::vect(file.path(root, "Regions", "BAM_BCR_NationalModel_UnBuffered.shp")) |> 
  terra::aggregate(x = _)

basins_na <- terra::vect("data/raw_data/hydrobasins/hybas_na_lev06_v1c.shp")
basins_ar <- terra::vect("data/raw_data/hydrobasins/hybas_ar_lev06_v1c.shp")

basins <- 
  rbind(basins_na, basins_ar) |>
  terra::project(x = _, y = bam_boundary) |>
  terra::crop(x = _, y = bam_boundary) 

# identify subbasins within bam_template
# second column is the summary value
keep <- terra::extract(
  !is.na(bam_template), 
  basins,
  fun     = function(x) any(x),   # keep if any cell is non-NA
  touches = FALSE,
)[,2]   # count any intersecting cell
                       
basins_subset <- basins[which(keep)]
plot(basins_subset)

# save 
terra::writeVector(basins_subset, file.path(ia_dir, "hydrobasins_masked.gpkg"), overwrite=TRUE)

