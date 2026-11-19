# ---
# title: Impact Assessment: crop Level 6 HydroBasins to study area (Canada)
# author: Mannfred Boehm
# created: August 19, 2025
# ---

library(terra)

root <- "G:/Shared drives/BAM_NationalModels5"

# any V5 1 km layer serves as the grid template
bam_template <- terra::rast(file.path(root, "PredictionRasters", "Biomass", "SCANFI", "1km", "SCANFIBalsamFir_1km_2020.tif"))

bam_boundary <-
  terra::vect(file.path(root, "Regions", "BAM_BCR_NationalModel_UnBuffered.shp")) |>
  terra::aggregate()

basins_na <- terra::vect("data/raw_data/hydrobasins/hybas_na_lev06_v1c.shp")
basins_ar <- terra::vect("data/raw_data/hydrobasins/hybas_ar_lev06_v1c.shp")

basins <-
  rbind(basins_na, basins_ar) |>
  terra::project(bam_boundary) |>
  terra::crop(bam_boundary)

# keep subbasins whose cell centres include at least one non-NA template cell
keep <- terra::extract(!is.na(bam_template), basins, fun = function(x) any(x), touches = FALSE)[, 2]

basins_subset <- basins[which(keep)]

terra::writeVector(basins_subset, file.path("data", "raw_data", "hydrobasins_masked.gpkg"), overwrite = TRUE)
