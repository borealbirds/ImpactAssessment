# ---
# title: Impact Assessment: crop Level 6 HydroBasins to study area (Canada)
# author: Mannfred Boehm
# created: August 19, 2025
# ---

library(terra)
library(tidyverse)

root <- "G:/Shared drives/BAM_NationalModels5"
ia_dir <- file.path(root, "data", "Extras", "sandbox_data", "impactassessment_sandbox")


bam_boundary <- terra::vect(file.path(root, "Regions", "BAM_BCR_NationalModel_Buffered.shp"))

basins_na <- terra::vect("data/raw_data/hydrobasins/hybas_na_lev06_v1c.shp")
basins_ar <- terra::vect("data/raw_data/hydrobasins/hybas_ar_lev06_v1c.shp")

basins <- 
  rbind(basins_na, basins_ar) |> 
  terra::project(x=_, y=bam_boundary) |> 
  terra::crop(x=_, y=ext(bam_boundary)) |> 
  terra::mask(x=_, mask=ext(bam_boundary))

# save cropped/masked soil carbon and pH data 
terra::writeVector(basins, file.path(ia_dir, "hydrobasins_masked.gpkg"), overwrite=TRUE)



# -----------------------------------------------------------
# subset subbasins to those associated with mine and oil & gas patches
# all mining locations from 1990-2015 are a subset of those in 2020
# therefore, we'll use the 2020 subbasins for all years

# PART I
# create polygons for mines and oil&gas

# use %d as a placeholder to insert year_y (2020)
year <- 2020
mines_filepath <- file.path(ia_dir, sprintf("mincan_mines_%d_masked.tif", year))

# import mining layer for year_y
mines_y <- terra::rast(mines_filepath)

# convert mines patches to a multi-polygon (one polygon per mine)
# this is necessary for terra::relate() to find the associated subbasin
mines_polygons_y <- terra::as.polygons(mines_y, dissolve = TRUE, na.rm = TRUE)
names(mines_polygons_y) <- "patch_id"

# import oil & gas layer (Hirsh-Pearson)
oilgas <- terra::rast(file.path(ia_dir, "hirshpearson_oil_gas_masked.tif"))
oilgas_polygons <- terra::as.polygons(oilgas, dissolve = TRUE, na.rm = TRUE)
names(oilgas_polygons) <- "patch_id"

# merge mines and oilgas footprints
combined_poly <- rbind(mines_polygons_y[mines_polygons_y$patch_id == 1],
                       oilgas_polygons[oilgas_polygons$patch_id == 1])
              



# PART II 
# import subbasins multi-polygon
subbasins <- terra::vect(file.path(ia_dir, "hydrobasins_masked.gpkg"))

# identify subbasins that intersect mines or oilgas
hits1 <- terra::relate(subbasins, combined_poly, relation = "intersects")  
hits2 <- unique(which(hits1, arr.ind = TRUE)[, 1])                          
all_subbasins <- subbasins[hits2]                                       

terra::writeVector(all_subbasins, file.path(ia_dir, "hydrobasins_subset_2020.gpkg"), overwrite=TRUE)




# -------------------------------------------------------------------
# get sum of subbasin areas, mean area, standard deviation, and range
# compute area of each polygon (km^2)
areas_km2 <- terra::expanse(all_subbasins, unit = "km")

# summary stats
subbasin_stats <- tibble(
  total_subbasins = length(areas_km2),
  total_area = sum(areas_km2),
  mean_area  = mean(areas_km2),
  sd_area    = sd(areas_km2),
  min_area   = min(areas_km2),
  max_area   = max(areas_km2)
)
