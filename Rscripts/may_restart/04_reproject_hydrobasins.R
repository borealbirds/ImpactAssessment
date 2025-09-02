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
# subset hydrobasins to subbasins associated with mine patches
# all mining locations from 1990-2015 are a subset of those in 2020
# therefore, we'll use the 2020 subbasins for all years

# PART I
# use %d as a placeholder to insert year_y (2020)
year <- 2020
mines_filepath <-
  file.path(ia_dir, sprintf("mincan_mines_%d_masked.tif", year))

# import mining layer for year_y
mines_y <- terra::rast(mines_filepath)

# group mines connected by rook's case as a patch
# assign each patch a unique identifier, set non-mine cells to NA
patches_y <- terra::patches(mines_y, directions = 4, zeroAsNA = TRUE)

# convert mine patches to a multi-polygon (one polygon per patch)
mines_polygons_y <- terra::as.polygons(patches_y, dissolve = TRUE, na.rm = TRUE)
names(mines_polygons_y) <- "patch_id"



# PART II
# import hydrobasins multi-polygon
subbasins_filepath <- file.path(ia_dir, "hydrobasins_masked.gpkg")
subbasins <- terra::vect(subbasins_filepath)

# helper function: for a single mine patch polygon, find associated watershed
subbasin_search <- function(mine_poly_y) {
  
  mine_centroid <- terra::centroids(mine_poly_y)
  
  # identify corresponding sub-basin
  hit <- which(terra::relate(mine_centroid, subbasins, relation = "within"))
  
  return(subbasins[hit])
  
} # close helper function


# identify subbasin for each mine patch
# seq_len creates a vector of integers {1, 2,...n}, where n is the number of mine polygons,
# i is the index where the above integers iterate through
print("assigning sub-basins to mines")
subbasins_for_mines <- lapply(seq_len(nrow(mines_polygons_y)), function(i) {
  subbasin_search(mines_polygons_y[i, ])
})

# combine all mine-subbasin SpatVectors into one multi-polygon
all_subbasins <- do.call(rbind, subbasins_for_mines)

# save cropped/masked soil carbon and pH data 
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
