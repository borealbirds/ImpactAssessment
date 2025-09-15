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
# subset subbasins to those associated with mine patches
# all mining locations from 1990-2015 are a subset of those in 2020
# therefore, we'll use the 2020 subbasins for all years

# PART I
# create mining patches by rook's case

# use %d as a placeholder to insert year_y (2020)
year <- 2020
mines_filepath <-
  file.path(ia_dir, sprintf("mincan_mines_%d_masked.tif", year))

# import mining layer for year_y
mines_y <- terra::rast(mines_filepath)

# group mines connected by rook's case as a patch
# assign each patch a unique identifier, set non-mine cells to NA
patches_y <- terra::patches(mines_y, directions = 4, zeroAsNA = TRUE)




# PART II
# trim patches so that we don't backfill urban areas

# import Hirsh-Pearson population density layer
# Creston: 664, Saulte Ste. Marie: 324; Fort Mac 1304
popden <- 
  terra::rast(file.path(ia_dir, "population_density.tif")) |> 
  terra::project(x = _, y = mines_y, method = "near")

# generate urban mask (TRUE for >=400 persons/km^2)
# Stats Can threshold for urbanized is 400: https://www150.statcan.gc.ca/n1/en/catalogue/92-164-X
# Hirsh-Pearson transformed pop density by: threshold == 3.333 * log10(popden + 1) 
density_threshold <- 3.333 * log10(400 + 1) 
urban <- popden >= density_threshold
urban[urban == 0] <- NA  # keep only TRUE as 1, rest NA

# crop to 
urban_polygons <- 
  terra::crop(x = urban, y = patches_y) |> 
  terra::as.polygons(x = _, na.rm = TRUE, dissolve = TRUE)

# convert mine patches to a multi-polygon (one polygon per patch)
mines_polygons_y <- terra::as.polygons(patches_y, dissolve = TRUE, na.rm = TRUE)
names(mines_polygons_y) <- "patch_id"

# remove a buffer cell if it's an urban area (52,505 km^2 from 53,198 km^2)
mines_polygons_y <- if (nrow(urban_polygons) > 0) terra::erase(mines_polygons_y, urban_polygons) else mines_polygons_y
mines_raster_y <- terra::rasterize(mines_polygons_y, mines_y, field = 1)


# PART III
# import subbasins multi-polygon
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
