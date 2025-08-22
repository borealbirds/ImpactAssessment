# ---
# title: Impact Assessment: train models on watersheds surrounding mines
# author: Mannfred Boehm
# created: May 7, 2025
# ---


library(tidyverse)
library(terra)

# set root path
root <- "G:/Shared drives/BAM_NationalModels5"
bam_boundary <- terra::vect(file.path(root, "Regions", "BAM_BCR_NationalModel_UnBuffered.shp"))


# for every year (1990-2020)
# 1. import year-matched mining layer
# 2. attach the mine to the appropriate watershed_w
# 3. import V5 covariate stack, year-matched disturbance layer, soil layer
# 4. stack everything imported in step 3 and convert to a dataframe
# 5. mask stack to exclude pixels above HF threshold + mine buffer
# 6. train model on low HF areas

# define function for training models on low HF areas surrounding a mine
# this function trains models for continuous biotic features

# define some objects for `backfill_mines_cont` to reference
mines_dir   = file.path(root, "gis", "other_landscape_covariates"),
stacks_dir  = file.path(root, "gis", "stacks"),
hb_filename = "hydrobasins_masked.gpkg",
soils_filename = "isric_soil_covariates_masked.tif"





backfill_mines_cont <- function(year, root){
  
  
  ### 
  # import year-matched mining layer
  ###
  
  # use %d as a placeholder to insert year_y
  mines_filepath <-
    file.path(root, "gis", "other_landscape_covariates", sprintf("mincan_mines_%d_masked.tif", year))
  
  # import mining layer for year_y
  mines_y <- terra::rast(mines_filepath)
  
  # group mines connected by rook's case as a patch
  # assign each patch a unique identifier, set non-mine cells to NA
  patches_y <- terra::patches(mines_y, directions = 4, zeroAsNA = TRUE)
  
  # convert mine patches to polygons (one polygon per mine)
  mines_polygons_y <- terra::as.polygons(patches_y, dissolve = TRUE, na.rm = TRUE)
  names(mines_polygons_y) <- "patch_id"
  
  
  
  ### 
  # identify each mine's subbasin as the training area
  ###
  
  subbasins_filepath <- file.path(root, "gis", "other_landscape_covariates", "hydrobasins_masked.gpkg")
  subbasins <- terra::vect(subbasins_filepath)
  
  # helper function: for a single mine polygon, find associated watershed
  subbasin_search <- function(mine_poly_y) {
  
     mine_centroid <- terra::centroids(mine_poly_y)
     
     # identify corresponding sub-basin, or nearest sub
     hit <- which(terra::relate(mine_centroid, subbasins, relation = "within"))
     
     return(subbasins[hit])
     
  } # close helper function
  
  
  # identify subbasin for each mine patch
  # seq_len creates a vector of integers {1, 2,...N}, where N is the number of mine polygons,
  # i is the index where the above integers iterate through
  subbasins_for_mines <- lapply(seq_len(nrow(mines_polygons_y)), function(i) {
    subbasin_search(mines_polygons_y[i, ])
  })
  
  
  
  ### 
  # get sum of subbasin areas, mean area, standard deviation, and range
  ###
  
  # combine all mine-subbasin SpatVectors into one multi-polygon
  all_subbasins <- do.call(rbind, subbasins_for_mines)
  
  # compute area of each polygon (km^2)
  areas_km2 <- terra::expanse(all_subbasins, unit = "km")
  
  # summary stats
  subbasin_stats <- tibble(
    total_area = sum(areas_km2),
    mean_area  = mean(areas_km2),
    sd_area    = sd(areas_km2),
    min_area   = min(areas_km2),
    max_area   = max(areas_km2)
  )
 
  

  ### 
  # import V5 covariate stack (BCR x year), year-matched disturbance layer, soil layer
  ###
 
  # get BCR(s) per subbasin
  # many subbasins overlap with multiple BCR ids, so
  # we give each unique combination of BCRs returned by bam_get_bcr an identifier
  # this way we can stack and train models once per unique BCR combo 
  # (instead of repeating stacking/training for each subbasin)
  bcr_list <- lapply(seq_len(nrow(all_subbasins)), function(i) {
        BAMexploreR::bam_get_bcr(version = "v5", ext = all_subbasins[i])
    })
  
  # create a unique multi-bcr key by pasting together bcr numbers
  # (53 unique multi- or single bcr combinations)
  bcr_keys <- vapply(bcr_list, function(v) {
    if (length(v) == 0) "" else paste(v, collapse = "_")
  }, character(1))
  
  # associate a bcr key info (535) to each mine-subbasin (535)
  all_subbasins$bcr_keys <- bcr_keys
  
  # identify which mine-subbasins belong to which bcr group
  bcr_groups <- split(seq_len(nrow(all_subbasins)), all_subbasins$bcr_keys)
  
  # create a new environment to 
  # cache mosaicked stacks per combo_key
  stack_cache <- new.env(parent = emptyenv())
  
  # fetch the appropriate V5 covariate stack for year_y
  get_stack_for_key <- function(key, year) {
    cache_id <- paste(key, year, sep = "_")
    if (exists(cache_id, envir = stack_cache)) return(get(cache_id, envir = stack_cache))
    
    # parse "can10_can40" -> c("can10", "can40")
    bcr_keys_split <- paste(strsplit(key, "_", fixed = TRUE)[[1]], paste0(year, ".tif"), sep="_")
    
    # load each BCR × year stack, then mosaic
    stacks <- lapply(bcr_keys_split, function(r) terra::rast(file.path(root, "gis", "stacks", r)))
    
    
    # MOSAICNIG?????
    stacks_mos <- lapply(stacks, function(s) terra::mosaic(s))
    
    assign(cache_id, stk, envir = stack_cache)
    stk
  }
  
  
     
return(subbasin_stats)
     
} # close backfill_mines_cont()






