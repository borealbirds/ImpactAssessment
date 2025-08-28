# ---
# title: Impact Assessment: train models on watersheds surrounding mines
# author: Mannfred Boehm
# created: May 7, 2025
# ---

library(BAMexploreR)
library(tidyverse)
library(terra)

# set root path
root <- "G:/Shared drives/BAM_NationalModels5"
bam_boundary <- terra::vect(file.path(root, "Regions", "BAM_BCR_NationalModel_UnBuffered.shp"))

# store predictor metadata as a reference
predictor_metadata <- 
  dplyr::tibble(BAMexploreR::predictor_metadata) |>
  dplyr::filter(version == "v5") |> 
  dplyr::select(predictor, definition, predictor_class) |> 
  dplyr::mutate(across('predictor', str_replace, 'Year', 'year')) |> 
  dplyr::mutate(across('predictor', str_replace, 'Method','method'))
  
categorical_predictors <- c("ABoVE_1km", "method", "NLCD_1km","MODISLCC_1km", 
                            "MODISLCC_5x5", "SCANFI_1km", "VLCE_1km")


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


years <- seq(from = 1990, to = 2020, by = 5)


backfill_mines_cont <- 
  function(year, root, categorical_predictors){
  
  
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
     
     # identify corresponding sub-basin
     hit <- which(terra::relate(mine_centroid, subbasins, relation = "within"))
     
     return(subbasins[hit])
     
  } # close helper function
  
  
  # identify subbasin for each mine patch
  # seq_len creates a vector of integers {1, 2,...N}, where N is the number of mine polygons,
  # i is the index where the above integers iterate through
  print("assigning sub-basins to mines")
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
    total_subbasins = length(areas_km2),
    total_area = sum(areas_km2),
    mean_area  = mean(areas_km2),
    sd_area    = sd(areas_km2),
    min_area   = min(areas_km2),
    max_area   = max(areas_km2)
  )
 
  

  ### 
  # find BCR(s) that intersect with subbasin_s and import V5 covariate stacks,
  # year-matched disturbance layer, and soil layer
  ###
 
  # get BCR(s) per subbasin
  # some subbasins overlap with multiple BCR ids, so
  # we give each unique combination of BCRs an identifier
  # this way we can stack and train models once per unique BCR combo 
  # (instead of repeating stacking/training for each subbasin)
  print("identifying BCR-subbasin intersections")
  bcr_list <- lapply(seq_len(nrow(all_subbasins)), function(i) {
        BAMexploreR::bam_get_bcr(version = "v5", ext = all_subbasins[i])
    })
  
  # create a unique multi-bcr key by pasting together bcr numbers
  # (53 unique multi- or single bcr combinations)
  bcr_keys <- vapply(bcr_list, function(v) {
    if (length(v) == 0) "" else paste(v, collapse = "_")
    }, character(1))
  
  # associate bcr key info (535) to each mine-subbasin (535)
  all_subbasins$bcr_keys <- bcr_keys
  
  # identify which mine-subbasins belong to which bcr group
  # for each group we'll train a separate model 
  # there are 52 groups vs 535 subbasins, so we save time by not 
  # redundantly training models on the same BCRs/BCR groups 
  bcr_groups <- split(seq_len(nrow(all_subbasins)), all_subbasins$bcr_keys)
  
  # build one dissolved polygon per unique BCR key
  # returns a named list: combo_polys[["can40_can41"]] = SpatVector (1 row)
  combo_polys <- lapply(names(bcr_groups), function(key) {
    idx <- bcr_groups[[key]]
    poly <- terra::aggregate(all_subbasins[idx, ], "bcr_keys")
    # dissolve keeps the attribute column; filter to the single row for this key
    poly[poly$bcr_keys == key, ]
  })
  names(combo_polys) <- names(bcr_groups)
  
  
  # fetch the appropriate V5 covariate stack(s) for year_y
  # mosaic stacks if the subbasin intersects multiple BCRs
  create_stack_for_subbasin <- function(key, year) {
    
    cache_id <- paste(key, year, sep = "_")
    
    if (exists(cache_id, envir = stack_cache)) {
      return(get(cache_id, envir = stack_cache))
    }
    
    combo_poly <- combo_polys[[key]]
    
    # convert "can10_can40" to c("can10.tif", "can40.tif")
    bcr_keys_split <- paste(strsplit(key, "_", fixed = TRUE)[[1]], paste0(year, ".tif"), sep="_")
    
    # load each BCR × year stack and crop to relevant subbasins
    print("loading covariate stacks")
    stacks <- lapply(bcr_keys_split, function(r) {
      
      terra::rast(file.path(root, "gis", "stacks", r)) |> 
      terra::crop(x = _, y = combo_poly) |> 
      terra::mask(x = _, mask = combo_poly)
      
    })
    # before mosaicking, apply `as.factor` to categorical layers
    mark_categoricals <- function(stack_s, categorical_predictors) {
      keep <- intersect(categorical_predictors, names(stack_s))
      for (marked_layer in keep) stack_s[[marked_layer]] <- terra::as.factor(stack_s[[marked_layer]])
      stack_s
    }
    
    print("applying as.factor to categorical layers in covariate stacks")
    stacks <- lapply(X = stacks, FUN = mark_categoricals, categorical_predictors = categorical_predictors)
    
    # single-layer template grid for this BCR-group
    grid_tmpl <- terra::rast(
      extent = terra::align(ext(combo_poly), stacks[[1]][[1]], snap = "out"),
      resolution = terra::res(stacks[[1]]),
      crs = terra::crs(stacks[[1]])
    ) 
    

    
    # mosaic categorical layers across BCR stacks
    # for a given categorical predictor:
    # 1) identify which stacks it exists in
    # 2) identify the layers it exists in
    # 3) mosaic those layers
    
    
    
    # in the case of a sub-basin intersecting multiple BCRs
    # mosaic BCR x year stacks together to create a single training space
    print("mosaic-ing covariate layers for sub-basins that intersect multiple BCRs")
    if (length(stacks) == 1) {
      
      stacks_mos <- stacks[[1]]
      
    } else {
      
      # extract all covariate names from all stacks
      all_names <- unique(unlist(lapply(stacks, names)))
      
      # mosaic by layer
      mosaic_layers <- lapply(all_names, function(name) {
        
        # collect layer_l from each stack where present
        rlist <- lapply(stacks, function(stack_s) {
          if (name %in% names(stack_s)) stack_s[[name]] else NULL
        })
        
        # remove NULLs if they exist
        rlist <- Filter(Negate(is.null), rlist)
        
        # if only one stacks has this layer, return it as-is
        # still, align result to the common grid_tmpl
        if (length(rlist) == 1) {
          out <- rlist[[1]]
          if (name %in% categorical_predictors) {
            out <- terra::resample(out, grid_tmpl, method = "near")
          } else {
            out <- terra::resample(out, grid_tmpl, method = "bilinear")
          }
          return(out)
        }
        
        
        # align all contributors to the BCR-group grid
        if (name %in% categorical_predictors) {
          rlist <- lapply(rlist, function(r)
            if (!terra::compareGeom(r, grid_tmpl, stopOnError = FALSE))
              terra::resample(r, grid_tmpl, method = "near") else r)
        } else {
          rlist <- lapply(rlist, function(r)
            if (!terra::compareGeom(r, grid_tmpl, stopOnError = FALSE))
              terra::resample(r, grid_tmpl, method = "bilinear") else r)
        }
        
        
        # stack contributors for this layer (extents should already match after crop/mask)
        layer_stack <- do.call(c, rlist)
        
        # mosaic layers shared between stacks
        if (name %in% categorical_predictors) {
        
          out <- terra::app(layer_stack, fun = function(x) terra::modal(x, ties = "first", na.rm = TRUE))
         
          
         } else {
            
          out <- terra::app(layer_stack, fun = mean, na.rm = TRUE)
      
         }
        
        out 
        
        }) # finish splitting and mosaicing multi-BCR groups
      
      # combine mosaics into a single stack 
      stacks_mos <- do.call(c, mosaic_layers)
      names(stacks_mos) <- all_names
    } # finish mosaicking
  
assign(cache_id, stacks_mos, envir = stack_cache)
return(subbasin_stats)
return(stacks_mos)
     
} # close backfill_mines_cont()






