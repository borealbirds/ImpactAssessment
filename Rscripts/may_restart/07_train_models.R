# ---
# title: Impact Assessment: train models on watersheds surrounding mines
# author: Mannfred Boehm
# created: May 7, 2025
# ---

library(tidyverse)
library(terra)

# set root path
root <- "G:/Shared drives/BAM_NationalModels5"
ia_dir <- file.path(root, "data", "Extras", "sandbox_data", "impactassessment_sandbox")



# train models per subbasin per year
backfill_mines_cont <- function(
    year,
    categorical_predictors = c("ABoVE_1km","method","NLCD_1km","MODISLCC_1km",
                               "MODISLCC_5x5","SCANFI_1km","VLCE_1km"),
    soil_cache = NULL,   # pass an env() to reuse across years; if NULL a fresh one is used
    quiet = FALSE
){
  
  if (is.null(soil_cache)) soil_cache <- new.env(parent = emptyenv())

  
 ###
 # load pre-mosaiced covariate stack for year_y
 ###
 stack_y <- terra::rast(file.path(ia_dir, sprintf("covariates_mosaiced_%d.tif", year)))
    
 # ensure categoricals are factors
 cats_present <- intersect(categorical_predictors, names(stack_y))
 for (nm in cats_present) stack_y[[nm]] <- terra::as.factor(stack_y[[nm]])
    
 # identify and store HF layers
 hf_layers <- intersect(c("CanHF_1km", "CanHF_5x5"), names(stack_y))
 hf_stack <- if (length(hf_layers)) stack_y[[hf_layers]] else NULL

 
 
 ###
 # import soils stack and add to covariate stack
 ###
 
 soil_src  <- terra::rast(file.path(ia_dir, "isric_soil_covariates_masked.tif"))
 
 # identifies the geometry of the current year’s covariate stack
 geom_key <- paste(nrow(stack_y), ncol(stack_y),
                   as.character(terra::ext(stack_y)),
                   terra::crs(stack_y), sep="|")
 
 # on the first loop (year) we import and align the soil stack
 # on subsequent iterations (years) we save time by using `get()`
 if (exists(geom_key, envir = soil_cache)) {
   soil_aligned <- get(geom_key, envir = soil_cache)
 } else {
   soil_aligned <- terra::project(
     soil_src, stack_y, method = "bilinear",
     filename = tempfile(fileext = ".tif"), overwrite = TRUE
   ) |>
     terra::crop(x= _, y = stack_y) |>
     terra::mask(x= _, mask = stack_y[[1]])
   assign(geom_key, soil_aligned, envir = soil_cache)
 }

                       
 names(soil_aligned) <- paste0("SOIL_", names(soil_aligned))
 stack_y <- c(stack_y, soil_aligned)          
 message("successfully added soil stack to covariate stack for year ", year)                        
                          
                          
 ###
 # import mines and build polygons/patches
 ###
    
 mines_y <- terra::rast(file.path(ia_dir, sprintf("mincan_mines_%d_masked.tif", year)))
 message("successfully imported mines for year ", year)
    
 # ensure mines are in same CRS/grid as covariates (use nearest since it's categorical)
 if (!terra::compareGeom(mines_y, stack_y, stopOnError = FALSE)) {
      mines_y <- terra::project(mines_y, stack_y, method = "near")
    }
    
 # generate mine patches (rook’s case) and convert to polygons
 patches_y <- terra::patches(mines_y, directions = 4, zeroAsNA = TRUE)
 mines_polygons_y <- terra::as.polygons(patches_y, dissolve = TRUE, na.rm = TRUE)
 names(mines_polygons_y) <- "patch_id"
    
    
 ###
 # rasterize polygon footprint onto the covariate grid to build a reusable mask
 ###
    
 # mine_mask: 1 where mine exists, NA elsewhere (aligns to stack_y grid)
 mine_mask <- terra::rasterize(mines_polygons_y, stack_y[[1]], field = 1, background = NA_real_)
 message("successfully generated mine_mask for year ", year)
  
 # non_mine_mask: 1 for valid training cells (non-mine), NA for mine cells
 non_mine_mask <- terra::ifel(is.na(mine_mask), 1, NA) * (!is.na(stack_y[[1]]))
 message("successfully generated non_mine_mask for year ", year)
    
 return(list(
   non_mine_mask = non_mine_mask,
   mine_mask     = mine_mask
 ))
}


# store predictor metadata as a reference
# predictor_metadata <- 
#   dplyr::tibble(BAMexploreR::predictor_metadata) |>
#   dplyr::filter(version == "v5") |> 
#   dplyr::select(predictor, definition, predictor_class) |> 
#   dplyr::mutate(across('predictor', str_replace, 'Year', 'year')) |> 
#   dplyr::mutate(across('predictor', str_replace, 'Method','method'))

