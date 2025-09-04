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
    years,
    ia_dir,
    categorical_predictors = c("ABoVE_1km","method","NLCD_1km","MODISLCC_1km",
                               "MODISLCC_5x5","SCANFI_1km","VLCE_1km"),
    quiet = FALSE
) {
  out <- vector("list", length(years))
  names(out) <- as.character(years)
  
  for (i in seq_along(years)) {
    y <- years[i]
    
    ###
    # load pre-mosaiced covariate stack for year_y
    ###
    stack_y <- terra::rast(covariate_stack_path(ia_dir, y))
    
    # ensure categoricals are factors
    cats_present <- intersect(categorical_predictors, names(stack_y))
    for (nm in cats_present) stack_y[[nm]] <- terra::as.factor(stack_y[[nm]])
    
    # index HF layers
    hf_layers <- intersect(c("CanHF_1km", "CanHF_5x5"), names(stack_y))
    hf_stack <- if (length(hf_layers)) stack_y[[hf_layers]] else NULL
    
    
    ###
    # import mines and build polygons/patches
    ###
    
    mines_y <- terra::rast(file.path(ia_dir, sprintf("mincan_mines_%d_masked.tif", y)))
    
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
    mine_mask <- if (nrow(mines_polygons_y) > 0) {
      terra::rasterize(mines_polygons_y, stack_y[[1]], field = 1, background = NA_real_)
    } else {
      stack_y[[1]] * NA_real_
    }
    
    # non_mine_mask: 1 for valid training cells (non-mine), NA for mine cells
    non_mine_mask <- terra::ifel(is.na(mine_mask), 1, NA)
    
    # package result (stop here; low-HF + training will use these)
    out[[i]] <- list(
      year           = y,
      stack          = stack_y,             # full mosaiced covariate stack (already in subbasin extent)
      hf_layers      = hf_stack,          # optional HF layer(s) to threshold later
      mines_rast     = mines_y,           # year-matched mines raster
      mines_patches  = patches_y,         # patches (SpatRaster)
      mines_polygons = mines_polygons_y,  # polygons (SpatVector)
      non_mine_mask  = non_mine_mask,     # 1 = keep, NA = mine footprint
      notes          = "Ready for low-HF masking (AND non_mine_mask) + model training."
    )
    
    if (!quiet) {
      message("[", y, "] loaded ",
              basename(cov_path), " and mines (",
              nlyr(stack_y), " layers; mine polys: ", nrow(mines_polygons_y), ")")
    }
  }
  
  # drop NULL entries for missing years
  Filter(Negate(is.null), out)
}


# store predictor metadata as a reference
# predictor_metadata <- 
#   dplyr::tibble(BAMexploreR::predictor_metadata) |>
#   dplyr::filter(version == "v5") |> 
#   dplyr::select(predictor, definition, predictor_class) |> 
#   dplyr::mutate(across('predictor', str_replace, 'Year', 'year')) |> 
#   dplyr::mutate(across('predictor', str_replace, 'Method','method'))

