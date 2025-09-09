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
ia_dir <- file.path(root, "data", "Extras", "sandbox_data", "impactassessment_sandbox")




# ----------------------------------------
# define predictor and response variables
# store predictor metadata as a reference
predictor_metadata <-
  dplyr::tibble(BAMexploreR::predictor_metadata) |>
  dplyr::filter(version == "v5") |>
  dplyr::select(predictor, definition, predictor_class) |>
  dplyr::mutate(across('predictor', str_replace, 'Year', 'year')) |>
  dplyr::mutate(across('predictor', str_replace, 'Method','method'))

soil_covs <- tibble(predictor = c("cec_0-5cm_mean_1000", "cec_100-200cm_mean_1000",
                                  "cec_15-30cm_mean_1000", "cec_30-60cm_mean_1000",  
                                  "cec_5-15cm_mean_1000", "cec_60-100cm_mean_1000", 
                                  "soc_0-5cm_mean_1000",  "soc_100-200cm_mean_1000",
                                  "soc_15-30cm_mean_1000", "soc_30-60cm_mean_1000",  
                                  "soc_5-15cm_mean_1000", "soc_60-100cm_mean_1000"),
                    predictor_class = rep("Soil Properties", 12))


actually_biotic_what <- c("StandardDormancy_1km", "StandardGreenup_1km", "Peatland_5x5", "Peatland_1km")
actually_biotic_df <- tibble(predictor = actually_biotic_what,
                             predictor_class = c("Annual Climate", "Annual Climate", "Wetland", "Wetland"))

abiotic_vars <-
  predictor_metadata |> 
  dplyr::filter(predictor_class %in% c("Annual Climate", "Climate Normals", "Topography", "Wetland", "Disturbance")) |> 
  tibble::add_row(predictor = "CAfire", predictor_class ="Time Since Disturbance") |> 
  dplyr::filter(!(predictor %in% actually_biotic_what)) |> 
  dplyr::bind_rows(soil_covs)

biotic_vars <-
  predictor_metadata |> 
  dplyr::filter(!(predictor_class %in% c(abiotic_vars$predictor_class, "Time", "Method"))) |> 
  dplyr::bind_rows(actually_biotic_df)
  



# -------------------------------------
# train models per subbasin per year
backfill_mines_cont <- function(
    year,
    categorical_predictors = c("ABoVE_1km", "NLCD_1km","MODISLCC_1km",
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
    
 
 
 ###
 # train models per mining patch
 ###
 
 # for each patch:
 # 1.identify the subbasin polygon the patch falls into (using hydrobasins_subset_2020.gpkg)
 # 2.crop/mask the covariate stack to that subbasin.
 # 3.apply both: non_mine_mask (exclude mine footprint) and # 
 #            "low HF mask"(i.e. keep pixels where CanHF_1km <= threshold). # 
 # 4.extract training data as a data.frame (covariates as predictors, biotic layer as response).
 #  # 5.train a model
 


 
 # map each mine patch to a subbasin ID
 all_subbasins <- terra::vect(file.path(ia_dir, "hydrobasins_subset_2020.gpkg"))
 
 # centroids preserve row order, so they align with mines_polygons_y
 mine_cent <- terra::centroids(mines_polygons_y)
 if (!terra::same.crs(mine_cent, all_subbasins)) {
   mine_cent <- terra::project(mine_cent, all_subbasins)
 }
 
 # 2) logical matrix: rows = patches, cols = subbasins
 hits <- terra::relate(mine_cent, all_subbasins, relation = "within")
 
 # 3) map each patch to ONE subbasin (first TRUE per row), NA if none
 has_hit <- rowSums(hits) > 0
 sb_id <- rep(NA_integer_, nrow(hits))
 sb_id[has_hit] <- max.col(hits[has_hit, , drop = FALSE], ties.method = "first")
 
 # 4) tidy map (lengths match by construction)
 map_df <- data.frame(patch_id = mines_polygons_y$patch_id, sb_id = sb_id)
 
 #unique subbasins to train
 sb_ids <- sort(unique(map_df$sb_id))
 
 # set low HF to <1 
 # from Hirsh-Pearson: we found that 82% of Canada’s land areas had a 
 # HF < 1 and therefore were considered intact
 lowHF_mask <- 
   terra::lapp(stack_y[["CanHF_1km"]],
                           \(v) ifelse(!is.finite(v), NA, ifelse(v < 1, 1, NA))) |> 
   as.factor()
   
 
 # plot (exported at 1000 x 751)
 # ggplot() +
 #   geom_spatraster(data=lowHF_mask) +
 #   geom_spatvector(data=all_subbasins, fill=NA) +
 #   scale_fill_manual(
 #     values = c("1" = "#CC79A7"),   # colour for lowHF pixels
 #     na.value = NA,              # transparent for NA
 #     guide = "none"              # remove legend if you don’t need it
 #   ) +
 #   theme_minimal() +
 #   coord_sf(crs = crs(all_subbasins))
   
 # problem: many subbasins don't have any pixels with HF < 1
 # counts_df <- terra::extract(
 #   lowHF_mask,
 #   all_subbasins,
 #   fun   = function(x) sum(!is.na(x)),  # count non-NA cells
 #   ID    = TRUE
 # )
 # 
 # mean(counts_df$CanHF_1km)
 # # 6169.381
 # 
 # hist(counts_df$CanHF_1km, breaks=100)
 
 
 
}




