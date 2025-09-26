# ---
# title: Impact Assessment: train models per subbasin
# author: Mannfred Boehm
# created: May 7, 2025
# ---

library(BAMexploreR)
library(furrr)
library(terra)
library(tidyterra)
library(tidyverse)
library(xgboost)


# set root path
root <- "G:/Shared drives/BAM_NationalModels5"
ia_dir <- file.path(root, "data", "Extras", "sandbox_data", "impactassessment_sandbox")

bam_boundary <- terra::vect(file.path(root, "Regions", "BAM_BCR_NationalModel_UnBuffered.shp"))




# PART I: define model covariates ----------------------------------------

# define predictor and response variables
# store predictor metadata as a reference
predictor_metadata <-
  dplyr::tibble(BAMexploreR::predictor_metadata) |>
  dplyr::filter(version == "v5") |>
  dplyr::select(predictor, definition, predictor_class) |>
  dplyr::mutate(across('predictor', str_replace, 'Year', 'year')) |>
  dplyr::mutate(across('predictor', str_replace, 'Method','method'))

# define soil covariate names
soil_covs <- tibble(predictor = c("cec_0-5cm_mean_1000", "cec_100-200cm_mean_1000",
                                  "cec_15-30cm_mean_1000", "cec_30-60cm_mean_1000",  
                                  "cec_5-15cm_mean_1000", "cec_60-100cm_mean_1000", 
                                  "soc_0-5cm_mean_1000",  "soc_100-200cm_mean_1000",
                                  "soc_15-30cm_mean_1000", "soc_30-60cm_mean_1000",  
                                  "soc_5-15cm_mean_1000", "soc_60-100cm_mean_1000"),
                    predictor_class = rep("Soil Properties", 12))

# convert some abiotic variables as biotic variables
actually_biotic_what <- c("StandardDormancy_1km", "StandardGreenup_1km", "Peatland_5x5", "Peatland_1km")
actually_biotic_df <- tibble(predictor = actually_biotic_what,
                             predictor_class = c("Annual Climate", "Annual Climate", "Wetland", "Wetland"))
# index abiotic variables
abiotic_vars <-
  predictor_metadata |> 
  dplyr::filter(predictor_class %in% c("Annual Climate", "Climate Normals", "Topography", "Wetland", "Disturbance")) |> 
  tibble::add_row(predictor = "CAfire", predictor_class ="Time Since Disturbance") |> 
  dplyr::filter(!(predictor %in% actually_biotic_what)) |> 
  dplyr::bind_rows(soil_covs)

# index biotic variables
biotic_vars <-
  predictor_metadata |> 
  dplyr::filter(!(predictor_class %in% c(abiotic_vars$predictor_class, "Time", "Method"))) |> 
  dplyr::bind_rows(actually_biotic_df)

# re-order biotic variables 
neworder <- readRDS(file = file.path(ia_dir, "biotic_variable_hierarchy.rds"))
biotic_vars <- biotic_vars[match(neworder, biotic_vars$predictor), ]



# PART II: import ancillary data ----------------------------------------

# import subbasins polygon
all_subbasins_subset <- vect(file.path(ia_dir, "hydrobasins_masked_merged_subset.gpkg"))

# import industry footprint polygon
combined_poly <- terra::vect(file.path(ia_dir, "combined_industry_footprint.gpkg")) 
  
# import low hf layer
lowhf_mask <- terra::rast(file.path(ia_dir, "CanHF_1km_lessthan1.tif"))




# PART III: build model training function ----------------------------------------

# train models per subbasin per year
# soil_cache: pass an env() to reuse across years; if NULL a fresh one is used
source(file.path(getwd(), "Rscripts", "may_restart", "train_and_backfill_subbasin_s.R"))
plan(multisession, workers = 3)  


train_and_backfill_per_year <- function(
    year,
    all_subbasins_subset,
    categorical_responses = c("ABoVE_1km", "NLCD_1km","MODISLCC_1km",
                               "MODISLCC_5x5","SCANFI_1km","VLCE_1km"),
    soil_cache = NULL,   
    quiet = FALSE){
  
 if (is.null(soil_cache)) soil_cache <- new.env(parent = emptyenv())

 ###
 # import covariate stack, soil stack, and low HF layer
 ###
 
 # import pre-mosaiced covariate stack for year_y
 stack_y <- terra::rast(file.path(ia_dir, sprintf("covariates_mosaiced_%d.tif", year)))
 
 # ensure categoricals are factors
 cats_present <- intersect(categorical_responses, names(stack_y))
 for (nm in cats_present) stack_y[[nm]] <- terra::as.factor(stack_y[[nm]])
 message("cateogrical variables converted to factors")
 
 # import soils stack and add to covariate stack
 soil_src  <- terra::rast(file.path(ia_dir, "isric_soil_covariates_masked.tif"))
 
 # identify the geometry of the current year’s covariate stack
 geom_key <- paste(nrow(stack_y), ncol(stack_y), as.character(terra::ext(stack_y)), crs(stack_y), sep="|")
 
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

 stack_y <- c(stack_y, soil_aligned)          
 message("successfully added soil stack to covariate stack for year ", year)                        
                          
 
 # project low hf layer to current stack
 lowhf_mask_y <- terra::project(x=lowhf_mask, y=stack_y, method = "near")
 
 # save so that we pass directories, not objects, to furrr
 cache_dir <- file.path(ia_dir, "cache_year_stacks")
 dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
 
 stack_path <- file.path(cache_dir, sprintf("stack_y_%d.tif", year))
 lowhf_path <- file.path(cache_dir, sprintf("lowhf_mask_y_%d.tif", year))
 terra::writeRaster(stack_y, stack_path, overwrite = TRUE)
 terra::writeRaster(lowhf_mask_y, lowhf_path, overwrite = TRUE)
 
 
 # use existing on-disk files for vectors (already have gpkg paths)
 subbasins_path <- file.path(ia_dir, "hydrobasins_masked_merged_subset.gpkg")
 combined_poly_path <- file.path(ia_dir, "combined_industry_footprint.gpkg")
 
 
 ###
 # train and backfill models per subbasin
 ###
 
 # apply model fitting and backfilling over all subbasins
 res <- future_map(
   seq_len(nrow(all_subbasins_subset)),
   \(i) train_and_backfill_subbasin_s(
     subbasin_index        = i,
     year                  = year,
     stack_y_path          = stack_path,
     lowhf_mask_y_path     = lowhf_path,
     abiotic_vars          = abiotic_vars,
     biotic_vars           = biotic_vars,
     categorical_responses = categorical_responses,
     combined_poly_path    = combined_poly_path,    
     subbasins_path        = subbasins_path,
     ia_dir                = ia_dir,
     neworder              = neworder,
     quiet                 = quiet
   ),
   .options  = furrr::furrr_options(packages = c("terra","xgboost","dplyr")),
   .progress = TRUE)
 
 invisible(res)
 
} # close train_and_backfill_per_year()

train_and_backfill_per_year(year = 2020, all_subbasins_subset = all_subbasins_subset)



