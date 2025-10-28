# ---
# title: Impact Assessment: train models per subbasin and backfill industry footprints
# author: Mannfred Boehm
# created: August 7, 2025
# ---

library(BAMexploreR)
library(furrr)
library(terra)
library(tidyverse)
library(xgboost)


# set root path -------------------------------------------------
root <- "G:/Shared drives/BAM_NationalModels5"
ia_dir <- file.path(root, "data", "Extras", "sandbox_data", "impactassessment_sandbox")



# define model covariates ----------------------------------------

# define predictor and response variables
# store predictor metadata as a reference
predictor_metadata <-
  dplyr::tibble(BAMexploreR::predictor_metadata) |>
  dplyr::filter(version == "v5") |>
  dplyr::select(predictor, definition, predictor_class) |>
  dplyr::mutate(dplyr::across('predictor', stringr::str_replace, 'Year', 'year')) |>
  dplyr::mutate(dplyr::across('predictor', stringr::str_replace, 'Method','method'))

# define soil covariate names
soil_covs <- tibble::tibble(predictor = c("cec_0-5cm_mean_1000", "cec_100-200cm_mean_1000",
                                  "cec_15-30cm_mean_1000", "cec_30-60cm_mean_1000",  
                                  "cec_5-15cm_mean_1000", "cec_60-100cm_mean_1000", 
                                  "soc_0-5cm_mean_1000",  "soc_100-200cm_mean_1000",
                                  "soc_15-30cm_mean_1000", "soc_30-60cm_mean_1000",  
                                  "soc_5-15cm_mean_1000", "soc_60-100cm_mean_1000"),
                    predictor_class = rep("Soil Properties", 12))

# convert some abiotic variables to biotic variables
actually_biotic_what <- c("StandardDormancy_1km", "StandardGreenup_1km", "Peatland_5x5", "Peatland_1km")
actually_biotic_df <- tibble::tibble(predictor = actually_biotic_what, predictor_class = c("Annual Climate", "Annual Climate", "Wetland", "Wetland"))

# define abiotic variables (V5 abiotic + CAfire + soil properties)
abiotic_vars <-
  predictor_metadata |> 
  dplyr::filter(predictor_class %in% c("Annual Climate", "Climate Normals", "Topography", "Wetland", "Disturbance")) |> 
  tibble::add_row(predictor = "CAfire", predictor_class ="Time Since Disturbance") |> 
  dplyr::filter(!(predictor %in% actually_biotic_what)) |> 
  dplyr::bind_rows(soil_covs)

# define biotic variables
biotic_vars <-
  predictor_metadata |> 
  dplyr::filter(!(predictor_class %in% c(abiotic_vars$predictor_class, "Time", "Method"))) |> 
  dplyr::bind_rows(actually_biotic_df)

# re-order biotic variables 
neworder <- readRDS(file = file.path(ia_dir, "biotic_variable_hierarchy.rds"))
biotic_vars <- biotic_vars[match(neworder, biotic_vars$predictor), ]



# define data paths ----------------------------------------

# subbasins polygon
all_subbasins_subset_path <- file.path(ia_dir, "hydrobasins_masked_merged_subset.gpkg")

# high human footprint raster
highhf_mask_path <- file.path(ia_dir, "hirshpearson", "CanHF_1km_morethan1.tif")

# low human footprint raster
lowhf_mask_path <- file.path(ia_dir, "hirshpearson", "CanHF_1km_lessthan1.tif")

# covariate stack for year y
stack_y_path <- terra::rast(file.path(ia_dir, sprintf("covariates_mosaiced_%d.tif", year)))

# training and backfilling function (subbasin level)
source(file.path(getwd(), "Rscripts", "may_restart", "08_train_and_backfill_subbasin_s.R"))


# define model training function ----------------------------------------

# train models per subbasin per year
train_and_backfill_per_year <- function(
    year,
    all_subbasins_subset_path,
    lowhf_mask_path,
    highhf_mask_path,
    stack_y_path,
    categorical_responses,
    quiet = FALSE){
  
 ###
 # import covariate stack, low HF layer for training, high HF layer for backfilling
 ###
 
 # import pre-mosaiced covariate stack for year_y
 stack_y <- terra::rast(stack_y_path)
 
 # ensure categoricals are factors
 cats_present <- intersect(categorical_responses, names(stack_y))
 for (nm in cats_present) stack_y[[nm]] <- terra::as.factor(stack_y[[nm]])

 # import low hf layer and project to current stack
 lowhf_mask <- terra::rast(lowhf_mask_path)
 lowhf_mask <- terra::project(x=lowhf_mask_y, y=stack_y, method = "near")
 
 # import high hf layer and project to current stack
 highhf_mask <- terra::rast(highhf_mask_path)
 highhf_mask <- terra::project(x=highhf_mask_y, y=stack_y, method = "near")
 
 # import subbasin boundaries
 all_subbasins_subset <- terra::vect(all_subbasins_subset_path)
 
 
 ###
 # train and backfill models per subbasin
 ###
 
 # apply model fitting and backfilling over all subbasins
 res <- train_and_backfill_subbasin_s(
     subbasin_index        = i,
     year                  = year,
     stack_y               = stack_y,
     lowhf_mask            = lowhf_mask,
     highhf_mask           = highhf_mask,
     all_subbasins_subset  = all_subbasins_subset, 
     abiotic_vars          = abiotic_vars,
     biotic_vars           = biotic_vars,
     categorical_responses = c("ABoVE_1km", "NLCD_1km","MODISLCC_1km",
                               "MODISLCC_5x5","SCANFI_1km","VLCE_1km"),
     ia_dir                = ia_dir,
     neworder              = neworder,
     quiet                 = quiet)
 
 invisible(res)
 
} # close train_and_backfill_per_year()



# train models and backfill biotic features for year y -----------------------------
train_and_backfill_per_year(year = 2020, 
                            all_subbasins_subset_path = all_subbasins_subset_path,
                            lowhf_mask_path = lowhf_mask_path,
                            highhf_mask_path = highhf_mask_path,
                            stack_y_path = stack_y_path,
                            categorical_responses = c("ABoVE_1km", "NLCD_1km","MODISLCC_1km",
                                                      "MODISLCC_5x5","SCANFI_1km","VLCE_1km"))

