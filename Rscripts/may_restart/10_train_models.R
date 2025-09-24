# ---
# title: Impact Assessment: train models per subbasin
# author: Mannfred Boehm
# created: May 7, 2025
# ---

library(BAMexploreR)
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


# for testing:
year <- 2020
subbasin_index <- 1
b<-biotic_cols[1]

# train models per subbasin per year
# soil_cache: pass an env() to reuse across years; if NULL a fresh one is used
mirai::daemons(3)
backfill_mines_cont <- function(
    year,
    all_subbasins_subset = all_subbasins_subset,
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
                          
 
 # project low hf layer
 lowhf_mask <- terra::project(x=lowhf_mask, y=stack_y, method = "near")
 
 ###
 # train models per subbasin
 ###
 # 1.crop/mask the covariate stack to subbasin s
 # 2. apply the low HF mask
 # 3. extract training data as a data.frame (covariates as predictors, biotic layer as response)
 # 4.train a model
 
 train_subbasin_s <- function(
    subbasin_index, 
    year, stack_y, 
    lowhf_mask, 
    abiotic_vars = abiotic_vars, 
    biotic_vars  = biotic_vars) {
   
   # isolate subbasin s
   subbasin_s <- all_subbasins_subset[subbasin_index]
   
   # crop covariate stack to subbasin
   cov_s <- 
     stack_y |>
     terra::crop(x = _, y = subbasin_s) |>
     terra::mask(x = _, mask = subbasin_s)
   
   # align low-HF to cov_s
   lowhf_mask_s <- terra::crop(lowhf_mask, cov_s, snap = "near")
   
   # for training: select pixels in covariate stack with low hf
   cov_train_s <- terra::mask(x = cov_s, mask = lowhf_mask_s)
   
   # convert covariate stack to a dataframe for modelling
   df_train <- 
     terra::as.data.frame(cov_train_s, xy=TRUE) |> 
     dplyr::as_tibble() |> 
     dplyr::rename(easting=x, northing=y) |> 
     dplyr::mutate(
       easting  = as.numeric(scale(easting)),
       northing = as.numeric(scale(northing)))
   
   # define predictors and responses
   abiotic_cols <- intersect(names(df_train), abiotic_vars$predictor)
   biotic_cols <- intersect(names(df_train), biotic_vars$predictor)
   biotic_cols <- na.omit(biotic_cols[match(neworder, biotic_cols)])
   
   # set output directory for current year
   out_dir <- file.path(ia_dir, "xgboost_models",
                        sprintf("year=%d", year),
                        sprintf("subbasin=%s", subbasin_index))
   
   dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
   
   # parallelize model fitting per biotic feature within subbasin_s
   purrr::map(.x = biotic_cols, .f = purrr::in_parallel(\(b, df_train, abiotic_cols, out_dir, categorical_responses){
   
     # check if continuous or categorical
     if (!(b %in% categorical_responses)){
         
       # keep only the rows where the response is observed
       idx <- which(!is.na(df_train[[b]]))
       
       # log transform the response to prevent negative values in predictions
       y <- log(as.numeric(df_train[[b]][idx]) + 1)
       
        # predictors (abiotic + coords, already scaled in df_train)
       X <- dplyr::select(df_train[idx, , drop = FALSE], all_of(abiotic_cols), easting, northing)
       X <- X[, colSums(!is.na(X)) > 0, drop = FALSE] # drop columns with all NAs
         
       # create DMatrix from low HF landscape
       dtrain_b <- xgboost::xgb.DMatrix(data = as.matrix(X), label = y, missing = NA)
       params_b <- xgboost::xgb.params(objective = "reg:squarederror", eval_metric = "rmse", max_depth = 3, eta = 0.2)
       
       # the `test_rmse_mean` metrics in $evaluation_log represent the average performance 
       # on the 20% holdout in each of the 5 folds
       cv_b <-xgboost::xgb.cv(params = params_b, data = dtrain_b, nrounds=5000, nfold=5, early_stopping_rounds = 50, verbose=FALSE)
       
       # fit model
       model_b <- xgboost::xgb.train(params = params_b, data = dtrain_b, nrounds = cv_b$early_stop$best_iteration, verbose = 0)
       
       # save model
       saveRDS(list(model = model_b, cv_log = cv_b$evaluation_log), file.path(out_dir, sprintf("model_%s_xgb.rds", b)))
       
     } else { # if b is categorical, use the following protocol
       
       # keep only the rows where the response is observed
       idx <- which(!is.na(df_train[[b]]))
       
       # multiclass xgboost labels need to start at 0
       # response as factor with compact levels, then 0..K-1 labels
       y_fac <- droplevels(df_train[[b]][idx])
       K     <- nlevels(y_fac)
       y <- as.integer(y_fac) - 1
       
       # predictors (abiotic + coords), drop all-NA cols and align rows
       X <- dplyr::select(df_train[idx, , drop = FALSE], dplyr::all_of(abiotic_cols), easting, northing)
       X <- X[, colSums(!is.na(X)) > 0, drop = FALSE]
      
       dtrain_b <- xgboost::xgb.DMatrix(data = X, label = y, missing = NA)
       params_b <- xgboost::xgb.params(objective = "multi:softprob", eval_metric = "mlogloss", max_depth = 3, eta = 0.2, num_class = K)
       cv_b <-xgboost::xgb.cv(params = params_b, data = dtrain_b, nrounds=5000, nfold=5, early_stopping_rounds = 50, verbose=FALSE)
       model_b <- xgboost::xgb.train(params = params_b, data = dtrain_b, nrounds = cv_b$early_stop$best_iteration, verbose = 0)
       
       saveRDS(list(model = model_b, cv_log = cv_b$evaluation_log), file.path(out_dir, sprintf("model_%s_xgb.rds", b)))
       
      } # close if/else (continuous vs categorical training)
    } # close training function inside of in_parallel
    
    # give arguments to purrr::map
    ),  df_train = df_train, abiotic_cols = abiotic_cols, out_dir = out_dir, categorical_responses = categorical_responses) # close purrr::map for training over all biotic features
  
  } # close train_subbasin_s()
     
 # apply model fitting over all subbasins
 lapply(X = seq_len(nrow(all_subbasins_subset)), FUN = train_subbasin_s,
        year = year, stack_y = stack_y, lowhf_mask = lowhf_mask,
        abiotic_vars = abiotic_vars, biotic_vars = biotic_vars)
 
 
 
 ###
 # backfill the industry footprints in subbasin_s
 ###
 
 # for backfilling: find industry pixels in subbasin_s rasterize onto cov_s grid
 industry_s <- terra::crop(combined_poly, subbasin_s)
 industry_mask_s <- terra::rasterize(industry_s, cov_s[[1]], field = 1, background = NA, touches = TRUE)     
 
 # indices to backfill = cells covered by footprint polygons
 backfill_idx <- which(!is.na(terra::values(industry_mask_s)))
 
 df_backfill <- 
   terra::as.data.frame(cov_s, xy = TRUE, na.rm = FALSE) |>
   dplyr::as_tibble() |>
   dplyr::rename(easting = x, northing = y) |>
   dplyr::mutate(
     easting  = as.numeric(scale(easting)),
     northing = as.numeric(scale(northing))) |> 
   dplyr::slice(backfill_idx)
} # close backfill_mines_cont() function

# use model to predict biotic feature b in industry footprint pixels (i.e. backfilling)
# also, inverse log the predictions to convert back to original scale
predictions_b <- exp(predict(model_b, newdata = )) - 1 

# after all training is done:
mirai::daemons(0)

