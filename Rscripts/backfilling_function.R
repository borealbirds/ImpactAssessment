# ---
# title: National Models 5.0 - create models to predict continuous biotic features from abiotic landscape
# author: Mannfred Boehm
# created: March 19, 2025
# ---


#1. run lines 1-132 in "test_backfill_SCANFI_bcr14.R"---- 
# this will import covariate rasters, and generate a dataframe `df_bcr14_2020`from 
# the predictor covariate stack "can14_2020.tif"
# We also set the threshold for low/high human footprint by defining `q50`





# ------------------------------------------------------------


# import backfilled (bf) landclass variables as predictors
filepaths <- (c(file.path(getwd(),"data","derived_data","backfilled_rasters","BCR14_MODISLCC_1km.tif"),
                file.path(getwd(),"data","derived_data","backfilled_rasters","BCR14_SCANFI_1km.tif"),
                file.path(getwd(),"data","derived_data","backfilled_rasters","BCR14_VLCE_1km.tif")))

backfilled_landcover <- terra::rast(filepaths)
names(backfilled_landcover) <- c("MODISLCC_1km_bf", "SCANFI_1km_bf", "VLCE_1km_bf")

# convert backfilled landclass variables to a dataframe 
df_backfilled_landcover <-
  backfilled_landcover |> 
  terra::as.data.frame(xy=TRUE) |>
  tibble::as_tibble() |> 
  dplyr::rename(lon=x, lat=y) 

# replace disturbed landcover values with backfilled values 
# (as "abiotic"predictors when backfilling continuous vegetation variables)
df_bcr14_2020 <- dplyr::left_join(df_bcr14_2020, df_backfilled_landcover, by=c("lon","lat"))

# define continuous vegetation features to backfill: 
# remove 5x5 since I can re-create these variables by scaling up from 1km after prediciton
# exclude var_class == "Landcover" as they are now predictors
biotic_vars <-
  nice_var_names |> 
  dplyr::filter(var_class %in% c("Greenup", "Biomass", "Wetland")) |> 
  dplyr::filter(!(var %in% c("WetOccur_1km", "WetOccur_5x5", "WetRecur_1km", "WetSeason_1km"))) |> # keep peatland but discard other water variables
  dplyr::filter(!grepl("5x5", var)) |> 
  dplyr::pull(var) |> 
  unique()

# drop biotic (bird) predictors that were thinned by VIF (see line 300 in "04.Stratify.R")
biotic_vars_thinned <- biotic_vars[which(biotic_vars %in% colnames(df_bcr14_2020))]

# re-order biotic (bird) predictors 
neworder <- c("Peatland_1km", "SCANFI_1km", "SCANFIclosure_1km", "SCANFIbiomass_1km", "SCANFIheight_1km", 
              "SCANFIprcC_1km", "SCANFIprcD_1km", "SCANFIBalsamFir_1km", "SCANFIBlackSpruce_1km",
              "SCANFIDouglasFir_1km", "SCANFIJackPine_1km", "SCANFILodgepolePine_1km", "SCANFIPonderosaPine_1km",
              "SCANFITamarack_1km", "SCANFIWhiteRedPine_1km", "StandardDormancy_1km", "StandardGreenup_1km")

biotic_vars_thinned <- biotic_vars_thinned[order(match(biotic_vars_thinned, neworder))]


# check the range of each of each continuous variable for negative values
# result: none of the covariates (that made it through VIF) have negative values
for (v in 1:length(biotic_vars_thinned)) {
  if (biotic_vars_thinned[v] %in% colnames(df_bcr14_2020)) {
    range_v <- range(df_bcr14_2020[,biotic_vars_thinned[v]], na.rm=T)
    covariate <- biotic_vars_thinned[v]
    df <- data.frame(covariate=covariate, range=range_v, index=v)
    print(df)
  } else NULL
}

# set predictor variables
# include landcover covariates as predictors
predictor_vars <- c(abiotic_vars, "lon", "lat", "SCANFI_1km_bf", "VLCE_1km_bf", "MODISLCC_1km_bf")


rm(nice_var_names, biotic_vars); gc()


# define a function for predicting "restored" biotic landscape features----
backfill_landscape <- function(i){
  
  # stage a unique subset of the BCR by removing NAs for biotic feature[i]
  df_bcr14_2020_i <- tidyr::drop_na(df_bcr14_2020, biotic_vars_thinned[i]) 
  
  # identify pixels with "high" and "low" human footprint
  # for biotic feature[i]
  CanHF_1km_present <- dplyr::filter(df_bcr14_2020_i, CanHF_1km > q50) 
  CanHF_1km_absent <- dplyr::filter(df_bcr14_2020_i, CanHF_1km <= q50)
  
  # log transform covariates to prevent negative values
  label_i <- log(as.numeric(CanHF_1km_absent[[biotic_vars_thinned[i]]]) + 1)
  
  # create DMatrix from low HF landscape
  # reminder: `predictor_vars` is `c(abiotic_vars, "lon", "lat", "SCANFI_1km_bf", "VLCE_1km_bf", "MODISLCC_1km_bf")`
  dtrain_i <- xgboost::xgb.DMatrix(data = as.matrix(CanHF_1km_absent[,predictor_vars]), label = label_i)
  
  params_i <- xgboost::xgb.params(objective = "reg:squarederror", eval_metric = "rmse", max_depth = 3, eta = 0.2)
  
  # the `test_rmse_mean` metrics in $evaluation_log represent the average performance 
  # on the 20% holdout in each of the 5 folds
  cv_i <-xgboost::xgb.cv(params = params_i, data = dtrain_i, nrounds=10000, nfold=5, early_stopping_rounds = 10, verbose=FALSE)
  
  model_i <- xgboost::xgb.train(params = params_i, data = dtrain_i, nrounds = cv_i$early_stop$best_iteration, verbose = 0)
  
  #log model info
  backfill_datalog[[i]] <- c(cv_i, model_i)
 
  # predict biotic feature[i] at "high" HF areas (i.e. backfilling)
  # also, inverse log the predictions to convert back to original scale
  predictions_i <- exp(predict(model_i, newdata = CanHF_1km_present[,abiotic_vars])) - 1 
  
  # then, union the landscape by adding back the locations with low HF which 
  # weren't subject to backfilling (via `add_row`)
  new_landscape_i <- 
    dplyr::tibble(!!biotic_vars_thinned[i] := predictions_i, 
                  lon = CanHF_1km_present$lon, 
                  lat = CanHF_1km_present$lat) |> 
    dplyr::add_row(CanHF_1km_absent[,c(biotic_vars_thinned[i], "lon", "lat")]) 
  
  
  # rasterize predictions
  # eventually we'll stack for `terra::predict()` of bird densities
  prediction_raster_i <- 
    terra::vect(new_landscape_i, geom = c("lon", "lat"), crs = crs(empty_raster)) |>
    terra::rasterize(y = empty_raster, field = biotic_vars_thinned[i])
  
  plot(prediction_raster_i, main=paste("backfilled", biotic_vars_thinned[i]))
  lines(bcr14_boundary, col="black", lwd=1)
  
  terra::writeRaster(prediction_raster_i, file.path(getwd(), "data", "derived_data", "backfilled_rasters", paste0("BCR14_", biotic_vars_thinned[i], ".tif")), overwrite=TRUE)
  
  
  print(paste("* created backfilled raster for: ", biotic_vars_thinned[i], " *"))
  return(prediction_raster_i)
  
}



#11. apply backfilling to all biotic covariates----

# define empty raster for placing backfilled features into
empty_raster <- terra::rast(extent=ext(stack_bcr14_2020), crs=crs(stack_bcr14_2020), resolution=1000)

# define empty data log
backfill_datalog <- list()

backfilled_rasters <- lapply(X=seq_along(biotic_vars_thinned), FUN=backfill_landscape)

# combine into one raster stack
backfilled_rasters2 <- do.call(c, backfilled_rasters)
names(backfilled_rasters2) <- biotic_vars_thinned

# set water bodies to zero
test <- terra::cover(backfilled_rasters2$SCANFIBalsamFir_1km, stack_bcr14_2020$SCANFIBalsamFir_1km)


saveRDS(backfilled_rasters2, file="C:/Users/mannf/Proton Drive/mannfredboehm/My files/Drive/boreal_avian_modelling_project/ImpactAssessment/data/derived_data/rds_files/backfilled_rasters.rds")
backfilled_rasters2 <- readRDS(file="C:/Users/mannf/Proton Drive/mannfredboehm/My files/Drive/boreal_avian_modelling_project/ImpactAssessment/data/derived_data/rds_files/backfilled_rasters.rds")




# then:
# 1. repeat predictions for all 18 biotic variables
# 2. rasterize predictions and stack
# 3. import bird models (`b.i`) for CAWA at year 2020
# 4. use terra::predict with new *biotic* stack, og *abiotic* stack, and bird model (`b.i`) to estimate  new bird densities

