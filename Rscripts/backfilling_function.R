

# ------------------------------------------------------------

# define empty raster for placing backfilled features into
empty_raster <- 
  terra::rast(extent=ext(stack_bcr14_2020), crs=crs(stack_bcr14_2020), resolution=1000)


# define a function for predicting "restored" biotic landscape features----
backfill_landscape <- function(i){
  
  # drop biotic predictors that were thinned by VIF (see line 300 in "04.Stratify.R")
  biotic_vars_thinned <- biotic_vars[which(biotic_vars %in% colnames(df_bcr14_2020))]
  
  # stage a unique spatial subset of the BCR by removing NAs for biotic feature[i]
  df_bcr14_2020_i <- tidyr::drop_na(df_bcr14_2020, biotic_vars_thinned[i]) 
  
  # identify pixels with "high" and "low" human footprint
  # for biotic feature[i]
  CanHF_1km_present <- dplyr::filter(df_bcr14_2020_i, CanHF_1km > q50) 
  CanHF_1km_absent <- dplyr::filter(df_bcr14_2020_i, CanHF_1km <= q50)
  
  # model biotic ~ abiotic using boosted regression trees
  xgbm_formula <- reformulate(termlabels = c(abiotic_vars, "lon", "lat"), 
                              response = biotic_vars_thinned[i])
  
  model_i <- 
    gbm::gbm(formula = xgbm_formula,
             data = CanHF_1km_absent,
             distribution="multinomial",
             n.trees = 2000,
             interaction.depth = 3,            
             shrinkage = 0.01,
             cv.folds = 5)
  
  # predict biotic feature[i] at "high" HF areas (i.e. backfilling)
  predictions_i <- gbm::predict.gbm(object=model_i, newdata=CanHF_1km_present, type="response") 
  
  # then, union the landscape by adding back the locations with low HF which 
  # weren't subject to backfilling (via `add_row`)
  new_landscape <- 
    dplyr::tibble(!!biotic_vars_thinned[i] := predictions_i, 
                  lon = CanHF_1km_present$lon, 
                  lat = CanHF_1km_present$lat) |> 
    dplyr::add_row(CanHF_1km_absent[,c(biotic_vars_thinned[i], "lon", "lat")]) 
  
  
  # rasterize predictions
  # eventually we'll stack for `terra::predict()` of bird densities
  predictions_raster_i <- 
    terra::vect(new_landscape, geom = c("lon", "lat"), crs = crs(empty_raster)) |>
    terra::rasterize(y = empty_raster, field = biotic_vars_thinned[i], fun=mean)
  
  print(paste("* created backfilled raster for: ", biotic_vars[i], " *"))
  return(predictions_raster_i)
  
}




#11. apply backfilling to all biotic covariates----
tmpcl <- clusterExport(cl, c("backfill_landscape", "abiotic_vars", "biotic_vars", "crs", "empty_raster"))

# this runs pretty quick on local for a single BCR x species, so don't need to test on cluster
backfilled_rasters <- lapply(X=seq_along(biotic_vars), FUN=backfill_landscape)

saveRDS(backfilled_rasters, file="C:/Users/mannf/Proton Drive/mannfredboehm/My files/Drive/boreal_avian_modelling_project/ImpactAssessment/data/derived_data/rds_files/backfilled_rasters.rds")
backfilled_rasters <- readRDS(file="C:/Users/mannf/Proton Drive/mannfredboehm/My files/Drive/boreal_avian_modelling_project/ImpactAssessment/data/derived_data/rds_files/backfilled_rasters.rds")


# [1] "* created backfilled raster for:  SCANFIWhiteRedPine_5x5  *"
# Error in gbm.fit(x = x, y = y, offset = offset, distribution = distribution,  : 
#                    The data set is too small or the subsampling rate is too large: `nTrain * bag.fraction <= 2 * n.minobsinnode + 1`
#                  In addition: Warning messages:
#                    1: In doTryCatch(return(expr), name, parentenv, handler) :
#                    display list redraw incomplete

# stack rasters
stack_from_list <- function(raster_list){
  
  raster_list[i]
  
}

backfilled_stack <- lapply(X=seq_along(backfilled_rasters), FUN=stack_from_list)
# then:
# 1. repeat predictions for all 18 biotic variables
# 2. rasterize predictions and stack
# 3. import bird models (`b.i`) for CAWA at year 2020
# 4. use terra::predict with new *biotic* stack, og *abiotic* stack, and bird model (`b.i`) to estimate  new bird densities

