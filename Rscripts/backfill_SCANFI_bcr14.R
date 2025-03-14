# ---
# title: National Models 5.0 - create models to predict biotic SCANFI classes from abiotic landscape
# author: Mannfred Boehm
# created: March 11, 2025
# ---


#1. run lines 1-161 in "test_backfill_SCANFI_bcr14.R"---- 
# this will import covariate rasters, and generate low and high HF dataframes 
# "CanHF_1km_absent" and "CanHF_1km_present"


# subset low HF dataset to only those with abiotic predictors and the biotic response
# include the response "SCANFI_1km" so that it can used for labelling 
# this is the data for training the model
predictor_vars <- c(abiotic_vars, "lon", "lat", "SCANFI_1km_rock", "SCANFI_1km_water")
CanHF_1km_absent_abiotic <- CanHF_1km_absent[,c(predictor_vars, "SCANFI_1km")]





#2. two-stage modelling----
# first, split SCANFI_1km classes between tree and non-tree vegetation
# all values for SCANFI_1km_rock and SCANFI_1km_water should be zero (checked by dplyr::count, and they are)
training_data_stage1 <- 
  CanHF_1km_absent_abiotic |> 
  dplyr::mutate(SCANFI_1km_broadclass = case_when(
    SCANFI_1km %in% c(1, 2, 4) ~ 0,  # bryoid, herbs, shrub (non-trees) 
    SCANFI_1km %in% c(5, 6, 7) ~ 1,  # broadleaf, conifer, mixed (trees)
    SCANFI_1km %in% c(3,8) ~ 2)) |>    # rock and water (determined)
  dplyr::mutate(SCANFI_1km_broadclass = factor(SCANFI_1km_broadclass, levels=c(0,1,2)))


# create DMatrix (a data structure for better speed when using xgboost)
# reminder: `predictor_vars` is `c(abiotic_vars, "lon", "lat", "SCANFI_1km_rock", "SCANFI_1km_water")`
dtrain_stage1 <- xgboost::xgb.DMatrix(data = as.matrix(training_data_stage1[,predictor_vars]), label = as.numeric(training_data_stage1$SCANFI_1km_broadclass) - 1)
params_stage1 <- xgboost::xgb.params(objective = "multi:softmax", eval_metric = "mlogloss", max_depth = 3, eta = 0.05, num_class = length(levels(training_data_stage1$SCANFI_1km_broadclass)))

# train stage 1 model (SCANFI non-trees vs trees vs abiotic)
cv_stage1 <-xgboost::xgb.cv(params = params_stage1, data = dtrain_stage1, nrounds=1200, nfold=5, early_stopping_rounds = 20)
model_stage1 <- xgboost::xgb.train(params = params_stage1, data = dtrain_stage1, nrounds = cv_stage1$early_stop$best_iteration, verbose = 0)




#3. build two multi-class models to predict (1) non-tree subclasses (bryoid, herb, shrub)----
# and (2) tree subclasses (broadleaf, conifer, mixed)
# note: the rock and water entries will be perfectly "predicted" since they are abiotic 

# filter data to non-trees entries
# SCANFI_1km is now a 4-factor column for non-treed classes (bryoids, herbs, low HF rock, shrubs)
# again, SCANFI_1km_water should all be zeros (checked by dplyr::count, and they are)
training_data_stage2_notrees <- 
  training_data_stage1 |> 
  dplyr::filter(SCANFI_1km_broadclass == 0) |> 
  dplyr::mutate(SCANFI_1km = factor(SCANFI_1km, levels = c(1, 2, 4)))

dtrain_non_tree <- xgb.DMatrix(
  data = as.matrix(training_data_stage2_notrees[, predictor_vars]), # exclude SCANFI_1km and SCANFI_1km_tree
  label = as.numeric(training_data_stage2_notrees$SCANFI_1km) - 1) # as.numeric() starts the classes at 1 but xgboost labels need to start at 0

params_non_tree <- 
  xgboost::xgb.params(objective = "multi:softmax", 
                      num_class = length(levels(training_data_stage2_notrees$SCANFI_1km)),
                      eval_metric = "mlogloss",
                      max_depth = 3,
                      eta = 0.05)

cv_notree <- xgboost::xgb.cv(params = params_non_tree, data = dtrain_non_tree, nrounds=1000, nfold=5, early_stopping_rounds = 20)
model_stage2_notree <- xgb.train(params = params_non_tree, data = dtrain_non_tree, nrounds = cv_notree$early_stop$best_iteration, verbose = 0)

# filter data to yes-trees entries
# SCANFI_1km is now a 3-factor column for yes-treed classes (broadlead, conifer, mixed)
# again, SCANFI_1km_rock and SCANFI_1km_water should all be zeros (checked by dplyr::count, and they are)
training_data_stage2_yestrees <- 
  training_data_stage1 |> 
  dplyr::filter(SCANFI_1km_broadclass == 1) |> 
  dplyr::mutate(SCANFI_1km = factor(SCANFI_1km, levels = c(5, 6, 7)))

dtrain_yes_tree <- xgb.DMatrix(
  data = as.matrix(training_data_stage2_yestrees[, predictor_vars]),
  label = as.numeric(training_data_stage2_yestrees$SCANFI_1km) - 1)

params_yes_tree <- 
  xgboost::xgb.params(objective = "multi:softmax",
                      num_class = length(levels(training_data_stage2_yestrees$SCANFI_1km)),
                      eval_metric = "mlogloss",
                      max_depth = 3,
                      eta = 0.05)

cv_yestree <- xgboost::xgb.cv(params = params_yes_tree, data = dtrain_yes_tree, nrounds=4000, nfold=5, early_stopping_rounds = 20)
model_stage2_yestree <- xgboost::xgb.train(params = params_yes_tree, data = dtrain_yes_tree, nrounds = cv_yestree$early_stop$best_iteration, verbose = 0)




#4. use two-stage models to backfill in high HF areas----

# stage datasets for backfilling (areas with high human footprint)
CanHF_1km_present_abiotic <- CanHF_1km_present[,c(predictor_vars, "SCANFI_1km")]

# even though we already have SCANFI classifications, in disturbed areas that original label might not 
# represent natural vegetation. By recoding, we know which pixels should be backfilled 
# (for non-water areas) and which should remain unchanged (water).
backfill_data_stage1 <- 
  CanHF_1km_present_abiotic |> 
  dplyr::mutate(SCANFI_1km_broadclass = case_when(
    SCANFI_1km %in% c(1, 2, 4) ~ 0,  # bryoid, herbs, shrub (non-trees) + rock (likely disturbed)
    SCANFI_1km %in% c(5, 6, 7) ~ 1,  # broadleaf, conifer, mixed (trees)
    SCANFI_1km %in% c(3, 8) ~ 2)) |>    # water (abiotic)
  dplyr::mutate(SCANFI_1km_broadclass = factor(SCANFI_1km_broadclass, levels=c(0,1,2)))

# reminder: predictor_vars <- c(abiotic_vars, "lon", "lat", "SCANFI_1km_rock", "SCANFI_1km_water")
# including the `label` in the DMatrix is practically useless unless we use them for comparitive purposes later on..
dbackfill <- xgboost::xgb.DMatrix(data = as.matrix(backfill_data_stage1[,predictor_vars]), label=as.numeric(backfill_data_stage1$SCANFI_1km_broadclass) - 1)

# stage 1: estimate broad SCANFI classes in high HF areas
backfill_data_stage1$SCANFI_1km_broadclass_predicted  <- predict(model_stage1, newdata = dbackfill)



# stage 2A:
# for non-treed predictions:
# filter high HF data for finer subclass predictions
backfill_data_stage2A <- dplyr::filter(backfill_data_stage1, SCANFI_1km_broadclass_predicted == 0)
dbackfill_stage2A <- xgboost::xgb.DMatrix(data = as.matrix(backfill_data_stage2A[, predictor_vars]))
pred_stage2A <- predict(model_stage2_notree, newdata = dbackfill_stage2A)

# convert numeric predictions (0, 1, 2) to original SCANFI classes using the factor mapping from training
unique(pred_stage2A) # [1] 2 3 1 no bryoids predicted..
pred_labels_stage2A <- factor(pred_stage2A, levels = c(0,1,2), labels = c("1", "2", "4"))
backfill_data_stage2A <- dplyr::mutate(backfill_data_stage2A, SCANFI_backfill = pred_labels_stage2A)

# stage 2B:
# for treed predictions:
# filter high HF data for finer subclass predictions
backfill_data_stage2B <- dplyr::filter(backfill_data_stage1, SCANFI_1km_broadclass_predicted == 1)
dbackfill_stage2B <- xgboost::xgb.DMatrix(data = as.matrix(backfill_data_stage2B[, predictor_vars]))
pred_stage2B <- predict(model_stage2_yestree, newdata = dbackfill_stage2B)

# convert numeric predictions (0, 1, 2) to original SCANFI classes using the factor mapping from training
unique(pred_stage2B)
pred_labels_stage2B <- factor(pred_stage2B, levels = c(0,1,2), labels = c("5", "6", "7"))
backfill_data_stage2B <- dplyr::mutate(backfill_data_stage2B, SCANFI_backfill = pred_labels_stage2B)

# stage 2C: "predict" (i.e. retain) abiotic SCANFI classes for water (but not rock)
backfill_water_stage2C <- 
  backfill_data_stage1 |> 
  dplyr::filter(SCANFI_1km == 8) |>   # only keep water pixels
  dplyr::mutate(SCANFI_backfill = "8")

# combine backfill predictions into single data frame
prediction_stage2 <- dplyr::bind_rows(backfill_data_stage2A, backfill_data_stage2B, backfill_water_stage2C)




#5. spatialize backfilled locations (high HF)----

# define an empty raster for placing backfilled features into
empty_raster <- terra::rast(extent=ext(stack_bcr14_2020), crs=crs(stack_bcr14_2020), resolution=1000)

prediction_stage2$SCANFI_backfill <- as.numeric(prediction_stage2$SCANFI_backfill)

prediction_stage2_raster <-
  terra::vect(prediction_stage2, geom = c("lon", "lat"), crs = crs(stack_bcr14_2020)) |>
  terra::rasterize(x=_, y = empty_raster, field = "SCANFI_backfill")

scanfi_cats <- data.frame(
  ID = c(1, 2, 3, 4, 5, 6, 7, 8),
  category = c("bryoid", "herb", "rock", "shrub", "broadleaf", "conifer", "mixed", "water"))

levels(prediction_stage2_raster) <- scanfi_cats

# fill in low HF areas (NAs in prediction raster because we only backfilled high HF)
# with low HF vegetation
#  `cover()` fills any NA cells in the first raster with values from the second
prediction_stage2_raster <- terra::cover(prediction_stage2_raster, stack_bcr14_2020$SCANFI_1km)
terra::writeRaster(prediction_stage2_raster, file.path(getwd(), "data", "derived_data", "backfilled_rasters", "BCR14_SCANFI_1km.tif"), overwrite=TRUE)

#6. visualize backfilling procedure----

my_colours <- c(
  "#006644",  # ID 1: bryoid
  "#006644",  # ID 2: herb
  "#D55E00",  # ID 3: rock
  "#006644",  # ID 4: shrub
  "#009E73",  # ID 5: broadleaf (treed, darker shade)
  "#009E73",  # ID 6: conifer (treed, darker shade)
  "#009E73",  # ID 7: mixed (treed, darker shade)
  "#0072B2"   # ID 8: water
)

# plot pre-backfilled landscape for "SCANFI_1km"
levels(stack_bcr14_2020$SCANFI_1km) <- scanfi_cats
terra::plot(stack_bcr14_2020$SCANFI_1km, col = my_colours)

# plot post-backfilled landscape for "SCANFI_1km"
terra::plot(prediction_stage2_raster, col = my_colours)


lines(bcr14_boundary, col="black", lwd=1)



#7. plot low vs high human footprint areas----
my_colours <- colorRampPalette(c("#0072B2", "#009E73", "#F0E442", "#D55E00"))(100)

# re-define data frames so that it's not dropping NAs from any particular covariate (e.g. SCANFI)
# i.e. use df_bcr14_2020 instead of df_bcr14_2020_i
CanHF_1km_present <- dplyr::filter(df_bcr14_2020, CanHF_1km > q50)
CanHF_1km_absent <- dplyr::filter(df_bcr14_2020, CanHF_1km <= q50) 

 

terra::vect(CanHF_1km_absent, geom = c("lon", "lat"), crs = crs(empty_raster)) |>
  terra::rasterize(y = empty_raster, field = "CanHF_1km") |> 
  terra::plot(col= colorRampPalette(c("#56B4E9", "#0072B2"))(1000))

terra::vect(CanHF_1km_present, geom = c("lon", "lat"), crs = crs(empty_raster)) |>
  terra::rasterize(y = empty_raster, field = "CanHF_1km") |> 
  terra::plot(col=my_colours)

lines(bcr14_boundary, col="black", lwd=1)






