# ---
# title: National Models 5.0 - create models to predict biotic SCANFI classes from abiotic landscape
# author: Mannfred Boehm
# created: March 11, 2025
# ---


#1. run lines 1-157 in "test_backfill_SCANFI_bcr14.R"---- 
# this will import covariate rasters, and generate low and high HF dataframes 
# "CanHF_1km_absent" and "CanHF_1km_present"


# subset low HF dataset to abiotic and biotic predictors (excluding MODIS and VLCE) 
# in theory, urban areas have been filtered out so SCANFI "rock" is natural rock
biotic_vars_thinned <- biotic_vars_thinned[!biotic_vars_thinned %in% c("SCANFI_1km", "MODISLCC_1km", "VLCE_1km")]
predictor_vars <- c(abiotic_vars, "lon", "lat", biotic_vars_thinned) 
CanHF_1km_absent_abiotic <- CanHF_1km_absent[,c(predictor_vars, "SCANFI_1km")]



#2. two-stage modelling----
# first, split SCANFI_1km classes between tree and non-tree vegetation
training_data_stage1 <- 
  CanHF_1km_absent_abiotic |> 
  dplyr::filter(SCANFI_1km != "8") |> # remove 2336 water areas (we will replace them later with `retain_water_stage2D`)
  dplyr::mutate(SCANFI_1km_broadclass = case_when(
    SCANFI_1km %in% c(1, 2, 4) ~ 0,  # bryoid, herbs, shrub (non-trees) 
    SCANFI_1km %in% c(5, 6, 7) ~ 1,  # broadleaf, conifer, mixed (trees)
    SCANFI_1km == 3 ~ 2)) |> # low HF ("natural") rock areas
  dplyr::mutate(SCANFI_1km_broadclass = factor(SCANFI_1km_broadclass, levels=c(0,1,2)))

count(training_data_stage1, SCANFI_1km_broadclass)
# SCANFI_1km_broadclass      n
# <fct>                    <int>
# 1 0                       5255
# 2 1                     120646
# 3 2                        527


# create DMatrix (a data structure for better speed when using xgboost)
# reminder: `predictor_vars` is c(abiotic_vars, "lon", "lat", biotic_vars_thinned) 
dtrain_stage1 <- xgboost::xgb.DMatrix(data = as.matrix(training_data_stage1[,predictor_vars]), label = as.numeric(training_data_stage1$SCANFI_1km_broadclass) - 1)
params_stage1 <- xgboost::xgb.params(objective = "multi:softmax", eval_metric = "mlogloss", max_depth = 3, eta = 0.05, num_class = length(levels(training_data_stage1$SCANFI_1km_broadclass)))

# train stage 1 model 
cv_stage1 <-xgboost::xgb.cv(params = params_stage1, data = dtrain_stage1, nrounds=1000, nfold=5, early_stopping_rounds = 20, verbose=FALSE)
model_stage1 <- xgboost::xgb.train(params = params_stage1, data = dtrain_stage1, nrounds = cv_stage1$early_stop$best_iteration, verbose = 0)




#3. build two multi-class models to learn to split SCANFI broadclasses----
# split broad class 0 into bryoid, herb, shrub
# split broad class 1 into broadleaf, conifer, mixed
# note: no model is built for broad class 3 because if model_stage1 predicts rock, 
# then it's rock (i.e. can't be split further)

# filter data to broad class 0 (non tree vegetation)
# SCANFI_1km is now a 3-factor column for non-treed classes (bryoids, herbs, shrubs)
# unique(training_data_stage2A$SCANFI_1km) [1] 4 2 1 
training_data_stage2A <- 
  training_data_stage1 |> # water areas already removed
  dplyr::filter(SCANFI_1km_broadclass == 0) |> 
  dplyr::mutate(SCANFI_1km = factor(SCANFI_1km, levels = c(1, 2, 4)))

count(training_data_stage2A, SCANFI_1km)
# A tibble: 3 × 2
# SCANFI_1km     n
#1 1            154 bryoid
#2 2           1017 herb
#3 4           4084 shrub

dtrain_non_tree <- xgb.DMatrix(
  data = as.matrix(training_data_stage2A[, predictor_vars]), # exclude SCANFI_1km and SCANFI_1km_tree
  label = as.numeric(training_data_stage2A$SCANFI_1km) - 1) # as.numeric() starts the classes at 1 but xgboost labels need to start at 0

params_non_tree <- 
  xgboost::xgb.params(objective = "multi:softmax", 
                      num_class = length(levels(training_data_stage2A$SCANFI_1km)),
                      eval_metric = "mlogloss",
                      max_depth = 3,
                      eta = 0.05)

cv_notree <- xgboost::xgb.cv(params = params_non_tree, data = dtrain_non_tree, nrounds=1000, nfold=5, early_stopping_rounds = 20, verbose=FALSE)
model_stage2_notree <- xgb.train(params = params_non_tree, data = dtrain_non_tree, nrounds = cv_notree$early_stop$best_iteration, verbose = 0)

# filter data for broad class 1
# SCANFI_1km is now a 3-factor column for yes-treed classes (broadlead, conifer, mixed)
training_data_stage2B <- 
  training_data_stage1 |> 
  dplyr::filter(SCANFI_1km_broadclass == 1) |> 
  dplyr::mutate(SCANFI_1km = factor(SCANFI_1km, levels = c(5, 6, 7)))

count(training_data_stage2B, SCANFI_1km)
#SCANFI_1km     n
#1 5          24593 broadleaf
#2 6          58274 conifer
#3 7          37779 mixed

dtrain_yes_tree <- xgb.DMatrix(
  data = as.matrix(training_data_stage2B[, predictor_vars]),
  label = as.numeric(training_data_stage2B$SCANFI_1km) - 1)

params_yes_tree <- 
  xgboost::xgb.params(objective = "multi:softmax",
                      num_class = length(levels(training_data_stage2B$SCANFI_1km)),
                      eval_metric = "mlogloss",
                      max_depth = 3,
                      eta = 0.05)

cv_yestree <- xgboost::xgb.cv(params = params_yes_tree, data = dtrain_yes_tree, nrounds=4000, nfold=5, early_stopping_rounds = 20, verbose=FALSE)
model_stage2_yestree <- xgboost::xgb.train(params = params_yes_tree, data = dtrain_yes_tree, nrounds = cv_yestree$early_stop$best_iteration, verbose = 0)




#4. use two-stage models to backfill in high HF areas----

# first, incorporate backfilled continuous variables at high HF locations
# (they will help to improve predictions of SCANFI classes)
df_backfilled <-
  terra::rast("C:/Users/mannf/Proton Drive/mannfredboehm/My files/Drive/boreal_avian_modelling_project/ImpactAssessment/data/derived_data/backfilled_rasters/BCR14_backfilled_continuous.tif") |> 
  terra::as.data.frame(xy=TRUE) |>
  tibble::as_tibble() |> 
  dplyr::rename(lon=x, lat=y) 


# isolate high HF locations and their abiotic features 
CanHF_1km_present_abiotic <- CanHF_1km_present[,c(abiotic_vars, "SCANFI_1km", "lon", "lat")]

# add in (backfilled) biotic features to high HF areas
backfill_data_stage1 <- dplyr::left_join(CanHF_1km_present_abiotic, df_backfilled, by=c("lon", "lat"))

# reminder: predictor_vars <- c(abiotic_vars, "lon", "lat", biotic_vars_thinned) 
# including the `label` in the DMatrix is practically useless unless we use them for comparative purposes later on..
dbackfill <- xgboost::xgb.DMatrix(data = as.matrix(backfill_data_stage1[,predictor_vars]))

# stage 1: estimate broad SCANFI classes in high HF areas
backfill_data_stage1$SCANFI_1km_broadclass_predicted  <- predict(model_stage1, newdata = dbackfill)
# count(backfill_data_stage1, SCANFI_1km_broadclass_predicted)
#  SCANFI_1km_broadclass_predicted      n
#1                               0   8613
#2                               1 120152
#3                               2     17

# stage 2A:
# for non-treed predictions:
# filter high HF data for finer subclass predictions
backfill_data_stage2A <- dplyr::filter(backfill_data_stage1, SCANFI_1km_broadclass_predicted == 0)
dbackfill_stage2A <- xgboost::xgb.DMatrix(data = as.matrix(backfill_data_stage2A[, predictor_vars]))
pred_stage2A <- predict(model_stage2_notree, newdata = dbackfill_stage2A)

# convert numeric predictions (0, 1, 2) to original SCANFI classes using the factor mapping from training
pred_labels_stage2A <- factor(pred_stage2A, levels = c(0,1,2), labels = c("1", "2", "4"))
backfill_data_stage2A <- dplyr::mutate(backfill_data_stage2A, SCANFI_backfill = pred_labels_stage2A)
# count(backfill_data_stage2A, SCANFI_backfill)
# SCANFI_backfill       n
# 2 2                8360
# 3 4                 253

# stage 2B:
# for treed predictions:
# filter high HF data for finer subclass predictions
backfill_data_stage2B <- dplyr::filter(backfill_data_stage1, SCANFI_1km_broadclass_predicted == 1)
dbackfill_stage2B <- xgboost::xgb.DMatrix(data = as.matrix(backfill_data_stage2B[, predictor_vars]))
pred_stage2B <- predict(model_stage2_yestree, newdata = dbackfill_stage2B)

# convert numeric predictions (0, 1, 2) to original SCANFI classes using the factor mapping from training
pred_labels_stage2B <- factor(pred_stage2B, levels = c(0,1,2), labels = c("5", "6", "7"))
backfill_data_stage2B <- dplyr::mutate(backfill_data_stage2B, SCANFI_backfill = pred_labels_stage2B)
count(backfill_data_stage2B, SCANFI_backfill)
# SCANFI_backfill     n
# 1 5               93670
# 2 6               6661
# 3 7               19821

# stage 2C: areas predicted as broad class 2 (equivalent to SCANFI_1km == 3 or "rock")
backfill_rock_stage2C <- 
  backfill_data_stage1 |> 
  dplyr::filter(SCANFI_1km_broadclass_predicted == 2) |>  
  dplyr::mutate(SCANFI_backfill = "3")

retain_water_stage2D <-
  CanHF_1km_present |> 
  dplyr::filter(SCANFI_1km == 8) |>   # only keep water pixels
  dplyr::mutate(SCANFI_backfill = as.character(SCANFI_1km))

# combine backfill predictions into single data frame
prediction_stage2 <- dplyr::bind_rows(backfill_data_stage2A, backfill_data_stage2B, backfill_rock_stage2C, retain_water_stage2D)
prediction_stage2$SCANFI_backfill <- as.numeric(prediction_stage2$SCANFI_backfill)
sort(unique(prediction_stage2$SCANFI_backfill)) #  2 3 4 5 6 7 8



#5. spatialize backfilled locations (high HF)----

# define an empty raster for placing backfilled features into
empty_raster <- terra::rast(extent=ext(stack_bcr14_2020), crs=crs(stack_bcr14_2020), resolution=1000)

prediction_stage2_raster <-
  terra::vect(prediction_stage2, geom = c("lon", "lat"), crs = crs(stack_bcr14_2020)) |>
  terra::rasterize(x=_, y = empty_raster, field = "SCANFI_backfill")

scanfi_cats <- data.frame(
  ID = c(1, 2, 3, 4, 5, 6, 7, 8),
  category = c("bryoid", "herb", "rock", "shrub", "broadleaf", "conifer", "mixed", "water"))

levels(prediction_stage2_raster) <- scanfi_cats

# `cover()` fills any NA cells in the first raster with values from the second
# so, fill in low HF areas (NAs in backfilled raster because we only backfilled high HF)
# with low HF vegetation areas in `stack_bcr14_2020$SCANFI_1km`
# also, fill in water areas, because they were not backfilled
prediction_stage2_raster <- terra::cover(prediction_stage2_raster, stack_bcr14_2020$SCANFI_1km)




#6. visualize backfilling procedure----

my_colours <- c(
  "#bae4b3",  # ID 1: bryoid
  "#bae4b3",  # ID 2: herb
  "#E69F00",  # ID 3: rock
  "#bae4b3",  # ID 4: shrub
  "#7FFF00",  # ID 5: broadleaf (treed, darker shade)
  "#006644",  # ID 6: conifer (treed, darker shade)
  "#9ACD32",  # ID 7: mixed (treed, darker shade)
  "#0072B2"   # ID 8: water
)

# plot pre-backfilled landscape for "SCANFI_1km"
levels(stack_bcr14_2020$SCANFI_1km) <- scanfi_cats
terra::plot(stack_bcr14_2020$SCANFI_1km, col = my_colours)
lines(bcr14_boundary, col="black", lwd=1)

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

# if it's good, write it
terra::writeRaster(prediction_stage2_raster, file.path(getwd(), "data", "derived_data", "backfilled_rasters", "BCR14_SCANFI_1km.tif"), overwrite=TRUE)





