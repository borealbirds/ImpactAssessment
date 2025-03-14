# ---
# title: National Models 5.0 - create models to predict biotic SCANFI classes from abiotic landscape
# author: Mannfred Boehm
# created: March 11, 2025
# ---


#1. run lines 1-134 in "test_backfill_SCANFI_bcr14.R"---- 
# this will import covariate rasters, and generate a dataframe `df_bcr14_2020`from 
# the predictor covariate stack "can14_2020.tif"
# We also set the threshold for low/high human footprint by defining `q50`




#4. test whether MODISLCC_1km land cover classes are accurately predicted---- 
# will eventually need to do the same for SCANFI_1km, VLCE_1km

# drop biotic predictors that were thinned by VIF (see line 300 in "04.Stratify.R")
biotic_vars_thinned <- biotic_vars[which(biotic_vars %in% colnames(df_bcr14_2020))]

# stage a unique spatial subset of the BCR by removing NAs for MODISLCC_1km
# also: treat water, barren, cropland MODIS classes as known features of the abiotic landscape
# i.e. we aren't trying to predict them. They are predictors.
# note: don't convert MODISLCC to a factor because xgb.DMatrix is expecting numerics
# note: urban classes will mostly be filtered out when we create a low CANHF_1km dataset
df_bcr14_2020_i <- 
  df_bcr14_2020 |> 
  tidyr::drop_na(MODISLCC_1km) |> 
  dplyr::mutate(MODISLCC_1km_crop = ifelse(MODISLCC_1km == 12, yes=1, no=0)) |> # create "cropland" predictor  (should perfectly predict MODIS class 12)
  dplyr::mutate(MODISLCC_1km_urban = ifelse(MODISLCC_1km == 13, yes=1, no=0)) |>  # create "urban" predictor (should perfectly predict MODIS class 13)
  dplyr::mutate(MODISLCC_1km_barren = ifelse(MODISLCC_1km == 16, yes=1, no=0)) |>  # create "barren" predictor (should perfectly predict MODIS class 16)
  dplyr::mutate(MODISLCC_1km_water = ifelse(MODISLCC_1km == 17, yes=1, no=0))  # create "water" predictor (should perfectly predict MODIS class 17)

# does every MODIS class appear in BCR14? 
dplyr::count(df_bcr14_2020_i, MODISLCC_1km)

# identify pixels with "high" and "low" human footprint 
CanHF_1km_present <- dplyr::filter(df_bcr14_2020_i, CanHF_1km > q50)
CanHF_1km_absent <- dplyr::filter(df_bcr14_2020_i, CanHF_1km <= q50) 

# subset low HF dataset to only those with abiotic predictors and the biotic response
# include the response "MODISLCC" so that it can used for labelling 
# note: # cropland is reduced from 14840 to 212 (98.5%), urban is reduced 99.5%
predictor_vars <- c(abiotic_vars, "lon", "lat", "MODISLCC_1km_crop", "MODISLCC_1km_urban", "MODISLCC_1km_barren", "MODISLCC_1km_water")
CanHF_1km_absent_abiotic <- CanHF_1km_absent[,c(predictor_vars, "MODISLCC_1km")] 



#5. two-stage modelling----
# stage 1: create broad MODISLCC_1km classes in low HF areas

# classes 15 (snow/ice) and 18 (unclassified) are not in BCR14
dplyr::count(CanHF_1km_absent_abiotic, MODISLCC_1km)

training_data_stage1 <- 
  CanHF_1km_absent_abiotic |> 
  dplyr::filter(!MODISLCC_1km %in% c(12, 13, 14)) |> # remove human disturbances (classes 12, 13, 14)
  dplyr::mutate(MODISLCC_1km_broadclass = case_when(
    MODISLCC_1km %in% c(1:5) ~ 0,  # trees
    MODISLCC_1km %in% c(6:10) ~ 1,  # non-tree vegetation
    MODISLCC_1km == 11 ~ 2, # wetland 
    MODISLCC_1km %in% c(16:17) ~ 3)) |> # barren and water
  dplyr::mutate(MODISLCC_1km_broadclass = factor(MODISLCC_1km_broadclass, levels=c(0,1,2,3)))


# create DMatrix (a data structure for better speed when using xgboost)
# reminder: `predictor_vars` is `c(abiotic_vars, "lon", "lat", "MODISLCC_1km_crop", "MODISLCC_1km_urban", "MODISLCC_1km_barren", "MODISLCC_1km_water")`
dtrain_stage1 <- xgboost::xgb.DMatrix(data = as.matrix(training_data_stage1[,predictor_vars]), label = as.numeric(training_data_stage1$MODISLCC_1km_broadclass) - 1)
params_stage1 <- xgboost::xgb.params(objective = "multi:softmax", eval_metric = "mlogloss", max_depth = 3, eta = 0.1, num_class = length(levels(training_data_stage1$MODISLCC_1km_broadclass)))

# train stage 1 model to classify broad MODIS classes in low HF areas
# cv_stage1$early_stop$best_iteration 3139
cv_stage1 <-xgboost::xgb.cv(params = params_stage1, data = dtrain_stage1, nrounds=5000, nfold=5, early_stopping_rounds = 20)
model_stage1 <- xgboost::xgb.train(params = params_stage1, data = dtrain_stage1, nrounds = cv_stage1$early_stop$best_iteration, verbose = 0)





#6. build two multi-class models to predict real MODIS classes from 
# broad classes 0 (tree class) and 1 (non-tree vegetation class)
# NOTE: we won't try to split the other broad classes (wetlands, HF, water, unclassified)
# we'll instead just port them over directly as "determined"

# train a model to split broad class 0 into 5 MODIS classes (tree types)
# filter for broadclass==0, then use `unique()` to check what classes actually got modelled before using `factor()`
training_data_stage2_class0 <- 
  training_data_stage1 |> 
  dplyr::filter(MODISLCC_1km_broadclass == 0) |> 
  dplyr::mutate(MODISLCC_1km = factor(MODISLCC_1km, levels = c(1, 2, 3, 4, 5)))

dtrain_class0 <- xgb.DMatrix(
  data = as.matrix(training_data_stage2_class0[, predictor_vars]), 
  label = as.numeric(training_data_stage2_class0$MODISLCC_1km) - 1) # as.numeric() starts the classes at 1 but xgboost labels need to start at 0

params_class0 <- 
  xgboost::xgb.params(objective = "multi:softmax", 
                      num_class = length(levels(training_data_stage2_class0$MODISLCC_1km)),
                      eval_metric = "mlogloss",
                      max_depth = 3,
                      eta = 0.1)

# cv_class0$early_stop$best_iteration 2291
cv_class0 <- xgboost::xgb.cv(params = params_class0, data = dtrain_class0, nrounds=3000, nfold=5, early_stopping_rounds = 20)
model_class0 <- xgb.train(params = params_class0, data = dtrain_class0, nrounds = cv_class0$early_stop$best_iteration, verbose = 0)


# train a model to split broad class 1 into 5 MODIS classes (non-tree vegetation)
training_data_stage2_class1 <- 
  training_data_stage1 |> 
  dplyr::filter(MODISLCC_1km_broadclass == 1) |> 
  dplyr::mutate(MODISLCC_1km = factor(MODISLCC_1km, levels = c(6,7,8,9,10))) 

dtrain_class1 <- xgb.DMatrix(
  data = as.matrix(training_data_stage2_class1[, predictor_vars]), 
  label = as.numeric(training_data_stage2_class1$MODISLCC_1km) - 1) # as.numeric() starts the classes at 1 but xgboost labels need to start at 0

params_class1 <- 
  xgboost::xgb.params(objective = "multi:softmax", 
                      num_class = length(levels(training_data_stage2_class1$MODISLCC_1km)),
                      eval_metric = "mlogloss",
                      max_depth = 3,
                      eta = 0.1)

cv_class1 <- xgboost::xgb.cv(params = params_class1, data = dtrain_class1, nrounds=1000, nfold=5, early_stopping_rounds = 20)
model_class1 <- xgb.train(params = params_class1, data = dtrain_class1, nrounds = cv_class1$early_stop$best_iteration, verbose = 0)




#4. use two-stage models to backfill in high HF areas----

# stage datasets for backfilling (areas with high human footprint)
CanHF_1km_present_abiotic <- CanHF_1km_present[,c(predictor_vars, "MODISLCC_1km")]

# even though we already have MODISLCC classifications, in disturbed areas that original label might not 
# represent natural vegetation. By recoding, we know which pixels should be backfilled 
# (for non-water areas) and which should remain unchanged (water).
backfill_data_stage1 <- 
  CanHF_1km_present_abiotic |> 
  dplyr::mutate(MODISLCC_1km_broadclass = case_when(
    MODISLCC_1km %in% c(1:5) ~ 0,  # trees
    MODISLCC_1km %in% c(6:10) ~ 1,  # non-tree vegetation
    MODISLCC_1km == 11 ~ 2, # wetland 
    MODISLCC_1km %in% c(12:14) ~ 3, # human footprints 
    MODISLCC_1km %in% c(16:17) ~ 4)) |> # barren and water
  dplyr::mutate(MODISLCC_1km_broadclass = factor(MODISLCC_1km_broadclass, levels=c(0,1,2,3,4)))

# reminder: predictor_vars <- c(abiotic_vars, "lon", "lat", "MODISLCC_1km_rock", "MODISLCC_1km_water")
# including the `label` in the DMatrix is practically useless unless we use them for comparitive purposes later on..
dbackfill <- xgboost::xgb.DMatrix(data = as.matrix(backfill_data_stage1[,predictor_vars]), label=as.numeric(backfill_data_stage1$MODISLCC_1km_broadclass) - 1)

# stage 1: estimate broad SCANFI classes in high HF areas
backfill_data_stage1$MODISLCC_1km_broadclass_predicted  <- predict(model_stage1, newdata = dbackfill)


# stage 2A:
# for broad class 0 predictions (trees):
# filter high HF data for finer subclass predictions (tree type)
backfill_data_stage2A <- dplyr::filter(backfill_data_stage1, MODISLCC_1km_broadclass_predicted == 0)
dbackfill_stage2A <- xgboost::xgb.DMatrix(data = as.matrix(backfill_data_stage2A[, predictor_vars]))
pred_stage2A <- predict(model_class0, newdata = dbackfill_stage2A)

# convert numeric predictions (0, 1, 2) to original MODIS classes using the factor mapping from training
unique(pred_stage2A) # [1] 0 4 3 no evegreen broadleaf or deciduous needle predicted
pred_labels_stage2A <- factor(pred_stage2A, levels = c(0,1,2,3,4), labels = c("1","2","3","4","5"))
backfill_data_stage2A <- dplyr::mutate(backfill_data_stage2A, MODIS_backfill = pred_labels_stage2A)


# stage 2B:
# for broad class 1 predictions (non-tree vegetation):
# filter high HF data for finer subclass predictions (vegetation type)
backfill_data_stage2B <- dplyr::filter(backfill_data_stage1, MODISLCC_1km_broadclass_predicted == 1)
dbackfill_stage2B <- xgboost::xgb.DMatrix(data = as.matrix(backfill_data_stage2B[, predictor_vars]))
pred_stage2B <- predict(model_class1, newdata = dbackfill_stage2B)

# convert numeric predictions (0, 1, 2) to original MODIS classes using the factor mapping from training
unique(pred_stage2B) # [1] 2 4 3 1 no closed shrubs predicted
pred_labels_stage2B <- factor(pred_stage2B, levels = c(0,1,2,3,4), labels = c("6","7","8","9","10"))
backfill_data_stage2B <- dplyr::mutate(backfill_data_stage2B, MODIS_backfill = pred_labels_stage2B)


# stage 2C and 2D: "predict" (i.e. retain) MODIS classes for water and wetlands, but not crops, urban, or barren
backfill_water_stage2C <- 
  backfill_data_stage1 |> 
  dplyr::filter(MODISLCC_1km == 17) |>   # only keep water pixels
  dplyr::mutate(MODIS_backfill = "17")

backfill_water_stage2D <- 
  backfill_data_stage1 |> 
  dplyr::filter(MODISLCC_1km == 11) |>   # only keep wetland pixels
  dplyr::mutate(MODIS_backfill = "11")

# combine backfill predictions into single data frame
prediction_stage2 <- dplyr::bind_rows(backfill_data_stage2A, backfill_data_stage2B, backfill_water_stage2C,backfill_water_stage2D)
prediction_stage2$MODIS_backfill <- as.numeric(prediction_stage2$MODIS_backfill)


#5. spatialize backfilled locations (high HF)----

# define an empty raster for placing backfilled features into
empty_raster <- terra::rast(extent=ext(stack_bcr14_2020), crs=crs(stack_bcr14_2020), resolution=1000)

# build backfilled raster
prediction_stage2_raster <-
  terra::vect(prediction_stage2, geom = c("lon", "lat"), crs = crs(stack_bcr14_2020)) |>
  terra::rasterize(x=_, y = empty_raster, field = "MODIS_backfill")

# fill in low HF areas (NAs in prediction raster because we only backfilled high HF)
# with low HF vegetation
#  `cover()` fills any NA cells in the first raster with values from the second
prediction_stage2_raster <- terra::cover(prediction_stage2_raster, stack_bcr14_2020$MODISLCC_1km)


# define classes for backfilled layer
modis_cats_backfill <- data.frame(
  ID = sort(unique(prediction_stage2$MODIS_backfill)), #  1  4  5  7  8  9 10 11 17 
  category = c("eg needle", "dec broadleaf", "mixed forest",
               "open shrub", "woody savanna", "savanna", "grassland", "wetland", "water"))

levels(prediction_stage2_raster) <- modis_cats_backfill

# define classes for original layer
modis_cats_original <- data.frame(
    ID = sort(unique(df_bcr14_2020_i$MODISLCC_1km)), # 1  2  3  4  5  7  8  9 10 11 12 13 14 15 16 17
    category = c("eg needle", "eg broadleaf", "dec needle", "dec broadleaf", "mixed forest",
              "open shrub", "woody savanna", "savanna", "grassland", "wetland",
             "cropland", "urban", "cropland/natural", "snow/ice", "barren", "water"))

levels(stack_bcr14_2020$MODISLCC_1km) <- modis_cats_original




#6. visualize backfilling procedure----

sort(unique(prediction_stage2$MODIS_backfill))
# 1  4  5  7  8  9 10 11 17

my_colours_backfill <- c(
  "#009E73",  # 1: evergreen needle
  "#009E73",  # 4: deciduous broadleaf
  "#009E73",  # 5: mixed forest
  "#006644",  # 7: open shrubs
  "#006644",  # 8: woody savanna
  "#006644",  # 9: savanna
  "#006644",  #10: grassland
  "#0072B2",  #11: wetland
  "#0072B2")   #17: water




sort(unique(df_bcr14_2020_i$MODISLCC_1km))
# 1  2  3  4  5  7  8  9 10 11 12 13 14 15 16 17

my_colours_original <- c(
  "#009E73",  # 1: evergreen needle
  "#009E73",  # 2: evergreen broadleaf
  "#009E73",  # 3: deciduous needle
  "#009E73",  # 4: deciduous broadleaf
  "#009E73",  # 5: mixed forest
  "#006644",  # 7: open shrubs
  "#006644",  # 8: woody savanna
  "#006644",  # 9: savanna
  "#006644",  #10: grassland
  "#0072B2",  #11: wetland
  "#D55E00",  #12: cropland
  "#D55E00",  #13: urban
  "#D55E00",  #14: cropland/natural
  "snow",     #15: snow/ice
  "grey",     #16: barren
  "#0072B2")   #17: water


  
# plot pre-backfilled landscape for "SCANFI_1km"
terra::plot(stack_bcr14_2020$MODISLCC_1km, col = my_colours_original)
lines(bcr14_boundary, col="black", lwd=1)

# plot post-backfilled landscape for "SCANFI_1km"
terra::plot(prediction_stage2_raster, col = my_colours_backfill)
lines(bcr14_boundary, col="black", lwd=1)




# save it if it's good
terra::writeRaster(prediction_stage2_raster, file.path(getwd(), "data", "derived_data", "backfilled_rasters", "BCR14_MODISLCC_1km.tif"), overwrite=TRUE)
