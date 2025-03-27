# ---
# title: National Models 5.0 - create models to predict biotic VLCE classes from abiotic landscape
# author: Mannfred Boehm
# created: March 11, 2025
# ---


#1. run lines 1-134 in "test_backfill_SCANFI_bcr14.R"---- 
# this will import covariate rasters, and generate a dataframe `df_bcr14_2020`from 
# the predictor covariate stack "can14_2020.tif"
# We also set the threshold for low/high human footprint by defining `q50`


# stage a unique spatial subset of the BCR by removing NAs for VLCE_1km
# we'll treat water (20), low HF barren (33), and low HF rock (32) VLCE classes 
# as known features of the abiotic landscape
df_bcr14_2020_i <-  tidyr::drop_na(df_bcr14_2020, VLCE_1km)

# does every VLCE class appear in BCR14? 
dplyr::count(df_bcr14_2020_i, VLCE_1km)

# identify pixels with "high" and "low" human footprint 
CanHF_1km_present <- dplyr::filter(df_bcr14_2020_i, CanHF_1km > q50)
CanHF_1km_absent <- dplyr::filter(df_bcr14_2020_i, CanHF_1km <= q50) 


# drop biotic predictors that were thinned by VIF (see line 300 in "04.Stratify.R")
biotic_vars_thinned <- biotic_vars[which(biotic_vars %in% colnames(df_bcr14_2020))]

# subset low HF dataset to abiotic and biotic predictors (excluding SCANFI and MODIS) 
biotic_vars_thinned <- biotic_vars_thinned[!biotic_vars_thinned %in% c("SCANFI_1km", "MODISLCC_1km", "VLCE_1km")]


# subset low HF dataset to abiotic predictors and the biotic response
# include the response "VLCE_1km" so that it can used for labelling 
predictor_vars <- c(abiotic_vars, "lon", "lat", biotic_vars_thinned)
CanHF_1km_absent_abiotic <- CanHF_1km_absent[,c(predictor_vars, "VLCE_1km")] 




#2. two-stage modelling----
# create broad VLCE_1km classes
# remove water (20) and low HF rock (32) and barren (33) areas:
# they are static so we will directly copy them over later. 
# also remove class 0 (no change): most "no change since 1987" areas also have high HF
# and even filtering out CanHF_1km areas leaves a significant "no change" signal in the Quebec-Windsor corridor
# which leads to most predictions being classed as "no change"
training_data_stage1 <- 
  CanHF_1km_absent_abiotic |> 
  dplyr::filter(!(VLCE_1km %in% c(0,20,32,33))) |>  
  dplyr::mutate(VLCE_1km_broadclass = case_when(
    VLCE_1km %in% c(80,81) ~ 0,  # wetland vs wetland-tree
    VLCE_1km %in% c(40,50,100) ~ 1, # bryoid, shrub, herb
    VLCE_1km %in% c(210:230) ~ 2)) |>  # trees 
  dplyr::mutate(VLCE_1km_broadclass = factor(VLCE_1km_broadclass, levels=c(0,1,2)))

count(training_data_stage1, VLCE_1km_broadclass)


# create DMatrix (a data structure for better speed when using xgboost)
# reminder: `predictor_vars` is `c(abiotic_vars, "lon", "lat")
dtrain_stage1 <- xgboost::xgb.DMatrix(data = as.matrix(training_data_stage1[,predictor_vars]), label = as.numeric(training_data_stage1$VLCE_1km_broadclass) - 1)
params_stage1 <- xgboost::xgb.params(objective = "multi:softmax", eval_metric = "mlogloss", max_depth = 3, eta = 0.1, num_class = length(levels(training_data_stage1$VLCE_1km_broadclass)))


# train stage 1 model (VLCE no change vs "natural" rocks vs vegetation)
# cv_stage1$early_stop$best_iteration 513 
cv_stage1 <-xgboost::xgb.cv(params = params_stage1, data = dtrain_stage1, nrounds=500, nfold=5, early_stopping_rounds = 20, verbose=FALSE)
model_stage1 <- xgboost::xgb.train(params = params_stage1, data = dtrain_stage1, nrounds = cv_stage1$early_stop$best_iteration, verbose = 0)



#3. build three multi-class models to predict real VLCE classes from----
# broad class 0 (wetland vs wetland-tree), 1 (bryoid vs herb vs shrub) and 2 (tree types). 

# train a model to split broad class 0 into 2 classes (wetland vs wetland-tree)
training_data_stage2_class0 <- 
  training_data_stage1 |> 
  dplyr::filter(VLCE_1km_broadclass == 0) |> 
  dplyr::mutate(VLCE_1km = factor(VLCE_1km, levels = c(80, 81)))

dtrain_class0 <- xgb.DMatrix(
  data = as.matrix(training_data_stage2_class0[, predictor_vars]), 
  label = as.numeric(training_data_stage2_class0$VLCE_1km) - 1) # as.numeric() starts the classes at 1 but xgboost labels need to start at 0

params_class0 <- 
  xgboost::xgb.params(objective = "multi:softmax", 
                      num_class = length(levels(training_data_stage2_class0$VLCE_1km)),
                      eval_metric = "mlogloss",
                      max_depth = 3,
                      eta = 0.05)

cv_class0 <- xgboost::xgb.cv(params = params_class0, data = dtrain_class0, nrounds=500, nfold=5, early_stopping_rounds = 20, verbose=FALSE)
model_class0 <- xgb.train(params = params_class0, data = dtrain_class0, nrounds = cv_class0$early_stop$best_iteration, verbose = 0)


# train a model to split broad class 1 into 3 classes (bryoid, shrub, herb)
training_data_stage2_class1 <- 
  training_data_stage1 |> 
  dplyr::filter(VLCE_1km_broadclass == 1) |> 
  dplyr::mutate(VLCE_1km = factor(VLCE_1km, levels = c(40, 50, 100)))

dtrain_class1 <- xgb.DMatrix(
  data = as.matrix(training_data_stage2_class1[, predictor_vars]), 
  label = as.numeric(training_data_stage2_class1$VLCE_1km) - 1) # as.numeric() starts the classes at 1 but xgboost labels need to start at 0

params_class1 <- 
  xgboost::xgb.params(objective = "multi:softmax", 
                      num_class = length(levels(training_data_stage2_class1$VLCE_1km)),
                      eval_metric = "mlogloss",
                      max_depth = 3,
                      eta = 0.05)

cv_class1 <- xgboost::xgb.cv(params = params_class1, data = dtrain_class1, nrounds=500, nfold=5, early_stopping_rounds = 20, verbose=FALSE)
model_class1 <- xgb.train(params = params_class1, data = dtrain_class1, nrounds = cv_class1$early_stop$best_iteration, verbose = 0)


# train a model to split broad class 2 (broadleaf vs conifer vs mixed)
training_data_stage2_class2 <- 
  training_data_stage1 |> 
  dplyr::filter(VLCE_1km_broadclass == 2) |> 
  dplyr::mutate(VLCE_1km = factor(VLCE_1km, levels = c(210, 220, 230)))

dtrain_class2 <- xgb.DMatrix(
  data = as.matrix(training_data_stage2_class2[, predictor_vars]), 
  label = as.numeric(training_data_stage2_class2$VLCE_1km) - 1) # as.numeric() starts the classes at 1 but xgboost labels need to start at 0

params_class2 <- 
  xgboost::xgb.params(objective = "multi:softmax", 
                      num_class = length(levels(training_data_stage2_class2$VLCE_1km)),
                      eval_metric = "mlogloss",
                      max_depth = 3,
                      eta = 0.05)

# cv_class2$early_stop$best_iteration 226
cv_class2 <- xgboost::xgb.cv(params = params_class2, data = dtrain_class2, nrounds=2000, nfold=5, early_stopping_rounds = 20, verbose=FALSE)
model_class2 <- xgb.train(params = params_class2, data = dtrain_class2, nrounds = cv_class2$early_stop$best_iteration, verbose = 0)



#4. use two-stage models to backfill in high HF areas----

# first, incorporate backfilled continuous variables at high HF locations
# (they will help to improve predictions of VLCE classes)
df_backfilled <-
  terra::rast("C:/Users/mannf/Proton Drive/mannfredboehm/My files/Drive/boreal_avian_modelling_project/ImpactAssessment/data/derived_data/backfilled_rasters/BCR14_backfilled_continuous.tif") |> 
  terra::as.data.frame(xy=TRUE) |>
  tibble::as_tibble() |> 
  dplyr::rename(lon=x, lat=y) 

# isolate high HF locations and their abiotic features 
CanHF_1km_present_abiotic <- CanHF_1km_present[,c(abiotic_vars, "VLCE_1km", "lon", "lat")]
nrow(CanHF_1km_present_abiotic) #128840

# remove water, wetland, and wetland-trees because we don't want to backfill there
# keep "no change" areas (VLCE class 0) because we *do* want to backfill vegetation there
backfill_data_stage1 <- 
  CanHF_1km_present_abiotic |> 
  dplyr::filter(!(VLCE_1km %in% c(20,80,81))) |> 
  dplyr::left_join(df_backfilled, by=c("lon", "lat")) # add in (backfilled) biotic features to high HF areas

nrow(backfill_data_stage1) # should be 125041 which is 128840 - (water + wetland + wetland-tree)

# reminder: predictor_vars <- c(abiotic_vars, "lon", "lat")
# including the `label` in the DMatrix is practically useless unless we use them for comparative purposes later on..
dbackfill <- xgboost::xgb.DMatrix(data = as.matrix(backfill_data_stage1[,predictor_vars]))

# stage 1: estimate broad VLCE classes in high HF areas
backfill_data_stage1$VLCE_1km_broadclass_predicted  <- predict(model_stage1, newdata = dbackfill)
count(backfill_data_stage1, VLCE_1km_broadclass_predicted)
# VLCE_1km_broadclass_predicted     n
# <dbl> <int>
# 1                                 0 1743 wetlands
# 2                                 1 3697 bryoid, shrub, herb
# 3                                 2 119601 trees

# stage 2A:
# predict VLCE class for predicted broad class 0 (wetland and wetland-tree):
backfill_data_stage2A <- dplyr::filter(backfill_data_stage1, VLCE_1km_broadclass_predicted == 0)
dbackfill_stage2A <- xgboost::xgb.DMatrix(data = as.matrix(backfill_data_stage2A[, predictor_vars]))
pred_stage2A <- predict(model_class0, newdata = dbackfill_stage2A)


# convert numeric predictions (0, 1, 2) to original VLCE classes using the factor mapping from training
unique(pred_stage2A) # [1] 0 only wetland predicted (no wetland-tree)
pred_labels_stage2A <- factor(pred_stage2A, levels = c(0,1), labels = c("80","81"))
backfill_data_stage2A <- dplyr::mutate(backfill_data_stage2A, VLCE_backfill = pred_labels_stage2A)




# stage 2B:
# predict VLCE class for predicted broad class 1 (bryoid, shrub, herb):
backfill_data_stage2B <- dplyr::filter(backfill_data_stage1, VLCE_1km_broadclass_predicted == 1)
dbackfill_stage2B <- xgboost::xgb.DMatrix(data = as.matrix(backfill_data_stage2B[, predictor_vars]))
pred_stage2B <- predict(model_class1, newdata = dbackfill_stage2B)


# convert numeric predictions (0, 1, 2) to original VLCE classes using the factor mapping from training
unique(pred_stage2B) # [1] 1 2  no bryoids predicted
pred_labels_stage2B <- factor(pred_stage2B, levels = c(0,1,2), labels = c("40","50","100"))
backfill_data_stage2B <- dplyr::mutate(backfill_data_stage2B, VLCE_backfill = pred_labels_stage2B)



# stage 2C:
# predict VLCE class for predicted broad class 2 (conifer, broadleaf, mixed):
backfill_data_stage2C <- dplyr::filter(backfill_data_stage1, VLCE_1km_broadclass_predicted == 2)
dbackfill_stage2C <- xgboost::xgb.DMatrix(data = as.matrix(backfill_data_stage2C[, predictor_vars]))
pred_stage2C <- predict(model_class2, newdata = dbackfill_stage2C)


# convert numeric predictions (0, 1, 2) to original VLCE classes using the factor mapping from training
unique(pred_stage2C) # [1] 0 2 1 all tree types are predicted
pred_labels_stage2C <- factor(pred_stage2C, levels = c(0,1,2), labels = c("210", "220", "230"))
backfill_data_stage2C <- dplyr::mutate(backfill_data_stage2C, VLCE_backfill = pred_labels_stage2C)



# stage 2D: 
# retain VLCE classes for water and wetland from areas with high HF
retain_aqueous_stage2D <- 
  CanHF_1km_present |> 
  dplyr::filter(VLCE_1km %in% c(20,80,81)) |> 
  dplyr::mutate(VLCE_backfill = as.character(VLCE_1km))


# combine backfill predictions into single data frame
# number of predictions and carry-overs so should be equal to nrow(CanHF_1km_present)
prediction_stage2 <- dplyr::bind_rows(backfill_data_stage2A, backfill_data_stage2B, backfill_data_stage2C, retain_aqueous_stage2D)
prediction_stage2$VLCE_backfill <- as.numeric(prediction_stage2$VLCE_backfill)
nrow(CanHF_1km_present_abiotic) #128840
nrow(prediction_stage2) #128840



#5. spatialize backfilled locations (high HF)----

# define an empty raster for placing backfilled features into
empty_raster <- terra::rast(extent=ext(stack_bcr14_2020), crs=crs(stack_bcr14_2020), resolution=1000)

# build backfilled raster
prediction_stage2_raster <-
  terra::vect(prediction_stage2, geom = c("lon", "lat"), crs = crs(stack_bcr14_2020)) |>
  terra::rasterize(x=_, y = empty_raster, field = "VLCE_backfill")

length(values(prediction_stage2_raster, na.rm=T)) #128840
length(values(stack_bcr14_2020$VLCE_1km, na.rm=T)) #274264


#  `cover()` fills any NA cells in the first raster with values from the second
# i.e. fill in NAs in high HF backfill layer with low HF pixels from original layer
prediction_stage2_raster <- terra::cover(prediction_stage2_raster, stack_bcr14_2020$VLCE_1km)
length(values(prediction_stage2_raster, na.rm=T)) #274264



# define classes
vlce_cats <- data.frame(
  ID = sort(unique(df_bcr14_2020_i$VLCE_1km)),  # 0  20  32  33  40  50  80  81 100 210 220 230
  category = c("no change", "water", "rock", "barren", "bryoid", "shrub",
               "wetland", "wetland-tree", "herb", "conifer", "broadleaf", "mixed forest"))

levels(prediction_stage2_raster) <- vlce_cats
levels(stack_bcr14_2020$VLCE_1km) <- vlce_cats



#6. visualize backfilling procedure----

my_colours <- c(
  "#CC79A7",  # 0: no change
  "#0072B2",  #20: water
  "#D55E00",  #32: rock
  "#D55E00",  #33: barren
  "#bae4b3",  #40: bryoid
  "#bae4b3",  #50: shrub
  "#56B4E9",  #80: wetland
  "#009E73",  #81: wetland-tree
  "#bae4b3",  #100: herb
  "#006644",  # 3: conifer
  "#4DAC26",  # 4: broadleaf
  "#74c476")  # 5: mixed forest
 
# plot pre-backfilled landscape for "VLCE_1km"
terra::plot(stack_bcr14_2020$VLCE_1km, col = my_colours)
lines(bcr14_boundary, col="black", lwd=1)

# plot post-backfilled landscape for "VLCE_1km"
terra::plot(prediction_stage2_raster, col = my_colours)
lines(bcr14_boundary, col="black", lwd=1)  


# save it if it's good
terra::writeRaster(prediction_stage2_raster, file.path(getwd(), "data", "derived_data", "backfilled_rasters", "BCR14_VLCE_1km.tif"), overwrite=TRUE)

