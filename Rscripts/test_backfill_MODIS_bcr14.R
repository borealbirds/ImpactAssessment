# ---
# title: National Models 5.0 - create models to predict biotic MODIS classes from abiotic landscape
# author: Mannfred Boehm
# created: March 11, 2025
# ---



#1. run lines 1-134 in "test_backfill_SCANFI_bcr14.R"---- 
# this will import covariate rasters, and generate a dataframe `df_bcr14_2020`from 
# the predictor covariate stack "can14_2020.tif"
# We also set the threshold for low/high human footprint by defining `q50`

#2. inspect MODIS LCCs
plot(as.factor(stack_bcr14_2020$MODISLCC_1km))

# MODIS Landcover Classes:
# 1 = evergreen needleleaf    10 = grassland
# 2 = evergreen broadleaf     11 = permanent wetland
# 3 = deciduous needleleaf    12 = cropland
# 4 = deciduous broadleaf     13 = urban
# 5 = mixed forest            14 = cropland/natural mosaic
# 6 = closed shrubland        15 = permanent snow/ice
# 7 = open shrubland          16 = barren
# 8 = woody savanna           17 = water
# 9 = savanna                 18 = unclassified


#4. test whether MODISLCC_1km land cover classes are accurately predicted---- 
# will eventually need to do the same for SCANFI_1km, VLCE_1km

# drop biotic predictors that were thinned by VIF (see line 300 in "04.Stratify.R")
biotic_vars_thinned <- biotic_vars[which(biotic_vars %in% colnames(df_bcr14_2020))]

# stage a unique spatial subset of the BCR by removing NAs for MODISLCC_1km
# also: treat water, barren, cropland MODIS classes as known features of the abiotic landscape
# i.e. we aren't trying to predict them. They are predictors.
# note: don't convert SCANFI_1km to a factor because xgb.DMatrix is expecting numerics
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
# in theory, this removes many locations with "rock" or "herb" SCANFI classes
# that are auto-correlated with urban areas and agriculture..but need to check.. 
CanHF_1km_present <- dplyr::filter(df_bcr14_2020_i, CanHF_1km > q50)
CanHF_1km_absent <- dplyr::filter(df_bcr14_2020_i, CanHF_1km <= q50) 

# subset low HF dataset to only those with abiotic predictors and the biotic response
# include the response "SCANFI_1km" so that it can used for labelling 
# note: # cropland is reduced from 14840 to 212 (98.5%)
predictor_vars <- c(abiotic_vars, "lon", "lat", "MODISLCC_1km_crop", "MODISLCC_1km_urban", "MODISLCC_1km_barren", "MODISLCC_1km_water")
CanHF_1km_absent_abiotic <- CanHF_1km_absent[,c(predictor_vars, "MODISLCC_1km")] 

# split data into training (80%) and hold-out (20%)
set.seed(123)
n <- nrow(CanHF_1km_absent_abiotic)
training_index <- sample(seq_len(n), size = round(0.80 * n))

training_data <- CanHF_1km_absent_abiotic[training_index, c(predictor_vars, "MODISLCC_1km")] 
holdout_data <- CanHF_1km_absent_abiotic[-training_index, c(predictor_vars, "MODISLCC_1km")]




#5. two-stage modelling----
# create broad SCANFI_1km classes
# all values for SCANFI_1km_rock and SCANFI_1km_water should be zero (checked by dplyr::count, and they are)

# classes 15 (snow/ice) and 18 (unclassified) are not in training data
dplyr::count(training_data, MODISLCC_1km)

training_data_stage1 <- 
  training_data |> 
  dplyr::mutate(MODISLCC_1km_broadclass = case_when(
    MODISLCC_1km %in% c(1:5) ~ 0,  # trees
    MODISLCC_1km %in% c(6:10) ~ 1,  # non-tree vegetation
    MODISLCC_1km == 11 ~ 2, # wetland 
    MODISLCC_1km %in% c(12:14) ~ 3, # human footprints 
    MODISLCC_1km %in% c(16:17) ~ 4)) |> # barren and water
  dplyr::mutate(MODISLCC_1km_broadclass = factor(MODISLCC_1km_broadclass, levels=c(0,1,2,3,4)))

# create DMatrix (a data structure for better speed when using xgboost)
# reminder: `predictor_vars` is `c(abiotic_vars, "lon", "lat", "SCANFI_1km_rock", "SCANFI_1km_water")`
dtrain_stage1 <- xgboost::xgb.DMatrix(data = as.matrix(training_data_stage1[,predictor_vars]), label = as.numeric(training_data_stage1$MODISLCC_1km_broadclass) - 1)

# set xgboost parameters
# train a model where each leaf of the boosted trees is associated with one of the SCANFI classes.
# then, compute "softmax" probabilities and return the class with the highest probability as the prediction.
# "softmax" is a function that generates a probability distribution over all possible classes 
# (compared to using a logisitic function for binary classification)
params_stage1 <- xgboost::xgb.params(objective = "multi:softmax", eval_metric = "mlogloss", max_depth = 3, eta = 0.1, num_class = length(levels(training_data_stage1$MODISLCC_1km_broadclass)))


# train stage 1 model (SCANFI non-trees vs trees vs abiotic)
# cv_stage1$early_stop$best_iteration 2182
cv_stage1 <-xgboost::xgb.cv(params = params_stage1, data = dtrain_stage1, nrounds=2500, nfold=5, early_stopping_rounds = 20)
model_stage1 <- xgboost::xgb.train(params = params_stage1, data = dtrain_stage1, nrounds = cv_stage1$early_stop$best_iteration, verbose = 0)



#6. build two multi-class models to predict real MODIS classes from 
# broad classes 0 (tree class) and 1 (non-tree vegetation class)
# NOTE: we won't try to split the other broad classes (wetlands, HF, water, unclassified)
# we'll instead just port them over directly as "determined"

# train a model to split broad class 0 into 5 MODIS classes (tree types)
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

cv_class0 <- xgboost::xgb.cv(params = params_class0, data = dtrain_class0, nrounds=2000, nfold=5, early_stopping_rounds = 20)
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

cv_class1 <- xgboost::xgb.cv(params = params_class1, data = dtrain_class1, nrounds=500, nfold=5, early_stopping_rounds = 20)
model_class1 <- xgb.train(params = params_class1, data = dtrain_class1, nrounds = cv_class1$early_stop$best_iteration, verbose = 0)




#7. now test the hierarchical models on the holdout data----
# use model_stage1 to predict broad MODIS classes
# if the prediction is broad class 0 (trees) pass the observation to model_class0 to get the specific tree type. 
# if the prediction is broad class 1 (non-tree vegetation) use model_class1 to determine whether it's shrub, grassland, etc.

# stage holdout data
# reminder: `predictor_vars` is `c(abiotic_vars, "lon", "lat", "SCANFI_1km_rock", "SCANFI_1km_water")`
# by filtering for vegetation, there are no rock or water classes to predict
holdout_data_stage1 <- 
  holdout_data |> 
  dplyr::mutate(MODISLCC_1km_broadclass = case_when(
    MODISLCC_1km %in% c(1:5) ~ 0,  # trees
    MODISLCC_1km %in% c(6:10) ~ 1,  # non-tree vegetation
    MODISLCC_1km == 11 ~ 2, # wetland 
    MODISLCC_1km %in% c(12:14) ~ 3, # human footprints 
    MODISLCC_1km %in% c(16:17) ~ 4)) |> # barren and water
  dplyr::mutate(MODISLCC_1km_broadclass = factor(MODISLCC_1km_broadclass, levels=c(0,1,2,3,4)))

dholdout <- xgboost::xgb.DMatrix(data = as.matrix(holdout_data_stage1[,predictor_vars]), label=as.numeric(holdout_data_stage1$MODISLCC_1km_broadclass) - 1)

# stage 1: estimate broad SCANFI classes from holdout data
holdout_data_stage1$MODISLCC_1km_broadclass_predicted <- predict(model_stage1, newdata = dholdout)


# stage 2A:
# for broad class 0 predictions:
# filter holdout data for further subclass predictions
holdout_class0 <- dplyr::filter(holdout_data_stage1, MODISLCC_1km_broadclass_predicted == 0)
dholdout_class0 <- xgboost::xgb.DMatrix(data = as.matrix(holdout_class0[, predictor_vars]))
pred_class0 <- predict(model_class0, newdata = dholdout_class0)

# convert numeric predictions (0, 1, 2) to original SCANFI classes using the factor mapping from training
unique(pred_class0) # 0 4 3
class0_labels <- factor(pred_class0, levels = c(0,3,4), labels = c("1", "4", "5"))
holdout_class0 <- dplyr::mutate(holdout_class0, MODISLCC_1km_prediction = class0_labels)


# stage 2B:
# for broad class 1 predictions:
# filter holdout data for further subclass predictions
holdout_class1 <- dplyr::filter(holdout_data_stage1, MODISLCC_1km_broadclass_predicted == 1)
dholdout_class1 <- xgboost::xgb.DMatrix(data = as.matrix(holdout_class1[, predictor_vars]))
pred_class1 <- predict(model_class1, newdata = dholdout_class1)

# convert numeric predictions (0, 1, 2) to original SCANFI classes using the factor mapping from training
unique(pred_class1) # 2 1 4 3
class1_labels <- factor(pred_class1, levels = c(1,2,3,4), labels = c("7","8","9","10"))
holdout_class1 <- dplyr::mutate(holdout_class1, MODISLCC_1km_prediction = class1_labels)


# stage 2C:
# stage determined MODIS classes for `bind_rows`
# filter for wetland, cropland, urban, crop/natural,snow/ice, barren, water, unclassified
holdout_determined <- 
  holdout_data_stage1 |> 
  dplyr::filter(MODISLCC_1km %in% c(11:18)) |> 
  dplyr::mutate(MODISLCC_1km_prediction = as.character(MODISLCC_1km))


# combine predictions into single data frame
prediction_stage2 <- dplyr::bind_rows(holdout_class0, holdout_class1, holdout_determined)




#8. evaluate model performance on holdout data----
confusion_matrix_broadclass <- table(actual = holdout_data_stage1$MODISLCC_1km_broadclass, predicted = holdout_data_stage1$MODISLCC_1km_broadclass_predicted)
confusion_matrix_broadclass
sum(diag(confusion_matrix_broadclass)) / sum(confusion_matrix_broadclass) # 0.8404
          #predicted
#actual     0     1     2     3     4
#0      19527   834     3     4     0
#1       3232  2015     3     4     0
#2          6     5     8     0     0
#3         16     5     0    45     0
#4          0     0     0     0    61


confusion_matrix_fineclass <- 
  table(actual = factor(prediction_stage2$MODISLCC_1km, levels = 1:18),
        predicted = factor(prediction_stage2$MODISLCC_1km_prediction, levels = 1:18))

confusion_matrix_fineclass
sum(diag(confusion_matrix_fineclass)) / sum(confusion_matrix_fineclass) # 0.7637

#       predicted
#actual 1     2     3     4     5     6     7     8     9    10    11    12    13    14    15    16    17    18
#1    803     0     0     3   503     0     0   172     1     0     0     0     0     0     0     0     0     0
#2      0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
#3      0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
#4      1     0     0  1123   879     0     0    73     0     0     0     0     0     0     0     0     0     0
#5    160     0     0   351 15704     0     0   584     3     1     0     0     0     0     0     0     0     0
#6      0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
#7      3     0     0     0     4     0     4    14     0     0     0     0     0     0     0     0     0     0
#8    226     0     0   212  2721     0     3  1891     5     0     0     0     0     0     0     0     0     0
#9      1     0     0     1    14     0     0    24     9     1     0     0     0     0     0     0     0     0
#10     3     0     0     2    45     0     1    48     2    13     0     0     0     0     0     0     0     0
#11     2     0     0     0     4     0     0     4     0     1    19     0     0     0     0     0     0     0
#12     0     0     0     0     1     0     0     0     0     0     0    35     0     0     0     0     0     0
#13     0     0     0     0     1     0     0     0     0     0     0     0     4     0     0     0     0     0
#14     0     0     0     8     6     0     0     3     2     0     0     0     0    27     0     0     0     0
#15     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
#16     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
#17     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0    61     0
#18     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
