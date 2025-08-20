# ---
# title: National Models 5.0 - create models to predict biotic VLCE classes from abiotic landscape
# author: Mannfred Boehm
# created: March 18, 2025
# ---



#1. run lines 1-134 in "test_backfill_SCANFI_bcr14.R"---- 
# this will import covariate rasters, and generate a dataframe `df_bcr14_2020`from 
# the predictor covariate stack "can14_2020.tif"
# We also set the threshold for low/high human footprint by defining `q50`

#2. inspect Virtual Land Cover Engine classes----
# from Hermosilla et al: https://drive.google.com/drive/u/1/folders/1KyUkTvc8VpiVhWX_GoCI6gNNWJSG7kZ9
plot(as.factor(stack_bcr14_2020$VLCE_1km))

# 0 = no change
# 20 = water
# 31 = snow_ice
# 32 = rock_rubble
# 33 = exposed_barren_land
# 40 = bryoids
# 50 = shrubs
# 80 = wetland
# 81 = wetland-treed
# 100 = herbs
# 210 = coniferous
# 220 = broadleaf
# 230 = mixedwood



#3. test whether VLCE_1km land cover classes are accurately predicted---- 
# will eventually need to do the same for SCANFI_1km, VLCE_1km

# drop biotic predictors that were thinned by VIF (see line 300 in "04.Stratify.R")
biotic_vars_thinned <- biotic_vars[which(biotic_vars %in% colnames(df_bcr14_2020))]

# stage a unique spatial subset of the BCR by removing NAs for VLCE_1km
# we'll treat water, barren, cropland VLCE classes as known features of the abiotic landscape
df_bcr14_2020_i <-  tidyr::drop_na(df_bcr14_2020, VLCE_1km)

# does every VLCE class appear in BCR14? 
dplyr::count(df_bcr14_2020_i, VLCE_1km)

# identify pixels with "high" and "low" human footprint 
CanHF_1km_present <- dplyr::filter(df_bcr14_2020_i, CanHF_1km > q50)
CanHF_1km_absent <- dplyr::filter(df_bcr14_2020_i, CanHF_1km <= q50) 

# this removes most locations classed as "barren" (17009 to 3943, 76.8%) or 
# "no change" (26519 to 1727, 93.5%)
dplyr::count(df_bcr14_2020_i, VLCE_1km)
dplyr::count(CanHF_1km_absent, VLCE_1km)


# subset low HF dataset to abiotic predictors and the biotic response
# include the response "VLCE_1km" so that it can used for labelling 
predictor_vars <- c(abiotic_vars, "lon", "lat")
CanHF_1km_absent_abiotic <- CanHF_1km_absent[,c(predictor_vars, "VLCE_1km")] 

# split data into training (80%) and hold-out (20%)
set.seed(123)
n <- nrow(CanHF_1km_absent_abiotic)
training_index <- sample(seq_len(n), size = round(0.80 * n))

training_data <- CanHF_1km_absent_abiotic[training_index, c(predictor_vars, "VLCE_1km")] 
holdout_data <- CanHF_1km_absent_abiotic[-training_index, c(predictor_vars, "VLCE_1km")]

# abiotic: class 31 (snow/ice) is not in dataset
# biotic:  classes 2 (evergreen broadleaf) and 6 (closed shrub) are not in the training or holdout data
dplyr::count(training_data, VLCE_1km)
dplyr::count(holdout_data, VLCE_1km)



#4. two-stage modelling----
# create broad VLCE_1km classes

training_data_stage1 <- 
  training_data |> 
  dplyr::filter(VLCE_1km != 20) |>  # remove water areas (they are static so we will directly copy them over later)
  dplyr::mutate(VLCE_1km_broadclass = case_when(
    VLCE_1km == 0 ~ 0,  # no change
    VLCE_1km %in% c(32:33) ~ 1,  # "natural" rock and barren
    VLCE_1km %in% c(40, 50, 80, 81, 100) ~ 2, #non-tree vegetation
    VLCE_1km %in% c(210:230) ~ 3)) |>  # trees 
  dplyr::mutate(VLCE_1km_broadclass = factor(VLCE_1km_broadclass, levels=c(0,1,2,3)))

count(training_data_stage1, VLCE_1km_broadclass)
# VLCE_1km_broadclass     n
# <fct>               <int>
# 1 0                    1365
# 2 1                    3281
# 3 2                    6824
# 4 3                   89850

# create DMatrix (a data structure for better speed when using xgboost)
# reminder: `predictor_vars` is `c(abiotic_vars, "lon", "lat")
dtrain_stage1 <- xgboost::xgb.DMatrix(data = as.matrix(training_data_stage1[,predictor_vars]), label = as.numeric(training_data_stage1$VLCE_1km_broadclass) - 1)
params_stage1 <- xgboost::xgb.params(objective = "multi:softmax", eval_metric = "mlogloss", max_depth = 3, eta = 0.1, num_class = length(levels(training_data_stage1$VLCE_1km_broadclass)))


# train stage 1 model (VLCE no change vs "natural" rocks vs vegetation)
# cv_stage1$early_stop$best_iteration 513 
cv_stage1 <-xgboost::xgb.cv(params = params_stage1, data = dtrain_stage1, nrounds=1000, nfold=5, early_stopping_rounds = 20)
model_stage1 <- xgboost::xgb.train(params = params_stage1, data = dtrain_stage1, nrounds = cv_stage1$early_stop$best_iteration, verbose = 0)



#5. build three multi-class models to predict real VLCE classes from 
# broad classes 1 (rock/barren), 2 (wetland vs non-tree vegetation classes) 
# and 3 (tree type). 
# NOTE: we won't try to split the (low) human footprint broad class 3
# we'll instead just port them over directly as "determined".

# train a model to split broad class 1 into 2 classes (rock or barren)
training_data_stage2_class1 <- 
  training_data_stage1 |> 
  dplyr::filter(VLCE_1km_broadclass == 1) |> 
  dplyr::mutate(VLCE_1km = factor(VLCE_1km, levels = c(32, 33)))

dtrain_class1 <- xgb.DMatrix(
  data = as.matrix(training_data_stage2_class1[, predictor_vars]), 
  label = as.numeric(training_data_stage2_class1$VLCE_1km) - 1) # as.numeric() starts the classes at 1 but xgboost labels need to start at 0

params_class1 <- 
  xgboost::xgb.params(objective = "multi:softmax", 
                      num_class = length(levels(training_data_stage2_class1$VLCE_1km)),
                      eval_metric = "mlogloss",
                      max_depth = 3,
                      eta = 0.05)

cv_class1 <- xgboost::xgb.cv(params = params_class1, data = dtrain_class1, nrounds=200, nfold=5, early_stopping_rounds = 20)
model_class1 <- xgb.train(params = params_class1, data = dtrain_class1, nrounds = cv_class1$early_stop$best_iteration, verbose = 0)


# train a model to split broad class 2 into 5 classes (bryoid, shrub, herb, wetland, wetland-tree)
training_data_stage2_class2 <- 
  training_data_stage1 |> 
  dplyr::filter(VLCE_1km_broadclass == 2) |> 
  dplyr::mutate(VLCE_1km = factor(VLCE_1km, levels = c(40, 50, 80, 81, 100)))

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
cv_class2 <- xgboost::xgb.cv(params = params_class2, data = dtrain_class2, nrounds=1000, nfold=5, early_stopping_rounds = 20)
model_class2 <- xgb.train(params = params_class2, data = dtrain_class2, nrounds = cv_class2$early_stop$best_iteration, verbose = 0)


# train a model to split broad class 3 into 3 classes (conifer, broadleaf, mixed)
training_data_stage2_class3 <- 
  training_data_stage1 |> 
  dplyr::filter(VLCE_1km_broadclass == 3) |> 
  dplyr::mutate(VLCE_1km = factor(VLCE_1km, levels = c(210, 220, 230)))

dtrain_class3 <- xgb.DMatrix(
  data = as.matrix(training_data_stage2_class3[, predictor_vars]), 
  label = as.numeric(training_data_stage2_class3$VLCE_1km) - 1) # as.numeric() starts the classes at 1 but xgboost labels need to start at 0

params_class3 <- 
  xgboost::xgb.params(objective = "multi:softmax", 
                      num_class = length(levels(training_data_stage2_class3$VLCE_1km)),
                      eval_metric = "mlogloss",
                      max_depth = 3,
                      eta = 0.05)

# cv_class3$early_stop$best_iteration 
cv_class3 <- xgboost::xgb.cv(params = params_class3, data = dtrain_class3, nrounds=2000, nfold=5, early_stopping_rounds = 20)
model_class3 <- xgb.train(params = params_class3, data = dtrain_class3, nrounds = cv_class3$early_stop$best_iteration, verbose = 0)




#6. now test the hierarchical models on the holdout data----
# use model_stage1 to predict broad VLCE classes
# if the prediction is broad class 0 (no change) there is no further subdividing
# if the prediction is broad class 1 (rock/barren), pass the observation to model_class1 to get the specific VLCE class. 
# if the prediction is broad class 2 (non-tree vegetation) use model_class2 to determine whether it's shrub, herb, etc.
# if the prediction is broad class 3 (trees) use model_class3 to determine if it's broadleaf, conifer, or mixed 

# stage holdout data
# reminder: `predictor_vars` is `c(abiotic_vars, "lon", "lat", "SCANFI_1km_rock", "SCANFI_1km_water")`
# by filtering for vegetation, there are no rock or water classes to predict
holdout_data_stage1 <- 
  holdout_data |> 
  dplyr::filter(VLCE_1km != 20) |>  # remove water areas (they are static so we will directly copy them over later)
  dplyr::mutate(VLCE_1km_broadclass = case_when(
    VLCE_1km == 0 ~ 0,  # no change
    VLCE_1km %in% c(32:33) ~ 1,  # "natural" rock and barren
    VLCE_1km %in% c(40, 50, 80, 81, 100) ~ 2, #non-tree vegetation
    VLCE_1km %in% c(210:230) ~ 3)) |>  # trees 
  dplyr::mutate(VLCE_1km_broadclass = factor(VLCE_1km_broadclass, levels=c(0,1,2,3)))

dplyr::count(holdout_data_stage1, VLCE_1km_broadclass)
# VLCE_1km_broadclass     n
# <fct>                     <int>
# 1 0                       362
# 2 1                       751
# 3 2                       1746
# 4 3                       22436

dholdout <- xgboost::xgb.DMatrix(data = as.matrix(holdout_data_stage1[,predictor_vars]), label=as.numeric(holdout_data_stage1$VLCE_1km_broadclass) - 1)

# stage 1: estimate broad VLCE classes from holdout data
holdout_data_stage1$VLCE_1km_broadclass_predicted <- predict(model_stage1, newdata = dholdout)
unique(holdout_data_stage1$VLCE_1km_broadclass_predicted) # 3 1 2 0 all broad classes were predicted


# stage 2A:
# for broad class 1 predictions:
# filter holdout data for further subclass predictions
holdout_class1 <- dplyr::filter(holdout_data_stage1, VLCE_1km_broadclass_predicted == 1)
dholdout_class1 <- xgboost::xgb.DMatrix(data = as.matrix(holdout_class1[, predictor_vars]))
pred_class1 <- predict(model_class1, newdata = dholdout_class1)

# convert numeric predictions (0, 1) to original VLCE classes using the factor mapping from training
unique(pred_class1) # 1 only barren predicted
class1_labels <- factor(pred_class1, levels = c(10, labels = c("33")))
holdout_class1 <- dplyr::mutate(holdout_class1, VLCE_1km_prediction = class1_labels)


# stage 2B:
# for broad class 1 predictions:
# filter holdout data for further subclass predictions
holdout_class2 <- dplyr::filter(holdout_data_stage1, VLCE_1km_broadclass_predicted == 2)
dholdout_class2 <- xgboost::xgb.DMatrix(data = as.matrix(holdout_class2[, predictor_vars]))
pred_class2 <- predict(model_class2, newdata = dholdout_class2)

# convert numeric predictions (0, 1) to original VLCE classes using the factor mapping from training
unique(pred_class2) # 2 1  only barren predicted
class1_labels <- factor(pred_class1, levels = c(1,2), labels = c("50","80"))
holdout_class1 <- dplyr::mutate(holdout_class1, VLCE_1km_prediction = class1_labels)


# stage 2C:
# for broad class 1 predictions:
# filter holdout data for further subclass predictions
holdout_class3 <- dplyr::filter(holdout_data_stage1, VLCE_1km_broadclass_predicted == 3)
dholdout_class3 <- xgboost::xgb.DMatrix(data = as.matrix(holdout_class3[, predictor_vars]))
pred_class3 <- predict(model_class3, newdata = dholdout_class3)

# convert numeric predictions (0, 1) to original VLCE classes using the factor mapping from training
unique(pred_class3) # 0 2 1 all three tree types predicted
class3_labels <- factor(pred_class3, levels = c(0,1,2), labels = c("210","220","230"))
holdout_class3 <- dplyr::mutate(holdout_class3, VLCE_1km_prediction = class3_labels)


# stage 2D:
# stage broad class 0 predictions as VLCE class predictions
holdout_class0 <-
  holdout_data_stage1 |> 
  dplyr::filter(VLCE_1km_broadclass_predicted == 0) |> 
  dplyr::mutate(VLCE_1km_prediction = as.character(VLCE_1km_broadclass_predicted))


# stage 2E:
# stage determined VLCE classes for `bind_rows`
# filter for water
holdout_determined <- 
  holdout_data|> 
  dplyr::filter(VLCE_1km == 20) |>
  dplyr::mutate(VLCE_1km_prediction = as.character(VLCE_1km))

count(holdout_determined, VLCE_1km_prediction)

# combine predictions into single data frame
prediction_stage2 <- dplyr::bind_rows(holdout_class0, holdout_class1, holdout_class2, holdout_class3, holdout_determined)



#7. evaluate model performance on holdout data----
confusion_matrix_broadclass <- table(actual = holdout_data_stage1$VLCE_1km_broadclass, predicted = holdout_data_stage1$VLCE_1km_broadclass_predicted)
confusion_matrix_broadclass
sum(diag(confusion_matrix_broadclass)) / sum(confusion_matrix_broadclass) # 0.8403

#      predicted
#actual 0     1     2     3
# 0   343     0     0    19
# 1     1     8     9   733
# 2     1     4    72  1669
# 3    10     1    40 22385

vlce_levels <- c("0","20","31","32","33","40","50","80","81","100","210","220","230")
confusion_matrix_fineclass <- 
  table(actual = factor(prediction_stage2$VLCE_1km, levels = vlce_levels),
        predicted = factor(prediction_stage2$VLCE_1km_prediction, levels = vlce_levels))

# NOTE: wetlands are not successfully predicted, so might need to treat as determined
confusion_matrix_fineclass
sum(diag(confusion_matrix_fineclass)) / sum(confusion_matrix_fineclass) # 0.562

#       predicted
#actual  0   20   31   32   33   40   50   80   81  100  210  220  230
# 0    343    0    0    0    0    0    0    0    0    0    4   10    5
# 20     0  473    0    0    0    0    0    0    0    0    0    0    0
# 31     0    0    0    0    0    0    0    0    0    0    0    0    0
# 32     0    0    0    0    0    0    0    0    0    0    9    0    7
# 33     1    0    0    0    0    0    8    0    0    0  368   48  300
# 40     0    0    0    0    0    0    0    0    0    0    7    0    2
# 50     0    0    0    0    0    0    0    0    0    0  673  176  408
# 80     0    0    0    0    0    0    3    0    0    0  110    5   52
# 81     0    0    0    0    0    0    0    0    0    0   16    0    3
# 100    1    0    0    0    0    0    1    0    0    0   73   76   67
# 210    2    0    0    0    0    0    1    0    0    0 7170  590 2550
# 220    3    0    0    0    0    0    0    0    0    0  753 2022 1086
# 230    5    0    0    0    0    0    0    0    0    0 3088  712 4408
