# ---
# title: Impact Assessment: create full covariate stacks for each BCR (for training)
# author: Mannfred Boehm
# created: May 7, 2025
# ---


library(tidyverse)
library(terra)

# set root path
root <- "G:/Shared drives/BAM_NationalModels5"

# import soil covariates
soil_covariates <- terra::rast(file.path(root, "gis", "other_landscape_covariates", "isric_soil_covariates_masked.tif"))

# import time-since-disturbance
CAfire <- terra::rast(file.path(root, "gis", "other_landscape_covariates", "CA_Forest_Fire_1985-2020_masked.tif"))


#2. load covariate stack_i and process----
# A) import V5 stack_i
# A) add soil and time-since-disturbance layers to stack_i
# B) convert stack_i into dataframe_i
# C) populate a list where every element is a dataframe_i from stack_i

# NEXT SCRIPT(s)
# C) create abiotic_vars and biotic_vars indices
# D) split dataframe_i into low_HF and high_HF
# E) train model to predict biotic_var_i on low_HF dataset
# F) identify sector to backfill in high_HF dataset
# G) backfill sector

# create file index where stacks are located (264 stacks found)
stack_directories <- list.files(file.path(root, "gis", "stacks"), pattern = "*\\.tif$", full.names = TRUE, recursive = FALSE)

# define function for lapply
stack_and_enframe <- function(stack_i) {
  
  # import stack_i from index
  stack_i <- terra::rast(stack_directories[i])
  
  # choose a random layer from stack_i to represent resolution, boundaries, etc.
  ref_layer <- stack_i$year
  
  # align extra covariates to stack_i
  # then, crop and mask extra covariates to the current BCR
  CAfire_i <- 
    terra::resample(x=CAfire, y=ref_layer, method="near") |> 
    terra::crop(x=_, y=ref_layer) |> 
    terra::mask(x=_, mask=ref_layer)
 
 
  soil_covariates_i <- 
    terra::resample(x=soil_covariates, y=stack_i) |> 
    terra::crop(x=_, y=stack_i) |> 
    terra::mask(x=_, mask=stack_i) 
    

  # add time since disturbance layer and soil layers to covariate stack
  stack_i_extra <- c(stack_i, CAfire_i, soil_covariates_i)
  
  # convert to raster stack to a dataframe
  stack_i_df <-
    terra::as.data.frame(x = stack_i_extra, xy = TRUE) |>
    tibble::as_tibble() |> 
    dplyr::rename(lon=x, lat=y) 

}

# NEED NAMES FOR KEEPING TRACK OF BCRS AND YEARS
# for every BCR x year, stack extra covariates to V5 covariates and turn into a dataframe
list_of_cov_dfs <- lapply(X = stack_directories, FUN = stack_and_enframe)


#3. stage datasets needed for modelling biotic landscape----

# import variable classes to 
# A) index all possible covariates used as predictors and 
# B) help separate biotic from abiotic covariates
nice_var_names <- 
  readr::read_csv(file.path(root, "data", "Extras", "sandbox_data", "impactassessment_sandbox", "nice_var_names_v5.csv")) |> 
  tidyr::drop_na() |> 
  unique()

# remove Peatland as they are considered biotic,
# remove *_5x5s since they're redundant with *_1km
# for hli3cl see line 67 in "08.CalculateExtrapolation.R"
abiotic_vars <-
  nice_var_names |> 
  dplyr::filter(var_class %in% c("Annual Climate", "Climate Normals", "Topography", "Wetland")) |> 
  tibble::add_row(var = "CAfire", var_class ="Time Since Disturbance") |> 
  tibble::add_row(var = "soil_carbon", var_class = "Soils") |> 
  tibble::add_row(var = "soil_ph", var_class = "Soils") |> 
  dplyr::filter(!(var %in% c("hli3cl_1km", "Peatland_1km"))) |> 
  dplyr::filter(!grepl("5x5", var)) |>  
  dplyr::pull(var) |> 
  unique()

# remove 5x5 since I can re-create these variables by
# scaling up from 1km after prediciton
# exclude var_class == "Landcover" in `filter()`
biotic_vars <-
  nice_var_names |> 
  dplyr::filter(var_class %in% c("Greenup", "Biomass", "Wetland", "Landcover")) |> 
  dplyr::filter(!(var %in% c("WetOccur_1km", "WetOccur_5x5", "WetRecur_1km", "WetSeason_1km"))) |> # keep peatland but discard other water variables
  dplyr::filter(!grepl("5x5", var)) |> 
  dplyr::pull(var) |> 
  unique()


# convert covariate stack to a dataframe 
# rows are XYZ, columns are XYZ
df_bcr14_2020 <-
  stack_bcr14_2020 |> 
  terra::as.data.frame(xy=TRUE) |>
  tibble::as_tibble() |> 
  dplyr::rename(lon=x, lat=y) 
  
# define threshold for "low" human footprint
q50 <- quantile(df_bcr14_2020$CanHF_1km, probs = 0.50, na.rm = TRUE)

# how does human footprint vary across this BCR?
# hist(df_bcr14_2020$CanHF_1km, main="CanHF_1km in BCR14")
# abline(v=q50, col="darkred", lwd=2)
# abline(v=mean(na.omit(df_bcr14_2020$CanHF_1km)), col="skyblue", lwd=2, lty="dashed")
       
quantile(na.omit(df_bcr14_2020$CanHF_1km))
#0%        25%        50%        75%       100% 
#0.0000000  0.9322898  5.2230940  9.6365070 54.5055199 





#4. test whether SCANFI_1km land cover classes are accurately predicted---- 
# will eventually need to do the same for VLCE_1km, MODISLCC_1km

# drop biotic predictors that were thinned by VIF (see line 300 in "04.Stratify.R")
biotic_vars_thinned <- biotic_vars[which(biotic_vars %in% colnames(df_bcr14_2020))]

# stage a unique spatial subset of the BCR by removing NAs for SCANFI_1km
# note: don't convert SCANFI_1km to a factor because xgb.DMatrix is expecting numerics
df_bcr14_2020_i <-  tidyr::drop_na(df_bcr14_2020, SCANFI_1km)
 

# identify pixels with "high" and "low" human footprint 
# in theory, this removes many locations with "rock" or "herb" SCANFI classes
# that are auto-correlated with urban areas and agriculture..but need to check.. 
CanHF_1km_present <- dplyr::filter(df_bcr14_2020_i, CanHF_1km > q50)
CanHF_1km_absent <- dplyr::filter(df_bcr14_2020_i, CanHF_1km <= q50) 

rm(nice_var_names); gc()

# subset low HF dataset to only those with abiotic predictors and the biotic response
# include the response "SCANFI_1km" so that it can used for labelling 
# in theory, urban areas have been filtered out so SCANFI "rock" is natural rock
predictor_vars <- c(abiotic_vars, "lon", "lat") 
CanHF_1km_absent_abiotic <- CanHF_1km_absent[,c(predictor_vars, "SCANFI_1km")]

# split data into training (80%) and hold-out (20%)
set.seed(123)
n <- nrow(CanHF_1km_absent_abiotic)
training_index <- sample(seq_len(n), size = round(0.80 * n))

training_data <- CanHF_1km_absent_abiotic[training_index, c(predictor_vars, "SCANFI_1km")] 
holdout_data <- CanHF_1km_absent_abiotic[-training_index, c(predictor_vars, "SCANFI_1km")]




#5. two-stage modelling----
# first, split SCANFI_1km classes between tree and non-tree vegetation
training_data_stage1 <- 
  CanHF_1km_absent_abiotic |> 
  dplyr::filter(SCANFI_1km != "8") |> # remove water areas (they are static so we will directly copy them over later)
  dplyr::mutate(SCANFI_1km_broadclass = case_when(
    SCANFI_1km %in% c(1, 2, 4) ~ 0,  # bryoid, herbs, shrub (non-trees) 
    SCANFI_1km %in% c(5, 6, 7) ~ 1,  # broadleaf, conifer, mixed (trees)
    SCANFI_1km == 3 ~ 2)) |> # low HF ("natural") rock areas
  dplyr::mutate(SCANFI_1km_broadclass = factor(SCANFI_1km_broadclass, levels=c(0,1,2)))

count(training_data_stage1, SCANFI_1km_broadclass)
# SCANFI_1km_broadclass      n
# <fct>                  <int>
#  0                       5255
#  1                     120646
#  2                        527

# create DMatrix (a data structure for better speed when using xgboost)
# reminder: `predictor_vars` is `c(abiotic_vars, "lon", "lat"), i.e. excludes `SCANFI_1km`
dtrain_stage1 <- xgboost::xgb.DMatrix(data = as.matrix(training_data_stage1[,predictor_vars]), label = as.numeric(training_data_stage1$SCANFI_1km_broadclass) - 1)

# set xgboost parameters
# train a model where each leaf of the boosted trees is associated with one of the SCANFI classes.
# then, compute "softmax" probabilities and return the class with the highest probability as the prediction.
# "softmax" is a function that generates a probability distribution over all possible classes 
# (compared to using a logisitic function for binary classification)
params_stage1 <- xgboost::xgb.params(objective = "multi:softmax", eval_metric = "mlogloss", max_depth = 3, eta = 0.05, num_class = length(levels(training_data_stage1$SCANFI_1km_broadclass)))

# train stage 1 model (SCANFI non-trees vs trees vs abiotic)
cv_stage1 <-xgboost::xgb.cv(params = params_stage1, data = dtrain_stage1, nrounds=1000, nfold=5, early_stopping_rounds = 20)
model_stage1 <- xgboost::xgb.train(params = params_stage1, data = dtrain_stage1, nrounds = cv_stage1$early_stop$best_iteration, verbose = 0)



#6. build two multi-class models to learn to split SCANFI broadclasses----
# split broad class 0 into bryoid, herb, shrub
# split broad class 1 into broadleaf, conifer, mixed
# note: no model is built for broad class 3 because if model_stage1 predicts rock, 
# then it's rock (i.e. can't be split further)

# filter data to broad class 0 (non tree vegetation)
# SCANFI_1km is now a 3-factor column for non-treed classes (bryoids, herbs, shrubs)
# unique(training_data_stage2A$SCANFI_1km) [1] 4 1 2 
training_data_stage2A <- 
  training_data_stage1 |> # water areas already removed
  dplyr::filter(SCANFI_1km_broadclass == 0) |> 
  dplyr::mutate(SCANFI_1km = factor(SCANFI_1km, levels = c(1, 2, 4)))

count(training_data_stage2A, SCANFI_1km)
# SCANFI_1km     n
# <fct>        <int>
# 1 1            154
# 2 2           1017
# 3 4           4084


dtrain_non_tree <- xgb.DMatrix(
  data = as.matrix(training_data_stage2A[, predictor_vars]), # exclude SCANFI_1km and SCANFI_1km_tree
  label = as.numeric(training_data_stage2A$SCANFI_1km) - 1) # as.numeric() starts the classes at 1 but xgboost labels need to start at 0

params_non_tree <- 
  xgboost::xgb.params(objective = "multi:softmax", 
                      num_class = length(levels(training_data_stage2A$SCANFI_1km)),
                      eval_metric = "mlogloss",
                      max_depth = 3,
                      eta = 0.05)

cv_notree <- xgboost::xgb.cv(params = params_non_tree, data = dtrain_non_tree, nrounds=1000, nfold=5, early_stopping_rounds = 20)
model_stage2_notree <- xgb.train(params = params_non_tree, data = dtrain_non_tree, nrounds = cv_notree$early_stop$best_iteration, verbose = 0)



# filter data for broad class 1
# SCANFI_1km is now a 3-factor column for yes-treed classes (broadlead, conifer, mixed)
# unique(training_data_stage2B$SCANFI_1km) [1] 6 7 5
training_data_stage2B <- 
  training_data_stage1 |> 
  dplyr::filter(SCANFI_1km_broadclass == 1) |> 
  dplyr::mutate(SCANFI_1km = factor(SCANFI_1km, levels = c(5, 6, 7)))

# count(training_data_stage2B, SCANFI_1km)
# SCANFI_1km     n
# <fct>       <int>
# 1 5          24593
# 2 6          58274
# 3 7          37779

dtrain_yes_tree <- xgb.DMatrix(
  data = as.matrix(training_data_stage2B[, predictor_vars]),
  label = as.numeric(training_data_stage2B$SCANFI_1km) - 1)

params_yes_tree <- 
  xgboost::xgb.params(objective = "multi:softmax",
                      num_class = length(levels(training_data_stage2B$SCANFI_1km)),
                      eval_metric = "mlogloss",
                      max_depth = 3,
                      eta = 0.05)

cv_yestree <- xgboost::xgb.cv(params = params_yes_tree, data = dtrain_yes_tree, nrounds=2000, nfold=5, early_stopping_rounds = 20)
model_stage2_yestree <- xgboost::xgb.train(params = params_yes_tree, data = dtrain_yes_tree, nrounds = cv_yestree$early_stop$best_iteration, verbose = 0)



#7. now test the 2-stage models on the holdout data----
# use model_stage1 to predict whether a point is treed or non-treed:
# if the prediction is non-treed, pass the observation to model_non_tree to get the specific class (bryoid, herb, or shrub).
# if the prediction is treed, use model_yes_tree to determine whether it's broadleaf, conifer, or mixed.

# stage holdout data
# reminder: `predictor_vars` is `c(abiotic_vars, "lon", "lat", "SCANFI_1km_rock", "SCANFI_1km_water")`
holdout_data_stage1 <- 
  holdout_data |> 
  dplyr::filter(SCANFI_1km != "8") |> # remove 2336 water areas (we can use terra::cover to replace them later)
  dplyr::mutate(SCANFI_1km_broadclass = case_when(
    SCANFI_1km %in% c(1, 2, 4) ~ 0,  # bryoid, herbs, shrub (non-trees) 
    SCANFI_1km %in% c(5, 6, 7) ~ 1,  # broadleaf, conifer, mixed (trees)
    SCANFI_1km == 3 ~ 2)) |> # low HF ("natural") rock areas
  dplyr::mutate(SCANFI_1km_broadclass = factor(SCANFI_1km_broadclass, levels=c(0,1,2)))

count(holdout_data_stage1, SCANFI_1km_broadclass)
# SCANFI_1km_broadclass     n
# <fct>                 <int>
# 1 0                      1071
# 2 1                     24122
# 3 2                       112

dholdout <- xgboost::xgb.DMatrix(data = as.matrix(holdout_data_stage1[,predictor_vars]), label=as.numeric(holdout_data_stage1$SCANFI_1km_broadclass) - 1)

# stage 1: estimate broad SCANFI classes from holdout data
# all three broad classes were predicted at least once:
# unique(holdout_data_stage1$SCANFI_1km_broadclass_predicted) 2 1 0
holdout_data_stage1$SCANFI_1km_broadclass_predicted <- predict(model_stage1, newdata = dholdout)
 
# stage 2A:
# filter holdout data for further subclass predictions
holdout_non_tree <- dplyr::filter(holdout_data_stage1, SCANFI_1km_broadclass_predicted == 0)
dholdout_non_tree <- xgboost::xgb.DMatrix(data = as.matrix(holdout_non_tree[, predictor_vars]))
pred_non_tree <- predict(model_stage2_notree, newdata = dholdout_non_tree) # unique(pred_non_tree) 2 1 (no bryoids predicted)


# convert numeric predictions (0, 1, 2) to original SCANFI classes 
non_tree_labels <- factor(pred_non_tree, levels = c(0,1,2), labels = c("1", "2", "4"))
holdout_non_tree <- dplyr::mutate(holdout_non_tree, SCANFI_1km_prediction = non_tree_labels)
# unique(holdout_non_tree$SCANFI_1km_prediction) # 4 2 confirm no bryoids predicted


# stage 2B:
# for areas predicted as broad class 1 (treed):
# filter holdout data for further subclass predictions
holdout_yes_tree <- dplyr::filter(holdout_data_stage1, SCANFI_1km_broadclass_predicted == 1)
dholdout_yes_tree <- xgboost::xgb.DMatrix(data = as.matrix(holdout_yes_tree[, predictor_vars]))
pred_yes_tree <- predict(model_stage2_yestree, newdata = dholdout_yes_tree) # unique(pred_yes_tree) 1 0 2 (all three tree types predicted)

# convert numeric predictions (0, 1, 2) to original SCANFI classes 
tree_labels <- factor(pred_yes_tree, levels = c(0,1,2), labels = c("5", "6", "7"))
holdout_yes_tree <- dplyr::mutate(holdout_yes_tree, SCANFI_1km_prediction = tree_labels) 
# unique(holdout_yes_tree$SCANFI_1km_prediction) # 6 5 7 confirm all tree types predicted


# stage 2C: areas predicted as broad class 2 (equivalent to SCANFI_1km == 3 == "rock")
holdout_rock <- 
  holdout_data_stage1 |> 
  dplyr::filter(SCANFI_1km_broadclass_predicted == 2) |>  
  dplyr::mutate(SCANFI_1km_prediction = "3")

# stage 2D: retain water areas and port over into holdout dataset
holdout_water <-
  holdout_data |> 
  dplyr::filter(SCANFI_1km == 8) |> 
  dplyr::mutate(SCANFI_1km_prediction = "8")

# combine predictions into single data frame
prediction_stage2 <- dplyr::bind_rows(holdout_non_tree, holdout_yes_tree, holdout_rock, holdout_water)




#8. evaluate model performance on holdout data----
confusion_matrix_broadclass <- table(actual = holdout_data_stage1$SCANFI_1km_broadclass, predicted = holdout_data_stage1$SCANFI_1km_broadclass_predicted)
confusion_matrix_broadclass
sum(confusion_matrix_broadclass)
#       predicted
# actual0     1     2
# 0    28  1043     0
# 1     9 24112     1
# 2     1   106     5

confusion_matrix_fineclass <- 
  table(actual = factor(prediction_stage2$SCANFI_1km, levels = 1:8),
        predicted = factor(prediction_stage2$SCANFI_1km_prediction, levels = 1:8))

confusion_matrix_fineclass
sum(diag(confusion_matrix_fineclass)) / sum(confusion_matrix_fineclass)

     # predicted
#actual1    2    3    4    5    6    7    8
# 1    0    0    0    2    1   35    3    0
# 2    0    4    0    0   52   68   76    0
# 3    0    0    5    1   12   80   14    0
# 4    0    1    0   21   80  530  198    0
# 5    0    0    0    0 2186 1692 1109    0
# 6    0    0    0    5  667 9478 1473    0
# 7    0    1    1    3  927 3766 2814    0
# 8    0    0    0    0    0    0    0  448


