# ---
# title: National Models 5.0 - create models to predict biotic MODIS classes from abiotic landscape
# author: Mannfred Boehm
# created: March 11, 2025
# ---


#1. run lines 1-134 in "test_backfill_SCANFI_bcr14.R"---- 
# this will import covariate rasters, and generate a dataframe `df_bcr14_2020`from 
# the predictor covariate stack "can14_2020.tif"
# We also set the threshold for low/high human footprint by defining `q50`


#2. from "test_backfill_MODIS_bcr14.R" we know that wetlands aren't accuratley predicted---- 
# however, now we have biotic features as predictors...

# stage a unique spatial subset of the BCR by removing NAs for MODISLCC_1km
# we'll treat water, low HF barren, and low HF cropland MODIS classes as known features of the abiotic landscape
df_bcr14_2020_i <-  tidyr::drop_na(df_bcr14_2020, MODISLCC_1km)

# does every MODIS class appear in BCR14? 
dplyr::count(df_bcr14_2020_i, MODISLCC_1km)

# identify pixels with "high" and "low" human footprint 
CanHF_1km_present <- dplyr::filter(df_bcr14_2020_i, CanHF_1km > q50)
CanHF_1km_absent <- dplyr::filter(df_bcr14_2020_i, CanHF_1km <= q50) 

# drop biotic predictors that were thinned by VIF (see line 300 in "04.Stratify.R")
biotic_vars_thinned <- biotic_vars[which(biotic_vars %in% colnames(df_bcr14_2020))]

# subset low HF dataset to abiotic and biotic predictors (excluding SCANFI and VLCE) 
biotic_vars_thinned <- biotic_vars_thinned[!biotic_vars_thinned %in% c("SCANFI_1km", "MODISLCC_1km", "VLCE_1km")]

predictor_vars <- c(abiotic_vars, "lon", "lat", biotic_vars_thinned)
CanHF_1km_absent_abiotic <- CanHF_1km_absent[,c(predictor_vars, "MODISLCC_1km")] 



#5. two-stage modelling----
# stage 1: create broad MODISLCC_1km classes in low HF areas


# remove human disturbance classes (12, 13, 14, 16) because any remaining human footprints are 
# very minor (less than the median) and we aren't interested in predicting human disturbances,
# we just want to train the model to predict natural vegetation
# we will retain wetlands, even though we don't have high confidence that they can be predicted
training_data_stage1 <- 
  CanHF_1km_absent_abiotic |> 
  dplyr::filter(!MODISLCC_1km %in% c(12:14, 16:17)) |> 
  dplyr::mutate(MODISLCC_1km_broadclass = case_when(
    MODISLCC_1km %in% c(1:5) ~ 0,  # trees
    MODISLCC_1km %in% c(6:10) ~ 1,  # non-tree vegetation
    MODISLCC_1km == 11 ~ 2)) |>  # wetland 
  dplyr::mutate(MODISLCC_1km_broadclass = factor(MODISLCC_1km_broadclass, levels=c(0,1,2)))


# create DMatrix (a data structure for better speed when using xgboost)
# reminder: `predictor_vars` is `c(abiotic_vars, "lon", "lat", "MODISLCC_1km_crop", "MODISLCC_1km_urban", "MODISLCC_1km_barren", "MODISLCC_1km_water")`
dtrain_stage1 <- xgboost::xgb.DMatrix(data = as.matrix(training_data_stage1[,predictor_vars]), label = as.numeric(training_data_stage1$MODISLCC_1km_broadclass) - 1)
params_stage1 <- xgboost::xgb.params(objective = "multi:softmax", eval_metric = "mlogloss", max_depth = 3, eta = 0.1, num_class = length(levels(training_data_stage1$MODISLCC_1km_broadclass)))

# train stage 1 model to classify broad MODIS classes in low HF areas
# cv_stage1$early_stop$best_iteration 3480
cv_stage1 <-xgboost::xgb.cv(params = params_stage1, data = dtrain_stage1, nrounds=2000, nfold=5, early_stopping_rounds = 20, verbose=FALSE)
model_stage1 <- xgboost::xgb.train(params = params_stage1, data = dtrain_stage1, nrounds = cv_stage1$early_stop$best_iteration, verbose = 0)





#6. build two multi-class models to predict real MODIS classes from 
# broad classes 0 (tree class) and 1 (non-tree vegetation class)
# NOTE: we won't try to split the other broad classes (wetlands, HF, water, unclassified)
# we'll instead just port them over directly as "determined"

# train a model to split broad class 0 into 5 MODIS classes (tree types)
# filter for broadclass==0, then use `unique()` to check what MODISLCC_1km classes actually got modelled before using `factor()`
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

# cv_class0$early_stop$best_iteration 2161
cv_class0 <- xgboost::xgb.cv(params = params_class0, data = dtrain_class0, nrounds=2000, nfold=5, early_stopping_rounds = 20, verbose=FALSE)
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

cv_class1 <- xgboost::xgb.cv(params = params_class1, data = dtrain_class1, nrounds=500, nfold=5, early_stopping_rounds = 20, verbose=FALSE)
model_class1 <- xgb.train(params = params_class1, data = dtrain_class1, nrounds = cv_class1$early_stop$best_iteration, verbose = 0)




#4. use two-stage models to backfill in high HF areas----

# first, incorporate backfilled continuous variables at high HF locations
# (they will help to improve predictions of MODIS classes)
df_backfilled <-
  terra::rast("C:/Users/mannf/Proton Drive/mannfredboehm/My files/Drive/boreal_avian_modelling_project/ImpactAssessment/data/derived_data/backfilled_rasters/BCR14_backfilled_continuous.tif") |> 
  terra::as.data.frame(xy=TRUE) |>
  tibble::as_tibble() |> 
  dplyr::rename(lon=x, lat=y) 

# isolate high HF locations and their abiotic features 
CanHF_1km_present_abiotic <- CanHF_1km_present[,c(abiotic_vars, "MODISLCC_1km", "lon", "lat")]
nrow(CanHF_1km_present_abiotic) #128840

# remove water and wetland areas because we don't want to backfill there
backfill_data_stage1 <- 
  CanHF_1km_present_abiotic |> 
  dplyr::filter(!(MODISLCC_1km %in% c(11, 17))) |>  # 540+1630
  dplyr::left_join(df_backfilled, by=c("lon", "lat")) # add in (backfilled) biotic features to high HF areas


nrow(backfill_data_stage1) # should be 126670 which is 128840 - (540+1630) 


# reminder: predictor_vars <- c(abiotic_vars, "lon", "lat")
# including the `label` in the DMatrix is practically useless unless we use them for comparative purposes later on..
dbackfill <- xgboost::xgb.DMatrix(data = as.matrix(backfill_data_stage1[,predictor_vars]))

# stage 1: estimate broad MODIS classes in high HF areas
backfill_data_stage1$MODISLCC_1km_broadclass_predicted  <- predict(model_stage1, newdata = dbackfill)
count(backfill_data_stage1, MODISLCC_1km_broadclass_predicted)
# MODISLCC_1km_broadclass_predicted     n
# <dbl>                              <int>
#                               0  97128 predicted trees
#                               1  29196 predicted non-tree veg
#                               2   346   predicted wetland

# stage 2A:
# for non-treed predictions:
# filter high HF data for finer subclass predictions
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


# stage 2C :
# # if predicted as broad class wetland, assign as MODIS class "wetland"
backfill_wetland_stage2C <- 
  backfill_data_stage1 |> 
  dplyr::filter(MODISLCC_1km_broadclass_predicted  == "2") |>    
  dplyr::mutate(MODIS_backfill = "11")


# stage 2D: 
# retain MODIS classes for water and wetland from areas with high HF
retain_aqueous_stage2D <- 
  CanHF_1km_present |> 
  dplyr::filter(MODISLCC_1km %in% c(11,17)) |>   # only keep water pixels
  dplyr::mutate(MODIS_backfill = as.character(MODISLCC_1km))



# combine backfill predictions into single data frame
# number of predictions and carry-overs so should be equal to nrow(CanHF_1km_present)
prediction_stage2 <- dplyr::bind_rows(backfill_data_stage2A, backfill_data_stage2B, backfill_wetland_stage2C, retain_aqueous_stage2D)
prediction_stage2$MODIS_backfill <- as.numeric(prediction_stage2$MODIS_backfill)
nrow(CanHF_1km_present_abiotic) #128840
nrow(prediction_stage2) #128840 



#5. spatialize backfilled locations (high HF)----

# define an empty raster for placing backfilled features into
empty_raster <- terra::rast(extent=ext(stack_bcr14_2020), crs=crs(stack_bcr14_2020), resolution=1000)

# build backfilled raster
prediction_stage2_raster <-
  terra::vect(prediction_stage2, geom = c("lon", "lat"), crs = crs(stack_bcr14_2020)) |>
  terra::rasterize(x=_, y = empty_raster, field = "MODIS_backfill")

length(values(prediction_stage2_raster, na.rm=T)) #128840
length(values(stack_bcr14_2020$MODISLCC_1km, na.rm=T)) #274264


#  `cover()` fills any NA cells in the first raster with values from the second
prediction_stage2_raster <- terra::cover(prediction_stage2_raster, stack_bcr14_2020$MODISLCC_1km)
length(values(prediction_stage2_raster, na.rm=T)) #274264


# define classes
modis_cats <- data.frame(
    ID = sort(unique(df_bcr14_2020_i$MODISLCC_1km)), # 1  2  3  4  5  7  8  9 10 11 12 13 14 15 16 17
    category = c("eg needle", "eg broadleaf", "dec needle", "dec broadleaf", "mixed forest",
              "open shrub", "woody savanna", "savanna", "grassland", "wetland",
             "cropland", "urban", "cropland/natural", "snow/ice", "barren", "water"))

levels(prediction_stage2_raster) <- modis_cats
levels(stack_bcr14_2020$MODISLCC_1km) <- modis_cats




#6. visualize backfilling procedure----

my_colours <- c(
  "#006644",  # 1: evergreen needle
  "#006644",  # 2: evergreen broadleaf
  "#4DAC26",  # 3: deciduous needle
  "#4DAC26",  # 4: deciduous broadleaf
  "#74c476",  # 5: mixed forest
  "#bae4b3",  # 7: open shrubs
  "#bae4b3",  # 8: woody savanna
  "#bae4b3",  # 9: savanna
  "#bae4b3",  #10: grassland
  "#0072B2",  #11: wetland
  "#D55E00",  #12: cropland
  "#D55E00",  #13: urban
  "#D55E00",  #14: cropland/natural
  "snow",     #15: snow/ice
  "gray",     #16: barren
  "#0072B2")   #17: water


  
# plot pre-backfilled landscape for "MODISLCC_1km"
terra::plot(stack_bcr14_2020$MODISLCC_1km, col = my_colours)
lines(bcr14_boundary, col="black", lwd=1)

# plot post-backfilled landscape for "MODISLCC_1km"
terra::plot(prediction_stage2_raster, col = my_colours)
lines(bcr14_boundary, col="black", lwd=1)


# save it if it's good
terra::writeRaster(prediction_stage2_raster, file.path(getwd(), "data", "derived_data", "backfilled_rasters", "BCR14_MODISLCC_1km.tif"), overwrite=TRUE)
