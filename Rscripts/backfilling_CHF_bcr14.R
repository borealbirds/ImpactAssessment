# ---
# title: National Models 5.0 - testing kriging assumptions for covariate layers
# author: Mannfred Boehm
# created: January 15, 2025
# ---



# import stratified covariate values as data frame (code copied from "08.CalculateExtrapolation.R" in V5 pipeline)
# subset `id` column to BCR (or some region) of interest (see visit.i$id that loads when `b.i` is loaded)
# append `id` column with latlong info



#1. Load packages----
print("* Loading packages on master *")

library(tidyverse)
library(terra)
library(xgboost)
library(parallel)



#2. Determine if testing and on local or cluster----
test <- TRUE
cc <- FALSE



#3. Set nodes for local vs cluster----
if(cc){ nodes <- 18}
if(!cc | test){ nodes <- 4}



#4. Create and register clusters----
print("* Creating clusters *")
cl <- makePSOCKcluster(nodes, type="PSOCK")



#5. Set root path----
print("* Setting root file path *")
if(cc){root <- "/scratch/mannfred/NationalModels"}
if(!cc){root <- "G:/Shared drives/BAM_NationalModels5"}

tmpcl <- clusterExport(cl, c("root"))



#6. Load packages on clusters----
print("* Loading packages on workers *")
tmpcl <- clusterEvalQ(cl, library(gbm))
tmpcl <- clusterEvalQ(cl, library(tidyverse))
tmpcl <- clusterEvalQ(cl, library(terra))


#7. Load covariate stack(s)----
print("* Loading covariate stack *")
stack_bcr14_2020 <- terra::rast(file.path(root, "gis", "stacks", "can14_2020.tif"))
plot(stack_bcr14_2020$CCNL_1km)


# import time of disturbance layer (CAfire) and crop to BCR14
# note that Hermosilla et al (2016) is in NAD_1983_Lambert_Conformal_Conic
CAfire <- 
  terra::rast(file.path(root, "gis", "disturbancetime", "CA_Forest_Fire_1985-2020.tif")) |> 
  terra::project(x=_, y=stack_bcr14_2020, method="near") |>  # reproject to match covariate stack
  terra::crop(x=_, y=stack_bcr14_2020) # crop to BCR14

names(CAfire) <- "CAfire"

# convert from time *of* disturbance to time *since* disturbance
# add 1 in denominator to avoid dividing by zero when time of disturbance is 2020
# output: values closer to 1 are more recently disturbed
values(CAfire) <- ifelse(test = CAfire[] == 0, yes = 0, no = 1 / ((max(values(CAfire)) - CAfire[]) + 1))

# add time since disturbance layer to covariate stack
stack_bcr14_2020 <- c(stack_bcr14_2020, CAfire)



# inspect distribution of disturbance times (after removing cells with "no change" aka "zero")
hist(values(CAfire)[which(values(CAfire) > 0)], main="dist. of disturbances")
plot(CAfire)

# overlay BCR14 boundary to sanity check
bcr14_boundary <- 
  terra::vect(file.path(root, "Regions", "BAM_BCR_NationalModel_Unbuffered.shp")) |> 
  terra::project(x=_, y=stack_bcr14_2020) |> 
  terra::crop(x=_, y=stack_bcr14_2020)

lines(bcr14_boundary, col="black", lwd=1)




#8. stage datasets needed for modelling biotic landscape----

# import variable classes to help separate biotic from abiotic
nice_var_names <- 
  readr::read_csv(file.path(root, "covariates_label.csv")) |> 
  dplyr::select(Label, Category) |> 
  dplyr::rename(var = Label, var_class = Category) |> 
  tidyr::drop_na() |> 
  unique()

# remove Peatland as they are considered biotic,
# remove WetOccur_5x5 since it's redundant with WetOccur_1km
# for hli3cl see line 67 in "08.CalculateExtrapolation.R"
abiotic_vars <-
  nice_var_names |> 
  dplyr::filter(var_class %in% c("Annual Climate", "Climate Normals", "Topography", "Wetland")) |> 
  tibble::add_row(var = "CAfire", var_class ="Time Since Disturbance") |> 
  dplyr::filter(!(var %in% c("hli3cl_1km", "Peatland_1km", "Peatland_5x5", "WetOccur_5x5"))) |> 
  dplyr::pull(var) |> 
  unique()

# remove 5x5 since I can re-create these variables by
# scaling up from 1km after prediciton
biotic_vars <-
  nice_var_names |> 
  dplyr::filter(var_class %in% c("Landcover", "Greenup", "Biomass", "LCC_MODIS", "Wetland")) |> 
  dplyr::filter(!(var %in% c("WetOccur_1km", "WetOccur_5x5", "WetRecur_1km", "WetSeason_1km"))) |> # keep peatland but discard other water variables
  dplyr::filter(!grepl("5x5", var)) |> 
  dplyr::pull(var) |> 
  unique()

# convert covariate stack to a dataframe (values will be used in regression)
# convert landcover classes from continuous to categorical
covs_df <-
  stack_bcr14_2020 |> 
  terra::as.data.frame(xy=TRUE) |>
  tibble::as_tibble() |> 
  dplyr::rename(lon=x, lat=y) |> 
  dplyr::mutate(SCANFI_1km = as.factor(SCANFI_1km),
                MODISLCC_1km = as.factor(MODISLCC_1km),
                MODISLCC_5x5 = as.factor(MODISLCC_5x5),
                VLCE_1km = as.factor(VLCE_1km))
  

# define threshold for "low" human footprint
q75 <- quantile(covs_df$CanHF_1km, probs = 0.75, na.rm = TRUE)

# how does human footprint vary across this BCR?
hist(covs_df$CanHF_1km, main="CanHF_1km in BCR14")
abline(v=q75, col="darkred", lwd=2)
abline(v=mean(na.omit(covs_df$CanHF_1km)), col="skyblue", lwd=2, lty="dashed")
       
quantile(na.omit(covs_df$CanHF_1km))
#0%        25%        50%        75%       100% 
#0.0000000  0.9322898  5.2230940  9.6365070 54.5055199 





#9. test whether SCANFI_1km land cover classes are accurately predicted---- 
# especially since time since disturbance is now an abiotic predictor

# drop biotic predictors that were thinned by VIF (see line 300 in "04.Stratify.R")
biotic_vars_thinned <- biotic_vars[which(biotic_vars %in% colnames(covs_df))]

# stage a unique spatial subset of the BCR by removing NAs for SCANFI_1km
i <- 17 #for SCANFI_1km
covs_df_i <- 
  tidyr::drop_na(covs_df, biotic_vars_thinned[i]) |>
  dplyr::mutate(label = as.numeric(SCANFI_1km) - 1)  # xgboost labels start at 0

# identify pixels with "high" and "low" human footprint
# for SCANFI_1km
CanHF_1km_present <- dplyr::filter(covs_df_i, CanHF_1km > q75) 
CanHF_1km_absent <- dplyr::filter(covs_df_i, CanHF_1km <= q75)

# subset low HF dataset to only those with abiotic predictors and the biotic response
#exclude "SCANFI_1km" since it will obviously correlate with "label"
CanHF_1km_absent_abiotic <- CanHF_1km_absent[,c(abiotic_vars, "lon", "lat", "label")]

# split data into training and hold-out
set.seed(123)
n <- nrow(CanHF_1km_absent_abiotic)
training_index <- sample(seq_len(n), size = round(0.8 * n))
predictor_vars <- c(abiotic_vars, "lon", "lat")
training_data <- CanHF_1km_absent_abiotic[training_index, c(predictor_vars, "label")] 
holdout_data <- CanHF_1km_absent_abiotic[-training_index, c(predictor_vars, "label")]

# create DMatrix (a data structure for better speed when using xgboost)
dtrain <- xgboost::xgb.DMatrix(data = as.matrix(training_data[,predictor_vars]), label = training_data$label)
dholdout <- xgboost::xgb.DMatrix(data = as.matrix(holdout_data[,predictor_vars]), label = holdout_data$label)

# set xgboost parameters
params <- 
  list(objective = "multi:softmax",  # direct multiclass prediction
       eval_metric = "mlogloss",
       num_class = length(levels(as.factor(training_data$label))),
       max_depth = 3,  
       eta = 0.1) # learning rate


# run cross-validation to tune the number of rounds
cv <- 
  xgboost::xgb.cv(params = params, 
                  data = dtrain, 
                  nrounds = 1000, 
                  nfold = 5, 
                  early_stopping_rounds = 10, 
                  verbose = 1)

#cv$early_stop
#$best_iteration
#[1] 927

#$best_score
#test-mlogloss 
#1.177003 

#$stopped_by_max_rounds
#[1] TRUE


# fit a model
final_model <- xgb.train(
  params = params,
  data = dtrain,
  nrounds = cv$early_stop$best_iteration,
  verbose = 1)

pred_holdout <- predict(final_model, dholdout)
conf_matrix <- table(predicted = pred_holdout, actual = holdout_data$label)
accuracy <- sum(diag(conf_matrix)) / sum(conf_matrix)
# [1] 0.5086198









# define empty raster for placing backfilled features into
empty_raster <- 
  terra::rast(extent=ext(stack_bcr14_2020), crs=crs(stack_bcr14_2020), resolution=1000)


#9. define a function for predicting "restored" biotic landscape features----
backfill_landscape <- function(i){
  
  # drop biotic predictors that were thinned by VIF (see line 300 in "04.Stratify.R")
  biotic_vars_thinned <- biotic_vars[which(biotic_vars %in% colnames(covs_df))]
  
  # stage a unique spatial subset of the BCR by removing NAs for biotic feature[i]
  covs_df_i <- tidyr::drop_na(covs_df, biotic_vars_thinned[i]) 
    
  # identify pixels with "high" and "low" human footprint
  # for biotic feature[i]
  CanHF_1km_present <- dplyr::filter(covs_df_i, CanHF_1km > q75) 
  CanHF_1km_absent <- dplyr::filter(covs_df_i, CanHF_1km <= q75)
  
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



# ------------------------------------------------------------
# for visualization purposes only: plot pre- and 
# post-backfilled landscape for "SCANFIBalsamFir_1km"
predictions_raster_i_smoothed <- terra::focal(predictions_raster_i, w = matrix(1,3,3), fun = mean, na.rm = TRUE) #3x3 matrix filled with 1s (each cell is given equal weight)

terra::plot(predictions_raster_i)
terra::plot(stack_bcr14_2020$SCANFIBalsamFir_1km)

lines(bcr14_boundary, col="black", lwd=1)

# for visualization purposes only: plot low vs high human footprint areas
# CanHF_present and CanHF_absent are made inside of the function, 
# so these maps are for "SCANFIBalsamFir_1km"
my_colours <- colorRampPalette(c("#0072B2", "#009E73", "#F0E442", "#D55E00"))(100)

terra::vect(CanHF_1km_absent, geom = c("lon", "lat"), crs = crs(empty_raster)) |>
terra::rasterize(y = empty_raster, field = "CanHF_1km") |> 
terra::plot(col= colorRampPalette(c("#56B4E9", "#0072B2"))(1000))

terra::vect(CanHF_1km_present, geom = c("lon", "lat"), crs = crs(empty_raster)) |>
terra::rasterize(y = empty_raster, field = "CanHF_1km") |> 
terra::plot(col=my_colours)

lines(bcr14_boundary, col="black", lwd=1)
# ------------------------------------------------------------




covs_df$SCANFI_1km |> unique()








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





#13. import bird models (`b.i`) for CAWA at year 2020----

# access gbm objects and append spp/bcr/boot info
root <- "G:/Shared drives/BAM_NationalModels5/output/bootstraps"
gbm_objs <- list.files(file.path(root), pattern = "*\\.R", full.names = TRUE, recursive = TRUE)


# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# HOW TO FILTER FOR YEAR????
# extract the species (FLBC), BCR, and bootstrap replicate from `gbm_objs`
sample_id <- 
  gbm_objs |> 
  basename() |> 
  stringr::str_split_fixed(pattern="_", n=3) |> 
  gsub("\\.R", "", x = _) |>
  tibble::as_tibble() |> 
  magrittr::set_colnames(c("spp", "bcr", "boot")) |> 
  dplyr::filter(spp == "CAWA")

gbm_data <- tibble(file_path = gbm_objs, spp = sample_id$spp, bcr = sample_id$bcr, boot = sample_id$boot) 



# write function that `terra::predict`s bird densities using V5 `gbm` models for CAWA 
repredict_with_backfill <- function(gbm_object, backfilled_stack) {

    
  # attempt to load the GBM object
  # "file_path" is a column name in `gbm_data`
  load(gbm_object[["file_path"]])
    
  # validate gbm object
  if (!exists("b.i") || is.null(b.i$n.trees) || b.i$n.trees <= 0) {
    warning(paste("Invalid GBM model in file:", gbm_object["file_path"]))
    return(NULL)
  } else {
    message("Successfully loaded: ", gbm_object["file_path"])
  }
  
  # predict density for a given spp x bcr x year
  s
  object_i <- backfilled_stack
  prediction_i <- terra::predict(object=SPATRASTER, model=b.i, )
  return(prediction_i)
}





#XYZ. check that covariates in low footprint dataset have----
# comparable range to original dataset
# this is important because we are modelling biotic ~ abiotic
# assuming that the low-HF dataset is representative of the entire BCR
subset_var_range_vs_og <- list()
for (i in seq_len(ncol(cov_clean_bcr14))){
  
  # omitting biotic covariates because we assume they will be affected by HF
  # and therefore datasets with low-HF are by definition going to have narrower
  # ranges of biotic variables
  if (colnames(CanHF_1km_absent)[i] %in% abiotic_vars) { 
    
    # identify a covariate of interest
    current_var <- colnames(CanHF_1km_absent)[i]
    
    # get range of current covariate for 
    original_quant <- quantile(na.omit(cov_clean_bcr14[[current_var]]))
    absent_quant <- quantile(na.omit(CanHF_1km_absent[[current_var]]))
    
    subset_var_range_vs_og[[i]] <- cbind(original_quant, absent_quant)
    names(subset_var_range_vs_og)[i] <- current_var
    
  }
  
  else(NULL)
  
}


#XYZ. check for spatial autocorrelation in residuals----


# extract residuals
CanHF_1km_absent_abiotic_mf$residuals <- stats::residuals(prcD_model)[,"Estimate"]
range(CanHF_1km_absent_abiotic_mf$residuals)

# check for autocorrelation in residuals
# define extent for rasterizing (xmin, xmax, ymin, ymax)
extent <- terra::ext(c(min(CanHF_1km_absent_abiotic_mf$lon),
                       max(CanHF_1km_absent_abiotic_mf$lon),
                       min(CanHF_1km_absent_abiotic_mf$lat),
                       max(CanHF_1km_absent_abiotic_mf$lat)))

# rasterize the set of points with low HF with associated residual variation
empty_raster <- terra::rast(extent=extent, resolution=0.01)
CanHF_1km_absent_abiotic_vect <- terra::vect(CanHF_1km_absent_abiotic_mf, geom=c("lon", "lat"), crs=crs)
CanHF_1km_absent_abiotic_rast <- terra::rasterize(x=CanHF_1km_absent_abiotic_vect, y=empty_raster, field="residuals") 

# test for spatial autocorrelation in residuals
autocor_test <- terra::autocor(CanHF_1km_absent_abiotic_rast, global=FALSE, method="moran")
moran_values <- terra::values(autocor_test)
summary(moran_values)
hist(moran_values, main = "Histogram of local Moran's I values", xlab = "Local Moran's I")











# OLD CODE



# `visit` is loaded with "04_NM5.0_data_stratify.R"
# filter for year(s) of interest to use as a `year` index for point locations in `cov_clean`
# saveRDS(visit, file="C:/Users/mannf/Proton Drive/mannfredboehm/My files/Drive/boreal_avian_modelling_project/ImpactAssessment/data/derived_data/rds_files/visit.rds")
visit <- readRDS(file="C:/Users/mannf/Proton Drive/mannfredboehm/My files/Drive/boreal_avian_modelling_project/ImpactAssessment/data/derived_data/rds_files/visit.rds")

visit_clean <-
  visit |> 
  dplyr::select(id, year) |> 
  dplyr::filter(year %in% 2020:2022)

# a single location ID may show up in as many as 5 BCRs
# this is because each BCR is buffered by 100km, creating BCR overlap
bcrlist |>
  dplyr::select(-id) |>
  (\(data) rowSums(data))() |>
  max()

cov_clean <- 
  cov |> 
  dplyr::select(where(is.numeric)) |>
  dplyr::left_join(dplyr::select(visit, c(lat, lon, id)), by = "id") |> # need lat-long for testing kriging assumptions
  dplyr::left_join(bcrlist, by = "id") |> # append bcr info
  dplyr::rowwise() |> # operate the following on each row
  dplyr::mutate(
    bcr1 = names(pick(can3:usa41423))[which(c_across(can3:usa41423))[1]], 
    bcr2 = names(pick(can3:usa41423))[which(c_across(can3:usa41423))[2]],
    bcr3 = names(pick(can3:usa41423))[which(c_across(can3:usa41423))[3]], 
    bcr4 = names(pick(can3:usa41423))[which(c_across(can3:usa41423))[4]], 
    bcr5 = names(pick(can3:usa41423))[which(c_across(can3:usa41423))[5]]) |>
  dplyr::ungroup()

# saveRDS(cov_clean, file="C:/Users/mannf/Proton Drive/mannfredboehm/My files/Drive/boreal_avian_modelling_project/ImpactAssessment/data/derived_data/rds_files/cov_clean.rds")
cov_clean <- 
  readRDS(file="C:/Users/mannf/Proton Drive/mannfredboehm/My files/Drive/boreal_avian_modelling_project/ImpactAssessment/data/derived_data/rds_files/cov_clean.rds") |> 
  as_tibble()

# focus dataset on BCR 14 and year >= 2020
cov_clean_bcr14 <-
  cov_clean |> 
  dplyr::filter(if_any(bcr1:bcr5, function(x) x == "can14")) |> 
  dplyr::select(-c(can3:bcr5)) |> # drop bcr columns after bcr of interest is selected
  tidyr::drop_na(lat, lon, CanHF_1km) |> 
  dplyr::filter(id %in% visit_clean$id)

# use bayesian spatial regression to model----
# the relationship between the abiotic landscape and biotic drivers of bird occurrence

# create dataset with only one biotic variable for prediction
# using `CanHF_1km_absent` so that model is informed by
# areas with low human footprint
# scale prcD so that it can be described by a beta dist.
CanHF_1km_absent_abiotic <-
  CanHF_1km_absent |>
  dplyr::select(any_of(c(abiotic_vars, "lat", "lon"))) |> 
  dplyr::mutate(SCANFIprcD_1km = CanHF_1km_absent$SCANFIprcD_1km/100) 

# check that response variable is between 0 and 1 (for using beta dist. downstream)
range(CanHF_1km_absent_abiotic$SCANFIprcD_1km)
# [1] 0.0002159509 0.9920376587

# scale and center abiotic predictors to make it easier to make an informed prior
# e.g. for the intercept, we can now guess at the average outcome, rather than the outcome at zero 
# remember to ensure that the response variable has no exact 0s or 1s, 
# as the beta distribution is defined on the open interval (0,1)
scaled_abiotic_vars <- scale(CanHF_1km_absent_abiotic[,abiotic_vars])
center_values <- attr(scaled_abiotic_vars, "scaled:center")
scale_values <- attr(scaled_abiotic_vars, "scaled:scale")

CanHF_1km_absent_abiotic[,abiotic_vars] <- as_tibble(scaled_abiotic_vars)




# define model
prcD_formula <- SCANFIprcD_1km ~ . - lat - lon + s(lat, lon)

# define priors
priors <- c(
  prior(normal(0, 1), class = "b"),        # Priors for coefficients on the logit scale for μ
  prior(normal(0, 1), class = "Intercept"),  # Prior for the intercept on the logit scale
  prior(exponential(1), class = "phi")       # Prior for the precision parameter
)

# fit a Bayesian spatial regression model
prcD_model <- 
  brms::brm(
    formula = prcD_formula,
    prior = NULL,
    data = CanHF_1km_absent_abiotic,
    family = brms::Beta(link = "logit"),  # ensures that the model’s predictions for the mean are constrained to (0, 1)
    cores = 4,            
    chains = 4,           
    iter = 1500,         
    control = list(adapt_delta = 0.80))  # for exploring parameter space. larger # for smaller step sizes.

# saveRDS(prcD_model, file="C:/Users/mannf/Proton Drive/mannfredboehm/My files/Drive/boreal_avian_modelling_project/ImpactAssessment/data/derived_data/rds_files/prcD_model.rds")
prcD_model<- readRDS(file="C:/Users/mannf/Proton Drive/mannfredboehm/My files/Drive/boreal_avian_modelling_project/ImpactAssessment/data/derived_data/rds_files/prcD_model.rds")

# `posterior_linpred` draws from the posterior distribution of parameter values for a single fitted model
# each draw contains a set of parameter values: e.g. regression coefficients (one per predictor), the intercept, a precision parameter ϕ in a beta regression, etc.
# for each draw, the set of parameter values (e.g. regression coefficients) is multiplied  by the corresponding predictor values from the new data.
# as a result, the average predicted values of prcD range from 0.033 to 0.956
# however, there are 505 NAs which means something went wrong with the prediction step (likely extrapolation) 

pred_mu <- brms::posterior_linpred(prcD_model, newdata = CanHF_1km_present, transform = TRUE)
summary(apply(pred_mu, 2, mean)) 

# posterior predictive check
brms::pp_check(prcD_model)

# k-fold cross validation
kfold_result <- brms::kfold(prcD_model, K = 10)



# use model to predict SCANFIprcD_1km for pixels with human footprint >q10
# using the model's posterior predictions
# the resulting matrix has 3000 rows and 43974 columns because for every new data point (each of the 48153 pixels where backfilling occurred), 
# the model generated 3000 predictive samples (4 chains with 750 post-warmup draws). 
prcD_predictions <- brms::posterior_predict(object = prcD_model, newdata = CanHF_1km_present[ ,c(abiotic_vars, "lon", "lat")], summary = FALSE)
# saveRDS(prcD_predictions, file="C:/Users/mannf/Proton Drive/mannfredboehm/My files/Drive/boreal_avian_modelling_project/ImpactAssessment/data/derived_data/rds_files/prcD_predictions.rds") 








#XYZ. use weighted averages of surrounding pixels (i.e. kriging)---- 
# to...

# A semivariogram quantifies the spatial autocorrelation of a variable by
# analyzing how its similarity (or dissimilarity) changes with distance.
# It measures the semi-variance, which is calculated as half the average 
# squared difference between data points separated by a given distance (h).

# generate variogram from pixels on the low-disturbance landscape
# (i.e. omit high human footprint pixels for creating variogram)
prcD_residuals_variogram <- gstat::variogram(object=residuals ~ 1, data=CanHF_1km_absent_abiotic_mf, locations= ~lon + lat)
plot(prcD_residuals_variogram)

# fit a model to the variogram
prcD_residuals_variogram_model <- gstat::fit.variogram(prcD_residuals_variogram, vgm(model = "Sph"))

          

#13. use variogram model to predict prcD at human footprint locations

# spatialize tables for kriging
# format coordinates for `gstat::krige`
CanHF_1km_absent_sf <- sf::st_as_sf(CanHF_1km_absent_abiotic_mf, coords = c("lon", "lat"))
CanHF_1km_present_sf <- sf::st_as_sf(CanHF_1km_present, coords = c("lon", "lat"))

# apply kriging
prcD_kriging_result <- 
  gstat::krige(
    formula = residuals ~ 1,
    locations = CanHF_1km_absent_sf,
    newdata = CanHF_1km_present_sf,
    model = prcD_residuals_variogram_model)
 



