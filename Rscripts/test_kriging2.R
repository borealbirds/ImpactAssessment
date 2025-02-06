# ---
# title: National Models 5.0 - testing kriging assumptions for covariate layers
# author: Mannfred Boehm
# created: January 15, 2025
# ---

# A semivariogram quantifies the spatial autocorrelation of a variable by
# analyzing how its similarity (or dissimilarity) changes with distance.
# It measures the semi-variance, which is calculated as half the average 
# squared difference between data points separated by a given distance (h).

# import covariate values as data frame (code copied from "08.CalculateExtrapolation.R" in V5 pipeline)
# subset `id` column to BCR (or some region) of interest (see visit.i$id that loads when `b.i` is loaded)
# append `id` column with latlong info
# run gstat::variogram


#1. Load packages----
print("* Loading packages on master *")
library(brms)
library(gstat)
library(tidyverse)
library(terra)
library(parallel)



#2. Determine if testing and on local or cluster----
test <- TRUE
cc <- FALSE



#3. Set nodes for local vs cluster----
if(cc){ nodes <- 48}
if(!cc | test){ nodes <- 1}



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
tmpcl <- clusterEvalQ(cl, library(gstat))
tmpcl <- clusterEvalQ(cl, library(gbm))
tmpcl <- clusterEvalQ(cl, library(tidyverse))
tmpcl <- clusterEvalQ(cl, library(terra))



#7. Load stratified data (raster data in table form) ----
print("* Loading data on master *")

load(file.path(root, "data", "04_NM5.0_data_stratify.R"))
rm(bird, covlist, offsets, gridlist, birdlist)
gc()



#8. create dataframe of continuous covariates with lat-long and BCR info----
# partly copied from "08.CalculateExtrapolation.R" in V5 pipeline
# which had `cov_clean[names(cov_clean)!="hli3cl_1km"]` but 
# I don't see that `hli3cl_1km` is a variable in `cov`?

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
cov_clean <- readRDS(file="C:/Users/mannf/Proton Drive/mannfredboehm/My files/Drive/boreal_avian_modelling_project/ImpactAssessment/data/derived_data/rds_files/cov_clean.rds")


# filter cov_clean to BCR 14 for testing purposes (smallest spatial extent)
cov_clean_bcr14 <-
  cov_clean |> 
  dplyr::filter(if_any(bcr1:bcr5, function(x) x == "can14")) |> 
  dplyr::select(-c(can3:bcr5)) |> # drop bcr columns after bcr of interest is selected
  tidyr::drop_na(lat, lon, CanHF_1km, SCANFIprcD_1km) # remove NAs from variables of interest
  # dplyr::sample_n(1000) # subset for quicker testing 

# how does human footprint vary across this BCR?
hist(cov_clean_bcr14$CanHF_1km)
quantile(na.omit(cov_clean_bcr14$CanHF_1km))
# 0%       25%       50%       75%      100% 
# 0.000000  8.432904 11.590043 16.560811 52.190685 

hist(cov_clean_bcr14$CanHF_5x5)
quantile(na.omit(cov_clean_bcr14$CanHF_5x5))
# 0%       25%       50%       75%      100% 
# 0.000000  6.151992  9.552637 13.576068 41.631847 


#9. Set crs and rasterize HF locations----
#NAD83(NSRS2007)/Conus Albers projection (epsg:5072)
crs <- "+proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs +type=crs"

# identify pixels with human footprint
CanHF_1km_present <- dplyr::filter(cov_clean_bcr14, CanHF_1km > 0) 
plot(terra::vect(CanHF_1km_present, geom = c("lon", "lat"), crs = crs))


# identify pixels with low human footprint
q10 <- quantile(cov_clean_bcr14$CanHF_1km, probs = 0.10, na.rm = TRUE)
CanHF_1km_absent <- dplyr::filter(cov_clean_bcr14, CanHF_1km <= q10)
plot(terra::vect(CanHF_1km_absent, geom = c("lon", "lat"), crs = crs))



#10. Model relationship between biotic and abiotic covariates----
# this is part 1 of 2 of regression kriging:
# regression kriging combines a regression model 
# (to explain the deterministic variation based on covariates)
# with kriging (to capture the spatially autocorrelated residuals).

# import variable classes to help separate biotic from abiotic
nice_var_names <- 
  readr::read_csv("C:/Users/mannf/Proton Drive/mannfredboehm/My files/Drive/boreal_avian_modelling_project/ImpactAssessment/data/raw_data/nice_var_names_v5.csv") |> 
  dplyr::filter(var != "hli3cl_1km")

abiotic_vars <-
  nice_var_names |> 
  dplyr::filter(var_class %in% c("Annual Climate", "Climate Normals", "Topography")) |> 
  dplyr::pull(var) 

#biotic_vars <-
 #nice_var_names |> 
  #dplyr::filter(var_class %in% c("Wetland", "Landcover", "Greenup", "Biomass")) |> 
  #dplyr::pull(var) 

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
# remember to ensure that your response variable has no exact 0s or 1s, 
# as the beta distribution is defined on the open interval (0,1)
CanHF_1km_absent_abiotic[abiotic_vars] <- as_tibble(scale(CanHF_1km_absent_abiotic[abiotic_vars]))


# check that covariates in low footprint dataset have
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




# 11. use bayesian spatial regression to model----
# the relationship between the abiotic landscape and biotic drivers of bird occurrence

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

# posterior predictive check
brms::pp_check(prcD_model)

# perform 10-fold cross validation
kfold_result <- brms::kfold(prcD_model, K = 10)

# saveRDS(prcD_model, file="C:/Users/mannf/Proton Drive/mannfredboehm/My files/Drive/boreal_avian_modelling_project/ImpactAssessment/data/derived_data/rds_files/prcD_model.rds")
prcD_model<- readRDS(file="C:/Users/mannf/Proton Drive/mannfredboehm/My files/Drive/boreal_avian_modelling_project/ImpactAssessment/data/derived_data/rds_files/prcD_model.rds")

# use model to predict SCANFIprcD_1km for pixels with human footprint using the model's posterior predictions
# the resulting matrix has 3000 rows and 48153 columns because for every new data point (each of the 48153 pixels where backfilling occurred), 
# the model generated 3000 predictive samples (4 chains with 750 post-warmup draws). 
prcD_predictions <- brms::posterior_predict(object = prcD_model, newdata = CanHF_1km_present, summary = FALSE)
# saveRDS(prcD_predictions, file="C:/Users/mannf/Proton Drive/mannfredboehm/My files/Drive/boreal_avian_modelling_project/ImpactAssessment/data/derived_data/rds_files/prcD_predictions.rds") 




#12. check for spatial autocorrelation in residuals----

# isolate the data rows (complete cases) actually used in model construction
CanHF_1km_absent_abiotic_mf <- CanHF_1km_absent_abiotic[rownames(model.frame(prcD_model)),]

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











#12. use weighted averages of surrounding pixels (i.e. kriging)---- 
# to...



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
 



