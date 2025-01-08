# ---
# title: Impact Assessment - make predictions with zeroed human footprint
# author: Mannfred M.A. Boehm
# created: November 26, 2024
# ---

# NOTE: this is a modified version of Elly's V5 prediction script:
#https://github.com/borealbirds/LandbirdModelsV5/blob/main/analysis/07.Predict/07.Predict.R

# 5 disturbance variables in "NationalModels_V5_VariableList.xlsx"
# CanHF_1km
# CanHF_5x5
# CCNL_1km
# GlobalHF_1km
# GlobalHF_5x5



#1. Load packages----
print("* Loading packages on master *")
library(tidyverse)
library(gbm)
library(here)
library(parallel)
library(Matrix)
library(terra)


#2. Determine if testing and on local or cluster----
test <- TRUE
cc <- FALSE

#3. Set nodes for local vs cluster----
if(cc){ nodes <- 32}
if(!cc | test){ nodes <- 4}

#4. Create and register clusters----
print("* Creating clusters *")
cl <- makePSOCKcluster(nodes, type="PSOCK")

#5. Set root path----
print("* Setting root file path *")
if(cc){root <- "/scratch/ecknight/NationalModels"}
if(!cc){root <- "G:/Shared drives/BAM_NationalModels5"}

tmpcl <- clusterExport(cl, c("root"))

#6. Load packages on clusters----
print("* Loading packages on workers *")
tmpcl <- clusterEvalQ(cl, library(gbm))
tmpcl <- clusterEvalQ(cl, library(tidyverse))
tmpcl <- clusterEvalQ(cl, library(Matrix))
tmpcl <- clusterEvalQ(cl, library(terra))





#WRITE PREDICTION FUNCTION##########
brt_predict <- function(i){
  
  #1. Get model settings---
  bcr.i <- model_index$bcr[i]
  spp.i <- model_index$spp[i]
  boot.i <- model_index$boot[i]
  year.i <- model_index$year[i]
  
  #2. Load model----
  load.i <- try(load(file.path(root, "output", "bootstraps", paste0(spp.i, "_", bcr.i, "_", boot.i, ".R"))))
  if(inherits(load.i, "try-error")){ return(NULL) }
  
  #3. Load raster stack----
  stack.i <- try(rast(file.path(root, "gis", "stacks", paste0(bcr.i, "_", year.i, ".tif"))))
  if(inherits(stack.i, "try-error")){ return(NULL) }
  stack.i$meth.i <- stack.i$method
  
  # List of disturbance variables to zero out
  disturbance_vars <- c("CanHF_1km", "CanHF_5x5", "CCNL_1km", "GlobalHF_1km", "GlobalHF_5x5")
  
  # Modify relevant layers in the stack to simulate no human disturbance
  
  for(var.i in disturbance_vars) {
    if(var.i %in% names(stack.i)) {
      stack.i[[var.i]] <- stack.i[[var.i]] * 0  # Set all values in the disturbance layer to zero
    } else {
      warning(paste(var.i, "not found in raster stack"))
    }
  }
  
  
  #4. Predict----
  pred.i <- try(terra::predict(model=b.i, object=stack.i, type="response"))
  
  #7. Save----
  if(!(file.exists(file.path(here(), "data", "derived_data", "predictions", spp.i)))){
    dir.create(file.path(here(), "data", "derived_data", "predictions", spp.i))
  }
  if(inherits(pred.i, "SpatRaster")){
    writeRaster(pred.i, file=file.path(here(), "data", "derived_data", "predictions", spp.i, paste0(spp.i, "_", bcr.i, "_", boot.i, "_", year.i, ".tiff")), overwrite=TRUE)
  }
  
}



#8. Export objects to cluster----
print("* Loading function on workers *")

tmpcl <- clusterExport(cl, c("brt_predict"))




#RUN MODELS#########

#1. Set desired years----
years <- seq(1985, 2020, 5)



#2. Create index of bootstrapped models----
booted <- 
  data.frame(path = list.files(file.path(root, "output", "bootstraps"), pattern="*.R", full.names=TRUE),
             file = list.files(file.path(root, "output", "bootstraps"), pattern="*.R")) |> 
  tidyr::separate(file, into=c("spp", "bcr", "boot"), sep="_", remove=FALSE) |> 
  dplyr::mutate(boot = as.numeric(str_sub(boot, -100, -3)))


#3. Create to do list----
# expand_grid finds all combinations of factor variables
# DEFINE BCR AND SPECIES HERE!
# the following species are in my "priority species list" AND
# are also in Erin's OSM report: LEYE, RCKI, PAWA, ALFL, WEWP, OSFL, CHSP

focal_spp <- c("LEYE", "RCKI", "PAWA", "ALFL", "WEWP", "OSFL", "CHSP")

model_index <- 
  booted |> 
  dplyr::select(bcr, spp, boot) |> 
  tidyr::expand_grid(year=years) |> 
  dplyr::arrange(spp, boot, year, bcr) |> 
  dplyr::filter(bcr="bcr61", spp %in% focal_spp)



#4. Run prediction function in parallel----
# output is a list of prediction rasters
print("* Making predictions *")
predictions <- parLapply(cl, X=1:nrow(model_index), fun=brt_predict)



#5. Close clusters----
print("* Shutting down clusters *")
stopCluster(cl)

if(cc){ q() }



#6. Generate population estimations from rasters

# for every element in `predictions` assign to a species/bcr/year identifier, 
# and take the mean population 

for (p in 1:length(predictions)){
  
  
  
}
pop_est_list <- lapply(X=predictions, FUN=terra::app, x=mod, fun="mean", na.rm=TRUE)
  


# reduce list to dataframe, with every row as a spp x bcr x year tuple 
# and a column for population estimate
# graph trend and plot trend on map (see Anna's CFS presentation)

