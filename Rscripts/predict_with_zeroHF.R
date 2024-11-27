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
  bcr.i <- loop$bcr[i]
  spp.i <- loop$spp[i]
  boot.i <- loop$boot[i]
  year.i <- loop$year[i]
  
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
  if(!(file.exists(file.path(root, "output", "predictions", spp.i)))){
    dir.create(file.path(root, "output", "predictions", spp.i))
  }
  if(inherits(pred.i, "SpatRaster")){
    writeRaster(pred.i, file=file.path(root, "output", "predictions", spp.i, paste0(spp.i, "_", bcr.i, "_", boot.i, "_", year.i, ".tiff")), overwrite=TRUE)
  }
  
}
