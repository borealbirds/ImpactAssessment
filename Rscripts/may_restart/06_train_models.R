# ---
# title: Impact Assessment: train models on watersheds surrounding mines
# author: Mannfred Boehm
# created: May 7, 2025
# ---


library(tidyverse)
library(terra)

# set root path
root <- "G:/Shared drives/BAM_NationalModels5"

# for every year (1990-2020)
# 1. import year-matched mining layer
# 2. attach the mine to the appropriate watershed_w
# 3. define training area as all watersheds touching watershed_w
# 4. import V5 covariate stack, year-matched disturbance layer, soil layer
# 5. stack everything imported in step 4 and convert to a dataframe
# 6. mask stack to exclude pixels above HF threshold + mine buffer
# 7. train model on low HF areas

# define function for training models on low HF areas surrounding a mine
# this function trains models for continuous biotic features
backfill_mines_cont <- function(year){
  
  mines_raster_file <- 
    list.files(
    path = file.path(root, "gis", "other_landscape_covariates"),
    pattern = "mincan_mines\\.tif$",  
    full.names = TRUE)
  
  
}






