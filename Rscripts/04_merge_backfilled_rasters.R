# ---
# title: National Models 5.0 - stage raster stacks for repredicting bird densities
# author: Mannfred Boehm
# created: March 25, 2025
# ---

library(tidyverse)
library(terra)



#1. aggregate 1x1 rasters into 5x5----

# set root path
root <- "G:/Shared drives/BAM_NationalModels5"

# import original full stack
stack_bcr14_2020 <- terra::rast(file.path(root, "gis", "stacks", "can14_2020.tif"))

# import backfilled continuous biotic features 
backfilled_biotic_stack <- terra::rast("C:/Users/mannf/Proton Drive/mannfredboehm/My files/Drive/boreal_avian_modelling_project/ImpactAssessment/data/derived_data/backfilled_rasters/BCR14_backfilled_continuous.tif")

# define a 5x5 moving window (all weights equal) to compute mean cell values
window <- matrix(1, nrow = 5, ncol = 5)

# apply `focal()` on the extended raster.
# using na.rm = TRUE means the mean is computed on the available (non-NA) cells.
fivebyfive_stacks <- terra::focal(backfilled_biotic_stack, w=window, na.policy="omit", fun=mean, na.rm=TRUE)

# rename aggregated layers
names(fivebyfive_stacks) <- gsub("1km", "5x5", names(fivebyfive_stacks))
  
# remove aggregated layers that don't exist in the reference stack
fivebyfive_stacks <- fivebyfive_stacks[[intersect(names(fivebyfive_stacks), names(stack_bcr14_2020))]]




#2. merge backfilled biotic features into one raster stack----

# import backfilled categorical biotic features
scanfi_raster <- terra::rast(file.path(getwd(), "data", "derived_data", "backfilled_rasters", "BCR14_SCANFI_1km.tif"))
modis_raster <- terra::rast(file.path(getwd(), "data", "derived_data", "backfilled_rasters", "BCR14_MODISLCC_1km.tif"))
vlce_raster <- terra::rast(file.path(getwd(), "data", "derived_data", "backfilled_rasters", "BCR14_VLCE_1km.tif"))

names(scanfi_raster) <- "SCANFI_1km"
names(modis_raster) <- "MODISLCC_1km"
names(vlce_raster) <- "VLCE_1km"

# aggregate MODIS (it's the one landclass covariate available in 5x5)
modis_fivebyfive <- terra::focal(modis_raster, w=window, na.policy="omit", fun="modal", na.rm=TRUE)
names(modis_fivebyfive) <- "MODISLCC_5x5"

# merge backfilled biotic features into one raster stack
backfilled_biotic_stack <- c(fivebyfive_stacks, scanfi_raster, modis_raster, 
                             vlce_raster, backfilled_biotic_stack,
                             modis_fivebyfive)



#3. define abiotic layers to create a full backfilled raster stack----

# define abiotic variables so that they can be ported over to the backfilled stack
# import variable classes to help separate biotic from abiotic
nice_var_names <- 
  readr::read_csv(file.path(root, "data", "Extras", "sandbox_data", "impactassessment_sandbox", "covariates_label_insert.csv")) |> 
  dplyr::select(`Covariate label`, Category) |> 
  dplyr::rename(var = `Covariate label`, var_class = Category) |> 
  tidyr::drop_na() |> 
  unique()


# define abiotic variables that we want to port over from original raster stack
abiotic_vars <-
  nice_var_names |> 
  dplyr::filter(!(var_class %in% c("Biomass", "Greenup", "Landcover"))) |> 
  dplyr::filter(!(var %in% c("hli3cl_1km","GlobalHF_1km","GlobalHF_5x5","usroad_1km","usroad_5x5", "Peatland_1km", "Peatland_5x5"))) |> 
  tibble::add_row(var = c("year", "method"), var_class = c("year", "method")) |> 
  dplyr::pull(var) |> 
  unique()

# check for covariates in x that aren't in y
setdiff(abiotic_vars, names(stack_bcr14_2020))
setdiff(names(stack_bcr14_2020), c(abiotic_vars,names(backfilled_biotic_stack)))

# isolate original abiotic rasters
abiotic_rasters <- stack_bcr14_2020[[abiotic_vars]]

# merge original abiotic layers with backfilled biotic layers
stack_bcr14_2020_bf <- c(abiotic_rasters, backfilled_biotic_stack)
setdiff(names(stack_bcr14_2020), names(stack_bcr14_2020_bf))

# set disturbances and roads to zero in backfilled stack
disturbance_layers <- c("CanHF_1km", "CanHF_5x5", "canroad_1km", "canroad_5x5", "CCNL_1km")

for (layer in disturbance_layers) {
  terra::values(stack_bcr14_2020_bf[[layer]]) <- 0
}

terra::writeRaster(stack_bcr14_2020_bf, file.path(getwd(), "data/derived_data/backfilled_rasters/can14_2020_backfilled.tif"), overwrite=TRUE)
