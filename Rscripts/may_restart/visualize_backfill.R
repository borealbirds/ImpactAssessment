# ---
# title: Impact Assessment: visualize and validate backfilled landscapes
# author: Mannfred Boehm
# created: December 1, 2025
# ---

library(terra)
library(tidyterra)
library(tidyverse)

root <- "G:/Shared drives/BAM_NationalModels5"
ia_dir <- file.path(root, "data", "Extras", "sandbox_data", "impactassessment_sandbox")

year <- 2020
stack_y <- rast(file.path(ia_dir, sprintf("covariates_mosaiced_%d.tif", year)))

# import subbasin boundaries and project to current stack
all_subbasins_subset <- terra::vect(file.path(ia_dir, "hydrobasins_masked_merged_subset.gpkg"))
all_subbasins_subset <- terra::project(x=all_subbasins_subset, y=stack_y)
sub <- all_subbasins_subset[57]

# pre backfilled landscape
pre57 <- mask(crop(stack_y, sub), sub)

# post backfilled landscape
post57 <- rast(file.path(ia_dir, "bart_models/2020/subbasin_57/subbasin_57_backfill.tif"))

# human footprint that was backfilled
highhf_mask <- terra::rast(file.path(ia_dir, "hirshpearson", "CanHF_1km_morethan1.tif"))
highhf_mask <- 
  terra::project(x=highhf_mask, y=stack_y, method = "near") |> 
  terra::crop(x=_, y=sub) |> 
  terra::mask(x=_, mask=sub)

lowhf_mask <- terra::rast(file.path(ia_dir, "hirshpearson", "CanHF_1km_lessthan1.tif"))
lowhf_mask <- 
  terra::project(x=lowhf_mask, y=stack_y, method = "near") |> 
  terra::crop(x=_, y=sub) |> 
  terra::mask(x=_, mask=sub)

# SCANFI
my_colours <- c(
  "#bae4b3",  # ID 1: bryoid
  "#bae4b3",  # ID 2: herb
  "#E69F00",  # ID 3: rock
  "#bae4b3",  # ID 4: shrub
  "#7FFF00",  # ID 5: broadleaf (treed, darker shade)
  "#006644",  # ID 6: conifer (treed, darker shade)
  "#9ACD32",  # ID 7: mixed (treed, darker shade)
  "#0072B2"   # ID 8: water
)

plot(pre57$SCANFI_1km, col=my_colours)
plot(all_subbasins_subset[57], add=TRUE )

plot(highhf_mask)
plot(all_subbasins_subset[57], add=TRUE )

plot(post57$SCANFI_1km, col=my_colours)

# plot human footprint
plot(lowhf_mask, col="#56B4E9")
plot(highhf_mask, col="#E69F00", add=TRUE)
plot(sub, add=TRUE)

# plot satellite
sat <- maptiles::get_tiles(
  x        = sub, # your AOI (SpatVector)
  provider = "Esri.WorldImagery",      # high-res imagery
  zoom     = 10,                       # bump to 11–12 if you need more detail
  crop     = TRUE
)

sat <- project(sat, stack_y, method = "near") |> 
  crop(x=_, y=sub) |> 
  mask(x=_, mask=sub) 

plot(sat)
# --------------------------------------------------------------------
# build confusion matrix (cm)

cm57 <- readRDS(file.path(ia_dir, "bart_models/2020/subbasin_57/subbasin_057_confusion.rds"))
all_levels <- sort(unique(c(cm57[[1]]$actual, cm57[[1]]$predicted)))

testscanfi <- table(factor(cm57[[1]]$actual,    levels = all_levels),
                   factor(cm57[[1]]$predicted, levels = all_levels))
testscanfi

metrics57 <- readRDS(file.path(ia_dir, "bart_models/2020/subbasin_57/subbasin_057_metrics.rds"))








# MODIS
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

# VLCE
my_colours <- c(
  "#CC79A7",  # 0: no change
  "#0072B2",  #20: water
  "#D55E00",  #32: rock
  "#D55E00",  #33: barren
  "#bae4b3",  #40: bryoid
  "#bae4b3",  #50: shrub
  "#56B4E9",  #80: wetland
  "#009E73",  #81: wetland-tree
  "#bae4b3",  #100: herb
  "#006644",  # 3: conifer
  "#4DAC26",  # 4: broadleaf
  "#74c476")  # 5: mixed forest


# TROUBLESHOOTING

diag_r <- crop(stack_y[[1]]
vals_idx <- rep(NA_real_, terra::ncell(stack_y[[1]]))
