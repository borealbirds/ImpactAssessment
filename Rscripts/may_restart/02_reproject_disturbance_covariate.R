# ---
# title: Impact Assessment: reproject time-since-disturbance to BAM crs
# author: Mannfred Boehm
# created: June 3, 2025
# ---

# reproject -> aggregate -> crop -> mask -> align (resample)

library(terra)
library(tidyverse)

# set root path
root <- "G:/Shared drives/BAM_NationalModels5"



# import time of disturbance layer (CAfire) and reproject (takes ~6 hours)
# note that Hermosilla et al (2016) is in NAD_1983_Lambert_Conformal_Conic
# downloaded from: https://opendata.nfis.org/mapserver/nfis-change_eng.html
CAfire <- terra::rast(file.path(root, "gis", "other_landscape_covariates", "CA_Forest_Fire_1985-2020.tif")) 

# group every 33 × 33 grid of 30m pixels into one 1000m x 1000m pixel (takes ~5 mins, but saves hours in reprojectiom time)
# `fact` tells terra how many raster cells to combine together when aggregating
fact <- round(1000 / res(CAfire)[1]) 
CAfire_agg <- terra::aggregate(CAfire, fact = fact, fun = "modal", na.rm = TRUE)

# now that we have a coarser raster, reprojection will be faster. then, crop and mask.
# first, create a template raster from the BAM boundary in EPSG:5072 (avoids GDAL errors)
bam_boundary <- terra::vect(file.path(root, "Regions", "BAM_BCR_NationalModel_Buffered.shp"))
bam_template <- terra::rast(ext(bam_boundary), resolution=1000, crs=crs(bam_boundary)) 

CAfire_mask <-
  terra::project(x=CAfire_agg, y=bam_template, method="near") |> # use "near" because years are categorical, not continuous
  terra::crop(x=_, y=bam_boundary) |> 
  terra::mask(x=_, mask=bam_boundary)


# convert from time *of* disturbance to time *since* disturbance
# add 1 in denominator to avoid dividing by zero when time of disturbance is 2020
# output: values closer to 1 are more recently disturbed
max_year <- terra::global(CAfire_mask, fun = "max", na.rm = TRUE)[1,1]
vals <- values(CAfire_mask)
values(CAfire_mask) <- ifelse(vals == 0, yes = 0, no = 1 / ((max_year - vals) + 1))


# save reprojected/cropped/masked time-since-disturbance layer
names(CAfire_mask) <- "CAfire"
terra::writeRaster(CAfire_mask, file.path(root, "gis", "other_landscape_covariates", "CA_Forest_Fire_1985-2020_masked.tif"), overwrite=TRUE)


