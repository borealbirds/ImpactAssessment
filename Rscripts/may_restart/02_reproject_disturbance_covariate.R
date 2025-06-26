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

# convert zeros and NaNs to NAs
vals <- values(CAfire_mask)
vals[vals == 0 | is.nan(vals)] <- NA
values(CAfire_mask) <- vals

# define a function that estimates "years since fire" from "year of fire" from any given year
find_year_since_fire <- 
  function(current_year) {
    ysf <- current_year - CAfire_mask
    ysf[ysf < 0] <- NA  # fire hasn't occurred yet
    
    # normalize to c(0,1): recent fires=1, older fires=0
    ysf_norm <- 1 / (ysf + 1)  # +1 avoids division by zero when ysf == 0
    names(ysf_norm) <- paste0("CAfire_", current_year)
    return(ysf_norm)
}

# generate one raster per year with "years since fire"
years <- seq(1990, 2020, by = 5)
CAfire_rasters <- purrr::map(.x = years, .f = find_year_since_fire)
names(CAfire_rasters) <- paste0("CAfire_", years)


# save reprojected/cropped/masked time-since-disturbance layer
purrr::iwalk(CAfire_rasters, ~ {
  terra::writeRaster(.x,
                     filename = file.path(root, "gis", "other_landscape_covariates", paste0(.y, "_masked.tif")),
                     overwrite = TRUE)})



