# ---
# title: Impact Assessment: conform oil and gas layer to BAM
# author: Mannfred Boehm
# created: September 15, 2025
# ---

library(terra)
library(tidyverse)


root <- "G:/Shared drives/BAM_NationalModels5"
ia_dir <- file.path(root, "data", "Extras", "sandbox_data", "impactassessment_sandbox")

# import Hirsh-Pearson oil and gas raster
oilgas <- terra::rast(file.path(ia_dir, "hirshpearson_oil_gas.tif")) 

# now that we have a coarser raster, reprojection will be faster. then, crop and mask.
# first, create a template raster from the BAM boundary in EPSG:5072 (avoids GDAL errors)
bam_boundary <- terra::vect(file.path(root, "Regions", "BAM_BCR_NationalModel_UnBuffered.shp"))
bam_template <- terra::rast(file.path(root, "PredictionRasters", "Biomass", "SCANFI", "1km", "SCANFIBalsamFir_1km_2020.tif"))

oilgas_mask <-
  terra::project(x=oilgas, y=bam_template, method="bilinear") |> # use "near" because years are categorical, not continuous
  terra::crop(x=_, y=bam_boundary) |> 
  terra::mask(x=_, mask=bam_boundary)

# convert 0s and NaNs to NA (terra::buffer needs this)
vals <- values(oilgas_mask)
vals[vals == 0 | is.nan(vals)] <- NA
values(oilgas_mask) <- vals

# create 5km buffer around oil and gas infrastructure
oilgas_buffered <- terra::buffer(x=oilgas_mask, width = 5*res(bam_template)[1]) 

# convert convert TRUE/FALSE to 1/0
oilgas_buffered <- terra::ifel(oilgas_buffered, 1, 0)

# remove buffers that overlap with urbanized areas
# import Hirsh-Pearson population density layer
# Creston: 664, Saulte Ste. Marie: 324; Fort Mac 1304
popden <- terra::rast(file.path(ia_dir, "hirshpearson_population_density.tif"))

# generate urban mask (TRUE for >=400 persons/km^2)
# Stats Can threshold for urbanized is 400: https://www150.statcan.gc.ca/n1/en/catalogue/92-164-X
# Hirsh-Pearson transformed pop density by: threshold == 3.333 * log10(popden + 1) 
density_threshold <- 3.333 * log10(400 + 1) 
urban <- popden >= density_threshold
urban[urban == 0] <- NA  # keep only TRUE as 1, rest NA

# now that the popden raster is simplified, project
urban_reproj <- terra::project(x = urban, y = oilgas_buffered, method = "near")
urban_reproj<- terra::ifel(urban_reproj, 1, 0)

# remove oil and gas buffer pixels that overlap with urban
# set oil&gas = 0 wherever urban == 1 (handles NA in urban by leaving oil/gas as-is)
oilgas_no_urban <- oilgas_buffered * ((urban_reproj != 1) | is.na(urban_reproj))

terra::writeRaster(oilgas_no_urban, filename = file.path(ia_dir, "hirshpearson_oil_gas_masked.tif"))



