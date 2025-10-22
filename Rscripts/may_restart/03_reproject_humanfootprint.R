# ---
# title: Impact Assessment: conform Hirsh-Pearson cumulative impact layer to BAM
# author: Mannfred Boehm
# created: October 20, 2025
# ---

library(terra)
library(tidyverse)


root <- "G:/Shared drives/BAM_NationalModels5"
ia_dir <- file.path(root, "data", "Extras", "sandbox_data", "impactassessment_sandbox")

# first, create a template raster from the BAM boundary in EPSG:5072 (avoids GDAL errors)
bam_boundary <- terra::vect(file.path(root, "Regions", "BAM_BCR_NationalModel_UnBuffered.shp"))
bam_template <- terra::rast(file.path(root, "PredictionRasters", "Biomass", "SCANFI", "1km", "SCANFIBalsamFir_1km_2020.tif"))

# import Hirsh-Pearson oil and gas raster
CanHF_1km <- 
  terra::rast(file.path(root, "CovariateRasters", "Disturbance", "cum_threat2020.02.18.tif")) |> 
  terra::project(x = _, y = bam_template) |>
  terra::resample(x = _, y = bam_template) |> 
  terra::crop(x = _, y = bam_template) |> 
  terra::mask(x = _, mask = bam_template)

# convert NaNs to NA (terra::buffer needs this)
vals <- values(CanHF_1km)
vals[is.nan(vals)] <- NA
values(CanHF_1km) <- vals

names(CanHF_1km) <- "CanHF_1km"
terra::writeRaster(CanHF_1km, filename = file.path(ia_dir, "hirshpearson", "CanHF_1km_masked.tif"))



# ----------------------------------------------
# create low and high human footprint layers

CanHF_1km <- terra::rast(file.path(ia_dir, "hirshpearson", "CanHF_1km_masked.tif"))

# set low HF to <1 
# from Hirsh-Pearson: "we found that 82% of Canada’s land areas had a 
# HF < 1 and therefore were considered intact"
lowhf_mask <- 
  terra::lapp(CanHF_1km, \(x) ifelse(!is.finite(x), NA, ifelse(x < 1, 1, NA))) |> 
  as.factor()

terra::writeRaster(lowhf_mask, file.path(ia_dir, "hirshpearson", "CanHF_1km_lessthan1.tif"), overwrite = TRUE)

highhf_mask <- 
  terra::lapp(CanHF_1km, \(x) ifelse(!is.finite(x), NA, ifelse(x >= 1, 1, NA))) |> 
  as.factor()

terra::writeRaster(highhf_mask, file.path(ia_dir, "hirshpearson", "CanHF_1km_morethan1.tif"), overwrite = TRUE)

