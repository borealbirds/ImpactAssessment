# ---
# title: Impact Assessment: conform Hirsh-Pearson mining layer to BAM
# author: Mannfred Boehm
# created: October 20, 2025
# ---

library(terra)
library(tidyverse)


root <- "G:/Shared drives/BAM_NationalModels5"
ia_dir <- file.path(root, "data", "Extras", "sandbox_data", "impactassessment_sandbox")

# import Hirsh-Pearson oil and gas raster
mines <- terra::rast(file.path(ia_dir, "hirshpearson", "hirshpearson_mines.tif")) 

# first, create a template raster from the BAM boundary in EPSG:5072 (avoids GDAL errors)
bam_boundary <- terra::vect(file.path(root, "Regions", "BAM_BCR_NationalModel_UnBuffered.shp"))
bam_template <- terra::rast(file.path(root, "PredictionRasters", "Biomass", "SCANFI", "1km", "SCANFIBalsamFir_1km_2020.tif"))

mines_mask <-
  terra::project(x=mines, y=bam_template, method="bilinear") |> 
  terra::crop(x=_, y=bam_boundary) |> 
  terra::mask(x=_, mask=bam_boundary)


