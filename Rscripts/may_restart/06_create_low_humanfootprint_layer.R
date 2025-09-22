# ---
# title: Impact Assessment: create low HF raster 


# author: Mannfred Boehm
# created: September 16, 2025
# ---

library(terra)

root <- "G:/Shared drives/BAM_NationalModels5"
ia_dir <- file.path(root, "data", "Extras", "sandbox_data", "impactassessment_sandbox")


# restrict analysis to BCRs that we have covariate data for
bam_template <- terra::rast(file.path(root, "PredictionRasters", "Biomass", "SCANFI", "1km", "SCANFIBalsamFir_1km_2020.tif"))


# import Hirsh-Pearson HF layer
CanHF_1km <- 
  terra::rast(file.path(root, "CovariateRasters", "Disturbance", "cum_threat2020.02.18.tif")) |> 
  terra::project(x = _, y = bam_template) |>
  terra::resample(x = _, y = bam_template) |> 
  terra::crop(x = _, y = bam_template) |> 
  terra::mask(x = _, mask = bam_template)

names(CanHF_1km) <- "CanHF_1km"


# set low HF to <1 
# from Hirsh-Pearson: "we found that 82% of Canada’s land areas had a 
# HF < 1 and therefore were considered intact"
lowhf_mask <- 
  terra::lapp(CanHF_1km, \(x) ifelse(!is.finite(x), NA, ifelse(x < 1, 1, NA))) |> 
  as.factor()


terra::writeRaster(lowhf_mask, file.path(ia_dir, "CanHF_1km_lessthan1.tif"), overwrite = TRUE)
