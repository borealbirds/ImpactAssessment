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



# ----------------------------------------------
# create "good sectors" high HF layer:
# CanHF >= 1 AND at least one target sector present (score > 0)
# excluded sectors: forestry_harvest, night_lights, population_density, nav_water
good_sectors <- c("built", "crop", "dam_and_associated_reservoir", "mines",
                  "oil_gas", "pasture", "rail", "roads")

good_sector_stack <- terra::rast(lapply(good_sectors, function(sec) {
  r <- terra::rast(file.path(ia_dir, "hirshpearson", paste0(sec, ".tif")))
  terra::project(r, CanHF_1km, method = "near")
}))

# sum sector scores; > 0 means at least one target sector is present
sector_sum <- terra::app(good_sector_stack, sum, na.rm = TRUE)

goodsectors_mask <-
  terra::ifel((CanHF_1km >= 1) & (sector_sum > 0), 1, NA) |>
  as.factor()

terra::writeRaster(goodsectors_mask,
                   file.path(ia_dir, "hirshpearson", "CanHF_1km_morethan1_goodsectors.tif"),
                   overwrite = TRUE)



# ----------------------------------------------
# per-sector connected-component (queen's case) patch summary
# mirrors the sector_mask logic in 15B / 12B

patch_summary <- do.call(rbind, lapply(good_sectors, function(sec) {

  sector_r  <- terra::project(
                 terra::rast(file.path(ia_dir, "hirshpearson", paste0(sec, ".tif"))),
                 CanHF_1km, method = "near")
  sector_mask <- terra::ifel((sector_r > 0) & (CanHF_1km >= 1), 1, NA)

  patches_r <- terra::patches(sector_mask, directions = 8, zeroAsNA = TRUE)

  area_r    <- terra::cellSize(sector_mask, unit = "km")
  patch_areas <- terra::zonal(area_r, patches_r, fun = "sum", na.rm = TRUE)

  data.frame(
    sector        = sec,
    n_patches     = nrow(patch_areas),
    total_area_km2 = sum(patch_areas[, 2]),
    mean_area_km2  = mean(patch_areas[, 2]),
    median_area_km2 = median(patch_areas[, 2]),
    min_area_km2   = min(patch_areas[, 2]),
    max_area_km2   = max(patch_areas[, 2])
  )
}))

print(patch_summary)

