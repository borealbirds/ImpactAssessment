# ---
# title: Impact Assessment: conform Hirsh-Pearson cumulative impact layer to BAM
# author: Mannfred Boehm
# created: October 20, 2025
# ---

# the Hirsh-Pearson (2022) cumulative footprint score, 300 m, averaged onto the V5 1 km grid,
# then split into training (HF < 1) and backfilled (HF >= 1) masks.
# input: a local copy of G:/.../CovariateRasters/Disturbance/cum_threat2020.02.18.tif

library(terra)
terraOptions(memfrac = 0.6)   # project() returns different values under memory pressure

root    <- "G:/Shared drives/BAM_NationalModels5"
hp_dir  <- file.path("data", "raw_data", "hirshpearson")

bam_template <- terra::rast(file.path(root, "PredictionRasters", "Biomass", "SCANFI", "1km", "SCANFIBalsamFir_1km_2020.tif"))

# each 1 km cell's AREA MEAN. never bilinear: from 300 m it spreads every cell over 2 x 2 km
# cells, which put a ring of footprint-free cells into the high-HF mask (until 2026-09-28)
CanHF_1km <-
  terra::rast(file.path(hp_dir, "raw_300m", "cum_threat2020.02.18.tif")) |>
  terra::project(bam_template, method = "average") |>
  terra::mask(bam_template)

names(CanHF_1km) <- "CanHF_1km"
terra::writeRaster(CanHF_1km, file.path(hp_dir, "CanHF_1km_masked.tif"), overwrite = TRUE)



# ----------------------------------------------
# create low and high human footprint layers

# re-read, so this section runs on its own without raw_300m/
CanHF_1km <- terra::rast(file.path(hp_dir, "CanHF_1km_masked.tif"))

# low HF is < 1. Hirsh-Pearson: "we found that 82% of Canada's land areas had a
# HF < 1 and therefore were considered intact"
lowhf_mask <-
  terra::lapp(CanHF_1km, \(x) ifelse(!is.finite(x), NA, ifelse(x < 1, 1, NA))) |>
  as.factor()

terra::writeRaster(lowhf_mask, file.path(hp_dir, "CanHF_1km_lessthan1.tif"), overwrite = TRUE)

highhf_mask <-
  terra::lapp(CanHF_1km, \(x) ifelse(!is.finite(x), NA, ifelse(x >= 1, 1, NA))) |>
  as.factor()

terra::writeRaster(highhf_mask, file.path(hp_dir, "CanHF_1km_morethan1.tif"), overwrite = TRUE)
