# ---
# title: Impact Assessment: train models per subbasin and backfill industry footprints
# author: Mannfred Boehm
# created: August 7, 2025
# ---

# one subbasin per task: BART models trained on the subbasin's low HF pixels predict its biotic
# covariates on its high HF pixels (08A) -> data/derived_data/bart_models/{year}/subbasin_{i}/
# on Fir: bash 07_train_and_backfill.sh 2020   (11 is chained after)
# locally: Rscript Rscripts/07_train_and_backfill.R <subbasin>
# a subbasin's error is caught, so the task still exits 0: grep the logs for "Error in subbasin"

library(BART)
library(terra)

ia_dir  <- if (nzchar(Sys.getenv("SLURM_JOB_ID"))) "/home/mannfred/scratch/impact_assessment" else getwd()
raw_dir <- file.path(ia_dir, "data", "raw_data")

subbasin_index <- as.integer(commandArgs(trailingOnly = TRUE)[1])
year           <- as.integer(Sys.getenv("YEARS", "2020"))   # the .sh submits one year per job
message("year ", year, ", subbasin ", subbasin_index)

for (f in c("08A_train_and_backfill_subbasin_s.R", "08B_deploy_gbart.R", "08B_deploy_mbart.R",
            "08B_deploy_pbart.R", "09_collect_metrics.R"))
  source(file.path(ia_dir, "Rscripts", f))



# -----------------------
# covariates

pm <- dplyr::filter(BAMexploreR::predictor_metadata, version == "v5")
peatland <- c("Peatland_1km", "Peatland_5x5")

soil_covs <- c("cec_0-5cm_mean_1000", "cec_100-200cm_mean_1000", "cec_15-30cm_mean_1000",
               "cec_30-60cm_mean_1000", "cec_5-15cm_mean_1000", "cec_60-100cm_mean_1000",
               "soc_0-5cm_mean_1000", "soc_100-200cm_mean_1000", "soc_15-30cm_mean_1000",
               "soc_30-60cm_mean_1000", "soc_5-15cm_mean_1000", "soc_60-100cm_mean_1000")

# predictors: V5's climate, topography, wetland and footprint covariates, plus CAfire and soil
abiotic_vars <- c(setdiff(pm$predictor[pm$predictor_class %in% c("Annual Climate", "Climate Normals",
                                                                 "Topography", "Wetland", "Disturbance")],
                          peatland),
                  "CAfire", soil_covs)

# V5's footprint covariates ("Disturbance"), set to 0 where the backfill is predicted (08A)
footprint_vars <- pm$predictor[pm$predictor_class == "Disturbance"]

# backfilled: V5's land cover and vegetation, plus peatland, in 00's order
biotic_vars <- intersect(readRDS(file.path(raw_dir, "biotic_variable_hierarchy.rds")),
                         c(pm$predictor[pm$predictor_class %in% c("Landcover", "Veg Species", "Veg Structure")],
                           peatland))

categorical_responses <- c("ABoVE_1km", "NLCD_1km", "MODISLCC_1km", "MODISLCC_5x5", "SCANFI_1km", "VLCE_1km")

# each land-cover layer's open-water class, left out of its training pixels so the backfill
# never predicts water (SCANFI 8, VLCE 20, MODIS IGBP 17). ABoVE and NLCD are not in the mosaic
water_classes <- list(SCANFI_1km = 8, VLCE_1km = 20, MODISLCC_1km = 17, MODISLCC_5x5 = 17)



# -----------------------
# import data

stack_y <- terra::rast(file.path(raw_dir, "covariates_mosaiced", sprintf("covariates_mosaiced_%d.tif", year)))

# training (HF < 1) and backfilled (HF >= 1) pixels. Hirsh-Pearson exists only for 2020, so
# another year needs its own dated masks, or HF_MASK_YEAR=2020 to use 2020's
mask_year <- as.integer(Sys.getenv("HF_MASK_YEAR", year))
hf_mask <- function(kind) {
  f <- file.path(raw_dir, "hirshpearson",
                 paste0("CanHF_1km_", kind, if (mask_year != 2020) paste0("_", mask_year), ".tif"))
  if (!file.exists(f)) stop("no footprint mask ", f, "; to use the 2020 masks, set HF_MASK_YEAR=2020")
  terra::project(terra::rast(f), stack_y, method = "near")
}
lowhf_mask  <- hf_mask("lessthan1")
highhf_mask <- hf_mask("morethan1")

all_subbasins_subset <-
  terra::vect(file.path(raw_dir, "hydrobasins_masked_merged_subset.gpkg")) |>
  terra::project(stack_y)



# -----------------------
# train and backfill

backfill_results <- tryCatch(
  train_and_backfill_subbasin_s(
    subbasin_index        = subbasin_index,
    year                  = year,
    stack_y               = stack_y,
    lowhf_mask            = lowhf_mask,
    highhf_mask           = highhf_mask,
    all_subbasins_subset  = all_subbasins_subset,
    abiotic_vars          = abiotic_vars,
    footprint_vars        = footprint_vars,
    biotic_vars           = biotic_vars,
    categorical_responses = categorical_responses,
    water_classes         = water_classes,
    ia_dir                = ia_dir),
  error = function(e) {
    message("Error in subbasin ", subbasin_index, ": ", conditionMessage(e))
    list(subbasin = subbasin_index, error = conditionMessage(e))
  })

print(backfill_results)
