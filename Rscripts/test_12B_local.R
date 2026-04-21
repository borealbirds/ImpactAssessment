# ---
# title: Local smoke-test for 12B/12C pipeline
# Generates synthetic data for one species x one small BCR, runs
# predict_species_bcr(), and validates the output structure.
# Run with: Rscript Rscripts/test_12B_local.R
# ---

suppressPackageStartupMessages({
  library(BAMexploreR)
  library(gbm)
  library(terra)
  library(tidyverse)
})


# ---- 1. Setup (mirrors 12B header exactly) ----

cc    <- FALSE
local <- TRUE

if (cc)            { ia_dir <- "/home/mannfred/scratch/impact_assessment" }
if (!cc && local)  { ia_dir <- getwd() }
if (!cc && !local) { ia_dir <- file.path("G:/Shared drives/BAM_NationalModels5", "data", "Extras", "sandbox_data", "impactassessment_sandbox") }

nm_root <- file.path(ia_dir, "test_data", "fake_nm_root")

bam_boundary         <- terra::vect(file.path(ia_dir, "data", "raw_data", "Regions", "BAM_BCR_NationalModel_Unbuffered.shp"))
all_subbasins_subset <- terra::vect(file.path(ia_dir, "data", "raw_data", "hydrobasins_masked_merged_subset.gpkg"))

bcr_subbasins_ref <- {
  hits <- terra::relate(all_subbasins_subset, bam_boundary, relation = "intersects")
  ij   <- which(hits, arr.ind = TRUE)
  tibble(
    sub_index = ij[, 1],
    HYBAS_ID  = all_subbasins_subset$first_HYBAS_ID[ij[, 1]],
    bcr_label = paste(bam_boundary$country[ij[, 2]], bam_boundary$subUnit[ij[, 2]], sep = "_"),
    bcr_code  = gsub("_", "", bcr_label)
  )
}

categorical_responses <- c("ABoVE_1km", "NLCD_1km", "MODISLCC_1km", "MODISLCC_5x5", "SCANFI_1km", "VLCE_1km")

predictor_metadata <-
  dplyr::tibble(BAMexploreR::predictor_metadata) |>
  dplyr::filter(version == "v5") |>
  dplyr::select(predictor, definition, predictor_class) |>
  dplyr::mutate(dplyr::across('predictor', stringr::str_replace, 'Year', 'year')) |>
  dplyr::mutate(dplyr::across('predictor', stringr::str_replace, 'Method', 'method'))

actually_biotic_what <- c("Peatland_5x5", "Peatland_1km")
actually_biotic_df   <- tibble::tibble(predictor = actually_biotic_what,
                                       predictor_class = c("Wetland", "Wetland"))

abiotic_vars <-
  predictor_metadata |>
  dplyr::filter(predictor_class %in% c("Annual Climate", "Climate Normals", "Topography",
                                        "Wetland", "Disturbance", "Time", "Method")) |>
  dplyr::filter(!(predictor %in% actually_biotic_what))

biotic_continuous_vars <-
  predictor_metadata |>
  dplyr::filter(!(predictor_class %in% abiotic_vars$predictor_class)) |>
  dplyr::bind_rows(actually_biotic_df) |>
  dplyr::filter(!predictor %in% categorical_responses) |>
  dplyr::pull(predictor)

disturbance_vars <-
  dplyr::tibble(BAMexploreR::predictor_metadata) |>
  dplyr::filter(version == "v5") |>
  dplyr::select(predictor, definition, predictor_class) |>
  dplyr::filter(predictor_class == "Disturbance")

industry_rasters <- list.files(file.path(ia_dir, "data", "raw_data", "hirshpearson"),
                                pattern = "\\.tif$", full.names = TRUE)
industry_names  <- sub(".tif", "", basename(industry_rasters))
industry_stack  <- setNames(lapply(industry_rasters, terra::rast), industry_names)

source(file.path(ia_dir, "Rscripts", "shapley_utils.R"))
source(file.path(ia_dir, "Rscripts", "12C_predict_species_bcr.R"))


# ---- 2. Test parameters ----

test_species      <- "BANS"
test_year         <- 2020
sectors           <- canonical_sectors()
test_coalition    <- c("roads")
test_coalition_id <- sectors_to_coalition_id(test_coalition, sectors)

# pick BCR with fewest subbasins to minimise compute time
test_bcr <- bcr_subbasins_ref |>
  dplyr::count(bcr_code) |>
  dplyr::arrange(n) |>
  dplyr::slice(1) |>
  dplyr::pull(bcr_code)

message("Test BCR: ", test_bcr, "  coalition_id: ", test_coalition_id)


# ---- 3. Create spatial template (coarse resolution, within test BCR) ----

bam_boundary$bcr_code_col <- gsub("_", "", paste(bam_boundary$country,
                                                   bam_boundary$subUnit, sep = "_"))
test_bcr_poly <- bam_boundary[bam_boundary$bcr_code_col == test_bcr, ]

hirsh_dir <- file.path(ia_dir, "data", "raw_data", "hirshpearson")
canHF_raw <- terra::rast(file.path(hirsh_dir, "CanHF_1km_morethan1.tif"))

test_bcr_proj  <- terra::project(test_bcr_poly, terra::crs(canHF_raw))
template_raw   <- terra::crop(canHF_raw, test_bcr_proj, mask = TRUE)
template       <- terra::aggregate(template_raw, fact = 5, fun = "max", na.rm = TRUE)
rm(canHF_raw, template_raw); gc()

n_cells <- terra::ncell(template)
message("Template: ", terra::nrow(template), " x ", terra::ncol(template),
        " (", n_cells, " cells, res=", round(terra::res(template)[1] / 1000, 1), " km)")


# ---- 4. Create output directories ----

dirs_needed <- c(
  file.path(nm_root, "output", "06_bootstraps", test_species),
  file.path(nm_root, "gis", "stacks"),
  file.path(ia_dir, "data", "derived_data", "bart_models_mosaics", test_year),
  file.path(ia_dir, "data", "derived_data", "predictions", test_species, test_bcr, test_year)
)
invisible(lapply(dirs_needed, dir.create, recursive = TRUE, showWarnings = FALSE))


# ---- 5. Synthesize GBM bootstrap models (b.list) ----

# Use all abiotics + 5 biotic continuous + all categorical.
# Full 46-biotic set causes 4,600 terra::values() calls (layer-by-layer on a
# 4,606-layer raster) which is prohibitively slow locally — 5 is enough to
# exercise every code path in predict_species_bcr().
n_biotic_test  <- 5L
biotic_subset  <- biotic_continuous_vars[seq_len(min(n_biotic_test, length(biotic_continuous_vars)))]
all_pred_names <- unique(c(abiotic_vars$predictor, biotic_subset, categorical_responses))
n_preds        <- length(all_pred_names)
set.seed(42)
fake_df      <- as.data.frame(matrix(rnorm(200 * n_preds), nrow = 200, ncol = n_preds,
                                      dimnames = list(NULL, all_pred_names)))
fake_df$y    <- rnorm(200)

message("Training synthetic GBM (", n_preds, " predictors, 10 trees)...")
m <- gbm::gbm(y ~ ., data = fake_df, n.trees = 10, distribution = "gaussian",
               interaction.depth = 1, verbose = FALSE)
attr(m, "bcr") <- test_bcr
b.list <- list(m, m)   # 2 bootstraps — sufficient to exercise loop logic

save(b.list, file = file.path(nm_root, "output", "06_bootstraps", test_species,
                               paste0("can_", test_bcr, ".Rdata")))
message("Saved fake b.list")


# ---- 6. Synthesize BCR covariate stack ----

var_names <- m$var.names
set.seed(43)
stack_obs_vals <- matrix(rnorm(n_cells * length(var_names)), nrow = n_cells)
stack_obs      <- terra::rast(template, nlyrs = length(var_names))
terra::values(stack_obs) <- stack_obs_vals
names(stack_obs) <- var_names

terra::writeRaster(stack_obs,
                   file.path(nm_root, "gis", "stacks",
                              paste0(test_bcr, "_", test_year, ".tif")),
                   overwrite = TRUE)
message("Saved fake BCR covariate stack (", length(var_names), " layers)")
rm(stack_obs_vals); gc()


# ---- 7. Synthesize BART backfill mosaic ----

biotic_cont_in_model <- intersect(setdiff(var_names, categorical_responses),
                                   biotic_continuous_vars)
cat_in_model         <- intersect(var_names, categorical_responses)
n_draws              <- 100L

draw_layer_names <- unlist(lapply(biotic_cont_in_model, function(v)
  sprintf("%s_draw_%03d", v, seq_len(n_draws))))
cat_layer_names  <- cat_in_model
all_bf_names     <- c(draw_layer_names, cat_layer_names)
n_bf_layers      <- length(all_bf_names)

message("Building BART mosaic: ", length(biotic_cont_in_model), " cont vars x 100 draws + ",
        length(cat_in_model), " cat vars = ", n_bf_layers, " layers")

set.seed(44)
bf_vals  <- matrix(abs(rnorm(n_cells * n_bf_layers, mean = 0.5)), nrow = n_cells)
stack_bf <- terra::rast(template, nlyrs = n_bf_layers)
terra::values(stack_bf) <- bf_vals
names(stack_bf) <- all_bf_names

terra::writeRaster(stack_bf,
                   file.path(ia_dir, "data", "derived_data", "bart_models_mosaics",
                              test_year, paste0(test_bcr, "_backfilled.tif")),
                   overwrite = TRUE)
message("Saved fake BART mosaic")
rm(bf_vals, stack_bf); gc()


# ---- 8. Synthesize observed bootstraps ----

set.seed(45)
obs_vals <- matrix(abs(rnorm(n_cells * 2, mean = 0.01, sd = 0.005)), nrow = n_cells)
obs_boot <- terra::rast(template, nlyrs = 2)
terra::values(obs_boot) <- obs_vals

terra::writeRaster(obs_boot,
                   file.path(ia_dir, "data", "derived_data", "predictions",
                              test_species, test_bcr, test_year, "observed_bootstraps.tif"),
                   overwrite = TRUE)
message("Saved fake observed bootstraps (2 layers)")
rm(obs_vals, obs_boot); gc()


# ---- 9. Run predict_species_bcr() ----

message("\n--- Running predict_species_bcr() ---")
result <- predict_species_bcr(
  species              = test_species,
  year                 = test_year,
  all_subbasins_subset = all_subbasins_subset,
  coalition            = test_coalition,
  coalition_id         = test_coalition_id,
  hirsh_dir            = hirsh_dir
)


# ---- 10. Validate output ----

message("\n--- Validating output ---")
expected_cols <- c("species", "subbasin", "bcr", "coalition_id",
                   "obs_total_mean", "obs_total_sd",
                   "obs_on_coalition_mean", "obs_on_coalition_sd",
                   "bf_on_coalition_mean", "bf_on_coalition_sd")

stopifnot(
  "result is a data.frame"      = is.data.frame(result),
  "result has rows"             = nrow(result) > 0,
  "expected columns present"   = all(expected_cols %in% names(result)),
  "species matches"             = all(result$species == test_species),
  "coalition_id matches"        = all(result$coalition_id == test_coalition_id),
  "obs_total_mean is finite"    = all(is.finite(result$obs_total_mean)),
  "bf_on_coalition_mean finite" = all(is.finite(result$bf_on_coalition_mean)),
  # NOTE: >= 0 not enforced — synthetic GBM on random data produces negative
  # predictions (gaussian, no log link); real models trained on bird data do not.
  "obs_total_mean is finite"    = all(is.finite(result$obs_total_mean)),
  "bf_on_coalition_mean finite2"= all(is.finite(result$bf_on_coalition_mean))
)

message("Result: ", nrow(result), " subbasin row(s) for BCR=", test_bcr)
message("ALL CHECKS PASSED")


# ---- 11. Cleanup ----

unlink(file.path(ia_dir, "test_data"), recursive = TRUE)
unlink(file.path(ia_dir, "data", "derived_data", "bart_models_mosaics", test_year,
                 paste0(test_bcr, "_backfilled.tif")))
unlink(file.path(ia_dir, "data", "derived_data", "predictions",
                 test_species, test_bcr), recursive = TRUE)
unlink(file.path(ia_dir, "data", "derived_data", "predictions_coalitions",
                 as.character(test_coalition_id), test_species, test_bcr),
       recursive = TRUE)

message("Test data cleaned up.")
