# Smoke test for Fix B (impute-and-flag in 08A) on the local machine.
# Mirrors 07_train_and_backfill.R setup (cc=FALSE, local=TRUE -> ia_dir=getwd(),
# uses the CAfire-FIXED local covariates_mosaiced_2020.tif), runs
# train_and_backfill_subbasin_s on a few subbasins, then reads the written
# subbasin_{i}_backfill.tif and reports the FINITE fraction of each continuous
# `_draw_*` layer over the high-HF backfill cells. Pre-fix this was ~1-3%
# (NaN cascade); post-fix it should approach ~100% (minus genuinely degenerate cov).
suppressMessages({library(BART); library(BAMexploreR); library(terra); library(tidyverse)})
terra::terraOptions(progress = 0)

cc <- FALSE; local <- TRUE
ia_dir <- getwd()
cat("ia_dir =", ia_dir, "\n")
subbasins <- as.integer(strsplit(Sys.getenv("SMOKE_SUBS", "1,3"), ",")[[1]])
cat("subbasins =", paste(subbasins, collapse=","), "\n")

# ---- predictor/response metadata (copied verbatim from 07) -------------------
predictor_metadata <- dplyr::tibble(BAMexploreR::predictor_metadata) |>
  dplyr::filter(version == "v5") |>
  dplyr::select(predictor, definition, predictor_class) |>
  dplyr::mutate(dplyr::across('predictor', stringr::str_replace, 'Year', 'year')) |>
  dplyr::mutate(dplyr::across('predictor', stringr::str_replace, 'Method','method'))
soil_covs <- tibble::tibble(predictor = c("cec_0-5cm_mean_1000","cec_100-200cm_mean_1000",
  "cec_15-30cm_mean_1000","cec_30-60cm_mean_1000","cec_5-15cm_mean_1000","cec_60-100cm_mean_1000",
  "soc_0-5cm_mean_1000","soc_100-200cm_mean_1000","soc_15-30cm_mean_1000","soc_30-60cm_mean_1000",
  "soc_5-15cm_mean_1000","soc_60-100cm_mean_1000"), predictor_class = rep("Soil Properties", 12))
actually_biotic_what <- c("Peatland_5x5", "Peatland_1km")
actually_biotic_df <- tibble::tibble(predictor = actually_biotic_what, predictor_class = c("Wetland","Wetland"))
abiotic_vars <- predictor_metadata |>
  dplyr::filter(predictor_class %in% c("Annual Climate","Climate Normals","Topography","Wetland","Disturbance")) |>
  tibble::add_row(predictor = "CAfire", predictor_class = "Time Since Disturbance") |>
  dplyr::filter(!(predictor %in% actually_biotic_what)) |>
  dplyr::bind_rows(soil_covs)
biotic_vars <- predictor_metadata |>
  dplyr::filter(!(predictor_class %in% c(abiotic_vars$predictor_class, "Time", "Method"))) |>
  dplyr::bind_rows(actually_biotic_df)
neworder <- readRDS(file.path(ia_dir, "data", "raw_data", "biotic_variable_hierarchy.rds"))
biotic_vars <- biotic_vars[match(neworder, biotic_vars$predictor), ]

make_logger <- function(logfile) {
  dir.create(dirname(logfile), recursive = TRUE, showWarnings = FALSE)
  function(fmt, ...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%T"), sprintf(fmt, ...)),
                         file = logfile, append = TRUE)
}
source(file.path(ia_dir, "Rscripts", "08A_train_and_backfill_subbasin_s.R"))

year <- 2020
stack_y <- terra::rast(file.path(ia_dir, "data", "raw_data", "covariates_mosaiced",
                                 sprintf("covariates_mosaiced_%d.tif", year)))
categorical_responses <- c("ABoVE_1km","NLCD_1km","MODISLCC_1km","MODISLCC_5x5","SCANFI_1km","VLCE_1km")
lowhf_mask  <- terra::project(terra::rast(file.path(ia_dir,"data","raw_data","hirshpearson","CanHF_1km_lessthan1.tif")), stack_y, method="near")
highhf_mask <- terra::project(terra::rast(file.path(ia_dir,"data","raw_data","hirshpearson","CanHF_1km_morethan1.tif")), stack_y, method="near")
all_subbasins_subset <- terra::project(terra::vect(file.path(ia_dir,"data","raw_data","hydrobasins_masked_merged_subset.gpkg")), stack_y)

for (si in subbasins) {
  cat("\n================ subbasin", si, "================\n")
  res <- tryCatch(
    train_and_backfill_subbasin_s(
      subbasin_index = si, year = year, stack_y = stack_y,
      lowhf_mask = lowhf_mask, highhf_mask = highhf_mask,
      abiotic_vars = abiotic_vars, biotic_vars = biotic_vars, ia_dir = ia_dir,
      quiet = FALSE, neworder = neworder, categorical_responses = categorical_responses,
      all_subbasins_subset = all_subbasins_subset, cc = cc),
    error = function(e) { cat("ERROR:", conditionMessage(e), "\n"); NULL })
  if (is.null(res)) next

  tif <- file.path(ia_dir, "data", "derived_data", "bart_models", year,
                   sprintf("subbasin_%s", si), sprintf("subbasin_%d_backfill.tif", si))
  if (!file.exists(tif)) { cat("no tif written (no layers created)\n"); next }
  r <- terra::rast(tif)
  draw1 <- grep("_draw_001$", names(r), value = TRUE)   # one draw per continuous biotic cov
  cat(sprintf("layers=%d  continuous biotic covs (via _draw_001)=%d\n", terra::nlyr(r), length(draw1)))
  for (ln in draw1) {
    v  <- terra::values(r[[ln]], mat = FALSE)
    bf <- v[!is.na(v) | is.nan(v)]                 # backfill cells = filled (NaN counts as filled-but-bad)
    # backfill cells are exactly the cells 08A wrote to; reconstruct as non-NA-template cells
    cells_written <- which(!is.na(terra::values(r[[1]], mat = FALSE)) | is.nan(terra::values(r[[1]], mat = FALSE)))
    vv <- v[cells_written]
    cat(sprintf("  %-28s  written=%d  finite=%.1f%%  NaN=%.1f%%\n",
                sub("_draw_001$","",ln), length(vv),
                100*mean(is.finite(vv)), 100*mean(is.nan(vv))))
  }
}
cat("\nSMOKE DONE\n")
