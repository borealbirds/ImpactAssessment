# ---
# title: Local regression test of 12F on CAWA can60 (real Fir inputs, 2 bootstraps)
# author: Mannfred Boehm
# ---
# Runs predict_species_all_coalitions() on the local mirror of Fir's CAWA can60 inputs in
# cluster_logs/localtest/ (gitignored; pulled from Fir 2026-09-24) and checks:
#   (1) all 255 tables and 9 arrays identical() to Fir's smoke 7 (cluster_logs/smoke7/), which
#       holds while 12F's predictions and masks are unchanged (a mask change, e.g. C2f, breaks it
#       by design);
#   (2) the direct part (phi_direct) summed over sectors and subbasins equals d0 - obs
#       recomputed independently per bootstrap.
# Run from the repo root: Rscript Rscripts/misc/verify_12f_can60_local.R   (~4 min, one core)
# Merged from the 2026-09-25 session harnesses (harness_setup.R + harness6.R).
suppressPackageStartupMessages({ library(BAMexploreR); library(gbm); library(terra); library(tidyverse) })

REPO <- "C:/Users/mannf/Drive/boreal_avian_modelling_project/ImpactAssessment"
LT   <- file.path(REPO, "cluster_logs/localtest")
Sys.setenv(TEST_BCR = "can60", TEST_N_BOOT = "2")

# ---- 12D's setup, paths pointed at the local mirror ----------------------------
nm_root <- file.path(LT, "nm"); ia_dir <- file.path(LT, "ia")
load(file.path(ia_dir, "data", "raw_data", "SpeciesPredictionTruncationValues.Rdata"))
bam_boundary         <- terra::vect(file.path(ia_dir, "data", "raw_data", "Regions", "BAM_BCR_NationalModel_Unbuffered.shp"))
all_subbasins_subset <- terra::vect(file.path(ia_dir, "data", "raw_data", "hydrobasins_masked_merged_subset.gpkg"))
bcr_subbasins_ref <- {
  hits <- terra::relate(all_subbasins_subset, bam_boundary, relation = "intersects")
  ij   <- which(hits, arr.ind = TRUE)
  tibble(sub_index = ij[, 1], HYBAS_ID = all_subbasins_subset$first_HYBAS_ID[ij[, 1]],
         bcr_label = paste(bam_boundary$country[ij[, 2]], bam_boundary$subUnit[ij[, 2]], sep = "_"),
         bcr_code  = gsub("_", "", bcr_label))
}
categorical_responses <- c("ABoVE_1km", "NLCD_1km", "MODISLCC_1km", "MODISLCC_5x5", "SCANFI_1km", "VLCE_1km")
predictor_metadata <- dplyr::tibble(BAMexploreR::predictor_metadata) |> dplyr::filter(version == "v5") |>
  dplyr::select(predictor, definition, predictor_class) |>
  dplyr::mutate(predictor = sub("Year", "year", predictor), predictor = sub("Method", "method", predictor))
actually_biotic_what <- c("Peatland_5x5", "Peatland_1km")
abiotic_vars <- predictor_metadata |>
  dplyr::filter(predictor_class %in% c("Annual Climate", "Climate Normals", "Topography", "Wetland", "Disturbance", "Time", "Method")) |>
  dplyr::filter(!(predictor %in% actually_biotic_what))
biotic_continuous_vars <- predictor_metadata |> dplyr::filter(!(predictor_class %in% abiotic_vars$predictor_class)) |>
  dplyr::bind_rows(tibble::tibble(predictor = actually_biotic_what, predictor_class = "Wetland")) |>
  dplyr::filter(!predictor %in% categorical_responses) |> dplyr::pull(predictor)
disturbance_vars <- dplyr::tibble(BAMexploreR::predictor_metadata) |> dplyr::filter(version == "v5") |>
  dplyr::select(predictor, definition, predictor_class) |> dplyr::filter(predictor_class == "Disturbance")
source(file.path(REPO, "Rscripts", "12E_shapley_utils.R"))
source(file.path(REPO, "Rscripts", "12G_gbm_tree_walk.R"))
Rcpp::sourceCpp(file.path(REPO, "Rscripts", "12G_gbm_tree_walk.cpp"))
year <- 2020; species <- "CAWA"
hirsh_dir <- file.path(ia_dir, "data", "raw_data", "hirshpearson")
sectors <- canonical_sectors()
target_ids <- c(sectors_to_coalition_id(sectors, sectors),
                vapply(sectors, function(s) sectors_to_coalition_id(s, sectors), numeric(1L)))

Sys.setenv(TEST_N_BOOT = "2")
source(file.path(REPO, "Rscripts", "12F_predict_species_all_coalitions.R"))
t <- system.time(res <- suppressWarnings(predict_species_all_coalitions(
  species, year = year, all_subbasins_subset = all_subbasins_subset, hirsh_dir = hirsh_dir,
  save_arrays_ids = coalition_array_ids(),
  rdata_files = file.path(nm_root, "output/06_bootstraps/CAWA/CAWA_can60.Rdata"),
  return_per_bcr = TRUE)))[["elapsed"]]
new <- res[[1]]; saveRDS(new, file.path(REPO, "cluster_logs/localtest/can60_2_12f.rds"))
ref <- readRDS(file.path(REPO, "cluster_logs/smoke7/CAWA_2020_can60.rds"))$result
cat(sprintf("== 12F 2 boots: %.0f s | tables identical: %d/255 | arrays identical: %d/9\n", t,
    sum(mapply(identical, ref$coalition_tables, new$coalition_tables)),
    sum(mapply(identical, ref$bcr_arrays, new$bcr_arrays))))
pd <- new$shapley_samples$phi_direct; cat("phi_direct dims:", dim(pd), "\n")
direct_N <- apply(pd, 3, sum)                                  # v_direct(N), BCR, per bootstrap

# independent d0 - obs over the full-coalition kept pixels
e <- new.env(); load(file.path(nm_root, "output/06_bootstraps/CAWA/CAWA_can60.Rdata"), envir = e); bl <- e$b.list[1:2]
load(file.path(ia_dir, "data/raw_data/SpeciesPredictionTruncationValues.Rdata"))
qsp <- q.out[q.out$spp == "CAWA", ]$densmax
q99 <- readRDS(file.path(ia_dir, "data/derived_data/predictions/CAWA/truncation_params.rds"))[["can60_2020"]]$q99
st  <- terra::rast(file.path(nm_root, "gis/stacks/can60_2020.tif"))
prj <- function(f) terra::values(terra::project(terra::rast(file.path(hirsh_dir, f)), st[[1]], method = "near"), mat = FALSE)
S   <- sapply(canonical_sectors(), function(s) { v <- prj(paste0(s, ".tif")) > 0; v[is.na(v)] <- FALSE; v })
hi  <- prj("CanHF_1km_morethan1.tif")
subs <- unique(bcr_subbasins_ref$sub_index[bcr_subbasins_ref$bcr_code == "can60"])
z   <- terra::values(terra::rasterize(all_subbasins_subset[subs, ], st[[1]], field = "first_HYBAS_ID"), mat = FALSE)
w   <- terra::values(terra::rast(file.path(ia_dir, "data/derived_data/predictions/CAWA/can60/2020/weight.tif")), mat = FALSE)
cells <- which(rowSums(S) > 0 & !is.na(hi) & hi >= 1 & !is.na(z))
X <- terra::values(st[[unique(unlist(lapply(bl, `[[`, "var.names")))]])[cells, , drop = FALSE]
X[is.nan(X)] <- NA_real_; X <- as.data.frame(X, check.names = FALSE)
cv <- intersect(bl[[1]]$var.names, categorical_responses); dv <- intersect(disturbance_vars$predictor, bl[[1]]$var.names)
ob <- terra::rast(file.path(ia_dir, "data/derived_data/predictions/CAWA/can60/2020/observed_bootstraps.tif"))
ind <- sapply(1:2, function(i) {
  Xd <- X; for (v in dv) Xd[[v]] <- 0; Xd <- as_model_factors(Xd, bl[[i]], cv)
  ww <- w[cells]; ok <- is.finite(ww) & ww != 0
  pd0 <- numeric(length(cells)); pd0[ok] <- gbm_predict_design(bl[[i]], gbm_design(bl[[i]], Xd[ok, , drop = FALSE]))
  pd0 <- pmin(pmin(pd0, qsp), q99) * ww * 100
  po  <- terra::values(ob[[i]], mat = FALSE)[cells] * ww * 100; po[is.na(po)] <- 0
  sum(pd0, na.rm = TRUE) - sum(po)
})
cat("footprint covariates in model:", paste(dv, collapse = ", "), "\n")
cat("v_direct(N) per bootstrap | 12F phi_direct:", round(direct_N, 3), "| independent:", round(ind, 3),
    "| max rel diff:", signif(max(abs(direct_N - ind) / pmax(abs(ind), 1)), 3), "\n")
