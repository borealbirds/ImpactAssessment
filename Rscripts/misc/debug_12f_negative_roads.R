# ---
# title: Test 3 - Does 12F's negative roads sign come from the road-pixel
#        definition, or from mean-vs-sampled BART draws (Jensen's inequality)?
# author: Mannfred Boehm
# ---
# Builds on 13B and 13C. Every isolated test so far gives the EXPECTED positive
# sign (zeroing disturbance at road pixels: +14.4%; BART-mean biotic raises the
# BRT prediction +20% vs observed), yet 12F's roads-only coalition returns a
# NEGATIVE impact for CAWA/OVEN. Two candidate causes remain:
#
#   (1) ROAD-PIXEL DEFINITION MISMATCH. 13B/13C define road pixels as
#       canroad_1km > 0 (the buffered road-density predictor in the BRT stack).
#       12C/12F define coalition pixels as roads.tif > 0 AND
#       CanHF_1km_morethan1.tif >= 1 (raw binary footprint in
#       data/raw_data/hirshpearson/). These are DIFFERENT pixel sets; 12F's
#       negative could live entirely in pixels 13B/13C never tested.
#
#   (2) MEAN-vs-SAMPLED DRAWS (Jensen's inequality). 13C Test C feeds the
#       per-pixel MEAN of 100 BART draws through the BRT once: BRT(E[draw]).
#       12C/12F instead SAMPLE one BART draw per scenario, predict, then average
#       the predictions: E[BRT(draw)]. The BRT is highly nonlinear, so
#       E[BRT(draw)] != BRT(E[draw]). If BART's posterior at road pixels has a
#       heavy low-vegetation tail, sampling can pull the mean prediction down.
#
# This script reproduces 12F's impact recipe LOCALLY for CAWA/can14, BRT
# bootstrap #1, and slices it two ways:
#
#   Test 1: impact (= bf - obs) per pixel-set partition (R-only, C-only, R&C),
#           holding biotic backfill fixed at the per-pixel mean draw.
#   Test 2: on the 12F coalition set, impact under BRT(E[draw]) (mean draw)
#           vs E[BRT(draw)] (sampled draws). Gap = Jensen effect.
#
# 12F impact recipe reproduced here (matches 12C_predict_species_bcr.R):
#   obs density = canonical observed_bootstraps.tif layer `boot_i`, x100.
#                 (already clamped by 12A steps 5-9; 12F's obs_on_coalition
#                  ingredient.)
#   bf  density = BRT on: disturbance predictors zeroed + categorical biotic
#                 swapped to BART backfill + continuous biotic swapped to BART
#                 backfill; pmin(pred, qsp) x100.
#   impact      = sum_pixels (bf - obs).  Negative => fewer birds predicted
#                 after removing the sector => 12F's backwards sign.
#
# NOTE ON TRUNCATED HANDOFF: the spec for this script was cut off mid-"Test 1".
# Decisions made here: (a) obs baseline uses the canonical observed_bootstraps
# (exactly 12F's obs_on_coalition source) rather than a local re-prediction;
# (b) Test 2's sampled estimate iterates all 100 draws deterministically -- the
# exact expectation 12C's with-replacement sampling converges to.
#
# Run from project root:  Rscript Rscripts/13D_road_definition_and_sampling.R
# ---

suppressPackageStartupMessages({
  library(terra)
  library(gbm)
  library(dplyr)
})

# ---- Execution context (LOCAL ONLY) -----------------------------------------

cc      <- FALSE
local   <- TRUE
ia_dir  <- getwd()
nm_root <- "G:/Shared drives/BAM_NationalModels5"

species  <- "CAWA"
bcr_code <- "can14"
year     <- 2020
boot_i   <- 1L

# ---- Paths ------------------------------------------------------------------

bcr_stack_path <- file.path(nm_root, "gis/stacks",
                            paste0(bcr_code, "_", year, ".tif"))
bf_mosaic_path <- file.path(ia_dir, "data/derived_data/bart_models_mosaics",
                            year, paste0(bcr_code, "_backfilled.tif"))
obs_boot_path  <- file.path(ia_dir, "data/derived_data/predictions",
                            species, bcr_code, year, "observed_bootstraps.tif")
hirsh_dir      <- file.path(ia_dir, "data/raw_data/hirshpearson")
roads_path     <- file.path(hirsh_dir, "roads.tif")
canhf_path     <- file.path(hirsh_dir, "CanHF_1km_morethan1.tif")

stopifnot(file.exists(bcr_stack_path))
stopifnot(file.exists(bf_mosaic_path))
stopifnot(file.exists(obs_boot_path))
stopifnot(file.exists(roads_path))
stopifnot(file.exists(canhf_path))

# ---- Species clamp ----------------------------------------------------------

load(file.path(ia_dir, "data/raw_data/SpeciesPredictionTruncationValues.Rdata"))
qsp <- q.out[q.out$spp == species, ]$q
message("species qsp clamp = ", qsp)

# ---- Load BRT bootstrap (single rep) ----------------------------------------

rdata_path <- list.files(file.path(nm_root, "output/06_bootstraps", species),
                         pattern = paste0(bcr_code, ".*\\.Rdata$"),
                         full.names = TRUE)[1]
stopifnot(length(rdata_path) == 1, file.exists(rdata_path))
message("loading BRT model: ", basename(rdata_path), "  (bootstrap #", boot_i, ")")

e <- new.env(parent = emptyenv())
load(rdata_path, envir = e)
model <- e$b.list[[boot_i]]
rm(e); gc()

model_vars <- model$var.names
cat_vars   <- model_vars[model$var.type > 0]

# gbm's $var.levels is an UNNAMED list indexed by predictor position; for factor
# predictors each element holds the training class values (e.g. VLCE_1km ->
# 20,33,50,80,100,210,220,230). Indexing it by name (model$var.levels[["x"]])
# silently returns NULL — the bug that makes 12C skip factor conversion. Look it
# up by position so prediction-frame factors carry the model's exact levels.
cat_levels <- setNames(
  lapply(cat_vars, function(v) as.character(model$var.levels[[match(v, model_vars)]])),
  cat_vars)

# Disturbance-class predictors in this model (zeroed in the bf recipe; matches
# 12C's `dist_shared` and 13B's `dist_in_model`).
dist_all <- dplyr::tibble(BAMexploreR::predictor_metadata) |>
  dplyr::filter(version == "v5", predictor_class == "Disturbance") |>
  dplyr::pull(predictor)
dist_in_model <- intersect(dist_all, model_vars)
message("model: ", length(model_vars), " predictors; ",
        length(cat_vars), " categorical")
message("disturbance predictors zeroed in bf recipe: ",
        paste(dist_in_model, collapse = ", "))
stopifnot("canroad_1km" %in% model_vars)

# ---- Load BCR observed stack ------------------------------------------------

stack_obs <- terra::rast(bcr_stack_path)
message("BCR observed stack: ", terra::nlyr(stack_obs), " layers, ",
        terra::ncell(stack_obs), " cells")
stopifnot("CanHF_1km" %in% names(stack_obs),
          "canroad_1km" %in% names(stack_obs))

# ---- Load BART backfill mosaic ----------------------------------------------

bf_full <- terra::rast(bf_mosaic_path)
message("BART backfill mosaic: ", terra::nlyr(bf_full), " layers")
if (!isTRUE(terra::compareGeom(bf_full, stack_obs, stopOnError = FALSE,
                               messages = FALSE))) {
  stop("BART mosaic grid does not match the BCR stack grid. 12C indexes the\n",
       "mosaic by cell directly; 13D does the same. Re-run 11_premosaic so the\n",
       "mosaic is built on the BCR prediction grid.")
}

# ---- Load canonical observed bootstraps -------------------------------------

obs_stack <- terra::rast(obs_boot_path)
message("observed_bootstraps.tif: ", terra::nlyr(obs_stack), " bootstrap layers")
if (!isTRUE(terra::compareGeom(obs_stack, stack_obs, stopOnError = FALSE,
                               messages = FALSE))) {
  stop("observed_bootstraps.tif grid does not match the BCR stack grid.")
}
stopifnot(boot_i <= terra::nlyr(obs_stack))

# ---- Biotic variable lists --------------------------------------------------

hierarchy <- readRDS(file.path(ia_dir, "data/raw_data/biotic_variable_hierarchy.rds"))
categorical_biotic <- c("VLCE_1km", "MODISLCC_1km", "MODISLCC_5x5",
                        "ABoVE_1km", "SCANFI_1km", "NLCD_1km")
continuous_biotic  <- setdiff(hierarchy, categorical_biotic)

# continuous biotic vars that the BRT uses AND that have BART draw layers
draw_lyr_all <- grep("_draw_[0-9]{3}$", names(bf_full), value = TRUE)
draw_covs_avail <- unique(sub("_draw_[0-9]{3}$", "", draw_lyr_all))
cont_swap <- intersect(intersect(continuous_biotic, draw_covs_avail), model_vars)
n_draws   <- length(grep(paste0("^", cont_swap[1], "_draw_[0-9]{3}$"),
                         names(bf_full)))
message("continuous biotic vars swapped (in BRT model & BART draws): ",
        length(cont_swap), "  | draws per var: ", n_draws)

# categorical biotic vars the BRT uses AND that we can backfill
cat_swap <- intersect(cat_vars, categorical_biotic)
message("categorical biotic vars swapped to BART backfill: ",
        if (length(cat_swap)) paste(cat_swap, collapse = ", ") else "(none)")

# ---- Define the two pixel sets ----------------------------------------------

# Set R (13B/13C): buffered road-density predictor > 0
road_v   <- terra::values(stack_obs[["canroad_1km"]], mat = FALSE)
road_idx <- which(!is.na(road_v) & road_v > 0)

# Set C (12C/12F): raw binary roads footprint > 0 AND CanHF_morethan1 >= 1.
# Project the hirshpearson rasters onto the BCR grid exactly as 12C does.
roads_r <- terra::project(terra::rast(roads_path), stack_obs[[1]], method = "near")
canhf_r <- terra::project(terra::rast(canhf_path), stack_obs[[1]], method = "near")
roads_fp_v <- terra::values(roads_r, mat = FALSE)
canhf_fp_v <- terra::values(canhf_r, mat = FALSE)
coal_idx <- which(!is.na(roads_fp_v) & roads_fp_v > 0 &
                  !is.na(canhf_fp_v) & canhf_fp_v >= 1)

if (length(road_idx) == 0) stop("no canroad_1km > 0 pixels in ", bcr_code)
if (length(coal_idx) == 0) stop("no 12F roads-coalition pixels in ", bcr_code)

cells  <- sort(union(road_idx, coal_idx))
in_R   <- cells %in% road_idx
in_C   <- cells %in% coal_idx
partition <- ifelse(in_R & in_C, "R_and_C",
              ifelse(in_R & !in_C, "R_only", "C_only"))

cat("\n", strrep("=", 78), "\n",
    "Test 1: road-pixel definition  --  ", species, " ", bcr_code, " ", year, "\n",
    strrep("=", 78), "\n", sep = "")
cat(sprintf("Set R  (canroad_1km > 0, the 13B/13C definition) : %7d pixels\n",
            length(road_idx)))
cat(sprintf("Set C  (roads.tif > 0 & CanHF_morethan1 >= 1, 12F): %7d pixels\n",
            length(coal_idx)))
cat(sprintf("  R and C (both)        : %7d   (%.1f%% of C)\n",
            sum(in_R & in_C), 100 * sum(in_R & in_C) / length(coal_idx)))
cat(sprintf("  R only (in 13B, not 12F): %7d\n", sum(in_R & !in_C)))
cat(sprintf("  C only (in 12F, NEVER tested by 13B/13C): %7d   (%.1f%% of C)\n",
            sum(!in_R & in_C), 100 * sum(!in_R & in_C) / length(coal_idx)))

# why are R pixels excluded from C? (diagnostic)
r_roadsfp <- roads_fp_v[road_idx]; r_canhf <- canhf_fp_v[road_idx]
cat(sprintf("\n  of Set R: %d have roads.tif > 0 ; %d have CanHF_morethan1 >= 1\n",
            sum(!is.na(r_roadsfp) & r_roadsfp > 0),
            sum(!is.na(r_canhf) & r_canhf >= 1)))

# ---- Extract observed predictor matrix at `cells` ---------------------------

X_full <- terra::values(stack_obs, dataframe = TRUE)

# Coerce categoricals to factors carrying the gbm model's exact training levels.
# dataframe = TRUE returns each categorical column as a factor of raster class
# labels; refactor onto cat_levels so gbm sees the levels it was trained on.
# Observed class values absent from the training set (e.g. VLCE 32/40/81) become
# NA and drop out via complete.cases.
for (v in cat_vars) {
  if (!v %in% names(X_full)) next
  X_full[[v]] <- factor(as.character(X_full[[v]]), levels = cat_levels[[v]])
}
X_obs_cells <- X_full[cells, model_vars, drop = FALSE]
rm(X_full); gc()

# ---- Observed density at `cells` (12F's obs_on_coalition ingredient) --------

obs_density <- terra::values(obs_stack[[boot_i]], mat = FALSE)[cells] * 100

# ---- BART backfill at `cells`: continuous draws -----------------------------
# Stored on log1p scale -> back-transform per draw (expm1), clean non-finite,
# clamp >= 0. Matches 12C lines 207-209 and 13C bart_mean_for_var().

message(Sys.time(), " | extracting ", length(cont_swap),
        " continuous biotic draw stacks at ", length(cells), " cells ...")
draws <- setNames(vector("list", length(cont_swap)), cont_swap)
for (v in cont_swap) {
  lyrs <- paste0(v, "_draw_", sprintf("%03d", seq_len(n_draws)))
  lyrs <- intersect(lyrs, names(bf_full))
  m <- terra::values(bf_full[[lyrs]], mat = TRUE)[cells, , drop = FALSE]
  m <- expm1(m)
  m[!is.finite(m)] <- 0
  m[m < 0] <- 0
  draws[[v]] <- m                      # [n_cells x n_draws]
}

# ---- BART backfill at `cells`: categorical codes ----------------------------
# 12C reads layer `v` (or `v_mean` fallback) as 1-based codes into model levels.

bf_cat_codes <- setNames(vector("list", length(cat_swap)), cat_swap)
for (v in cat_swap) {
  lyr <- if (v %in% names(bf_full)) v else paste0(v, "_mean")
  if (!lyr %in% names(bf_full)) {
    message("  WARNING: no backfilled layer for categorical var ", v,
            " — keeping observed values")
    bf_cat_codes[[v]] <- NULL
    next
  }
  bf_cat_codes[[v]] <- terra::values(bf_full[[lyr]], mat = FALSE)[cells]
}
cat_swap_use <- names(Filter(Negate(is.null), bf_cat_codes))

rm(bf_full); gc()

# ---- bf prediction recipe (matches 12C) -------------------------------------
# cont_vals: named list var -> numeric vector of length(cells), giving the
# continuous biotic values to inject. Returns birds/km^2 per cell, NA where
# the predictor row is incomplete.

predict_bf <- function(cont_vals) {
  X <- X_obs_cells
  for (dv in dist_in_model) X[[dv]] <- 0
  for (v in names(cont_vals)) if (v %in% names(X)) X[[v]] <- cont_vals[[v]]
  for (v in cat_swap_use) {
    if (!v %in% names(X)) next
    # Backfilled categorical rasters store actual class values (e.g. VLCE
    # 20,33,...,230), same domain as the observed layer. Refactor onto the
    # model's training levels; the mbart no-data code 0 falls out as NA.
    X[[v]] <- factor(as.character(bf_cat_codes[[v]]), levels = cat_levels[[v]])
  }
  ok <- complete.cases(X)
  p  <- rep(NA_real_, nrow(X))
  if (any(ok)) {
    p[ok] <- gbm::predict.gbm(model, X[ok, , drop = FALSE],
                              n.trees = model$n.trees, type = "response")
  }
  pmin(p, qsp) * 100
}

# ---- Test 1: impact per pixel-set partition (mean-draw biotic) ---------------

cont_meandraw <- lapply(draws, rowMeans)         # per-pixel mean over draws
bf_meandraw   <- predict_bf(cont_meandraw)
impact_md     <- bf_meandraw - obs_density       # 12F's per-pixel v(S)

ok1 <- is.finite(obs_density) & is.finite(bf_meandraw)
cat(sprintf("\ncomplete-case pixels (finite obs & bf): %d / %d\n",
            sum(ok1), length(cells)))

cat("\nimpact = bf(mean-draw) - obs   [birds/km^2; negative = 12F's backwards sign]\n")
cat(sprintf("%-22s %8s %12s %12s %12s %10s\n",
            "partition", "n", "sum_obs", "sum_bf", "sum_impact", "%neg_px"))
cat(strrep("-", 78), "\n")
part_summary <- list()
for (lab in c("R_and_C", "R_only", "C_only", "ALL_R", "ALL_C")) {
  sel <- switch(lab,
    R_and_C = ok1 & partition == "R_and_C",
    R_only  = ok1 & partition == "R_only",
    C_only  = ok1 & partition == "C_only",
    ALL_R   = ok1 & in_R,
    ALL_C   = ok1 & in_C)
  if (!any(sel)) { cat(sprintf("%-22s %8d\n", lab, 0)); next }
  so <- sum(obs_density[sel]); sb <- sum(bf_meandraw[sel])
  cat(sprintf("%-22s %8d %12.0f %12.0f %+12.0f %9.1f%%\n",
              lab, sum(sel), so, sb, sb - so,
              100 * mean(impact_md[sel] < 0)))
  part_summary[[lab]] <- c(n = sum(sel), sum_obs = so, sum_bf = sb,
                           sum_impact = sb - so)
}
cat(strrep("-", 78), "\n")
cat("Read: if C_only is strongly negative while R_* is positive, the negative\n",
    "12F roads sign lives in pixels 13B/13C never tested -> cause (1) confirmed.\n",
    sep = "")

# ---- Test 2: mean-draw vs sampled-draw on the 12F coalition set --------------

cat("\n", strrep("=", 78), "\n",
    "Test 2: BRT(E[draw]) vs E[BRT(draw)] on 12F coalition set C\n",
    "        (Jensen's inequality; both with disturbance zeroed + cat backfill)\n",
    strrep("=", 78), "\n", sep = "")

message(Sys.time(), " | predicting ", n_draws, " sampled-draw scenarios ...")
pred_draws <- matrix(NA_real_, nrow = length(cells), ncol = n_draws)
for (k in seq_len(n_draws)) {
  cont_k <- lapply(draws, function(m) m[, k])
  pred_draws[, k] <- predict_bf(cont_k)
}

bf_sampled <- rowMeans(pred_draws, na.rm = TRUE)   # E[BRT(draw)] per pixel

coal_pos <- which(in_C)
ok2 <- is.finite(obs_density) & is.finite(bf_meandraw) & is.finite(bf_sampled)
sel_C <- ok2 & in_C
nC <- sum(sel_C)

so_C <- sum(obs_density[sel_C])
sb_md_C   <- sum(bf_meandraw[sel_C])
sb_smp_C  <- sum(bf_sampled[sel_C])

cat(sprintf("coalition pixels analysed (finite all): %d / %d\n", nC, length(coal_idx)))
cat(sprintf("\n  observed                         sum=%12.0f  mean=%.4f\n",
            so_C, mean(obs_density[sel_C])))
cat(sprintf("  bf, mean draw   BRT(E[draw])     sum=%12.0f  mean=%.4f\n",
            sb_md_C, mean(bf_meandraw[sel_C])))
cat(sprintf("  bf, sampled     E[BRT(draw)]     sum=%12.0f  mean=%.4f\n",
            sb_smp_C, mean(bf_sampled[sel_C])))
cat(sprintf("\n  impact, mean draw  (sum bf - sum obs) = %+12.0f   (%+.1f%% of obs)\n",
            sb_md_C - so_C, 100 * (sb_md_C - so_C) / so_C))
cat(sprintf("  impact, sampled    (sum bf - sum obs) = %+12.0f   (%+.1f%% of obs)\n",
            sb_smp_C - so_C, 100 * (sb_smp_C - so_C) / so_C))
cat(sprintf("\n  Jensen gap  E[BRT(draw)] - BRT(E[draw]) = %+12.0f birds  (%+.1f%% of mean-draw bf)\n",
            sb_smp_C - sb_md_C, 100 * (sb_smp_C - sb_md_C) / sb_md_C))
cat(sprintf("  fraction of coalition pixels with sampled < mean-draw: %.1f%%\n",
            100 * mean(bf_sampled[sel_C] < bf_meandraw[sel_C])))

# per-pixel shape of the 100-draw prediction distribution (skew check)
row_mean_C <- rowMeans(pred_draws[sel_C, , drop = FALSE], na.rm = TRUE)
row_med_C  <- apply(pred_draws[sel_C, , drop = FALSE], 1, median, na.rm = TRUE)
cat(sprintf("\n  per-pixel 100-draw prediction distribution (coalition set):\n"))
cat(sprintf("    mean(rowMean)  = %.4f   mean(rowMedian) = %.4f\n",
            mean(row_mean_C), mean(row_med_C)))
cat(sprintf("    rowMean < rowMedian in %.1f%% of pixels  (left-skew => sampling pulls DOWN)\n",
            100 * mean(row_mean_C < row_med_C)))
cat("\nRead: if impact(sampled) is negative while impact(mean-draw) is positive,\n",
    "Jensen's inequality from BRT nonlinearity explains 12F's sign -> cause (2).\n",
    sep = "")

# ---- Verdict ----------------------------------------------------------------

cat("\n", strrep("=", 78), "\n", "VERDICT\n", strrep("=", 78), "\n", sep = "")
imp_C_only <- if (!is.null(part_summary[["C_only"]]))
  part_summary[["C_only"]]["sum_impact"] else NA_real_
cause1 <- isTRUE(imp_C_only < 0) &&
  isTRUE(part_summary[["ALL_R"]]["sum_impact"] >= 0)
cause2 <- isTRUE((sb_smp_C - so_C) < 0) && isTRUE((sb_md_C - so_C) >= 0)
cat(sprintf("  cause (1) road-pixel definition  : %s\n",
            if (cause1) "SUPPORTED (C_only negative, R positive)"
            else "not supported"))
cat(sprintf("  cause (2) mean-vs-sampled draws  : %s\n",
            if (cause2) "SUPPORTED (sampled flips negative, mean-draw positive)"
            else "not supported"))
if (!cause1 && !cause2)
  cat("  neither isolated cause flips the sign here — see the per-partition and\n",
      "  Jensen numbers above for the partial contribution of each.\n", sep = "")

# ---- Save outputs -----------------------------------------------------------

out_dir <- file.path(ia_dir, "output_tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

per_pixel <- data.frame(
  cell        = cells,
  in_R        = in_R,
  in_C        = in_C,
  partition   = partition,
  obs         = obs_density,
  bf_meandraw = bf_meandraw,
  bf_sampled  = bf_sampled,
  impact_meandraw = bf_meandraw - obs_density,
  impact_sampled  = bf_sampled  - obs_density
)
write.csv(per_pixel,
          file.path(out_dir, sprintf("13D_road_def_sampling_%s_%s_%d_boot%d.csv",
                                     species, bcr_code, year, boot_i)),
          row.names = FALSE)

saveRDS(list(per_pixel    = per_pixel,
             pred_draws_C = pred_draws[in_C, , drop = FALSE],
             coal_cells   = cells[in_C],
             part_summary = part_summary,
             test2 = list(n = nC, sum_obs = so_C,
                          sum_bf_meandraw = sb_md_C,
                          sum_bf_sampled  = sb_smp_C),
             road_idx = road_idx, coal_idx = coal_idx,
             dist_in_model = dist_in_model,
             cont_swap = cont_swap, cat_swap = cat_swap_use),
        file.path(out_dir, sprintf("13D_road_def_sampling_%s_%s_%d_boot%d.rds",
                                   species, bcr_code, year, boot_i)))

message("\nwrote results to ", out_dir)
message(Sys.time(), " | done.")
