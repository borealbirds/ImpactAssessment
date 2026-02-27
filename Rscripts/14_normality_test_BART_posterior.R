# ---
# title:   Test normality assumption of gbart() posterior for backfilled covariates
# author:  Mannfred Boehm
# purpose: 12B_predict_species_bcr.R constructs counterfactual covariate stacks as
#            mu_stack[[v]] + z * sd_stack[[v]]    (z = 0, ±1.96)
#          where mu and sd are derived from the gbart() posterior via the log-normal
#          formulas in 08B_deploy_gbart.R.  This is only accurate if the ±1.96*SD
#          interval on the original scale closely matches the empirical 2.5/97.5th
#          percentiles of the back-transformed posterior draws.
#
#          This script tests that assumption by fitting gbart() (with identical
#          arguments to 08B) on n_samples random (subbasin, covariate) pairs from
#          year 2020, then comparing:
#            (a) Shapiro-Wilk test on log1p-scale posterior draws (is the log-scale
#                posterior actually normal?)
#            (b) containment rate: does b_mean ± 1.96*b_sd contain the empirical
#                95% posterior interval?
#            (c) signed relative error of the normal CI bounds vs. empirical quantiles
#
# context: Alliance Canada cluster, single node, multi-core
# output:  data/derived_data/rds_files/normality_test_results.rds
# ---

# --- suggested SLURM header ---------------------------------------------------
# #!/bin/bash
# #SBATCH --job-name=normtest_gbart
# #SBATCH --time=06:00:00
# #SBATCH --mem=64G
# #SBATCH --cpus-per-task=16
# #SBATCH --output=logs/normtest_%j.log
# module load r/4.4.0
# Rscript Rscripts/normality_test_gbart_posterior.R
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(BART)
  library(BAMexploreR)
  library(terra)
  library(tidyverse)
  library(parallel)
})

# ---- settings ----------------------------------------------------------------

ia_dir    <- "/home/mannfred/scratch/impact_assessment"
year      <- 2020
n_samples <- 500L   # number of random (subbasin, covariate) pairs to test
n_sw_px   <- 30L    # pixels per run on which to apply Shapiro-Wilk
set.seed(42)

n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1"))
message(sprintf("[%s] starting — %d core(s), %d pairs", Sys.time(), n_cores, n_samples))

# ---- file paths (passed by value to workers; avoids terra fork issues) -------

stack_path     <- file.path(ia_dir, "data", "raw_data", "covariates_mosaiced",
                            sprintf("covariates_mosaiced_%d.tif", year))
lowhf_path     <- file.path(ia_dir, "data", "raw_data", "hirshpearson", "CanHF_1km_lessthan1.tif")
highhf_path    <- file.path(ia_dir, "data", "raw_data", "hirshpearson", "CanHF_1km_morethan1.tif")
subbasins_path <- file.path(ia_dir, "data", "raw_data", "hydrobasins_masked_merged_subset.gpkg")
hier_path      <- file.path(ia_dir, "data", "raw_data", "biotic_variable_hierarchy.rds")

# ---- predictor metadata (mirrors 07_train_and_backfill.R) -------------------

categorical_responses <- c("ABoVE_1km", "NLCD_1km", "MODISLCC_1km",
                           "MODISLCC_5x5", "SCANFI_1km", "VLCE_1km")

predictor_metadata <-
  dplyr::tibble(BAMexploreR::predictor_metadata) |>
  dplyr::filter(version == "v5") |>
  dplyr::select(predictor, definition, predictor_class) |>
  dplyr::mutate(dplyr::across("predictor", stringr::str_replace, "Year",   "year")) |>
  dplyr::mutate(dplyr::across("predictor", stringr::str_replace, "Method", "method"))

soil_covs <- tibble::tibble(
  predictor = c(
    "cec_0-5cm_mean_1000",   "cec_100-200cm_mean_1000", "cec_15-30cm_mean_1000",
    "cec_30-60cm_mean_1000", "cec_5-15cm_mean_1000",    "cec_60-100cm_mean_1000",
    "soc_0-5cm_mean_1000",   "soc_100-200cm_mean_1000", "soc_15-30cm_mean_1000",
    "soc_30-60cm_mean_1000", "soc_5-15cm_mean_1000",    "soc_60-100cm_mean_1000"),
  predictor_class = rep("Soil Properties", 12))

actually_biotic_what <- c("Peatland_5x5", "Peatland_1km")
actually_biotic_df   <- tibble::tibble(predictor      = actually_biotic_what,
                                       predictor_class = c("Wetland", "Wetland"))

abiotic_vars <-
  predictor_metadata |>
  dplyr::filter(predictor_class %in% c("Annual Climate", "Climate Normals",
                                        "Topography", "Wetland", "Disturbance")) |>
  tibble::add_row(predictor = "CAfire", predictor_class = "Time Since Disturbance") |>
  dplyr::filter(!(predictor %in% actually_biotic_what)) |>
  dplyr::bind_rows(soil_covs)

biotic_vars <-
  predictor_metadata |>
  dplyr::filter(!(predictor_class %in% c(abiotic_vars$predictor_class, "Time", "Method"))) |>
  dplyr::bind_rows(actually_biotic_df)

neworder    <- readRDS(hier_path)
biotic_vars <- biotic_vars[match(neworder, biotic_vars$predictor), ]

# continuous biotic covariates present in the covariate stack
stack_layers      <- names(terra::rast(stack_path))
biotic_continuous <- intersect(setdiff(neworder, categorical_responses), stack_layers)
message(sprintf("[%s] %d continuous biotic covariates in stack", Sys.time(), length(biotic_continuous)))

# ---- sample n_samples random (subbasin, covariate) pairs --------------------

n_subbasins <- nrow(terra::vect(subbasins_path))
pairs <- data.frame(
  sample_id = seq_len(n_samples), #n=500
  subbasin  = sample(seq_len(n_subbasins), n_samples, replace = TRUE),
  covariate = sample(biotic_continuous,    n_samples, replace = TRUE),
  stringsAsFactors = FALSE)

message(sprintf("[%s] sampled %d pairs; dispatching to workers", Sys.time(), n_samples))

# ---- per-pair worker function ------------------------------------------------
# Everything terra-related is loaded from disk inside this function so that
# terra's C++ external pointers are valid within the forked worker process.

test_one_pair <- function(pair,
                          stack_path, lowhf_path, highhf_path, subbasins_path,
                          abiotic_vars, biotic_vars, neworder,
                          categorical_responses, n_sw_px) {

  suppressPackageStartupMessages({ library(BART); library(terra); library(dplyr) })

  k       <- pair$sample_id
  sub_idx <- pair$subbasin
  b       <- pair$covariate

  stub <- list(sample_id = k, subbasin = sub_idx, covariate = b, status = NA_character_)

  tryCatch({

    # load and align spatial objects ------------------------------------------
    stack_y   <- terra::rast(stack_path)
    subbasins <- terra::vect(subbasins_path) |> terra::project(terra::crs(stack_y))
    lowhf     <- terra::rast(lowhf_path)
    highhf    <- terra::rast(highhf_path)

    poly  <- subbasins[sub_idx]
    cov_s <- terra::crop(stack_y, poly) |> terra::mask(poly)

    # resample masks to the subbasin grid (mirrors 08A_train_and_backfill_subbasin_s.R)
    lowhf_s  <- terra::resample(lowhf,  cov_s, method = "near") |> terra::mask(poly)
    highhf_s <- terra::resample(highhf, cov_s, method = "near") |> terra::mask(poly)

    # build training and backfill data frames ---------------------------------
    df_full     <- terra::as.data.frame(cov_s, xy = TRUE, na.rm = FALSE, cells = TRUE)
    df_train    <- terra::as.data.frame(terra::mask(cov_s, lowhf_s), xy = TRUE, na.rm = FALSE)
    backfill_idx <- which(terra::values(highhf_s) == 1)
    df_backfill  <- df_full[backfill_idx, , drop = FALSE]

    if (nrow(df_backfill) == 0) return(c(stub, list(status = "no_highhf_pixels")))

    # coordinate scaling (mirrors 08A) ----------------------------------------
    xy_ok     <- which(!is.na(df_full$x) & !is.na(df_full$y))
    xy_center <- colMeans(df_full[xy_ok, c("x","y")], na.rm = TRUE)
    xy_scale  <- apply(df_full[xy_ok, c("x","y")], 2, sd, na.rm = TRUE)

    df_train[, c("x","y")]    <- sweep(df_train[, c("x","y")],    2, xy_center, "-")
    df_train[, c("x","y")]    <- sweep(df_train[, c("x","y")],    2, xy_scale,  "/")
    df_backfill[, c("x","y")] <- sweep(df_backfill[, c("x","y")], 2, xy_center, "-")
    df_backfill[, c("x","y")] <- sweep(df_backfill[, c("x","y")], 2, xy_scale,  "/")

    # predictor selection (mirrors 08A) ----------------------------------------
    abiotic_cols <- intersect(names(df_train), abiotic_vars$predictor)
    biotic_cols  <- intersect(names(df_train), biotic_vars$predictor)
    biotic_cols  <- na.omit(biotic_cols[match(neworder, biotic_cols)])
    biotic_cols_cont_local <- setdiff(biotic_cols, categorical_responses)

    if (!(b %in% names(df_train)))                      return(c(stub, list(status = "covariate_absent")))
    idx <- which(!is.na(df_train[[b]]))
    if (length(idx) < 10)                               return(c(stub, list(status = "too_few_train")))
    if (length(unique(df_train[[b]][idx])) < 2)         return(c(stub, list(status = "invariant")))

    b_pos    <- match(b, biotic_cols_cont_local)
    b_before <- if (!is.na(b_pos) && b_pos > 1L)
                  biotic_cols_cont_local[seq_len(b_pos - 1L)] else character(0)
    predictors <- c(abiotic_cols, b_before, "x", "y")

    # build BART design matrices (mirrors 08A) --------------------------------
    df_tb <- df_train[idx, predictors, drop = FALSE]
    df_tb <- df_tb[, colSums(!is.na(df_tb)) > 0, drop = FALSE]
    sds   <- sapply(df_tb, function(x) sd(as.numeric(x), na.rm = TRUE))
    df_tb <- df_tb[, !is.na(sds) & sds > 0, drop = FALSE]

    df_bb <- df_backfill[, intersect(colnames(df_tb), names(df_backfill)), drop = FALSE]
    df_bb <- df_bb[, colSums(!is.na(df_bb)) > 0, drop = FALSE]
    df_tb <- df_tb[, colnames(df_bb), drop = FALSE]

    y <- as.numeric(df_train[[b]][idx])
    if (length(unique(log1p(y))) <= 2) return(c(stub, list(status = "degenerate_logy")))

    # 90/10 holdout split (mirrors 08B_deploy_gbart.R) ------------------------
    set.seed(abs(as.integer(sprintf("%d%03d", sub_idx, which(biotic_cols == b)))))
    ho    <- sample(seq_along(y), size = round(0.10 * length(y)))
    df_tb <- df_tb[-ho, , drop = FALSE]
    y_tr  <- y[-ho]

    if (length(unique(log1p(y_tr))) <= 2) return(c(stub, list(status = "degenerate_after_split")))

    # fit gbart — identical arguments to 08B_deploy_gbart.R ------------------
    fit <- BART::gbart(
      x.train = as.matrix(df_tb),
      y.train = log1p(y_tr),
      x.test  = as.matrix(df_bb),
      type    = "wbart",
      k       = 3,
      ntree   = 50L,
      ndpost  = 700L,
      nskip   = 300L,
      sparse  = TRUE,
      sigest  = sd(log1p(y_tr)),
      sigdf   = 3
    )

    # yhat.test is [ndpost x n_px]; each column is one backfill pixel
    n_px <- ncol(fit$yhat.test)

    # ---- log1p-scale diagnostics --------------------------------------------
    # The posterior at the log1p scale is what BART models directly.
    # If it is normal here, the log-normal back-transformation is valid.

    mu_log  <- fit$yhat.test.mean               # posterior mean    [n_px]
    var_log <- apply(fit$yhat.test, 2, var)     # posterior variance [n_px]
    sd_log  <- sqrt(var_log)

    # empirical 2.5 / 97.5 quantiles on log1p scale
    q025_log <- apply(fit$yhat.test, 2, quantile, probs = 0.025)
    q975_log <- apply(fit$yhat.test, 2, quantile, probs = 0.975)

    # normal-approximation CI on log1p scale  (mu ± 1.96*sigma)
    norm_lo_log <- mu_log - 1.96 * sd_log
    norm_hi_log <- mu_log + 1.96 * sd_log

    # Shapiro-Wilk on a random sample of pixels (n must be 3–5000; ndpost=700 is fine)
    sw_j  <- sample(seq_len(n_px), min(n_sw_px, n_px))
    sw_pv <- vapply(sw_j, function(j) {
      x <- fit$yhat.test[, j]
      if (length(unique(x)) < 3L) return(NA_real_)
      tryCatch(shapiro.test(x)$p.value, error = function(e) NA_real_)
    }, numeric(1))

    # ---- original-scale diagnostics ----------------------------------------
    # These mirror the exact back-transformation in 08B_deploy_gbart.R and
    # the CI structure used in 12B's make_counterfactual_stack().

    b_mean <- expm1(mu_log)
    b_sd   <- sqrt(expm1(var_log) * exp(2 * mu_log + var_log))   # log-normal SD

    post_orig <- expm1(fit$yhat.test)   # back-transformed posterior [ndpost x n_px]
    q025_orig <- apply(post_orig, 2, quantile, probs = 0.025)
    q975_orig <- apply(post_orig, 2, quantile, probs = 0.975)

    # normal-approximation CI on original scale  (b_mean ± 1.96 * b_sd)
    # this is what make_counterfactual_stack() produces with z = ±1.96
    norm_lo_orig <- b_mean - 1.96 * b_sd
    norm_hi_orig <- b_mean + 1.96 * b_sd

    # ---- summary statistics per run -----------------------------------------

    # containment: does the normal CI fully contain the empirical 95% interval?
    contains_log  <- mean(norm_lo_log  <= q025_log  & norm_hi_log  >= q975_log)
    contains_orig <- mean(norm_lo_orig <= q025_orig & norm_hi_orig >= q975_orig)

    # signed relative error: (norm_bound - empirical_bound) / |empirical_bound|
    # positive  → normal CI is wider/shifted (conservative)
    # negative  → normal CI is narrower/shifted (anti-conservative)
    eps <- 1e-6
    err_lo_log  <- mean((norm_lo_log  - q025_log)  / (abs(q025_log)  + eps))
    err_hi_log  <- mean((norm_hi_log  - q975_log)  / (abs(q975_log)  + eps))
    err_lo_orig <- mean((norm_lo_orig - q025_orig) / (abs(q025_orig) + eps))
    err_hi_orig <- mean((norm_hi_orig - q975_orig) / (abs(q975_orig) + eps))

    list(
      sample_id          = k,
      subbasin           = sub_idx,
      covariate          = b,
      status             = "ok",
      n_train            = nrow(df_tb),
      n_test             = n_px,
      # Shapiro-Wilk (log1p scale): are posterior draws normally distributed?
      sw_frac_reject_05  = mean(sw_pv < 0.05, na.rm = TRUE),
      sw_pval_median     = median(sw_pv, na.rm = TRUE),
      # log1p-scale CI agreement
      contains_log       = contains_log,
      err_lo_log         = err_lo_log,
      err_hi_log         = err_hi_log,
      # original-scale CI agreement (what 12B actually uses)
      contains_orig      = contains_orig,
      err_lo_orig        = err_lo_orig,
      err_hi_orig        = err_hi_orig,
      # posterior spread on log1p scale (context for how diffuse the posteriors are)
      mean_sd_log        = mean(sd_log)
    )

  }, error = function(e) c(stub, list(status = paste0("error: ", conditionMessage(e)))))
}

# ---- dispatch ----------------------------------------------------------------

results_raw <- parallel::mclapply(
  X                     = split(pairs, pairs$sample_id),
  FUN                   = test_one_pair,
  stack_path            = stack_path,
  lowhf_path            = lowhf_path,
  highhf_path           = highhf_path,
  subbasins_path        = subbasins_path,
  abiotic_vars          = abiotic_vars,
  biotic_vars           = biotic_vars,
  neworder              = neworder,
  categorical_responses = categorical_responses,
  n_sw_px               = n_sw_px,
  mc.cores              = n_cores,
  mc.preschedule        = FALSE   # dynamic scheduling; faster workers pick up the next task
)

results <- dplyr::bind_rows(lapply(results_raw, as.data.frame))

# ---- save --------------------------------------------------------------------

out_path <- file.path(ia_dir, "data", "derived_data", "rds_files", "normality_test_results.rds")
saveRDS(results, out_path)
message(sprintf("[%s] saved %d rows → %s", Sys.time(), nrow(results), out_path))

# ---- console summary ---------------------------------------------------------

ok   <- dplyr::filter(results, status == "ok")
skip <- dplyr::filter(results, status != "ok")

message(sprintf(
  "\n=== Results ===\nok: %d | skipped/errored: %d (of %d requested)\n",
  nrow(ok), nrow(skip), n_samples))

if (nrow(skip) > 0) {
  message("Skip/error breakdown:")
  print(dplyr::count(skip, status))
}

if (nrow(ok) > 0) {
  message(sprintf(
"
--- (a) Shapiro-Wilk on log1p-scale posterior draws ---
  Fraction of pixels rejecting normality (alpha=0.05):  %.3f  [median across pairs]
  Median Shapiro-Wilk p-value:                          %.4f  [median across pairs]
  (values close to 0 / high rejection rate → non-normal posterior at log1p scale)

--- (b) CI containment: normal approx CI ⊇ empirical 95%% interval ---
  log1p scale:    %.3f  (ideal = 1.00)
  original scale: %.3f  (ideal = 1.00; this is what 12B uses)

--- (c) Signed relative error  (norm_bound − empirical_bound) / |empirical_bound| ---
  log1p  lower bound (ideal ≈ 0): %+.4f   upper: %+.4f
  orig   lower bound (ideal ≈ 0): %+.4f   upper: %+.4f
  (positive → normal CI overshoots empirical; negative → undershoots)",
    median(ok$sw_frac_reject_05, na.rm = TRUE),
    median(ok$sw_pval_median,    na.rm = TRUE),
    mean(ok$contains_log,  na.rm = TRUE),
    mean(ok$contains_orig, na.rm = TRUE),
    mean(ok$err_lo_log,  na.rm = TRUE), mean(ok$err_hi_log,  na.rm = TRUE),
    mean(ok$err_lo_orig, na.rm = TRUE), mean(ok$err_hi_orig, na.rm = TRUE)
  ))
}
