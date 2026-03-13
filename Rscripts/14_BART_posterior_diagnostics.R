# ---
# title:   Test gbart() posterior diagnostics: normality assumption and draw adequacy
# author:  Mannfred Boehm
# purpose: Two diagnostic questions, answered by fitting gbart() (identical arguments
#          to 08B) on n_samples random (subbasin, covariate) pairs from year 2020:
#
#          (A) NORMALITY (legacy diagnostic, Gaussian sampling approach now retired):
#              Was the log1p-scale posterior actually Gaussian?  12B previously sampled
#              expm1(logmean + logsd * ε) (ε ~ N(0,1)).  Tests:
#              (a) Shapiro-Wilk test on log1p-scale draws
#              (b) containment rate: does expm1(logmean ± 1.96·logsd) contain the
#                  empirical 2.5/97.5 percentiles?
#              (c) signed relative error of the normal CI bounds
#
#          (B) DRAW ADEQUACY (new diagnostic, for current draw-based approach):
#              Does a random subsample of n_draws=50 from the 700-draw posterior
#              adequately represent the full posterior?  12B now resamples 1 of 50
#              stored draws per scenario rather than drawing from a Gaussian.
#              Tests (20 repeated subsamples of 50, to get stability):
#              (d) mean absolute relative error of q025/q50/q975 from 50 draws
#                  vs. the full 700-draw posterior (per pixel, averaged across pixels)
#              (e) containment rate: does the 50-draw 95% interval contain the
#                  700-draw 95% interval?
#              (f) SD of the above metrics across the 20 repeated subsamples
#
# context: Alliance Canada cluster, single node, multi-core
# output:  data/derived_data/rds_files/bart_posterior_diagnostics.rds
# ---

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

    post_orig <- expm1(fit$yhat.test)   # back-transformed posterior [ndpost x n_px]
    q025_orig <- apply(post_orig, 2, quantile, probs = 0.025)
    q975_orig <- apply(post_orig, 2, quantile, probs = 0.975)

    # back-transform log1p-scale CI bounds — matches what 12B samples:
    # expm1(logmean + logsd * ε), so ±1.96 bounds become expm1(logmean ± 1.96·logsd)
    norm_lo_orig <- expm1(mu_log - 1.96 * sd_log)
    norm_hi_orig <- expm1(mu_log + 1.96 * sd_log)

    # ---- (A) normality summary statistics ------------------------------------

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

    # ---- (B) draw adequacy: does n_draws=50 represent the full posterior? ----
    # Reference quantiles from all 700 draws (the "truth").
    n_draws_test <- 50L
    n_reps_da    <- 20L   # repeated subsamples to assess stability

    q025_full <- apply(fit$yhat.test, 2, quantile, probs = 0.025)
    q50_full  <- apply(fit$yhat.test, 2, quantile, probs = 0.50)
    q975_full <- apply(fit$yhat.test, 2, quantile, probs = 0.975)

    # For each rep: subsample n_draws_test rows, compute quantiles, compare to full
    rep_stats <- vapply(seq_len(n_reps_da), function(r) {
      idx50  <- sample(nrow(fit$yhat.test), n_draws_test)
      mat50  <- fit$yhat.test[idx50, , drop = FALSE]
      q025_s <- apply(mat50, 2, quantile, probs = 0.025)
      q50_s  <- apply(mat50, 2, quantile, probs = 0.50)
      q975_s <- apply(mat50, 2, quantile, probs = 0.975)

      # containment: 50-draw interval contains the full-posterior 95% interval
      contains50 <- mean(q025_s <= q025_full & q975_s >= q975_full)

      # mean absolute relative error per quantile (averaged across pixels)
      err_q025 <- mean(abs(q025_s - q025_full) / (abs(q025_full) + eps))
      err_q50  <- mean(abs(q50_s  - q50_full)  / (abs(q50_full)  + eps))
      err_q975 <- mean(abs(q975_s - q975_full) / (abs(q975_full) + eps))

      c(contains50 = contains50,
        err_q025   = err_q025,
        err_q50    = err_q50,
        err_q975   = err_q975)
    }, numeric(4))

    # mean and SD across the n_reps_da subsamples (SD reflects subsample variability)
    da_contains_mean <- mean(rep_stats["contains50", ])
    da_contains_sd   <- sd(rep_stats["contains50", ])
    da_err_q025_mean <- mean(rep_stats["err_q025",   ])
    da_err_q025_sd   <- sd(rep_stats["err_q025",   ])
    da_err_q50_mean  <- mean(rep_stats["err_q50",    ])
    da_err_q50_sd    <- sd(rep_stats["err_q50",    ])
    da_err_q975_mean <- mean(rep_stats["err_q975",   ])
    da_err_q975_sd   <- sd(rep_stats["err_q975",   ])

    list(
      sample_id          = k,
      subbasin           = sub_idx,
      covariate          = b,
      status             = "ok",
      n_train            = nrow(df_tb),
      n_test             = n_px,
      # (A) Shapiro-Wilk (log1p scale): are posterior draws normally distributed?
      sw_frac_reject_05  = mean(sw_pv < 0.05, na.rm = TRUE),
      sw_pval_median     = median(sw_pv, na.rm = TRUE),
      # (A) log1p-scale CI agreement
      contains_log       = contains_log,
      err_lo_log         = err_lo_log,
      err_hi_log         = err_hi_log,
      # (A) original-scale CI agreement (what 12B previously used)
      contains_orig      = contains_orig,
      err_lo_orig        = err_lo_orig,
      err_hi_orig        = err_hi_orig,
      # (A) posterior spread on log1p scale (context for how diffuse the posteriors are)
      mean_sd_log        = mean(sd_log),
      # (B) draw adequacy: 50-draw subsample vs. full 700-draw posterior
      # containment: fraction of pixels where 50-draw 95% interval ⊇ full 95% interval
      da_contains_mean   = da_contains_mean,
      da_contains_sd     = da_contains_sd,
      # mean absolute relative error of q025/q50/q975 (mean and SD across 20 subsamples)
      da_err_q025_mean   = da_err_q025_mean,
      da_err_q025_sd     = da_err_q025_sd,
      da_err_q50_mean    = da_err_q50_mean,
      da_err_q50_sd      = da_err_q50_sd,
      da_err_q975_mean   = da_err_q975_mean,
      da_err_q975_sd     = da_err_q975_sd
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

out_path <- file.path(ia_dir, "data", "derived_data", "rds_files", "bart_posterior_diagnostics.rds")
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
=== (A) NORMALITY (legacy Gaussian sampling — now retired) ===

--- (a) Shapiro-Wilk on log1p-scale posterior draws ---
  Fraction of pixels rejecting normality (alpha=0.05):  %.3f  [median across pairs]
  Median Shapiro-Wilk p-value:                          %.4f  [median across pairs]
  (values close to 0 / high rejection rate → non-normal posterior at log1p scale)

--- (b) CI containment: normal approx CI ⊇ empirical 95%% interval ---
  log1p scale:    %.3f  (ideal = 1.00)
  original scale: %.3f  (ideal = 1.00; expm1(logmean ± 1.96·logsd), as in old 12B)

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

  message(sprintf(
"
=== (B) DRAW ADEQUACY (50-draw subsample vs. full 700-draw posterior) ===
    Each metric is the mean across %d pairs; (SD) is SD across 20 subsamples per pair.

--- (d) Mean absolute relative error of quantiles (50 draws vs. 700 draws) ---
  q025:  %.4f  (SD %.4f)   (ideal = 0)
  q50:   %.4f  (SD %.4f)   (ideal = 0)
  q975:  %.4f  (SD %.4f)   (ideal = 0)

--- (e) Containment: 50-draw 95%% interval ⊇ full-posterior 95%% interval ---
  mean containment rate: %.3f  (SD %.4f)   (ideal = 1.00)
  (values < 1 mean some pixels' tails are underrepresented by the 50 stored draws)",
    nrow(ok),
    mean(ok$da_err_q025_mean, na.rm = TRUE), mean(ok$da_err_q025_sd, na.rm = TRUE),
    mean(ok$da_err_q50_mean,  na.rm = TRUE), mean(ok$da_err_q50_sd,  na.rm = TRUE),
    mean(ok$da_err_q975_mean, na.rm = TRUE), mean(ok$da_err_q975_sd, na.rm = TRUE),
    mean(ok$da_contains_mean, na.rm = TRUE), mean(ok$da_contains_sd, na.rm = TRUE)
  ))
}
