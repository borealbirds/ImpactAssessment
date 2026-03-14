# ---
# title:   BART posterior draw adequacy diagnostic
# author:  Mannfred Boehm
# purpose: Test whether n_draws=100 random draws subsampled from the 700-draw
#          gbart() posterior adequately represent the full posterior.
#          12B resamples 1 of 100 stored draws per counterfactual scenario; this
#          script verifies that 100 draws is enough to cover the posterior's shape
#          and tails regardless of their distribution.
#
#          For each of n_samples random (subbasin, covariate) pairs, fits gbart()
#          (identical arguments to 08B), then repeats 20 subsamples of 100 draws
#          and compares their quantiles to the full 700-draw posterior:
#            (a) mean absolute relative error of q025/q50/q975
#                (100-draw subsample vs. full 700-draw posterior, per pixel)
#            (b) SD of the above across the 20 repeated subsamples
#                (reflects how stable a single stored set of 100 draws is)
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
  sample_id = seq_len(n_samples),
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
                          categorical_responses) {

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

    # ---- draw adequacy: does n_draws=50 represent the full posterior? --------
    # Reference quantiles from all 700 draws (the "truth").
    n_draws_test <- 100L
    n_reps_da    <- 20L   # repeated subsamples to assess stability
    eps          <- 1e-6

    q025_full <- apply(fit$yhat.test, 2, quantile, probs = 0.025)
    q50_full  <- apply(fit$yhat.test, 2, quantile, probs = 0.50)
    q975_full <- apply(fit$yhat.test, 2, quantile, probs = 0.975)

    # For each rep: subsample n_draws_test rows, compute quantiles, compare to full
    rep_stats <- vapply(seq_len(n_reps_da), function(r) {
      idx100 <- sample(nrow(fit$yhat.test), n_draws_test)
      mat100 <- fit$yhat.test[idx100, , drop = FALSE]
      q025_s <- apply(mat100, 2, quantile, probs = 0.025)
      q50_s  <- apply(mat100, 2, quantile, probs = 0.50)
      q975_s <- apply(mat100, 2, quantile, probs = 0.975)

      # mean absolute relative error per quantile (averaged across pixels)
      err_q025 <- mean(abs(q025_s - q025_full) / (abs(q025_full) + eps))
      err_q50  <- mean(abs(q50_s  - q50_full)  / (abs(q50_full)  + eps))
      err_q975 <- mean(abs(q975_s - q975_full) / (abs(q975_full) + eps))

      c(err_q025 = err_q025,
        err_q50  = err_q50,
        err_q975 = err_q975)
    }, numeric(3))

    # mean and SD across the n_reps_da subsamples (SD reflects subsample variability)
    da_err_q025_mean <- mean(rep_stats["err_q025", ])
    da_err_q025_sd   <- sd(rep_stats["err_q025",   ])
    da_err_q50_mean  <- mean(rep_stats["err_q50",  ])
    da_err_q50_sd    <- sd(rep_stats["err_q50",    ])
    da_err_q975_mean <- mean(rep_stats["err_q975", ])
    da_err_q975_sd   <- sd(rep_stats["err_q975",   ])

    list(
      sample_id        = k,
      subbasin         = sub_idx,
      covariate        = b,
      status           = "ok",
      n_train          = nrow(df_tb),
      n_test           = n_px,
      # mean absolute relative error of q025/q50/q975 (mean and SD across 20 subsamples)
      da_err_q025_mean = da_err_q025_mean,
      da_err_q025_sd   = da_err_q025_sd,
      da_err_q50_mean  = da_err_q50_mean,
      da_err_q50_sd    = da_err_q50_sd,
      da_err_q975_mean = da_err_q975_mean,
      da_err_q975_sd   = da_err_q975_sd
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
=== DRAW ADEQUACY (100-draw subsample vs. full 700-draw posterior) ===
    Each metric is the mean across %d pairs; (SD) is SD across 20 subsamples per pair.

--- Mean absolute relative error of quantiles (100 draws vs. 700 draws) ---
  q025:  %.4f  (SD %.4f)   (ideal = 0)
  q50:   %.4f  (SD %.4f)   (ideal = 0)
  q975:  %.4f  (SD %.4f)   (ideal = 0)",
    nrow(ok),
    mean(ok$da_err_q025_mean, na.rm = TRUE), mean(ok$da_err_q025_sd, na.rm = TRUE),
    mean(ok$da_err_q50_mean,  na.rm = TRUE), mean(ok$da_err_q50_sd,  na.rm = TRUE),
    mean(ok$da_err_q975_mean, na.rm = TRUE), mean(ok$da_err_q975_sd, na.rm = TRUE)
  ))
}
