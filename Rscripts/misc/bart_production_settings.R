# ---
# title:   Are 08B's gbart() settings good enough? Burn-in, mixing, trees, draws, validation
# author:  Mannfred Boehm
# created: November 4, 2026
# ---

# bart_vs_xgboost_benchmark.R found slow mixing on one subbasin with 5 predictors. This repeats the
# check on real production fits: pairs of (subbasin, continuous covariate) that 07 fitted, with
# production's full predictor set and training pixels, sampled across production's range of
# training sizes and weighted by the covariate's importance in the bird models.
#
# Per pair, as 08A and 08B build it (prep_bart_design() is sourced from 08A):
#   chains    n_chains chains of nskip + n_kept sweeps at each ntree, on production's 90% split.
#             x.test = a sample of the real backfill pixels + the 10% holdout, so one set of
#             chains gives burn-in, R-hat, effective draws, every thinning and multi-chain option
#             (they are windows of these chains) and held-out skill.
#   psis      PSIS-LOO from one chain's first 700 kept draws: LOO without refitting.
#   footprint the subbasin's real footprint, shifted onto its training land, held out. This is the
#             geometry the backfill faces; a random holdout or LOO leaves every neighbour in.
#             Scored with observed earlier covariates, as 08B's holdout is, and with cascaded
#             ones, as a real backfill is.
#   size      the same holdout scored by a model on 90% of the 90%, for what 10% of the pixels buys.
#
# One deviation: the earlier covariates at backfill pixels come from a quick cascade (100 + 100
# sweeps, no blank model) rather than production's. They are x.test inputs only.
#
# Results go to work_dir, one file per fit, and a rerun skips what exists. Stops at DEADLINE.
# Summarised by bart_production_settings_summary.R.
#
# Run from the project root, one job at a time:
#   BPS_WORK=<dir> DEADLINE="2026-11-04 15:20:00" Rscript Rscripts/misc/bart_production_settings.R

suppressPackageStartupMessages({
  library(BART)
  library(terra)
  library(dplyr)
})

# 1. settings --------------------------------------------------------------------------------

ia_dir   <- getwd()
year     <- 2020L
work_dir <- Sys.getenv("BPS_WORK")
deadline <- as.POSIXct(Sys.getenv("DEADLINE", "2099-01-01 00:00:00"))
if (work_dir == "") stop("set BPS_WORK")
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)

n_chains  <- as.integer(Sys.getenv("N_CHAINS", "4"))
n_kept    <- as.integer(Sys.getenv("N_KEPT", "1700"))
arms      <- as.integer(strsplit(Sys.getenv("ARMS", "50,200,100"), ",")[[1]])
n_pairs   <- as.integer(Sys.getenv("N_PAIRS", "16"))
only_pair <- as.integer(Sys.getenv("ONLY_PAIR", NA))
nskip     <- 300L           # 08B
ndpost    <- 700L           # 08B
n_cap     <- 13000L         # gbart() holds about 5 copies of [kept x n]; 5.4 GB free here
n_mid_max <- 4000L          # the 100-tree arm runs only below this n, for time
n_hf      <- 200L           # backfill pixels in x.test
n_test    <- 400L           # holdout pixels scored
seed0     <- 20261104L

terraOptions(memfrac = 0.3, progress = 0)
say <- function(...) message(format(Sys.time(), "[%H:%M:%S] "), sprintf(...))
past_deadline <- function() Sys.time() > deadline
logp <- function(...) invisible(NULL)

source(file.path(ia_dir, "Rscripts", "08A_train_and_backfill_subbasin_s.R"))  # prep_bart_design()


# 2. 07's predictor definitions --------------------------------------------------------------

meta <- dplyr::filter(BAMexploreR::predictor_metadata, version == "v5") |>
  dplyr::select(predictor, predictor_class) |>
  dplyr::mutate(predictor = sub("Method", "method", sub("Year", "year", predictor)))
soil <- paste0(rep(c("cec", "soc"), each = 6), "_",
               c("0-5", "100-200", "15-30", "30-60", "5-15", "60-100"), "cm_mean_1000")
wetland_biotic <- c("Peatland_5x5", "Peatland_1km")
abiotic <- meta |>
  dplyr::filter(predictor_class %in% c("Annual Climate", "Climate Normals", "Topography",
                                       "Wetland", "Disturbance")) |>
  tibble::add_row(predictor = "CAfire", predictor_class = "Time Since Disturbance") |>
  dplyr::filter(!(predictor %in% wetland_biotic)) |>
  dplyr::bind_rows(tibble::tibble(predictor = soil, predictor_class = "Soil Properties"))
biotic <- c(setdiff(meta$predictor[!(meta$predictor_class %in% c(abiotic$predictor_class, "Time", "Method"))],
                    character(0)), wetland_biotic)
neworder <- readRDS(file.path(ia_dir, "data", "raw_data", "biotic_variable_hierarchy.rds"))
categorical <- c("ABoVE_1km", "NLCD_1km", "MODISLCC_1km", "MODISLCC_5x5", "SCANFI_1km", "VLCE_1km")


# 3. the pairs -------------------------------------------------------------------------------

# Production's pairs, n below n_cap, by training-size quartile of all pairs, one pair per
# subbasin, covariate weighted by its summed relative influence in V5's bird models. Strata take
# turns, so a run stopped early still spans the range.
sample_pairs <- function() {
  metrics <- as.data.frame(readRDS(file.path(ia_dir, "data", "derived_data", "rds_files", "model_metrics",
                                             year, "continuous_train_metrics.rds")))
  influence <- readRDS(file.path(ia_dir, "data", "derived_data", "rds_files", "bam_predictor_importance_v5.rds")) |>
    dplyr::group_by(predictor) |>
    dplyr::summarise(influence = sum(mean_rel_inf), .groups = "drop")
  breaks <- quantile(metrics$ntrain, c(0, .25, .5, .75, 1))
  breaks[5] <- n_cap
  candidates <- metrics |>
    dplyr::inner_join(influence, by = c("covariate" = "predictor")) |>
    dplyr::filter(ntrain <= n_cap) |>
    dplyr::mutate(stratum = cut(ntrain, breaks, labels = FALSE, include.lowest = TRUE))
  set.seed(seed0)
  picked <- list()
  used <- integer(0)
  for (i in seq_len(n_pairs)) {
    pool <- dplyr::filter(candidates, stratum == c(2L, 3L, 1L, 4L)[(i - 1L) %% 4L + 1L], !(subbasin %in% used))
    row <- pool[sample(nrow(pool), 1, prob = pool$influence), ]
    used <- c(used, row$subbasin)
    picked[[i]] <- row
  }
  dplyr::bind_rows(picked) |>
    dplyr::transmute(pair = dplyr::row_number(), subbasin, covariate, stratum, ntrain, r2_train = r2_median)
}

pairs <- sample_pairs()
print(as.data.frame(pairs))
write.csv(pairs, file.path(work_dir, "pairs.csv"), row.names = FALSE)


# 4. fitting ---------------------------------------------------------------------------------

# 08B's call
fit_gbart <- function(x, y, x_test, ntree, nskip, ndpost, seed) {
  set.seed(seed)
  t0 <- proc.time()[[3]]
  capture.output(
    fit <- suppressWarnings(
      BART::gbart(x.train = x, y.train = y, x.test = x_test, type = "wbart", k = 3, ntree = ntree,
                  ndpost = ndpost, nskip = nskip, sparse = TRUE, sigest = sd(y), sigdf = 3,
                  printevery = 1000000L))
  )
  fit$seconds <- proc.time()[[3]] - t0
  fit
}

# 08A's frames for one subbasin, footprint covariates at 0 on the backfill pixels. Keeps the
# metres grid position of every pixel (x_m, y_m) before 08A's coordinate scaling.
subbasin_frames <- function(s) {
  stack_y   <- rast(file.path(ia_dir, "data", "raw_data", "covariates_mosaiced",
                              sprintf("covariates_mosaiced_%d.tif", year)))
  subbasins <- project(vect(file.path(ia_dir, "data", "raw_data", "hydrobasins_masked_merged_subset.gpkg")),
                       crs(stack_y))
  poly  <- subbasins[s]
  cov_s <- mask(crop(stack_y, poly), poly)
  hf <- function(kind) {
    r <- rast(file.path(ia_dir, "data", "raw_data", "hirshpearson", sprintf("CanHF_1km_%s.tif", kind)))
    mask(project(crop(r, project(poly, crs(r))), cov_s, method = "near"), poly)
  }
  low  <- hf("lessthan1")
  high <- hf("morethan1")
  train_cells <- which(values(low, mat = FALSE) == 1)
  fill_cells  <- which(values(high, mat = FALSE) == 1)
  df <- as.data.frame(cov_s, xy = TRUE, na.rm = FALSE)
  stopifnot(nrow(df) == ncell(cov_s))
  xy_ok     <- which(!is.na(df$x) & !is.na(df$y))
  xy_center <- colMeans(df[xy_ok, c("x", "y")])
  xy_scale  <- apply(df[xy_ok, c("x", "y")], 2, sd)
  df$x_m <- df$x
  df$y_m <- df$y
  df[, c("x", "y")] <- sweep(sweep(df[, c("x", "y")], 2, xy_center, "-"), 2, xy_scale, "/")
  list(train = df[train_cells, ], fill = df[fill_cells, ], res = res(cov_s)[1])
}

# Training pixels under the footprint shifted by whole cells: offsets until 10% of the training
# pixels are covered, keeping any offset that holds the total to 30%. Scattered shapes stay
# footprint-shaped, however the copies land.
footprint_holdout <- function(train, fill, cell, seed) {
  set.seed(seed)
  key <- function(x, y) paste(round(x / cell), round(y / cell))
  train_key <- key(train$x_m, train$y_m)
  span_x <- diff(range(train$x_m))
  span_y <- diff(range(train$y_m))
  held <- rep(FALSE, nrow(train))
  for (try in seq_len(2000L)) {
    dx <- round(runif(1, -span_x, span_x) / cell) * cell
    dy <- round(runif(1, -span_y, span_y) / cell) * cell
    hit <- train_key %in% key(fill$x_m + dx, fill$y_m + dy)
    if (mean(held | hit) <= 0.30) held <- held | hit
    if (mean(held) >= 0.10) break
  }
  which(held)
}

# The earlier covariates in the hierarchy, predicted in order at `target` rows from `rows` of
# train, each model seeing the ones before it as the backfill does (08A's cascade, shortened).
quick_cascade <- function(train, target, b_before, abiotic_cols, seed) {
  for (j in seq_along(b_before)) {
    v <- b_before[j]
    idx <- which(!is.na(train[[v]]))
    if (length(idx) < 10L || length(unique(train[[v]][idx])) < 2L) {
      target[[v]] <- if (length(idx)) mean(train[[v]][idx]) else NA_real_
      next
    }
    design <- prep_bart_design(train[idx, c(abiotic_cols, b_before[seq_len(j - 1L)], "x", "y"), drop = FALSE],
                               target, logp)
    y <- log1p(train[[v]][idx])
    if (length(unique(y)) <= 2L) {
      target[[v]] <- mean(train[[v]][idx])
      next
    }
    fit <- fit_gbart(as.matrix(design$train), y, as.matrix(design$backfill), 50L, 100L, 100L, seed + j)
    target[[v]] <- colMeans(expm1(fit$yhat.test))
  }
  target
}

# PSIS-LOO of the first ndpost kept draws. Predictions on the original scale, as 08B scores them.
psis_loo <- function(fit, y_log, seconds_fit) {
  t0 <- proc.time()[[3]]
  keep <- seq_len(min(ndpost, nrow(fit$yhat.train)))
  f <- fit$yhat.train[keep, , drop = FALSE]
  sigma <- fit$sigma[nskip + keep]
  ll <- dnorm(matrix(y_log, length(keep), length(y_log), byrow = TRUE), f,
              matrix(sigma, length(keep), length(y_log)), log = TRUE)
  result <- suppressWarnings(loo::loo(ll, r_eff = rep(1, ncol(ll)), save_psis = TRUE, cores = 1))
  w <- weights(result$psis_object, log = FALSE, normalize = TRUE)
  list(pred_loo = colSums(w * expm1(f)), pred_insample = colMeans(expm1(f)),
       y = expm1(y_log), pareto_k = result$diagnostics$pareto_k,
       elpd = result$estimates["elpd_loo", ], seconds = proc.time()[[3]] - t0,
       cpo_log = if (!is.null(fit$CPO)) fit$CPO else NULL)
}

run_pair <- function(p) {
  id <- sprintf("pair_%02d", p$pair)
  path <- function(...) file.path(work_dir, paste0(id, "_", ..., ".rds"))
  do <- function(file, expr) {
    if (file.exists(file)) return(invisible(readRDS(file)))
    if (past_deadline()) stop("deadline")
    value <- expr
    saveRDS(value, file, compress = FALSE)
    invisible(value)
  }

  say("%s: subbasin %d, %s, production n = %d", id, p$subbasin, p$covariate, p$ntrain)
  frames <- subbasin_frames(p$subbasin)
  train  <- frames$train
  fill   <- frames$fill
  b      <- p$covariate
  footprint_cols <- intersect(names(fill), abiotic$predictor[abiotic$predictor_class == "Disturbance"])
  fill[, footprint_cols] <- 0

  abiotic_cols <- intersect(names(train), abiotic$predictor)
  biotic_cols  <- intersect(names(train), biotic)
  biotic_cols  <- as.character(na.omit(biotic_cols[match(neworder, biotic_cols)]))
  cont_cols    <- setdiff(biotic_cols, categorical)
  b_before     <- cont_cols[seq_len(match(b, cont_cols) - 1L)]
  predictors   <- c(abiotic_cols, b_before, "x", "y")

  idx <- which(!is.na(train[[b]]))
  set.seed(seed0 + p$pair)
  fill_sample <- fill[sample(nrow(fill), min(n_hf, nrow(fill))), , drop = FALSE]

  # held out in the shape of the footprint, then the quick cascade on what is left, predicting
  # both the backfill pixels and the held-out ones
  fp_held <- footprint_holdout(train[idx, ], fill, frames$res, seed0 + p$pair)
  fp_test <- if (length(fp_held) > n_test) sample(fp_held, n_test) else fp_held
  rest <- train[idx[-fp_held], ]
  cascaded <- do(path("cascade"), {
    t0 <- proc.time()[[3]]
    target <- rbind(fill_sample[, names(train)], train[idx[fp_test], ])
    target[, intersect(footprint_cols, names(target))] <- 0
    obs_fill <- fill_sample[[b]]
    out <- quick_cascade(rest, target, b_before, abiotic_cols, seed0 + 1000L * p$pair)
    list(target = out, obs_fill = obs_fill, seconds = proc.time()[[3]] - t0)
  })
  fill_cascaded <- cascaded$target[seq_len(nrow(fill_sample)), , drop = FALSE]
  fp_cascaded   <- cascaded$target[nrow(fill_sample) + seq_along(fp_test), , drop = FALSE]

  # production's design and 90/10 split
  design <- prep_bart_design(train[idx, predictors, drop = FALSE], fill_cascaded, logp)
  x_all  <- as.matrix(design$train)
  y_all  <- as.numeric(train[[b]][idx])
  set.seed(abs(as.integer(sprintf("%d%03d", p$subbasin, which(biotic_cols == b)))))
  holdout <- sample(seq_along(y_all), size = round(0.10 * length(y_all)))
  x_fit <- x_all[-holdout, , drop = FALSE]
  y_fit <- log1p(y_all[-holdout])
  scored <- if (length(holdout) > n_test) sample(holdout, n_test) else holdout
  x_test <- rbind(as.matrix(design$backfill), x_all[scored, , drop = FALSE])

  do(path("meta"), list(
    pair = p, n_fit = nrow(x_fit), n_predictors = ncol(x_fit), n_biotic_before = length(b_before),
    n_fill = nrow(fill), n_fill_sample = nrow(fill_sample), n_holdout = length(holdout),
    y_scored = y_all[scored], obs_fill = cascaded$obs_fill,
    n_fp_held = length(fp_held), cascade_seconds = cascaded$seconds,
    sd_y_log = sd(y_fit), mean_y = mean(y_all)))

  # chains
  for (ntree in arms[arms != 100L | nrow(x_fit) <= n_mid_max]) for (chain in seq_len(n_chains)) {
    do(path(sprintf("t%03d_c%d", ntree, chain)), {
      fit <- fit_gbart(x_fit, y_fit, x_test, ntree, nskip, n_kept, seed0 + 100000L * p$pair + 1000L * chain + ntree)
      say("%s: %d trees, chain %d: %.0f s, sigma %.3f", id, ntree, chain, fit$seconds, mean(fit$sigma[-seq_len(nskip)]))
      if (chain == 1L) do(path(sprintf("t%03d_psis", ntree)), psis_loo(fit, y_fit, fit$seconds))
      list(f = fit$yhat.test, sigma = fit$sigma, seconds = fit$seconds)
    })
  }

  # footprint-shaped holdout, at 08B's length
  rest_rows <- setdiff(seq_along(idx), fp_held)
  test_obs  <- train[idx[fp_test], ]
  test_obs[, intersect(footprint_cols, names(test_obs))] <- 0
  both <- rbind(test_obs[, names(train)], fp_cascaded[, names(train)])
  design_fp <- prep_bart_design(train[idx[rest_rows], predictors, drop = FALSE], both, logp)
  for (ntree in intersect(c(50L, 200L), arms)) {
    do(path(sprintf("fp_t%03d", ntree)), {
      fit <- fit_gbart(as.matrix(design_fp$train), log1p(y_all[rest_rows]), as.matrix(design_fp$backfill),
                       ntree, nskip, ndpost, seed0 + 7000L * p$pair + ntree)
      list(f = fit$yhat.test, sigma = fit$sigma[-seq_len(nskip)], y = y_all[fp_test],
           n_observed = length(fp_test), seconds = fit$seconds)
    })
  }

  # the same holdout from 90% of the fitting pixels, against the main chains' first 700 draws
  for (chain in 1:2) {
    do(path(sprintf("size_c%d", chain)), {
      set.seed(seed0 + 500L * p$pair + chain)
      sub <- sample(nrow(x_fit), round(0.9 * nrow(x_fit)))
      small <- fit_gbart(x_fit[sub, , drop = FALSE], y_fit[sub], x_test, 50L, nskip, ndpost, seed0 + 9000L * p$pair + chain)
      list(small = small$yhat.test, sigma_small = small$sigma[-seq_len(nskip)])
    })
  }
  say("%s done", id)
}


# 5. run -------------------------------------------------------------------------------------

todo <- if (is.na(only_pair)) seq_len(nrow(pairs)) else only_pair
for (i in todo) {
  outcome <- tryCatch(run_pair(pairs[i, ]), error = function(e) e)
  if (inherits(outcome, "error")) {
    say("pair %d: %s", i, conditionMessage(outcome))
    if (identical(conditionMessage(outcome), "deadline")) break
  }
}
say("finished")
