# ---
# title:   BART vs xgboost on a one-watershed backfilling problem
# author:  Mannfred Boehm
# created: November 3, 2026
# ---
#
# Worked example for a presentation: backfill one SCANFI covariate on the high-HF pixels of one
# small subbasin, from models trained on its low-HF pixels, with BART::gbart() and xgboost.
#
# What is held equal: pixels, the 5 abiotic predictors, the log1p response, the folds, one
# thread, and the deliverable, 100 draws per pixel (BART: 100 posterior draws; xgboost: 100
# bootstrap refits). BART keeps 08B's priors and is not tuned; xgboost is shown with its
# defaults and tuned (nested, so the tuning never sees the test pixels).
#
# BART's draws and the bootstrap's are not the same object. Both describe the mean function
# f(x), so neither covers a new pixel's value. For coverage and CRPS each draw gets noise on
# the log1p scale from the same Gaussian model: BART's sigma draw, or bootstrap b's out-of-bag
# residual SD.
#
# Convergence (section 6): gbart() starts the sparse prior at sweep nskip / 2, so a different
# burn-in is a different sampler. 4 long chains therefore run at nskip = 300, and every candidate
# burn-in is a later window of those same chains. Burn-in and mixing are separate questions, and
# mixing decides how far 100 draws from a chain are from 100 bootstrap refits, which are
# independent by construction. Sections 7-9 run BART as 08B does and as arm "BART long", with 10
# times the sweeps, to see whether any of this reaches the held-out scores.
#
# SCANFIprcC is OTHER conifers, not all conifers: the 8 named species + prcC + prcD sum to ~100.
#
# Stages are cached in out_dir: delete a stage's .rds to redo it. Run from the project root:
# Rscript Rscripts/misc/bart_vs_xgboost_benchmark.R

suppressPackageStartupMessages({
  library(BART)
  library(xgboost)
  library(ranger)
  library(posterior)
  library(terra)
  library(tidyverse)
  library(patchwork)
})

# 1. settings ------------------------------------------------------------------------------

ia_dir         <- getwd()
year           <- 2020
subbasin_index <- 388L   # smallest boreal (can61) subbasin with prcC SD >= 15, >= 1000 training and >= 150 backfill pixels
response       <- "SCANFIprcC_1km"
n_predictors   <- 5L
n_draws        <- 100L
seed           <- 20261103L

# 08B_deploy_gbart.R's call. Production keeps 700 draws and stores a random 100.
bart_args <- list(type = "wbart", k = 3, ntree = 50L, sparse = TRUE, sigdf = 3,
                  printevery = 1000000L)
ndpost    <- 700L
burn_in   <- 300L

# convergence
n_chains        <- 4L
chain_length    <- 20000L                    # kept sweeps after the 300 burn-in
warmup_sweeps   <- 1000L                     # the first sweeps, which gbart() does not return
candidate_burns <- c(300L, 1000L, 2000L, 3000L, 4000L, 5000L, 6000L, 7000L, 8000L)
chain_lengths   <- c(700L, 2000L, 5000L, 10000L)
n_rhat_pixels   <- 200L                      # random backfill pixels for pixel-level R-hat
rhat_max        <- 1.01                      # Vehtari et al. 2021
mixing_sweeps   <- 10000L                    # kept sweeps per chain in the settings comparison
burn_long       <- 5000L                     # "BART long": past sigma's drift (section 6), then
long_thin       <- 50L                       # 100 draws, one every 50 sweeps: 10,000 sweeps in all
mixing_settings <- tribble(~setting,                ~ntree, ~sparse,
                           "50 trees, sparse (08B)", 50L,   TRUE,
                           "50 trees, not sparse",   50L,   FALSE,
                           "200 trees, sparse",      200L,  TRUE)

# cross-validation: k folds, repeated with a new block grid or a new shuffle
n_folds   <- 5L
n_repeats <- c(spatial = 5L, random = 1L)   # random CV is the side comparison, and the slowest to tune
block_m   <- 10000                           # spatial blocks, metres (see section 3)

# xgboost: its defaults (it has no default nrounds in R; 100 is the Python API's), and the
# grid tuned by inner CV, which includes those defaults. Early stopping picks nrounds for each
# grid point.
xgb_default <- list(params = list(objective = "reg:squarederror", tree_method = "hist",
                                  eta = 0.3, max_depth = 6, min_child_weight = 1,
                                  subsample = 1, nthread = 1),
                    nrounds = 100L)
xgb_grid <- expand.grid(eta = c(0.03, 0.1, 0.3), max_depth = c(2L, 4L, 6L),
                        min_child_weight = c(1, 5, 20), subsample = c(0.7, 1))
xgb_max_rounds <- 1500L
xgb_patience   <- 30L
n_inner_folds  <- 3L

n_timing_reps <- 5L

window_km <- 30L   # fig 5 zooms to the square this wide that holds the most backfill pixels

out_dir <- file.path(ia_dir, "data", "derived_data", "rds_files", "bart_vs_xgboost")
fig_dir <- file.path(ia_dir, "output_figures", "bart_vs_xgboost")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

say <- function(...) message(format(Sys.time(), "[%H:%M:%S] "), sprintf(...))

# read a stage from out_dir, or run it and save it
cached <- function(name, expr) {
  path <- file.path(out_dir, paste0(name, ".rds"))
  if (file.exists(path)) {
    say("%s: cached", name)
    return(readRDS(path))
  }
  value <- expr
  saveRDS(value, path)
  value
}


# 2. one subbasin's pixels, as 08A crops them ------------------------------------------------

stack_y   <- rast(file.path(ia_dir, "data", "raw_data", "covariates_mosaiced",
                            sprintf("covariates_mosaiced_%d.tif", year)))
subbasins <- project(vect(file.path(ia_dir, "data", "raw_data",
                                    "hydrobasins_masked_merged_subset.gpkg")), crs(stack_y))
poly  <- subbasins[subbasin_index]
cov_s <- mask(crop(stack_y, poly), poly)

hf_mask <- function(kind) {
  path <- file.path(ia_dir, "data", "raw_data", "hirshpearson", sprintf("CanHF_1km_%s.tif", kind))
  mask(resample(rast(path), cov_s, method = "near"), poly)
}

# 07's abiotic_vars without V5's "Disturbance" class, which is the footprint itself (10C drops it
# too), and without 08A's x and y, so the 5 predictors are all a model gets
meta <- dplyr::filter(BAMexploreR::predictor_metadata, version == "v5")
candidates <- meta$predictor[meta$predictor_class %in%
                               c("Annual Climate", "Climate Normals", "Topography", "Wetland")]
candidates <- c(candidates, "CAfire", grep("^(cec|soc)_", names(stack_y), value = TRUE))
candidates <- intersect(setdiff(candidates, c("Peatland_1km", "Peatland_5x5")), names(stack_y))

pixels <- as.data.frame(c(cov_s[[c(response, candidates)]],
                          setNames(hf_mask("lessthan1"), "low"),
                          setNames(hf_mask("morethan1"), "high")),
                        xy = TRUE, na.rm = FALSE)

# xgboost takes NA natively and BART needs 08A's imputation, so both get complete rows: a
# candidate over 2% blank here is dropped, and the few pixels the rest leave blank go too
train_rows <- which(pixels$low == 1 & !is.na(pixels[[response]]))
fill_rows  <- which(pixels$high == 1)
na_share   <- map_dbl(candidates, function(v) mean(is.na(pixels[c(train_rows, fill_rows), v])))
varies     <- map_lgl(candidates, function(v) isTRUE(sd(pixels[train_rows, v], na.rm = TRUE) > 0))
usable     <- na_share <= 0.02 & varies
say("%d candidates, %d usable; dropped: %s", length(candidates), sum(usable),
    paste(candidates[!usable], collapse = ", "))
candidates <- candidates[usable]
complete   <- complete.cases(pixels[, candidates])
train_rows <- train_rows[complete[train_rows]]
fill_rows  <- fill_rows[complete[fill_rows]]

train <- pixels[train_rows, ]
fill  <- pixels[fill_rows, ]
y_pct <- train[[response]]
y_log <- log1p(y_pct)
say("subbasin %d: %d training pixels, %d backfill pixels; %s mean %.1f, SD %.1f, range %.1f-%.1f",
    subbasin_index, nrow(train), nrow(fill), response, mean(y_pct), sd(y_pct),
    min(y_pct), max(y_pct))


# 3. folds -----------------------------------------------------------------------------------

# Random folds leave a held-out pixel's neighbours in training, so they reward what a tree can
# do with smooth climate layers: find a place. Backfilling asks for land the models did not see,
# so spatial blocks are the headline. 10 km blocks hold test pixels a median 2.7 km (95th
# percentile 5.3 km) from the nearest training pixel; the real backfill pixels are 1.4 km (4.1 km)
# away. Blocks go to folds largest first, onto the emptiest fold, from a grid with a random offset.
make_folds <- function(scheme, rows, k, seed) {
  set.seed(seed)
  if (scheme == "random") return(sample(rep_len(seq_len(k), length(rows))))
  offset <- runif(2, 0, block_m)
  block  <- paste(floor((train$x[rows] + offset[1]) / block_m),
                  floor((train$y[rows] + offset[2]) / block_m))
  sizes  <- sort(table(block), decreasing = TRUE)
  load   <- numeric(k)
  fold_of_block <- setNames(integer(length(sizes)), names(sizes))
  for (b in names(sizes)) {
    f <- which.min(load + runif(k, 0, 1e-6))
    fold_of_block[b] <- f
    load[f] <- load[f] + sizes[[b]]
  }
  unname(fold_of_block[block])
}


# 4. the 5 predictors ------------------------------------------------------------------------

# Greedy forward selection on spatial-CV RMSE, using a third learner (random forest) so that
# neither contender picks its own features. Adding a predictor must improve held-out skill, so
# collinear climate layers do not crowd out the rest. Selection sees every training pixel: any
# optimism it brings is shared by both models.
predictors <- cached("predictors", {
  folds <- make_folds("spatial", seq_len(nrow(train)), n_folds, seed)
  cv_rmse <- function(vars) {
    sse <- 0
    for (f in seq_len(n_folds)) {
      test <- folds == f
      rf   <- ranger(x = train[!test, vars, drop = FALSE], y = y_log[!test], num.trees = 200,
                     num.threads = 4, seed = seed)
      sse  <- sse + sum((y_log[test] - predict(rf, train[test, vars, drop = FALSE])$predictions)^2)
    }
    sqrt(sse / nrow(train))
  }
  chosen <- tibble(predictor = character(), cv_rmse_log1p = numeric())
  for (i in seq_len(n_predictors)) {
    left   <- setdiff(candidates, chosen$predictor)
    rmses  <- map_dbl(left, function(v) cv_rmse(c(chosen$predictor, v)))
    chosen <- add_row(chosen, predictor = left[which.min(rmses)], cv_rmse_log1p = min(rmses))
    say("predictor %d: %s (spatial-CV RMSE %.3f on log1p)", i, left[which.min(rmses)], min(rmses))
  }
  chosen
})
write_csv(left_join(predictors, select(meta, predictor, definition, predictor_class),
                    by = "predictor"),
          file.path(out_dir, "predictors.csv"))

x_train <- as.matrix(train[, predictors$predictor])
x_fill  <- as.matrix(fill[, predictors$predictor])


# 5. the two methods -------------------------------------------------------------------------

# 08B's call. With keepevery > 1 gbart() warns while sizing its log-likelihood, which we don't use.
fit_bart <- function(x, y, x_test, nskip, ndpost, seed, keepevery = 1L,
                     ntree = bart_args$ntree, sparse = bart_args$sparse) {
  set.seed(seed)
  capture.output(
    fit <- suppressWarnings(
      BART::gbart(x.train = x, y.train = y, x.test = x_test,
                  type = bart_args$type, k = bart_args$k, ntree = ntree,
                  ndpost = ndpost, nskip = nskip, keepevery = keepevery, sparse = sparse,
                  sigest = sd(y), sigdf = bart_args$sigdf,
                  printevery = bart_args$printevery))
  )
  fit
}

# 100 draws of f at x_test, and each one's sigma. fit$sigma holds every sweep, burn-in included.
bart_draws <- function(x, y, x_test, nskip, ndpost, seed, keepevery = 1L) {
  fit  <- fit_bart(x, y, x_test, nskip, ndpost, seed, keepevery)
  keep <- sample(ndpost, n_draws)
  list(f = fit$yhat.test[keep, , drop = FALSE], sigma = fit$sigma[nskip + keep * keepevery])
}

# 100 refits on row bootstraps of the training pixels, each with its out-of-bag residual SD
xgb_draws <- function(x, y, x_test, config, seed) {
  set.seed(seed)
  n      <- nrow(x)
  d_test <- xgb.DMatrix(x_test, nthread = 1)
  f      <- matrix(NA_real_, n_draws, nrow(x_test))
  sigma  <- numeric(n_draws)
  for (b in seq_len(n_draws)) {
    rows <- sample.int(n, n, replace = TRUE)
    oob  <- setdiff(seq_len(n), rows)
    fit  <- xgb.train(config$params,
                      xgb.DMatrix(x[rows, , drop = FALSE], label = y[rows], nthread = 1),
                      nrounds = config$nrounds, verbose = 0)
    f[b, ]   <- predict(fit, d_test)
    sigma[b] <- sqrt(mean((y[oob] - predict(fit, x[oob, , drop = FALSE]))^2))
  }
  list(f = f, sigma = sigma)
}

# best of xgb_grid on inner folds of the training pixels it is given
tune_xgb <- function(x, y, inner_folds) {
  data  <- xgb.DMatrix(x, label = y, nthread = 1)
  folds <- unname(split(seq_along(inner_folds), inner_folds))
  grid  <- pmap_dfr(xgb_grid, function(eta, max_depth, min_child_weight, subsample) {
    params <- list(objective = "reg:squarederror", tree_method = "hist", eta = eta,
                   max_depth = max_depth, min_child_weight = min_child_weight,
                   subsample = subsample, nthread = 1)
    cv <- xgb.cv(params, data, nrounds = xgb_max_rounds, folds = folds, metrics = "rmse",
                 early_stopping_rounds = xgb_patience, maximize = FALSE, verbose = 0)
    tibble(eta, max_depth, min_child_weight, subsample, nrounds = cv$best_iteration,
           rmse = cv$evaluation_log$test_rmse_mean[cv$best_iteration])
  })
  best <- slice_min(grid, rmse, n = 1, with_ties = FALSE)
  list(params = list(objective = "reg:squarederror", tree_method = "hist", eta = best$eta,
                     max_depth = best$max_depth, min_child_weight = best$min_child_weight,
                     subsample = best$subsample, nthread = 1),
       nrounds = best$nrounds)
}

# log1p back to % cover. 12F clamps at 0 too; the cap at 100 only matters once noise is added.
to_pct <- function(z) pmin(pmax(expm1(z), 0), 100)

# one Gaussian noise draw per element, with row d's sigma
with_noise <- function(f, sigma) f + rnorm(length(f), 0, rep(sigma, times = ncol(f)))

# CRPS of each column's sample: E|X - y| - E|X - X'| / 2
crps_sample <- function(x, y) {
  m <- nrow(x)
  colMeans(abs(sweep(x, 2, y))) - colSums(apply(x, 2, sort) * (2 * seq_len(m) - m - 1)) / m^2
}

pixel_scores <- function(draws, y) {
  f   <- to_pct(draws$f)
  new <- to_pct(with_noise(draws$f, draws$sigma))
  f_q <- apply(f,   2, quantile, probs = c(0.05, 0.95))
  p_q <- apply(new, 2, quantile, probs = c(0.05, 0.95))
  tibble(y = y, pred = colMeans(f), f_lo = f_q[1, ], f_hi = f_q[2, ],
         p_lo = p_q[1, ], p_hi = p_q[2, ], crps = crps_sample(new, y))
}


# 6. is a burn-in of 300 enough? ------------------------------------------------------------

# Two questions. Burn-in: have the chains forgotten where they started? That shows as drift all
# 4 chains share. Mixing: once they have, how long does a chain take to get around the posterior?
# That shows as chains that disagree. R-hat fails on slow mixing whatever the burn-in, so burn-in
# is read from shared drift, and R-hat and ESS by chain length. Both are reported as effect sizes:
# at this ESS neither Monte Carlo standard errors nor a t test on 4 chain means are reliable.

chains <- cached("bart_chains", map(seq_len(n_chains), function(chain) {
  say("chain %d of %d: %d + %d sweeps", chain, n_chains, burn_in, chain_length)
  fit <- fit_bart(x_train, y_log, x_fill, burn_in, chain_length, seed + chain)
  list(sigma = fit$sigma, f = fit$yhat.test, splits = rowSums(fit$varcount))
}))

# gbart() returns f only after burn-in. To watch the burn-in itself, these chains keep sweep 2
# onward. They differ from production's only in starting the sparse prior at once.
warmup <- cached("bart_warmup", map_dfr(seq_len(n_chains), function(chain) {
  fit <- fit_bart(x_train, y_log, x_fill, 1L, warmup_sweeps, seed + 10L + chain)
  tibble(chain = chain, sweep = seq_len(warmup_sweeps) + 1L, sigma = fit$sigma[-1],
         footprint_pct = rowMeans(to_pct(fit$yhat.test)), total_splits = rowSums(fit$varcount))
}))

# [sweep after burn_in x chain] for each quantity. sigma and the footprint mean are what the
# draws feed downstream; total splits is tree size.
trace <- list(
  sigma         = sapply(chains, function(ch) ch$sigma[-seq_len(burn_in)]),
  footprint_pct = sapply(chains, function(ch) rowMeans(to_pct(ch$f))),
  total_splits  = sapply(chains, function(ch) ch$splits)
)
reference <- seq(chain_length / 2 + 1, chain_length)   # the chains' second half

# Burn-in b: the 700 draws production would keep after discarding b sweeps, against the
# reference, in the quantity's own units and in posterior SDs
burn_in_table <- map_dfr(candidate_burns, function(burn) {
  window <- (burn - burn_in) + seq_len(ndpost)
  imap_dfr(trace, function(m, quantity) {
    w <- m[window, ]
    r <- m[reference, ]
    tibble(burn_in = burn, quantity = quantity, mean_window = mean(w), mean_reference = mean(r),
           shift = mean(w) - mean(r), shift_in_sd = (mean(w) - mean(r)) / sd(r))
  })
})
write_csv(burn_in_table, file.path(out_dir, "convergence_burn_in.csv"))
print(as.data.frame(burn_in_table), digits = 3)

# the sweep at which each warm-up chain first enters the reference's central 90%
reference_band <- imap_dfr(trace, function(m, quantity) {
  q <- quantile(m[reference, ], c(0.05, 0.95))
  tibble(quantity = quantity, q05 = q[[1]], q95 = q[[2]])
})
warmup_entry <- warmup |>
  pivot_longer(c(sigma, footprint_pct, total_splits), names_to = "quantity") |>
  left_join(reference_band, by = "quantity") |>
  group_by(quantity, chain) |>
  summarise(first_inside = min(sweep[value >= q05 & value <= q95]), .groups = "drop_last") |>
  summarise(latest_chain = max(first_inside), .groups = "drop")
write_csv(warmup_entry, file.path(out_dir, "convergence_warmup_entry.csv"))
print(as.data.frame(warmup_entry))

# the first candidate from which sigma stays within 0.5 posterior SD of its long-run mean
sigma_settled <- burn_in_table |>
  filter(quantity == "sigma") |>
  arrange(desc(burn_in)) |>
  filter(cumall(abs(shift_in_sd) <= 0.5)) |>
  pull(burn_in) |>
  min()
say("sigma settled by sweep %d (production's burn-in: %d)", sigma_settled, burn_in)

# mixing: R-hat and ESS by chain length, from production's burn-in and a much later start. If
# burn-in were all that was wrong, the later start would cure it at equal length.
set.seed(seed)
rhat_pixels <- sample(ncol(chains[[1]]$f), n_rhat_pixels)

mixing_by_length <- cached("mixing_by_length", {
  map_dfr(c(burn_in, max(candidate_burns)), function(burn) {
    map_dfr(chain_lengths, function(len) {
      say("mixing: burn-in %d, %d kept sweeps", burn, len)
      rows   <- (burn - burn_in) + seq_len(len)
      global <- imap_dfr(trace, function(m, quantity) {
        tibble(quantity = quantity, rhat = rhat(m[rows, ]), ess_bulk = ess_bulk(m[rows, ]),
               ess_tail = ess_tail(m[rows, ]))
      })
      pixel_stats <- do.call(rbind, map(rhat_pixels, function(p) {
        w <- sapply(chains, function(ch) ch$f[rows, p])
        c(rhat = rhat(w), ess_bulk = ess_bulk(w), ess_tail = ess_tail(w))
      }))
      bind_rows(global,
                tibble(quantity = "pixel (median)", rhat = median(pixel_stats[, "rhat"]),
                       ess_bulk = median(pixel_stats[, "ess_bulk"]),
                       ess_tail = median(pixel_stats[, "ess_tail"])),
                tibble(quantity = "pixel (worst 5%)",
                       rhat = unname(quantile(pixel_stats[, "rhat"], 0.95)),
                       ess_bulk = unname(quantile(pixel_stats[, "ess_bulk"], 0.05)),
                       ess_tail = unname(quantile(pixel_stats[, "ess_tail"], 0.05)))) |>
        mutate(burn_in = burn, kept_sweeps = len, .before = 1)
    })
  })
})
write_csv(mixing_by_length, file.path(out_dir, "convergence_mixing.csv"))
print(as.data.frame(mixing_by_length), digits = 4)

# sweeps per effective draw: what a chain would run for 100 draws to act like 100 independent ones
longest <- filter(mixing_by_length, burn_in == max(candidate_burns),
                  kept_sweeps == max(chain_lengths), quantity == "pixel (median)")
sweeps_per_draw <- max(chain_lengths) * n_chains / longest$ess_bulk
say("%.0f sweeps per effective draw at the median pixel: 100 independent-like draws need %s sweeps",
    sweeps_per_draw, format(round(n_draws * sweeps_per_draw), big.mark = ","))

# What a stored sample of 100 draws gets right, per pixel on % cover, against the 4 chains'
# reference: production's (a random 100 of the 700 after burn-in 300), and the long arm's
pooled      <- do.call(rbind, map(chains, function(ch) to_pct(ch$f[reference, ])))
pooled_mean <- colMeans(pooled)
pooled_sd   <- apply(pooled, 2, sd)
rm(pooled)
set.seed(seed)
sample_quality <- map_dfr(seq_len(n_chains), function(chain) {
  arms <- list("BART"      = sort(sample(ndpost, n_draws)),
               "BART long" = (burn_long - burn_in) + long_thin * seq_len(n_draws))
  imap_dfr(arms, function(rows, arm) {
    kept  <- to_pct(chains[[chain]]$f[rows, ])
    shift <- abs(colMeans(kept) - pooled_mean)
    tibble(arm = arm, chain = chain, ess_of_100_pixel_median = median(apply(kept, 2, ess_basic)),
           mean_shift_pct_median = median(shift), mean_shift_pct_q95 = unname(quantile(shift, 0.95)),
           mean_shift_in_sd_median = median(shift / pooled_sd),
           sd_ratio_median = median(apply(kept, 2, sd) / pooled_sd))
  })
})
write_csv(sample_quality, file.path(out_dir, "convergence_sample_quality.csv"))
print(as.data.frame(sample_quality), digits = 3)

# Which of 08B's settings drives this? sigma_shift is the first 700 kept sweeps against the
# chain's second half, in posterior SDs: the same drift as above.
mixing_by_setting <- cached("mixing_by_setting", pmap_dfr(mixing_settings, function(setting, ntree, sparse) {
  runs <- map(seq_len(n_chains), function(chain) {
    say("%s: chain %d of %d", setting, chain, n_chains)
    seconds <- system.time(
      fit <- fit_bart(x_train, y_log, x_fill, burn_in, mixing_sweeps, seed + 20L + chain,
                      ntree = ntree, sparse = sparse))[["elapsed"]]
    list(seconds = seconds, f = fit$yhat.test[, rhat_pixels], sigma = fit$sigma[-seq_len(burn_in)])
  })
  stats <- do.call(rbind, map(seq_along(rhat_pixels), function(p) {
    w <- sapply(runs, function(r) r$f[, p])
    c(rhat = rhat(w), ess_bulk = ess_bulk(w))
  }))
  sigma   <- sapply(runs, function(r) r$sigma)
  late    <- sigma[seq(mixing_sweeps / 2 + 1, mixing_sweeps), ]
  seconds <- sum(map_dbl(runs, "seconds"))
  tibble(setting = setting, seconds_per_chain = seconds / n_chains,
         sigma_early = mean(sigma[seq_len(ndpost), ]), sigma_late = mean(late),
         sigma_shift_in_sd = (mean(sigma[seq_len(ndpost), ]) - mean(late)) / sd(late),
         rhat_pixel_median = median(stats[, "rhat"]), ess_pixel_median = median(stats[, "ess_bulk"]),
         sweeps_per_effective_draw = mixing_sweeps * n_chains / median(stats[, "ess_bulk"]),
         effective_draws_per_minute = 60 * median(stats[, "ess_bulk"]) / seconds)
}))
write_csv(mixing_by_setting, file.path(out_dir, "convergence_by_setting.csv"))
print(as.data.frame(mixing_by_setting), digits = 3)


# 7. cross-validation -----------------------------------------------------------------------

# Each repeat is cached alone and its seeds depend only on r, so raising n_repeats runs only the
# new repeats.
run_cv <- function(scheme) {
  map_dfr(seq_len(n_repeats[[scheme]]), function(r) cached(sprintf("cv_%s_rep%d", scheme, r), {
    folds <- make_folds(scheme, seq_len(nrow(train)), n_folds, seed + 100L * r)
    map_dfr(seq_len(n_folds), function(k) {
      say("cv %s: repeat %d of %d, fold %d of %d", scheme, r, n_repeats[[scheme]], k, n_folds)
      test  <- which(folds == k)
      fit   <- which(folds != k)
      s     <- seed + 1000L * r + k
      tuned <- tune_xgb(x_train[fit, ], y_log[fit], make_folds(scheme, fit, n_inner_folds, s))
      x_fit  <- x_train[fit, ]
      x_test <- x_train[test, , drop = FALSE]
      draws  <- list(
        "BART"            = bart_draws(x_fit, y_log[fit], x_test, burn_in, ndpost, s),
        "xgboost default" = xgb_draws(x_fit, y_log[fit], x_test, xgb_default, s),
        "xgboost tuned"   = xgb_draws(x_fit, y_log[fit], x_test, tuned, s))
      if (scheme == "spatial")
        draws[["BART long"]] <- bart_draws(x_fit, y_log[fit], x_test, burn_long, n_draws, s,
                                           keepevery = long_thin)
      imap_dfr(draws, function(d, method) {
        mutate(pixel_scores(d, y_pct[test]), method = method, scheme = scheme, rep = r, fold = k,
               pixel = test, y_fit_mean = mean(y_pct[fit]), nrounds = tuned$nrounds,
               eta = tuned$params$eta, max_depth = tuned$params$max_depth)
      })
    })
  }))
}

cv <- bind_rows(run_cv("spatial"), run_cv("random"))

# each repeat predicts every pixel once, so each repeat gives one score per method. R2 is skill
# over predicting the training mean, so it can go below 0.
cv_metrics <- cv |>
  group_by(scheme, method, rep) |>
  summarise(rmse = sqrt(mean((y - pred)^2)), mae = mean(abs(y - pred)), bias = mean(pred - y),
            r2 = 1 - sum((y - pred)^2) / sum((y - y_fit_mean)^2),
            crps = mean(crps),
            coverage_90 = mean(y >= p_lo & y <= p_hi), width_90 = mean(p_hi - p_lo),
            coverage_90_f_only = mean(y >= f_lo & y <= f_hi), width_90_f_only = mean(f_hi - f_lo),
            .groups = "drop")
cv_summary <- cv_metrics |>
  group_by(scheme, method) |>
  summarise(across(rmse:width_90_f_only, mean), .groups = "drop")
write_csv(cv_metrics, file.path(out_dir, "cv_metrics.csv"))
write_csv(cv_summary, file.path(out_dir, "cv_summary.csv"))
print(as.data.frame(cv_summary), digits = 3)
print(as.data.frame(count(distinct(cv, scheme, rep, fold, eta, max_depth, nrounds), scheme, eta, max_depth)))


# 8. compute time ---------------------------------------------------------------------------

# Every workflow does the real job, training on all the low-HF pixels and predicting the backfill
# ones, on one thread (BART::gbart is single-threaded, and so is production). Workflows take
# turns so that a busy moment on the laptop does not land on one of them.
tuned_full <- cached("xgb_tuned_full",
                     tune_xgb(x_train, y_log,
                              make_folds("spatial", seq_len(nrow(train)), n_inner_folds, seed)))

timing <- cached("timing", {
  workflows <- list(
    "BART long"                        = function(s) bart_draws(x_train, y_log, x_fill, burn_long, n_draws, s, keepevery = long_thin),
    "BART, keep 100 of 700 (as 08B)"   = function(s) bart_draws(x_train, y_log, x_fill, burn_in, ndpost, s),
    "BART, keep all 100"               = function(s) bart_draws(x_train, y_log, x_fill, burn_in, n_draws, s),
    "xgboost default, 100 refits"      = function(s) xgb_draws(x_train, y_log, x_fill, xgb_default, s),
    "xgboost tuned, 100 refits"        = function(s) xgb_draws(x_train, y_log, x_fill, tuned_full, s),
    "xgboost tuning alone"             = function(s) tune_xgb(x_train, y_log, make_folds("spatial", seq_len(nrow(train)), n_inner_folds, s)),
    "xgboost default, one fit"         = function(s) {
      predict(xgb.train(xgb_default$params, xgb.DMatrix(x_train, label = y_log, nthread = 1),
                        nrounds = xgb_default$nrounds, verbose = 0), x_fill)
    })
  map_dfr(seq_len(n_timing_reps), function(r) {
    imap_dfr(workflows, function(run, workflow) {
      gc(verbose = FALSE)
      seconds <- system.time(run(seed + r))[["elapsed"]]
      say("timing %d of %d: %s, %.1f s", r, n_timing_reps, workflow, seconds)
      tibble(rep = r, workflow = workflow, seconds = seconds)
    })
  })
})
timing_summary <- timing |>
  group_by(workflow) |>
  summarise(median = median(seconds), min = min(seconds), max = max(seconds), .groups = "drop") |>
  arrange(desc(median))
write_csv(timing_summary, file.path(out_dir, "timing_summary.csv"))
print(as.data.frame(timing_summary), digits = 3)


# 9. the backfill ---------------------------------------------------------------------------

backfill <- cached("backfill", list(
  "BART"            = bart_draws(x_train, y_log, x_fill, burn_in, ndpost, seed),
  "BART long"       = bart_draws(x_train, y_log, x_fill, burn_long, n_draws, seed, keepevery = long_thin),
  "xgboost default" = xgb_draws(x_train, y_log, x_fill, xgb_default, seed),
  "xgboost tuned"   = xgb_draws(x_train, y_log, x_fill, tuned_full, seed)))

backfill_summary <- imap_dfr(backfill, function(d, method) {
  f <- to_pct(d$f)
  footprint <- rowMeans(f)
  tibble(method = method, footprint_mean = mean(footprint),
         footprint_q05 = unname(quantile(footprint, 0.05)),
         footprint_q95 = unname(quantile(footprint, 0.95)),
         pixel_sd_median = median(apply(f, 2, sd)))
}) |>
  add_row(method = "observed", footprint_mean = mean(fill[[response]], na.rm = TRUE))
write_csv(backfill_summary, file.path(out_dir, "backfill_summary.csv"))
print(as.data.frame(backfill_summary), digits = 3)


# 10. figures, sized for a 16:9 slide --------------------------------------------------------

ink      <- "#0b0b0b"
ink_soft <- "#52514e"
surface  <- "#fcfcfb"
rule     <- "#e4e3de"
chain_colours  <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100")
method_colours <- c("BART" = "#2a78d6", "xgboost tuned" = "#eb6834",
                    "xgboost default" = "#1baf7a", "BART long" = "#eda100")
cover_ramp  <- c("#cde2fb", "#6da7ec", "#256abf", "#0d366b")
spread_ramp <- c("#fde6dc", "#f4a27f", "#eb6834", "#7a2c10")

theme_slide <- function(base_size = 15) {
  theme_minimal(base_size = base_size) +
    theme(plot.background  = element_rect(fill = surface, colour = NA),
          panel.background = element_rect(fill = surface, colour = NA),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = rule, linewidth = 0.3),
          text             = element_text(colour = ink),
          axis.text        = element_text(colour = ink_soft),
          strip.text       = element_text(colour = ink, face = "bold", hjust = 0),
          plot.title       = element_text(face = "bold"),
          plot.title.position = "plot",
          plot.subtitle    = element_text(colour = ink_soft),
          plot.caption     = element_text(colour = ink_soft, hjust = 0),
          plot.caption.position = "plot",
          legend.position  = "none")
}

save_slide <- function(plot, name, width = 13.33, height = 7.5) {
  ggsave(file.path(fig_dir, name), plot, width = width, height = height, dpi = 200, bg = surface)
}

# a caption under a multi-panel figure
figure_caption <- function(text) {
  plot_annotation(caption = text,
                  theme = theme(plot.background = element_rect(fill = surface, colour = NA),
                                plot.caption = element_text(colour = ink_soft, hjust = 0, size = 13)))
}

# [sweep x chain] matrix to long, as trailing means of k sweeps, every step-th sweep
rolling <- function(m, first_sweep, k, step) {
  smooth <- apply(m, 2, function(v) as.numeric(stats::filter(v, rep(1 / k, k), sides = 1)))
  rows   <- seq(k, nrow(m), by = step)
  map_dfr(seq_len(ncol(m)), function(i) {
    tibble(chain = factor(i), sweep = first_sweep + rows - 1L, value = smooth[rows, i])
  })
}

# burn-in: sigma zoomed in, where 300 looks fine, and zoomed out, where it does not
sigma_all <- sapply(chains, function(ch) ch$sigma)
sigma_ref <- mean(trace$sigma[reference, ])

first_sweeps <- map_dfr(seq_len(n_chains), function(i) {
  tibble(chain = factor(i), sweep = 1:1000, value = sigma_all[1:1000, i])
})
fig_zoom_in <- ggplot(first_sweeps, aes(sweep, value, colour = chain)) +
  geom_line(linewidth = 0.4) +
  geom_vline(xintercept = burn_in, colour = ink_soft, linetype = "dashed", linewidth = 0.6) +
  annotate("text", x = burn_in + 20, y = Inf, vjust = 1.5, hjust = 0, colour = ink_soft, size = 4.6,
           label = sprintf("burn-in %d", burn_in)) +
  scale_colour_manual(values = chain_colours) +
  scale_x_continuous(labels = scales::comma) +
  labs(x = "Sweep", y = "Residual SD (sigma, log1p scale)", title = "First 1,000 sweeps: looks settled by 300") +
  theme_slide()

fig_zoom_out <- ggplot(rolling(sigma_all, 1L, 500L, 25L), aes(sweep, value, colour = chain)) +
  geom_hline(yintercept = sigma_ref, colour = ink_soft, linewidth = 0.5) +
  geom_vline(xintercept = burn_in, colour = ink_soft, linetype = "dashed", linewidth = 0.6) +
  geom_line(linewidth = 0.6) +
  scale_colour_manual(values = chain_colours) +
  scale_x_continuous(labels = scales::comma) +
  labs(x = "Sweep", y = "Sigma, mean of the last 500 sweeps",
       title = sprintf("All 20,300 sweeps: falling until about %s", scales::comma(sigma_settled))) +
  theme_slide()
save_slide(fig_zoom_in + fig_zoom_out +
             figure_caption("Each colour is one chain. Solid line: mean of sweeps 10,301 to 20,300."),
           "fig1_burn_in.png", height = 5.7)

# mixing: the footprint mean wanders without a shared trend, and R-hat by chain length
fig_wander <- ggplot(rolling(trace$footprint_pct, burn_in + 1L, 500L, 25L), aes(sweep, value, colour = chain)) +
  geom_line(linewidth = 0.6) +
  scale_colour_manual(values = chain_colours) +
  scale_x_continuous(labels = scales::comma) +
  labs(x = "Sweep", y = "Mean cover on the footprint (%)", title = "Chains wander apart and back") +
  theme_slide()

start_labels <- sprintf("after burn-in %s", scales::comma(c(burn_in, max(candidate_burns))))
rhat_long <- mixing_by_length |>
  filter(quantity %in% c("pixel (median)", "pixel (worst 5%)")) |>
  mutate(start = factor(sprintf("after burn-in %s", scales::comma(burn_in)), levels = start_labels))

fig_rhat <- ggplot(rhat_long, aes(kept_sweeps, rhat, colour = start, linetype = quantity)) +
  geom_hline(yintercept = rhat_max, colour = ink_soft, linewidth = 0.5) +
  annotate("text", x = min(chain_lengths), y = rhat_max, vjust = -0.6, hjust = 0, colour = ink_soft,
           size = 4.2, label = sprintf("target %.2f", rhat_max)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.6) +
  scale_x_log10(breaks = chain_lengths, labels = scales::comma) +
  scale_colour_manual(values = chain_colours[1:2], name = NULL) +
  scale_linetype_manual(values = c("pixel (median)" = "solid", "pixel (worst 5%)" = "dashed"), name = NULL) +
  labs(x = "Sweeps kept per chain (4 chains)", y = "R-hat, backfill pixels",
       title = "R-hat by chain length") +
  theme_slide() +
  theme(legend.position = "top", legend.justification = "left", legend.box = "vertical",
        legend.margin = margin(0, 0, 0, 0))
save_slide(fig_wander + fig_rhat +
             figure_caption("Left: each colour is one chain, as the mean of its last 500 sweeps. Right: 200 random backfill pixels."),
           "fig2_mixing.png", height = 5.7)

# compute time: mean and SD of the runs. The tuned row adds each run's tuning and refits first.
tune_label  <- "xgboost tuned: grid search, then 100 refits"
time_labels <- c(
  "BART long"                      = sprintf("BART long: %s sweeps, keep every %dth",
                                             scales::comma(burn_long + n_draws * long_thin), long_thin),
  "BART, keep 100 of 700 (as 08B)" = "BART as 08B: 1,000 sweeps, keep 100 of 700",
  "BART, keep all 100"             = "BART shortest: 400 sweeps, keep all 100",
  "xgboost tuning alone"           = tune_label,
  "xgboost tuned, 100 refits"      = tune_label,
  "xgboost default, 100 refits"    = "xgboost default: 100 refits",
  "xgboost default, one fit"       = "xgboost default: one fit, no uncertainty")
time_points <- timing |>
  mutate(label  = time_labels[workflow],
         family = if_else(startsWith(workflow, "BART"), "BART", "xgboost")) |>
  group_by(label, family, rep) |>
  summarise(seconds = sum(seconds), .groups = "drop_last") |>
  summarise(mean = mean(seconds), sd = sd(seconds), .groups = "drop") |>
  arrange(desc(family), mean) |>
  mutate(label = fct_inorder(label),
         text  = sprintf(if_else(mean < 10, "%.2f ± %.2f s", "%.1f ± %.1f s"), mean, sd))
write_csv(time_points, file.path(out_dir, "timing_mean_sd.csv"))

fig_time <- ggplot(time_points, aes(mean, label, colour = family)) +
  geom_errorbar(aes(xmin = mean - sd, xmax = mean + sd), width = 0.3, linewidth = 0.7) +
  geom_point(size = 3) +
  geom_text(aes(x = mean + sd, label = text), hjust = -0.14, colour = ink, size = 4.8) +
  scale_colour_manual(values = c("BART" = "#2a78d6", "xgboost" = "#eb6834")) +
  scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0.01, 0.16))) +
  labs(x = sprintf("Seconds on one thread, mean ± SD of %d runs", n_timing_reps), y = NULL,
       title = sprintf("Wall time to backfill %d pixels from %s training pixels", nrow(fill), scales::comma(nrow(train)))) +
  theme_slide()
save_slide(fig_time, "fig3_compute_time.png", height = 5.2)

# cross-validation on the spatial folds, the headline (section 3); cv_metrics.csv has the random ones
metric_labels <- c(rmse = "RMSE (% cover)", crps = "CRPS (% cover)",
                   coverage_90 = "Coverage of 90% interval", width_90 = "Width of 90% interval (% cover)")
method_order  <- c("xgboost default", "xgboost tuned", "BART long", "BART")

cv_long <- cv_metrics |>
  filter(scheme == "spatial") |>
  pivot_longer(all_of(names(metric_labels)), names_to = "metric") |>
  mutate(metric = factor(metric_labels[metric], levels = metric_labels),
         method = factor(method, levels = method_order))
cv_means <- cv_long |>
  group_by(metric, method) |>
  summarise(value = mean(value), .groups = "drop")
nominal <- tibble(metric = factor(metric_labels["coverage_90"], levels = metric_labels), value = 0.9)

fig_cv <- ggplot(cv_long, aes(value, method, colour = method)) +
  geom_vline(data = nominal, aes(xintercept = value), colour = ink_soft, linetype = "dashed") +
  geom_point(size = 2.4, alpha = 0.55) +
  geom_point(data = cv_means, size = 5, shape = 124) +
  facet_wrap(vars(metric), nrow = 1, scales = "free_x") +
  scale_colour_manual(values = method_colours) +
  labs(x = NULL, y = NULL, title = sprintf("Held-out performance, %g km spatial blocks", block_m / 1000)) +
  theme_slide(14) +
  theme(panel.spacing.x = unit(1.4, "lines"))
save_slide(fig_cv, "fig4_cross_validation.png", height = 4.8)

# the backfilled landscape, laid out as figure_athabasca_osa_backfill.R: the subbasin for context,
# then the window. A window centred on a cell corner has its edges on cell edges.
map_xgb <- "xgboost tuned"
cell    <- res(cov_s)[1]
half    <- window_km * 1000 / 2
centres <- expand_grid(x = seq(min(fill$x), max(fill$x), by = cell) + cell / 2,
                       y = seq(min(fill$y), max(fill$y), by = cell) + cell / 2) |>
  mutate(held = map2_int(x, y, function(cx, cy) sum(abs(fill$x - cx) < half & abs(fill$y - cy) < half)))
best  <- slice_max(centres, held, n = 1, with_ties = FALSE)
win_x <- best$x + c(-half, half)
win_y <- best$y + c(-half, half)
in_window <- function(d) filter(d, between(x, win_x[1], win_x[2]), between(y, win_y[1], win_y[2]))
say("%d km window holds %d of %d backfill pixels", window_km, best$held, nrow(fill))

basin      <- as.vector(ext(poly))
context_r  <- 0.55 * max(basin[["xmax"]] - basin[["xmin"]], basin[["ymax"]] - basin[["ymin"]])
context_x  <- mean(basin[c("xmin", "xmax")]) + c(-context_r, context_r)
context_y  <- mean(basin[c("ymin", "ymax")]) + c(-context_r, context_r)

draw_summary <- function(d) {
  f <- to_pct(d$f)
  in_window(tibble(x = fill$x, y = fill$y, cover = colMeans(f), spread = apply(f, 2, sd)))
}
bart_fill <- draw_summary(backfill[["BART"]])
xgb_fill  <- draw_summary(backfill[[map_xgb]])

# what stays on the training pixels, with the backfill's cover on the rest
landscape <- function(on_fill) {
  bind_rows(in_window(tibble(x = train$x, y = train$y, cover = train[[response]])), on_fill)
}
observed_map <- landscape(in_window(tibble(x = fill$x, y = fill$y, cover = fill[[response]])))
bart_map     <- landscape(select(bart_fill, x, y, cover))
xgb_map      <- landscape(select(xgb_fill, x, y, cover))
footprint_map <- in_window(as.data.frame(cov_s[["CanHF_1km"]], xy = TRUE))

# every mapped cell, blank where nothing was backfilled
spread_map <- function(on_fill) {
  left_join(select(observed_map, x, y), select(on_fill, x, y, spread), by = c("x", "y"))
}

cover_limits  <- range(observed_map$cover, bart_map$cover, xgb_map$cover, na.rm = TRUE)
spread_limits <- c(0, max(bart_fill$spread, xgb_fill$spread))

map_theme <- function() {
  list(theme_void(base_size = 12),
       theme(plot.background = element_rect(fill = surface, colour = NA),
             plot.title = element_text(colour = ink, size = 15, hjust = 0.5, margin = margin(b = 6)),
             legend.position = "right", legend.title = element_text(colour = ink, size = 11),
             legend.text = element_text(colour = ink_soft, size = 10),
             plot.margin = margin(4, 6, 4, 6)))
}

map_coord <- function(x_range, y_range) {
  coord_sf(crs = "EPSG:5072", xlim = x_range, ylim = y_range, expand = FALSE, datum = NA)
}

map_panel <- function(data, value, title, fill_scale) {
  ggplot() +
    geom_raster(data = data, aes(x, y, fill = .data[[value]])) +
    fill_scale +
    map_coord(win_x, win_y) +
    labs(title = title) +
    map_theme()
}

scale_bar <- function(x_range, y_range, km) {
  span <- diff(x_range)
  x0   <- x_range[1] + 0.06 * span
  x1   <- x0 + km * 1000
  y    <- y_range[1] + 0.25 * span
  tick <- 0.015 * span
  list(annotate("segment", x = c(x0, x0, x1), xend = c(x1, x0, x1),
                y = c(y, y - tick, y - tick), yend = c(y, y + tick, y + tick),
                colour = "white", linewidth = 0.8),
       annotate("text", x = (x0 + x1) / 2, y = y - 0.035 * span, label = sprintf("%d km", km),
                colour = "white", size = 5, vjust = 1))
}

satellite_panel <- function(x_range, y_range, zoom, km, title, marker = NULL) {
  box   <- ext(x_range[1], x_range[2], y_range[1], y_range[2])
  tiles <- maptiles::get_tiles(as.polygons(box, crs = crs(cov_s)), provider = "Esri.WorldImagery",
                               zoom = zoom, crop = TRUE, cachedir = file.path(out_dir, "tiles"))
  ggplot() +
    tidyterra::geom_spatraster_rgb(data = crop(tiles, box)) +
    tidyterra::geom_spatvector(data = poly, fill = NA, colour = "white", linewidth = 0.4) +
    marker +
    scale_bar(x_range, y_range, km) +
    map_coord(x_range, y_range) +
    labs(title = title) +
    map_theme()
}

cover_scale <- function(guide) {
  scale_fill_viridis_c(option = "D", limits = cover_limits, na.value = NA, guide = guide,
                       name = "Other-conifer\ncover (%)")
}
spread_scale <- function(guide) {
  scale_fill_viridis_c(option = "C", limits = spread_limits, na.value = rule, guide = guide,
                       name = sprintf("SD of the %d\ndraws (%% cover)", n_draws))
}
footprint_scale <- scale_fill_viridis_c(option = "A", direction = -1, na.value = NA,
                                        name = "Human\nfootprint")

window_box <- annotate("rect", xmin = win_x[1], xmax = win_x[2], ymin = win_y[1], ymax = win_y[2],
                       fill = NA, colour = "#eda100", linewidth = 0.8)

fig_maps <- wrap_plots(
  satellite_panel(context_x, context_y, 10, 25, sprintf("Subbasin %d, northeast BC", subbasin_index),
                  marker = window_box),
  map_panel(observed_map, "cover", sprintf("Observed cover, %d", year), cover_scale("none")),
  map_panel(bart_map, "cover", "BART backfill", cover_scale("none")),
  map_panel(xgb_map, "cover", sprintf("%s backfill", map_xgb), cover_scale("colourbar")),
  satellite_panel(win_x, win_y, 12, 10, sprintf("%d km window", window_km)),
  map_panel(footprint_map, "CanHF_1km", sprintf("Human footprint score, %d", year), footprint_scale),
  map_panel(spread_map(bart_fill), "spread", "BART uncertainty", spread_scale("none")),
  map_panel(spread_map(xgb_fill), "spread", sprintf("%s uncertainty", map_xgb), spread_scale("colourbar")),
  ncol = 4) +
  figure_caption(sprintf("The window holds %d of the %d backfilled pixels, the most of any %d km square. Grey: training pixels, not backfilled.\n%s",
                         best$held, nrow(fill), window_km,
                         maptiles::get_credit("Esri.WorldImagery")))

save_slide(fig_maps, "fig5_backfill_maps.png")
say("figures in %s", fig_dir)
