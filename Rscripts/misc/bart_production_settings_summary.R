# ---
# title:   Summarise bart_production_settings.R
# author:  Mannfred Boehm
# created: November 4, 2026
# ---

# Reads the per-fit files in BPS_WORK and writes one CSV per question to
# data/derived_data/rds_files/bart_production_settings/. Row windows are of kept sweeps, so row 1
# is sweep nskip + 1.
#
# Effective draws come from BETWEEN chains. Two independent chains' means differ by what the
# stored draws fail to pin down, so ESS = posterior variance / (half the mean squared difference
# of chain means). Within-chain estimates cannot see a chain that has not yet left its mode.
#
# Run from the project root: BPS_WORK=<dir> Rscript Rscripts/misc/bart_production_settings_summary.R

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(posterior)
})

work_dir <- Sys.getenv("BPS_WORK")
out_dir  <- Sys.getenv("BPS_OUT", file.path("data", "derived_data", "rds_files", "bart_production_settings"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
nskip  <- 300L
ndpost <- 700L
n_rhat <- 100L     # backfill pixels for pixel-level R-hat

pair_ids <- sort(unique(sub("_meta\\.rds$", "", list.files(work_dir, "_meta\\.rds$"))))
rd <- function(id, what) {
  f <- file.path(work_dir, sprintf("%s_%s.rds", id, what))
  if (file.exists(f)) readRDS(f) else NULL
}

# interval coverage and CRPS of the predictive distribution: f plus noise, on the log1p scale
crps_sample <- function(x, y) {
  n <- length(x)
  xs <- sort(x)
  mean(abs(x - y)) - sum((2 * seq_len(n) - n - 1) * xs) / n^2
}
holdout_scores <- function(f, sigma, y) {
  pred <- colMeans(expm1(f))
  noisy <- f + rnorm(length(f), 0, rep(sigma, times = ncol(f)))
  z <- pmax(expm1(noisy), 0)
  lo <- apply(z, 2, quantile, 0.05)
  hi <- apply(z, 2, quantile, 0.95)
  flo <- apply(expm1(f), 2, quantile, 0.05)
  fhi <- apply(expm1(f), 2, quantile, 0.95)
  tibble(rmse = sqrt(mean((pred - y)^2)),
         r2 = 1 - sum((pred - y)^2) / sum((y - mean(y))^2),
         crps = mean(map_dbl(seq_along(y), function(i) crps_sample(z[, i], y[i]))),
         coverage90 = mean(y >= lo & y <= hi),
         width90 = mean(hi - lo),
         coverage90_f_only = mean(y >= flo & y <= fhi))
}

# ESS and SD of a stored set, against the arm's posterior. rows_by_chain: a list of row indices.
stored_set <- function(m, rows_by_chain, post_var) {
  means <- map(seq_along(m), function(c) colMeans(m[[c]][rows_by_chain[[c]], , drop = FALSE]))
  sds   <- map(seq_along(m), function(c) apply(m[[c]][rows_by_chain[[c]], , drop = FALSE], 2, sd))
  pairs <- combn(length(m), 2, simplify = FALSE)
  sq    <- Reduce(`+`, map(pairs, function(p) (means[[p[1]]] - means[[p[2]]])^2)) / length(pairs)
  ess   <- post_var / (sq / 2)
  list(ess = ess, sd_ratio = Reduce(`+`, sds) / length(sds) / sqrt(post_var), mean = means)
}

# The same for a set pooled over two chains: the three ways of splitting four chains in two
pooled_pair_ess <- function(means, post_var) {
  splits <- list(list(1:2, 3:4), list(c(1, 3), c(2, 4)), list(c(1, 4), 2:3))
  sq <- Reduce(`+`, map(splits, function(s) ((means[[s[[1]][1]]] + means[[s[[1]][2]]]) / 2 -
                                             (means[[s[[2]][1]]] + means[[s[[2]][2]]]) / 2)^2)) / 3
  post_var / (sq / 2)
}

set.seed(1)
arm_rows <- list()
config_rows <- list()
holdout_rows <- list()
psis_rows <- list()
fp_rows <- list()
size_rows <- list()
burn_rows <- list()

for (id in pair_ids) {
  meta <- rd(id, "meta")
  files <- list.files(work_dir, sprintf("^%s_t[0-9]+_c[0-9]+\\.rds$", id))
  arms <- sort(unique(as.integer(sub(".*_t([0-9]+)_c.*", "\\1", files))))
  nh <- meta$n_fill_sample
  y  <- meta$y_scored
  hf <- seq_len(nh)
  ho <- nh + seq_along(y)
  key <- tibble(pair = id, subbasin = meta$pair$subbasin, covariate = meta$pair$covariate,
                n_fit = meta$n_fit, p = meta$n_predictors)

  for (ntree in arms) {
    chains <- map(sort(files[grepl(sprintf("_t%03d_c", ntree), files)]),
                  function(f) readRDS(file.path(work_dir, f)))
    n_c <- length(chains)
    if (n_c < 2L) next
    K <- min(map_int(chains, function(ch) nrow(ch$f)))
    late <- if (K >= 1500L) 1001:K else (K %/% 2L + 1L):K
    prod <- seq_len(min(ndpost, K))
    f_hf  <- map(chains, function(ch) ch$f[seq_len(K), hf, drop = FALSE])
    agg   <- map(chains, function(ch) as.matrix(rowMeans(expm1(ch$f[seq_len(K), hf, drop = FALSE]))))
    sig   <- map(chains, function(ch) as.matrix(ch$sigma[nskip + seq_len(K)]))
    secs  <- mean(map_dbl(chains, "seconds"))
    sweep_s <- secs / (nskip + K)

    # the posterior, from every chain's late window, so between-chain spread counts
    pool <- function(m) do.call(rbind, map(m, function(x) x[late, , drop = FALSE]))
    var_hf  <- apply(pool(f_hf), 2, var)
    var_agg <- apply(pool(agg), 2, var)
    var_sig <- apply(pool(sig), 2, var)
    mean_hf <- colMeans(pool(f_hf))

    # R-hat and within-pooled ESS
    pix <- sort(sample(nh, min(n_rhat, nh)))
    stat <- function(m, j, rows, fun) fun(sapply(m, function(x) x[rows, j]))
    rhat_pix <- function(rows) map_dbl(pix, function(j) stat(f_hf, j, rows, rhat))
    ess_pix  <- map_dbl(pix, function(j) stat(f_hf, j, late, ess_bulk))

    # burn-in: sigma's approach. First sweep from which the mean of the chains, in 50-sweep
    # windows, stays within 0.5 posterior SD of the late mean for 100 sweeps
    trace <- rowMeans(sapply(chains, function(ch) ch$sigma[seq_len(nskip + K)]))
    ma <- stats::filter(trace, rep(1 / 50, 50), sides = 1)
    ref <- mean(pool(sig))
    ok <- abs(ma - ref) < 0.5 * sqrt(var_sig)
    runs <- rle(as.vector(ok))
    start <- cumsum(c(1, head(runs$lengths, -1)))
    entry <- start[which(runs$values %in% TRUE & runs$lengths >= 100)[1]]

    # the shift of a window from the late mean, signed, in posterior SDs and, for the aggregate,
    # as a share of the backfill's contrast with the observed footprint
    fill_mean <- mean(pool(agg))
    contrast  <- fill_mean - mean(meta$obs_fill, na.rm = TRUE)
    shift <- function(rows) {
      sig_shift <- mean(map_dbl(sig, function(x) mean(x[rows, ]))) - ref
      agg_shift <- mean(map_dbl(agg, function(x) mean(x[rows, ]))) - fill_mean
      pix_shift <- rowMeans(sapply(f_hf, function(x) colMeans(x[rows, , drop = FALSE]))) - mean_hf
      tibble(sigma_shift_sd = sig_shift / sqrt(var_sig), agg_shift_sd = agg_shift / sqrt(var_agg),
             agg_shift_pct_contrast = 100 * agg_shift / abs(contrast),
             pixel_abs_shift_sd = median(abs(pix_shift) / sqrt(var_hf)))
    }

    arm_rows[[length(arm_rows) + 1L]] <- bind_cols(key, tibble(
      ntree = ntree, n_chains = n_c, kept = K, seconds_per_sweep = sweep_s, seconds_08B_length = sweep_s * (nskip + ndpost),
      sigma_late = ref, sigma_rel_sd_y = ref / meta$sd_y_log,
      rhat_sigma_prod = stat(sig, 1, prod, rhat), rhat_sigma_late = stat(sig, 1, late, rhat),
      rhat_agg_prod = stat(agg, 1, prod, rhat), rhat_agg_late = stat(agg, 1, late, rhat),
      rhat_pixel_prod_median = median(rhat_pix(prod)), rhat_pixel_late_median = median(rhat_pix(late)),
      rhat_pixel_late_q95 = unname(quantile(rhat_pix(late), 0.95)),
      ess_pixel_late_median = median(ess_pix),
      sweeps_per_effective_draw = length(late) * n_c / median(ess_pix),
      entry_sweep_sigma = entry,
      fill_mean = fill_mean, obs_mean = mean(meta$obs_fill, na.rm = TRUE), contrast = contrast,
      n_fill = meta$n_fill))

    burn_rows[[length(burn_rows) + 1L]] <- bind_cols(key, tibble(ntree = ntree), bind_rows(
      mutate(shift(prod), window = "burn 300, 700 kept (08B)"),
      mutate(shift(intersect(301:1000, seq_len(K))), window = "sweeps 601-1300"),
      mutate(shift(intersect(701:1000, seq_len(K))), window = "sweeps 1001-1300")))

    # stored sets: 100 draws from one chain, as each option would store them
    options <- list(
      "burn 300, random 100 of 700 (08B)" = function() sort(sample(ndpost, 100L)),
      "burn 300, every 7th of 700"        = function() seq(7L, ndpost, by = 7L),
      "burn 1000, every 7th of 700"       = function() 700L + seq(7L, ndpost, by = 7L),
      "burn 300, every 17th of 1700"      = function() seq(17L, 1700L, by = 17L))
    for (opt in names(options)) {
      rows <- options[[opt]]()
      if (max(rows) > K) next
      by_chain <- rep(list(rows), n_c)
      if (grepl("random", opt)) by_chain <- map(seq_len(n_c), function(i) options[[opt]]())
      px <- stored_set(f_hf, by_chain, var_hf)
      ag <- stored_set(agg, by_chain, var_agg)
      config_rows[[length(config_rows) + 1L]] <- bind_cols(key, tibble(
        ntree = ntree, option = opt, sweeps = nskip + max(rows), seconds = sweep_s * (nskip + max(rows)),
        ess_pixel = median(px$ess), ess_agg = ag$ess, sd_ratio = median(px$sd_ratio),
        rmse_vs_post_sd = median(sqrt(rowMeans(sapply(px$mean, function(m) (m - mean_hf)^2)) / var_hf))))
    }

    if (n_c == 4L) {
      rows <- rep(list(seq(7L, ndpost, by = 7L)), n_c)
      px <- stored_set(f_hf, rows, var_hf)
      ag <- stored_set(agg, rows, var_agg)
      config_rows[[length(config_rows) + 1L]] <- bind_cols(key, tibble(
        ntree = ntree, option = "2 chains, each burn 300 + 700 kept", sweeps = 2L * (nskip + ndpost),
        seconds = sweep_s * 2 * (nskip + ndpost), ess_pixel = median(pooled_pair_ess(px$mean, var_hf)),
        ess_agg = pooled_pair_ess(ag$mean, var_agg), sd_ratio = NA_real_, rmse_vs_post_sd = NA_real_))
    }

    # held-out skill on 08B's 10% holdout, one chain at a time, by window
    for (w in c("08B window", "late window")) {
      rows <- if (w == "08B window") prod else late
      s <- map_dfr(chains, function(ch) holdout_scores(ch$f[rows, ho, drop = FALSE], ch$sigma[nskip + rows], y))
      holdout_rows[[length(holdout_rows) + 1L]] <- bind_cols(key, tibble(ntree = ntree, window = w), summarise(s, across(everything(), mean)))
    }
    for (k in c(2L, 4L)) {
      if (k > n_c) next
      sets <- if (k == 2L) list(1:2, 3:4) else list(1:4)
      s <- map_dfr(sets, function(set) {
        pooled <- do.call(rbind, map(chains[set], function(ch) ch$f[prod, ho, drop = FALSE]))
        holdout_scores(pooled, unlist(map(chains[set], function(ch) ch$sigma[nskip + prod])), y)
      })
      holdout_rows[[length(holdout_rows) + 1L]] <- bind_cols(key, tibble(ntree = ntree, window = sprintf("%d chains pooled, 08B length each", k)),
                                                             summarise(s, across(everything(), mean)))
    }

    # PSIS-LOO against this chain's real holdout
    ps <- rd(id, sprintf("t%03d_psis", ntree))
    if (!is.null(ps)) {
      real <- holdout_scores(chains[[1]]$f[prod, ho, drop = FALSE], chains[[1]]$sigma[nskip + prod], y)
      psis_rows[[length(psis_rows) + 1L]] <- bind_cols(key, tibble(
        ntree = ntree, rmse_insample = sqrt(mean((ps$pred_insample - ps$y)^2)),
        rmse_psis_loo = sqrt(mean((ps$pred_loo - ps$y)^2)), rmse_holdout = real$rmse,
        share_k_over_0.7 = mean(ps$pareto_k > 0.7), median_k = median(ps$pareto_k),
        psis_seconds = ps$seconds))
    }

    # footprint-shaped holdout at 08B's length
    fp <- rd(id, sprintf("fp_t%03d", ntree))
    if (!is.null(fp)) {
      n_o <- fp$n_observed
      real <- holdout_scores(chains[[1]]$f[prod, ho, drop = FALSE], chains[[1]]$sigma[nskip + prod], y)
      fp_rows[[length(fp_rows) + 1L]] <- bind_cols(key, tibble(ntree = ntree, n_footprint_test = n_o, n_held = meta$n_fp_held),
        rename_with(holdout_scores(fp$f[, seq_len(n_o), drop = FALSE], fp$sigma, fp$y), function(x) paste0(x, "_fp_observed")),
        rename_with(holdout_scores(fp$f[, n_o + seq_len(n_o), drop = FALSE], fp$sigma, fp$y), function(x) paste0(x, "_fp_cascaded")),
        rename_with(real, function(x) paste0(x, "_random")))
    }
  }

  # 90% of the fitting pixels against all of them, same holdout, same length
  chains50 <- map(sort(files[grepl("_t050_c", files)]), function(f) readRDS(file.path(work_dir, f)))
  for (chain in 1:2) {
    z <- rd(id, sprintf("size_c%d", chain))
    if (is.null(z) || length(chains50) < chain) next
    full <- chains50[[chain]]
    size_rows[[length(size_rows) + 1L]] <- bind_cols(key, tibble(
      chain = chain,
      rmse_90 = holdout_scores(z$small[, ho, drop = FALSE], z$sigma_small, y)$rmse,
      rmse_100 = holdout_scores(full$f[seq_len(ndpost), ho, drop = FALSE], full$sigma[nskip + seq_len(ndpost)], y)$rmse))
  }
  rm(chains50)
  message(id, " summarised")
}

results <- list(arms = bind_rows(arm_rows), burn_in = bind_rows(burn_rows), stored_draws = bind_rows(config_rows),
                holdout = bind_rows(holdout_rows), psis = bind_rows(psis_rows),
                footprint = bind_rows(fp_rows), size = bind_rows(size_rows))
iwalk(results, function(x, name) write.csv(x, file.path(out_dir, paste0(name, ".csv")), row.names = FALSE))

options(width = 200, pillar.width = 200)
show <- function(x, title) { cat("\n==", title, "\n"); print(as.data.frame(x), digits = 3) }

show(select(results$arms, pair, covariate, n_fit, p, ntree, seconds_08B_length, sigma_late, rhat_sigma_late, rhat_agg_late,
            rhat_pixel_prod_median, rhat_pixel_late_median, rhat_pixel_late_q95, sweeps_per_effective_draw, entry_sweep_sigma),
     "per pair and arm")
show(results$arms |> group_by(ntree) |> summarise(
       pairs = n(), rhat_sigma_late = median(rhat_sigma_late), rhat_agg_late = median(rhat_agg_late),
       rhat_pixel_prod = median(rhat_pixel_prod_median), rhat_pixel_late = median(rhat_pixel_late_median),
       rhat_pixel_q95 = median(rhat_pixel_late_q95), share_pairs_rhat_under_1.05 = mean(rhat_pixel_late_median < 1.05),
       sweeps_per_effective = median(sweeps_per_effective_draw), entry_sweep = median(entry_sweep_sigma, na.rm = TRUE),
       entry_max = max(entry_sweep_sigma, na.rm = TRUE), never_entered = sum(is.na(entry_sweep_sigma)),
       seconds_08B = median(seconds_08B_length)), "by arm (median over pairs)")
show(results$burn_in |> group_by(ntree, window) |> summarise(
       pairs = n(), sigma_shift_sd = median(sigma_shift_sd), agg_shift_sd = median(agg_shift_sd),
       agg_abs_shift_pct_contrast = median(abs(agg_shift_pct_contrast)), agg_shift_pct_contrast_signed = median(agg_shift_pct_contrast),
       pixel_abs_shift_sd = median(pixel_abs_shift_sd), .groups = "drop"), "burn-in windows vs late window")
show(results$stored_draws |> group_by(ntree, option) |> summarise(
       pairs = n(), sweeps = first(sweeps), seconds = median(seconds), ess_pixel = median(ess_pixel), ess_pixel_q10 = quantile(ess_pixel, .1),
       ess_agg = median(ess_agg), sd_ratio = median(sd_ratio), rmse_vs_post_sd = median(rmse_vs_post_sd), .groups = "drop"),
     "stored 100 draws by option")
show(results$holdout |> group_by(ntree, window) |> summarise(
       pairs = n(), across(c(rmse, r2, crps, coverage90, width90, coverage90_f_only), median), .groups = "drop"), "holdout")
if (nrow(results$psis)) show(results$psis, "PSIS-LOO")
if (nrow(results$footprint)) show(select(results$footprint, pair, covariate, n_fit, ntree, n_footprint_test, rmse_random, rmse_fp_observed, rmse_fp_cascaded,
                                         r2_random, r2_fp_observed, r2_fp_cascaded, coverage90_random, coverage90_fp_observed, coverage90_fp_cascaded), "footprint-shaped holdout")
if (nrow(results$size)) show(results$size |> group_by(pair, covariate, n_fit) |> summarise(rmse_90 = mean(rmse_90), rmse_100 = mean(rmse_100), .groups = "drop") |> mutate(change_pct = 100 * (rmse_100 / rmse_90 - 1)), "90% vs 100% of the fitting pixels")
