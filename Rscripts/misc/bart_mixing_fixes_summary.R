# ---
# title:   Summarise bart_mixing_fixes.R
# author:  Mannfred Boehm
# created: November 4, 2026
# ---

# One row per pair and variant, against bart_production_settings.R's own chains (t050, t200) on
# the same fits. Rows are kept sweeps: 1-700 is 08B's window, the last 700 the late window.
# Effective draws are between chains, as in bart_production_settings_summary.R.
#
# Run from the project root: BPS_WORK=<dir> MIX_WORK=<dir> Rscript Rscripts/misc/bart_mixing_fixes_summary.R

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(tidyr)
  library(posterior)
})

work_dir <- Sys.getenv("BPS_WORK")
mix_dir  <- Sys.getenv("MIX_WORK")
out_dir  <- Sys.getenv("BPS_OUT", file.path("data", "derived_data", "rds_files", "bart_production_settings"))
options(width = 200, pillar.width = 200)

crps_sample <- function(x, y) {
  n <- length(x)
  xs <- sort(x)
  mean(abs(x - y)) - sum((2 * seq_len(n) - n - 1) * xs) / n^2
}
holdout_scores <- function(f, sigma, y) {
  pred <- colMeans(expm1(f))
  z <- pmax(expm1(f + rnorm(length(f), 0, rep(sigma, times = ncol(f)))), 0)
  flo <- apply(expm1(f), 2, quantile, 0.05)
  fhi <- apply(expm1(f), 2, quantile, 0.95)
  c(rmse = sqrt(mean((pred - y)^2)),
    crps = mean(map_dbl(seq_along(y), function(i) crps_sample(z[, i], y[i]))),
    coverage_f = mean(y >= flo & y <= fhi))
}

# per-chain ESS: posterior variance over half the mean squared difference of chain means
between_ess <- function(m, rows, post_var) {
  means <- map(m, function(x) colMeans(x[rows, , drop = FALSE]))
  pairs <- combn(length(m), 2, simplify = FALSE)
  sq <- Reduce(`+`, map(pairs, function(p) (means[[p[1]]] - means[[p[2]]])^2)) / length(pairs)
  median(post_var / (sq / 2))
}

summarise_chains <- function(chains, meta, id, variant) {
  nh <- meta$n_fill_sample
  hf <- seq_len(nh)
  ho <- nh + seq_along(meta$y_scored)
  y  <- meta$y_scored
  K  <- min(map_int(chains, function(ch) nrow(ch$f)))
  rows08 <- seq_len(700L)
  rowsL  <- (K - 699L):K
  f_hf <- map(chains, function(ch) ch$f[seq_len(K), hf, drop = FALSE])
  pool <- do.call(rbind, map(f_hf, function(x) x[rowsL, , drop = FALSE]))
  post_var <- apply(pool, 2, var)
  pix <- seq(1L, nh, length.out = min(100L, nh))
  rh <- function(rows) map_dbl(pix, function(j) rhat(sapply(f_hf, function(x) x[rows, j])))
  agg <- map(chains, function(ch) rowMeans(expm1(ch$f[seq_len(K), hf, drop = FALSE])))
  sig_rows <- function(ch, rows) { s <- ch$sigma; n <- length(s) - nrow(ch$f); s[n + rows] }

  one <- map(chains, function(ch) holdout_scores(ch$f[rows08, ho, drop = FALSE], sig_rows(ch, rows08), y))
  all4 <- holdout_scores(do.call(rbind, map(chains, function(ch) ch$f[rows08, ho, drop = FALSE])),
                         unlist(map(chains, sig_rows, rows08)), y)
  use <- map(chains, "use")
  use_cor <- NA_real_
  if (!any(map_lgl(use, is.null))) {
    pairs <- combn(length(use), 2, simplify = FALSE)
    use_cor <- mean(map_dbl(pairs, function(p) cor(use[[p[1]]], use[[p[2]]], method = "spearman")))
  }
  tibble(pair = id, variant = variant,
         rhat_late = median(rh(rowsL)), rhat_08B = median(rh(rows08)),
         rhat_agg = rhat(sapply(agg, function(a) a[rowsL])),
         ess_per_chain = between_ess(f_hf, rows08, post_var),
         sd_ratio = median(sqrt(apply(f_hf[[1]][rows08, , drop = FALSE], 2, var) / post_var)),
         within_ess = median(map_dbl(pix, function(j) ess_bulk(sapply(f_hf, function(x) x[rowsL, j])))),
         rmse_1 = mean(map_dbl(one, "rmse")), crps_1 = mean(map_dbl(one, "crps")),
         cov_f_1 = mean(map_dbl(one, "coverage_f")),
         rmse_4 = all4[["rmse"]], crps_4 = all4[["crps"]], cov_f_4 = all4[["coverage_f"]],
         sigma_late = mean(map_dbl(chains, function(ch) mean(tail(ch$sigma, 700L)))),
         fill_mean = mean(map_dbl(agg, function(a) mean(a[rowsL]))),
         ess_full = between_ess(f_hf, seq_len(K), post_var),
         seconds_per_1000 = 1000 * mean(map_dbl(chains, "seconds")) /
           (if (is.null(chains[[1]]$sweeps)) 300 + K else chains[[1]]$sweeps),
         use_cor = use_cor)
}

set.seed(1)
rows <- list()
files <- list.files(mix_dir, "^pair_[0-9]+_.*_c[0-9]+[.]rds$")
all_keys <- unique(sub("_c[0-9]+[.]rds$", "", files))
curve_keys <- all_keys[grepl("_fp_", all_keys)]
keys <- setdiff(all_keys, curve_keys)
ids <- sort(unique(sub("^(pair_[0-9]+)_.*", "\\1", keys)))
for (id in ids) {
  meta <- readRDS(file.path(work_dir, paste0(id, "_meta.rds")))
  for (t in c(50L, 200L)) {
    fs <- list.files(work_dir, sprintf("^%s_t%03d_c[0-9][.]rds$", id, t), full.names = TRUE)
    rows[[length(rows) + 1L]] <- summarise_chains(map(fs, readRDS), meta, id, sprintf("08B %d trees", t))
  }
  for (key in keys[startsWith(keys, paste0(id, "_"))]) {
    variant <- sub(paste0("^", id, "_"), "", key)
    fs <- list.files(mix_dir, sprintf("^%s_c[0-9]+[.]rds$", key), full.names = TRUE)
    if (length(fs) < 2L) next
    rows[[length(rows) + 1L]] <- summarise_chains(map(fs, readRDS), meta, id, variant)
  }
}
res <- bind_rows(rows)
write.csv(res, file.path(out_dir, "mixing_fixes.csv"), row.names = FALSE)

cat("\n== per pair and variant\n")
print(as.data.frame(select(res, pair, variant, rhat_late, rhat_08B, rhat_agg, ess_per_chain, ess_full, sd_ratio, within_ess,
                           rmse_1, rmse_4, cov_f_4, sigma_late, fill_mean, seconds_per_1000, use_cor)), digits = 3)

ref <- filter(res, variant == "08B 50 trees") |> select(pair, rmse0 = rmse_1, crps0 = crps_1, rmse40 = rmse_4, s0 = seconds_per_1000, fill0 = fill_mean)
cat("\n== by variant: medians over pairs; ratios against 08B's 50 trees, one chain\n")
print(as.data.frame(res |> left_join(ref, by = "pair") |> group_by(variant) |> summarise(
  pairs = n(), rhat_late = median(rhat_late), rhat_worst = max(rhat_late), rhat_agg = median(rhat_agg),
  ess_per_chain = median(ess_per_chain), ess_full = median(ess_full), within_ess = median(within_ess), sd_ratio = median(sd_ratio),
  rmse_vs_base = median(rmse_1 / rmse0), rmse_worst = max(rmse_1 / rmse0), crps_vs_base = median(crps_1 / crps0),
  rmse4_vs_base4 = median(rmse_4 / rmse40), cov_f_1 = median(cov_f_1), cov_f_4 = median(cov_f_4),
  time_vs_base = median(seconds_per_1000 / s0), use_cor = median(use_cor)) |> arrange(rhat_late)), digits = 3)

# learning curves on the footprint-shaped holdout: 08B's burn-in, then 1,000 rows, one per 20 sweeps
windows <- list(`sweeps 301-1,000` = 1:35, `1,001-3,000` = 36:135, `3,001-10,000` = 136:485, `10,001-20,300` = 486:1000)
curve_rows <- list()
for (key in curve_keys) {
  id <- sub("^(pair_[0-9]+)_.*", "\\1", key)
  fp <- readRDS(file.path(mix_dir, paste0(id, "_design.rds")))$fp
  chains <- map(list.files(mix_dir, sprintf("^%s_c[0-9]+[.]rds$", key), full.names = TRUE), readRDS)
  n_o <- fp$n_observed
  y <- fp$y_test
  for (w in names(windows)) {
    f <- do.call(rbind, map(chains, function(ch) ch$f[windows[[w]], , drop = FALSE]))
    score <- function(cols) {
      z <- expm1(f[, cols, drop = FALSE])
      pred <- colMeans(z)
      c(rmse = sqrt(mean((pred - y)^2)), r2 = 1 - sum((pred - y)^2) / sum((y - mean(y))^2),
        coverage_f = mean(y >= apply(z, 2, quantile, 0.05) & y <= apply(z, 2, quantile, 0.95)))
    }
    o <- score(seq_len(n_o))
    k <- score(n_o + seq_len(n_o))
    curve_rows[[length(curve_rows) + 1L]] <- tibble(
      pair = id, sampler = sub("^pair_[0-9]+_fp_", "", key), window = w, n_chains = length(chains),
      rmse_observed = o[["rmse"]], r2_observed = o[["r2"]], cov_f_observed = o[["coverage_f"]],
      rmse_cascaded = k[["rmse"]], r2_cascaded = k[["r2"]], cov_f_cascaded = k[["coverage_f"]])
  }
}
if (length(curve_rows)) {
  curves <- bind_rows(curve_rows) |> mutate(window = factor(window, names(windows)))
  write.csv(curves, file.path(out_dir, "footprint_learning_curves.csv"), row.names = FALSE)
  cat("\n== footprint-shaped holdout by sweep window\n")
  print(as.data.frame(arrange(curves, pair, sampler, window)), digits = 3)
  cat("\n== RMSE against gbart at 50 trees in its first window, about 08B's own (median over pairs)\n")
  ref <- curves |> filter(sampler == "gbart50", window == names(windows)[1]) |>
    select(pair, ref_observed = rmse_observed, ref_cascaded = rmse_cascaded)
  print(as.data.frame(curves |> left_join(ref, by = "pair") |>
    mutate(ratio_observed = rmse_observed / ref_observed, ratio_cascaded = rmse_cascaded / ref_cascaded) |>
    group_by(sampler, window) |> summarise(
      pairs = n(), observed = median(ratio_observed), observed_worst = max(ratio_observed),
      cascaded = median(ratio_cascaded), cascaded_worst = max(ratio_cascaded),
      cov_f_observed = median(cov_f_observed), cov_f_cascaded = median(cov_f_cascaded), .groups = "drop")), digits = 3)
}
