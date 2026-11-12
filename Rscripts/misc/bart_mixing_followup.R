# ---
# title:   Why long gbart chains drift, and R-hat by chain length
# author:  Mannfred Boehm
# created: November 5, 2026
# ---

# Run after bart_mixing_fixes.R, from the project root, with MIX (its work dir), RUN (the BPS_WORK of
# bart_production_settings.R) and BPS_OUT set. Prints predictor use and sigma by sweep window, then
# writes chain_length.csv.

w <- Sys.getenv("MIX")
use_of <- function(p, v) {
  fs <- list.files(w, sprintf("^%s_%s_c[0-9][.]rds$", p, v), full.names = TRUE)
  m <- sapply(fs, function(f) { u <- readRDS(f)$use; u / sum(u) })
  list(mean = rowMeans(m), chains = m)
}
for (p in c("pair_03", "pair_07")) for (v in c("base", "gbart_long50", "dbarts_long200")) {
  u <- use_of(p, v)$mean
  if (is.null(names(u))) names(u) <- paste0("V", seq_along(u))
  top <- sort(u, decreasing = TRUE)[1:5]
  xy <- sum(u[names(u) %in% c("x", "y")])
  cat(sprintf("%s %-15s top-5 share %.2f; x+y share %.3f; predictors for half the splits %d: %s\n", p, v, sum(top), xy,
              which(cumsum(sort(u, decreasing = TRUE)) >= 0.5)[1], paste(sprintf("%s %.2f", names(top), top), collapse = ", ")))
}

cat("\nsigma by sweep window (burn 300, then 20,000 sweeps), mean over chains\n")
for (p in c("pair_01", "pair_03", "pair_07")) for (v in c("fp_gbart50", "fp_dbarts200")) {
  fs <- list.files(w, sprintf("^%s_%s_c[0-9][.]rds$", p, v), full.names = TRUE)
  s <- sapply(fs, function(f) { z <- readRDS(f)$sigma
    if (length(z) > 5000) tail(z, 20000L)[seq(20L, 20000L, 20L)] else tail(z, 1000L) })
  m <- rowMeans(s)
  win <- list(1:35, 36:135, 136:485, 486:750, 751:1000)
  cat(sprintf("%s %-13s %s\n", p, v, paste(sprintf("%.4f", sapply(win, function(r) mean(m[r]))), collapse = "  ")))
}

suppressPackageStartupMessages({library(dplyr); library(purrr); library(posterior)})
options(width = 200)
w <- Sys.getenv("MIX"); run <- Sys.getenv("RUN"); out <- Sys.getenv("BPS_OUT")
rows <- list()
for (v in c("dbarts_long200", "dbarts_long50")) for (id in sprintf("pair_%02d", c(3, 7, 1, 6))) {
  fs <- list.files(w, sprintf("^%s_%s_c[0-9][.]rds$", id, v), full.names = TRUE)
  if (length(fs) < 4L) next
  meta <- readRDS(file.path(run, paste0(id, "_meta.rds")))
  hf <- seq_len(meta$n_fill_sample)
  ch <- map(fs, function(f) readRDS(f)$f[, hf, drop = FALSE])
  post_var <- apply(do.call(rbind, ch), 2, var)
  pix <- seq(1L, length(hf), length.out = 100L)
  for (burn in c(0L, 100L)) for (L in c(250L, 350L, 500L, 700L, 900L)) {
    r <- burn + seq_len(L)
    if (max(r) > 1000L) next
    rh <- map_dbl(pix, function(j) rhat(sapply(ch, function(x) x[r, j])))
    means <- map(ch, function(x) colMeans(x[r, , drop = FALSE]))
    pr <- combn(4, 2, simplify = FALSE)
    sq <- Reduce(`+`, map(pr, function(p) (means[[p[1]]] - means[[p[2]]])^2)) / length(pr)
    rows[[length(rows) + 1L]] <- tibble(pair = id, n_fit = meta$n_fit, variant = v, burn_in = 3000L + 20L * burn, kept_sweeps = 20L * L,
                                        rhat = median(rh), rhat_q95 = unname(quantile(rh, .95)), ess_per_chain = median(post_var / (sq / 2)))
  }
}
x <- bind_rows(rows)
write.csv(x, file.path(out, "chain_length.csv"), row.names = FALSE)
print(as.data.frame(x |> filter(burn_in == 3000L, variant == "dbarts_long200")), digits = 3)
print(as.data.frame(x |> filter(burn_in == 5000L, variant == "dbarts_long200")), digits = 3)
print(as.data.frame(x |> filter(burn_in == 3000L, variant == "dbarts_long50", kept_sweeps %in% c(10000L, 18000L))), digits = 3)
