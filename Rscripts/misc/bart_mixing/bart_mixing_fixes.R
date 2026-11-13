# ---
# title:   What would make BART mix on production's fits?
# author:  Mannfred Boehm
# created: November 4, 2026
# ---

# bart_production_settings.R found pixel R-hat of 1.9 at 50 trees and 1.7 at 200. This tries the
# likely causes one at a time on the same fits: same training rows, backfill pixels, holdout,
# chain length and windows, so every variant is comparable with that run's chains. A variant
# that lowers R-hat by underfitting is no use, so _summary.R reports held-out skill beside it.
#
# Variants (50 trees unless named):
#   base         08B's gbart() call again, recording which predictors each chain uses
#   sparse_off   sparse = FALSE. Does the sparse prior lock chains onto different predictors?
#   decorrelated one predictor per cluster at |Spearman r| > 0.9
#   no_xy        without the scaled coordinates
#   n1000        1,000 random training pixels: are near-duplicate neighbours the cause?
#   dbarts       dbarts::bart2(), which adds change and swap tree moves to birth and death
#   gfr_warm     stochtree::bart(), each chain started from its own 10 grow-from-root sweeps
#   dbarts_long50, dbarts_long200, gbart_long50   chains of 23,000 sweeps
#   dbarts_dec200  dbarts_long200 on the decorrelated predictors
#
# Run from the project root, one job at a time. BPS_WORK holds bart_production_settings.R's files:
#   BPS_WORK=<dir> MIX_WORK=<dir> MINUTES=60 Rscript Rscripts/misc/bart_mixing/bart_mixing_fixes.R
# Optional: VARIANTS=a,b  PAIRS=3,7  N_CHAINS=4

lines <- readLines(file.path("Rscripts", "misc", "bart_production_settings.R"))
eval(parse(text = lines[seq_len(grep("^# 5[.] run", lines) - 1L)]))

mix_dir  <- Sys.getenv("MIX_WORK")
if (mix_dir == "") stop("set MIX_WORK")
dir.create(mix_dir, recursive = TRUE, showWarnings = FALSE)
stop_at  <- Sys.time() + 60 * as.numeric(Sys.getenv("MINUTES", "60"))
pair_ids <- as.integer(strsplit(Sys.getenv("PAIRS", "3,7,1,6"), ",")[[1]])
wanted   <- Sys.getenv("VARIANTS", "sparse_off,dbarts,gfr_warm,base,decorrelated,no_xy,n1000")
wanted   <- strsplit(wanted, ",")[[1]]


# 1. each pair's design, as run_pair() built it ----------------------------------------------

pair_design <- function(p) {
  file <- file.path(mix_dir, sprintf("pair_%02d_design.rds", p$pair))
  if (file.exists(file)) return(readRDS(file))
  id     <- sprintf("pair_%02d", p$pair)
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
  fp_held <- footprint_holdout(train[idx, ], fill, frames$res, seed0 + p$pair)
  fp_test <- if (length(fp_held) > n_test) sample(fp_held, n_test) else fp_held
  cascaded <- readRDS(file.path(work_dir, paste0(id, "_cascade.rds")))
  fill_cascaded <- cascaded$target[seq_len(nrow(fill_sample)), , drop = FALSE]
  fp_cascaded   <- cascaded$target[nrow(fill_sample) + seq_along(fp_test), , drop = FALSE]
  design <- prep_bart_design(train[idx, predictors, drop = FALSE], fill_cascaded, logp)
  x_all  <- as.matrix(design$train)
  y_all  <- as.numeric(train[[b]][idx])
  set.seed(abs(as.integer(sprintf("%d%03d", p$subbasin, which(biotic_cols == b)))))
  holdout <- sample(seq_along(y_all), size = round(0.10 * length(y_all)))
  scored  <- if (length(holdout) > n_test) sample(holdout, n_test) else holdout
  meta    <- readRDS(file.path(work_dir, paste0(id, "_meta.rds")))
  stopifnot(isTRUE(all.equal(y_all[scored], meta$y_scored)))

  # the footprint-shaped holdout's design
  rest_rows <- setdiff(seq_along(idx), fp_held)
  test_obs  <- train[idx[fp_test], ]
  test_obs[, intersect(footprint_cols, names(test_obs))] <- 0
  both <- rbind(test_obs[, names(train)], fp_cascaded[, names(train)])
  design_fp <- prep_bart_design(train[idx[rest_rows], predictors, drop = FALSE], both, logp)

  out <- list(x = x_all[-holdout, , drop = FALSE], y = log1p(y_all[-holdout]),
              x_test = rbind(as.matrix(design$backfill), x_all[scored, , drop = FALSE]),
              fp = list(x = as.matrix(design_fp$train), y = log1p(y_all[rest_rows]),
                        x_test = as.matrix(design_fp$backfill), y_test = y_all[fp_test],
                        n_observed = length(fp_test)))
  saveRDS(out, file, compress = FALSE)
  out
}


# 2. the samplers ----------------------------------------------------------------------------

late <- function(m) m[max(1L, nrow(m) - 699L):nrow(m), , drop = FALSE]

run_gbart <- function(x, y, x_test, seed, ntree = 50L, sparse = TRUE, k = 3, power = 2) {
  set.seed(seed)
  t0 <- proc.time()[[3]]
  capture.output(fit <- suppressWarnings(BART::gbart(
    x.train = x, y.train = y, x.test = x_test, type = "wbart", k = k, power = power, ntree = ntree,
    ndpost = n_kept, nskip = nskip, sparse = sparse, sigest = sd(y), sigdf = 3, printevery = 1000000L)))
  list(f = fit$yhat.test, sigma = fit$sigma, seconds = proc.time()[[3]] - t0,
       use = colMeans(late(fit$varcount)))
}

run_dbarts <- function(x, y, x_test, seed, ntree = 50L) {
  t0 <- proc.time()[[3]]
  fit <- dbarts::bart2(x, y, test = x_test, n.trees = ntree, n.samples = n_kept, n.burn = nskip,
                       n.chains = 1L, n.threads = 1L, k = 3, power = 2, base = 0.95, sigest = sd(y),
                       sigdf = 3, sigquant = 0.9, seed = seed, verbose = FALSE, keepTrees = FALSE,
                       keepTrainingFits = FALSE, combineChains = TRUE)
  list(f = fit$yhat.test, sigma = c(fit$first.sigma, fit$sigma), seconds = proc.time()[[3]] - t0,
       use = colMeans(late(fit$varcount)))
}

# seeds differ and num_chains = 1, so chains share no warm start and R-hat stays honest
run_gfr <- function(x, y, x_test, seed, ntree = 50L, num_gfr = 10L) {
  t0 <- proc.time()[[3]]
  fit <- stochtree::bart(X_train = x, y_train = y, X_test = x_test, num_gfr = num_gfr, num_burnin = 0L,
                         num_mcmc = n_kept,
                         general_params = list(num_chains = 1L, random_seed = seed, verbose = FALSE),
                         mean_forest_params = list(num_trees = ntree))
  list(f = t(fit$y_hat_test), sigma = c(rep(NA_real_, nskip), sqrt(fit$sigma2_global_samples)),
       seconds = proc.time()[[3]] - t0, use = NULL)
}

# Long chains: burn-in 3,000, then 20,000 sweeps storing every 20th, so the last 700 rows span
# 14,000 sweeps. Is R-hat high because 700 sweeps is about one autocorrelation time?
long <- list(burn = 3000L, kept = 1000L, thin = 20L)

# A learning curve on the footprint-shaped holdout: 08B's burn-in, then 20,000 sweeps. A random
# holdout's pixels sit beside training pixels, so it rewards fine local detail a backfill cannot use
curve <- list(burn = 300L, kept = 1000L, thin = 20L)

# bart2() counts n.samples and n.burn in sweeps and stores n.samples / n.thin rows, unlike gbart()'s
# ndpost, which counts stored rows
run_dbarts_long <- function(x, y, x_test, seed, ntree, spec = long) {
  t0 <- proc.time()[[3]]
  fit <- dbarts::bart2(x, y, test = x_test, n.trees = ntree, n.samples = spec$kept * spec$thin, n.burn = spec$burn,
                       n.thin = spec$thin, n.chains = 1L, n.threads = 1L, k = 3, power = 2, base = 0.95,
                       sigest = sd(y), sigdf = 3, sigquant = 0.9, seed = seed, verbose = FALSE,
                       keepTrees = FALSE, keepTrainingFits = FALSE, combineChains = TRUE)
  stopifnot(nrow(fit$yhat.test) == spec$kept)
  list(f = fit$yhat.test, sigma = c(fit$first.sigma, fit$sigma), seconds = proc.time()[[3]] - t0,
       use = colMeans(late(fit$varcount)), sweeps = spec$burn + spec$kept * spec$thin)
}

run_gbart_long <- function(x, y, x_test, seed, ntree, spec = long) {
  set.seed(seed)
  t0 <- proc.time()[[3]]
  capture.output(fit <- suppressWarnings(BART::gbart(
    x.train = x, y.train = y, x.test = x_test, type = "wbart", k = 3, ntree = ntree, ndpost = spec$kept,
    nskip = spec$burn, keepevery = spec$thin, sparse = TRUE, sigest = sd(y), sigdf = 3, printevery = 1000000L)))
  stopifnot(nrow(fit$yhat.test) == spec$kept)
  list(f = fit$yhat.test, sigma = fit$sigma, seconds = proc.time()[[3]] - t0,
       use = colMeans(late(fit$varcount)), sweeps = spec$burn + spec$kept * spec$thin)
}

thin_predictors <- function(x, r = 0.9) {
  x <- x[, apply(x, 2, sd) > 0, drop = FALSE]
  clusters <- cutree(hclust(as.dist(1 - abs(cor(x, method = "spearman"))), "complete"), h = 1 - r)
  colnames(x)[!duplicated(clusters)]
}

variants <- list(
  base         = function(d, s) run_gbart(d$x, d$y, d$x_test, s),
  sparse_off   = function(d, s) run_gbart(d$x, d$y, d$x_test, s, sparse = FALSE),
  decorrelated = function(d, s) { keep <- thin_predictors(d$x)
                                  run_gbart(d$x[, keep, drop = FALSE], d$y, d$x_test[, keep, drop = FALSE], s) },
  no_xy        = function(d, s) { keep <- setdiff(colnames(d$x), c("x", "y"))
                                  run_gbart(d$x[, keep, drop = FALSE], d$y, d$x_test[, keep, drop = FALSE], s) },
  n1000        = function(d, s) { set.seed(s %% 1000L); rows <- sample(nrow(d$x), min(1000L, nrow(d$x)))
                                  run_gbart(d$x[rows, , drop = FALSE], d$y[rows], d$x_test, s) },
  dbarts       = function(d, s) run_dbarts(d$x, d$y, d$x_test, s),
  gfr_warm     = function(d, s) run_gfr(d$x, d$y, d$x_test, s),
  dbarts_long50  = function(d, s) run_dbarts_long(d$x, d$y, d$x_test, s, 50L),
  dbarts_long200 = function(d, s) run_dbarts_long(d$x, d$y, d$x_test, s, 200L),
  dbarts_dec200  = function(d, s) { keep <- thin_predictors(d$x)
                                    run_dbarts_long(d$x[, keep, drop = FALSE], d$y, d$x_test[, keep, drop = FALSE], s, 200L) },
  gbart_long50   = function(d, s) run_gbart_long(d$x, d$y, d$x_test, s, 50L),
  fp_gbart50     = function(d, s) run_gbart_long(d$fp$x, d$fp$y, d$fp$x_test, s, 50L, curve),
  fp_dbarts50    = function(d, s) run_dbarts_long(d$fp$x, d$fp$y, d$fp$x_test, s, 50L, curve),
  fp_dbarts200   = function(d, s) run_dbarts_long(d$fp$x, d$fp$y, d$fp$x_test, s, 200L, curve)
)


# 3. run: every pair of a variant before the next variant ------------------------------------

pair_table <- read.csv(file.path(work_dir, "pairs.csv"))
designs <- lapply(pair_ids, function(i) { say("design for pair %d", i); pair_design(pair_table[i, ]) })
names(designs) <- pair_ids
for (v in wanted) for (i in pair_ids) for (chain in seq_len(n_chains)) {
  file <- file.path(mix_dir, sprintf("pair_%02d_%s_c%d.rds", i, v, chain))
  if (file.exists(file)) next
  if (Sys.time() > stop_at) { say("out of time at %s, pair %d", v, i); quit(save = "no", status = 0) }
  result <- tryCatch(variants[[v]](designs[[as.character(i)]], seed0 + 100000L * i + 1000L * chain + 50L),
                     error = function(e) e)
  if (inherits(result, "error")) { say("%s pair %d chain %d: %s", v, i, chain, conditionMessage(result)); next }
  saveRDS(result, file, compress = FALSE)
  say("%s pair %d chain %d: %.0f s, sigma %.4f", v, i, chain, result$seconds, mean(tail(result$sigma, 700L)))
}
say("finished")
