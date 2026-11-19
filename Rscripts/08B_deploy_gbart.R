# ---
# title: Impact Assessment: Gaussian BART for one continuous covariate (sourced by 07)
# author: Mannfred Boehm
# ---

# fits log1p(y) (so predictions stay positive) on 90% of the training pixels and checks it on
# the other 10%. returns the posterior mean backfill (raw scale, for 08A's cascade), 100
# posterior draws as layers <b>_draw_001..100 (log1p scale; what 12F uses) and the metrics
deploy_gbart <- function(b, y, design, seed, subbasin_index, year, logp) {

  y <- as.numeric(y)
  draw_names <- paste0(b, "_draw_", sprintf("%03d", seq_len(100L)))
  n_bf <- nrow(design$backfill)

  set.seed(seed)
  holdout_idx <- sample(seq_len(length(y)), size = round(0.10 * length(y)))
  df_holdout  <- design$train[holdout_idx, ]
  df_train    <- design$train[-holdout_idx, ]
  y_tr        <- y[-holdout_idx]

  # two values or fewer (subbasin 8, SCANFIBalsamFir_5x5): the training mean everywhere
  if (length(unique(log1p(y_tr))) <= 2) {
    logp("[%s] degenerate (≤2 unique values after log1p)", b)
    return(list(backfill = mean(y_tr),
                layers   = setNames(rep(list(rep(log1p(mean(y_tr)), n_bf)), 100L), draw_names)))
  }

  # one backfilled pixel: the training mean
  if (n_bf == 1) {
    logp("[%s] single backfill pixel: skipped gbart()", b)
    return(list(backfill = mean(y_tr),
                layers   = setNames(rep(list(mean(log1p(y_tr))), 100L), draw_names)))
  }

  fit <- BART::gbart(x.train = as.matrix(df_train),
                     y.train = log1p(y_tr),
                     x.test  = as.matrix(design$backfill),
                     type    = "wbart",
                     k       = 3,
                     ntree   = 50L,
                     ndpost  = 700L,
                     nskip   = 300L,
                     sparse  = TRUE,                  # the sampler favours informative predictors
                     sigest  = sd(log1p(y_tr)),
                     sigdf   = 3)
  logp("gbart fit to [%s]", b)

  # yhat.test is [draws x pixels], log1p scale: 100 of its rows are kept as the draws
  b_mean <- colMeans(expm1(fit$yhat.test))
  draws  <- fit$yhat.test[sample(nrow(fit$yhat.test), 100L), , drop = FALSE]

  top_var <- names(sort(colMeans(fit$varprob), decreasing = TRUE))[1]
  metrics_train <- collect_metrics_gbart(fit, y_tr, b, subbasin_index, year, top_var)

  # holdout, on the predictors the fit kept
  pred_holdout <- predict(fit, newdata = as.matrix(df_holdout[, colnames(fit$varprob), drop = FALSE]))
  metrics_holdout <- collect_holdout_metrics_gbart(y[holdout_idx], expm1(colMeans(pred_holdout)),
                                                   b, subbasin_index, year, top_var)

  list(backfill = b_mean,
       layers   = setNames(lapply(seq_len(100L), function(d) draws[d, ]), draw_names),
       metrics  = list(metrics_train, metrics_holdout))
}
