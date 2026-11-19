# ---
# title: Impact Assessment: BART training and holdout metrics, one row per model (sourced by 07)
# author: Mannfred Boehm
# ---

# gbart, training pixels. fit$yhat.train is [draws x pixels] on the log1p scale; y is raw
collect_metrics_gbart <- function(fit, y, covariate, subbasin, year, top_var) {

  yhat <- expm1(fit$yhat.train)
  mu   <- colMeans(yhat)

  # 95% posterior predictive interval: residual noise added on the log1p scale
  ndpost <- nrow(fit$yhat.train)
  ntrain <- ncol(fit$yhat.train)
  noise  <- matrix(rnorm(ndpost * ntrain, sd = rep(fit$sigma, each = ntrain)), nrow = ndpost, ncol = ntrain)
  yrep   <- expm1(fit$yhat.train + noise)
  coverage95 <- mean(y >= apply(yrep, 2, quantile, 0.025) & y <= apply(yrep, 2, quantile, 0.975))

  # Bayesian R^2 (Gelman), one per posterior draw
  ss_total <- sum((y - mean(y))^2)
  r2_post  <- apply(yhat, 1L, function(yhat_draw) 1 - sum((y - yhat_draw)^2) / ss_total)
  r2_lohi <- quantile(r2_post, c(.025, .975))

  data.frame(covariate  = covariate,
             subbasin   = subbasin,
             year       = year,
             ntrain     = length(y),
             rmse       = sqrt(mean((y - mu)^2)),
             mae        = mean(abs(y - mu)),
             r2_median  = median(r2_post),
             r2_lo      = unname(r2_lohi[1]),
             r2_hi      = unname(r2_lohi[2]),
             sigma_mean = mean(fit$sigma),   # log1p scale
             sigma_sd   = sd(fit$sigma),
             coverage95 = coverage95,
             top_var    = top_var,
             split      = "train")
}


# gbart, holdout pixels: observed vs posterior mean, both raw scale
collect_holdout_metrics_gbart <- function(y_obs, yhat_mean, covariate, subbasin, year, top_var) {

  resid <- y_obs - yhat_mean

  data.frame(covariate = covariate,
             subbasin  = subbasin,
             year      = year,
             n_holdout = length(y_obs),
             rmse      = sqrt(mean(resid^2)),
             mae       = mean(abs(resid)),
             r2        = 1 - sum(resid^2) / sum((y_obs - mean(y_obs))^2),
             top_var   = top_var,
             split     = "holdout")
}


# mbart, accuracy and mean entropy of the modal class. `pred` is predict(fit) on the pixels
# whose true classes are `actual`; present are the fit's classes, 1..K
mbart_metrics <- function(pred, actual, present) {
  prob <- matrix(pred$prob.test.mean, ncol = pred$K, byrow = TRUE)
  list(accuracy     = mean(present[max.col(prob, ties.method = "first")] == actual, na.rm = TRUE),
       mean_entropy = mean(-rowSums(prob * log(pmax(prob, 1e-12)))))
}

# mbart, training pixels (mbart2 does not keep its training predictions, so re-predict)
collect_metrics_mbart <- function(fit, X_train, y, covariate, subbasin, year, top_var, present) {

  pred <- predict(fit, newdata = X_train[, colnames(fit$varprob[[1]]), drop = FALSE])
  m    <- mbart_metrics(pred, y, present)

  data.frame(covariate    = covariate,
             subbasin     = subbasin,
             year         = year,
             n_obs        = length(y),
             accuracy     = m$accuracy,
             mean_entropy = m$mean_entropy,
             top_var      = top_var,
             split        = "train")
}

# mbart, holdout pixels
collect_holdout_metrics_mbart <- function(pred, actual_holdout_class, covariate, subbasin, year, top_var, present) {

  m <- mbart_metrics(pred, actual_holdout_class, present)

  data.frame(covariate    = covariate,
             subbasin     = subbasin,
             year         = year,
             n_holdout    = length(actual_holdout_class),
             accuracy     = m$accuracy,
             mean_entropy = m$mean_entropy,
             top_var      = top_var,
             split        = "holdout")
}
