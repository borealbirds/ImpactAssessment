# helper: collect metrics from a BART fit
# fit$yhat.train: matrix (ndpost x n_train) on log1p scale
# fit$sigma: length ndpost, residual sd on log1p scale
# y: observed values on original scale (not transformed)
collect_metrics_gbart <- function(fit, y, covariate, subbasin, year, top_var = NA_character_) {
  
  
  # Transform posterior draws to original scale
  yhat_draws_orig <- expm1(fit$yhat.train)   # ndpost x ntrain, ORIGINAL scale
  mu <- colMeans(yhat_draws_orig)            # posterior mean per observation on original scale
  var_f <- apply(yhat_draws_orig, 2L, var)
  
  # If you want a predictive se: we need the predictive variance on original scale.
  # We won't mix log-scale sigma directly here; instead compute empirical predictive variance:
  sig2_empirical <- mean(apply(yhat_draws_orig, 2, var))
  
  # For posterior predictive draws: add noise in log-space, then transform
  sigma_draws_log <- fit$sigma               # ndpost (log1p scale)
  ndpost <- nrow(fit$yhat.train)
  ntrain <- ncol(fit$yhat.train)
  
  # build yrep in log-space then transform:
  # yrep_log_draws = fit$yhat.train + noise(log-scale)
  noise_mat <- matrix(rnorm(ndpost * ntrain, sd = rep(sigma_draws_log, each = ntrain)), 
                      nrow = ndpost, ncol = ntrain)
  yrep_log_draws <- fit$yhat.train + noise_mat
  yrep_draws_orig <- expm1(yrep_log_draws)   # transform to original scale
  
  # compute 95% posterior predictive intervals on original scale
  pred_lower <- apply(yrep_draws_orig, 2, quantile, 0.025)
  pred_upper <- apply(yrep_draws_orig, 2, quantile, 0.975)
  
  # empirical coverage (original-scale)
  coverage95 <- mean(y >= pred_lower & y <= pred_upper)
  
  # residuals against posterior mean (original-scale)
  resid <- y - mu
  rmse <- sqrt(mean(resid^2))
  mae  <- mean(abs(resid))
  
  # Bayesian R^2 (Gelman) on ORIGINAL scale:
  ss_total <- sum((y - mean(y))^2)
  # apply over posterior draws (using yhat_draws_orig)
  R2_post <- apply(yhat_draws_orig, 1L, function(yhat_draw) {
    1 - sum((y - yhat_draw)^2) / ss_total
  })
  r2_med  <- median(R2_post)
  r2_lohi <- quantile(R2_post, c(.025, .975))
  
  data.frame(
    covariate   = covariate,
    subbasin    = subbasin,
    year        = year,
    ntrain      = length(y),
    rmse        = rmse,
    mae         = mae,
    r2_median   = r2_med,
    r2_lo       = unname(r2_lohi[1]),
    r2_hi       = unname(r2_lohi[2]),
    sigma_mean  = mean(fit$sigma),   # note: this is on log1p scale
    sigma_sd    = sd(fit$sigma),     # on log1p scale
    coverage95  = coverage95,
    top_var     = top_var,
    split       = "train"
  )
}
