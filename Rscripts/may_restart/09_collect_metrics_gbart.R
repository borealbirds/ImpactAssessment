# helper: collect metrics from a BART fit
collect_metrics_gbart <- function(fit, y, covariate, subbasin, year, top_var = NA_character_) {
  yhat_draws <- fit$yhat.train                # assumes ndpost x n
  mu    <- colMeans(yhat_draws)
  var_f <- apply(yhat_draws, 2L, var)
  sig2  <- mean(fit$sigma^2)
  se_pred <- sqrt(var_f + sig2)
  
  sigma_draws <- fit$sigma  # length ndpost
  
  # add gaussian noise to simulate the posterior predictive distribution
  yrep_draws <- fit$yhat.train + matrix(
    rnorm(length(fit$yhat.train), sd = rep(sigma_draws, each = ncol(fit$yhat.train))),
    nrow = nrow(fit$yhat.train))
  
  # compute 95% posterior predictive intervals
  pred_lower <- apply(yrep_draws, 2, quantile, 0.025)
  pred_upper <- apply(yrep_draws, 2, quantile, 0.975)
  
  # compute empirical coverage
  # the overlap between the posterior of y and the observed y
  # ~95 is desirable because it's a balance between accurate but not over-fit
  coverage95 <- mean(y >= pred_lower & y <= pred_upper)
  
  fhat_mean  <- colMeans(yhat_draws)
  resid      <- y - fhat_mean
  
  rmse <- sqrt(mean(resid^2))
  mae  <- mean(abs(resid))
  
  # Bayesian R^2 (Gelman): compute per draw, summarise
  ss_total <- sum( (y - mean(y))^2 )
  R2_post  <- apply(yhat_draws, 1L, function(mu) 1 - sum((y - mu)^2)/ss_total)
  r2_med   <- median(R2_post)
  r2_lohi  <- quantile(R2_post, c(.025, .975))
  
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
    sigma_mean  = mean(fit$sigma),
    sigma_sd    = stats::sd(fit$sigma),
    coverage95  = coverage95,
    top_var     = top_var
  )
}
