collect_holdout_metrics_gbart <- function(y_obs, yhat_mean,
                                          covariate, subbasin, year, top_var) {
  
  resid <- y_obs - yhat_mean
  
  rmse <- sqrt(mean(resid^2))
  mae  <- mean(abs(resid))
  
  r2    <- 1 - sum(resid^2) / sum( (y_obs - mean(y_obs))^2 )
  
  data.frame(
    covariate   = covariate,
    subbasin    = subbasin,
    year        = year,
    n_holdout   = length(y_obs),
    rmse        = rmse,
    mae         = mae,
    r2          = r2,
    top_var     = top_var,
    split       = "holdout"
  )
}
