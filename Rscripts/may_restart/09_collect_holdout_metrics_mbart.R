collect_holdout_metrics_mbart <- function(fit, pred, actual_holdout_class,
                                          covariate, subbasin, year, top_var, present) {
  
  # get matrix with most likely classes per holdout pixel (rows)
  # `pred` was generated using `fit`, 
  # so it's limited to the classes in `fit$cats`
  K_model   <- pred$K
  prob_ecol <- matrix(pred$prob.test.mean, ncol = K_model, byrow = TRUE)
  colnames(prob_ecol) <- present[fit$cats]
  
  # predictions
  yhat_holdout_class <- present[max.col(prob_ecol, ties.method = "first")]
  
  # accuracy
  accuracy <- mean(yhat_holdout_class == actual_holdout_class, na.rm = TRUE)
  
  # entropy
  eps <- 1e-12
  entropy <- -rowSums(prob_ecol * log(pmax(prob_ecol, eps)))
  
  # return metrics only (not predictions)
  data.frame(
    covariate    = covariate,
    subbasin     = subbasin,
    year         = year,
    n_holdout    = length(y_holdout),
    accuracy     = accuracy,
    mean_entropy = mean(entropy),
    top_var      = top_var,
    split        = "holdout"
  )
}
