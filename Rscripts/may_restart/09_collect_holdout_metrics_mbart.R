collect_holdout_metrics_mbart <- function(fit, pred, y_holdout,
                                          covariate, subbasin, year, top_var, present) {
  
  K_model    <- pred$K
  
  # pixel-major -> reshape into matrix (np × K_model)
  np <- length(pred$prob.test.mean) / K_model
  prob_mean_model <- matrix(
    pred$prob.test.mean,
    ncol = K_model,
    byrow = TRUE
  )
  
  # mbart columns already correspond to present[1:K′]
  prob_ecol <- prob_mean_model
  colnames(prob_ecol) <- present
  
  # predicted classes in ecological order
  pred_idx <- max.col(prob_ecol, ties.method = "first")
  pred_class <- present[pred_idx]
  
  # accuracy
  accuracy <- mean(pred_class == y_holdout, na.rm = TRUE)
  
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
