collect_metrics_mbart <- function(fit, X_train, y, 
                                  covariate, subbasin, year, top_var, present) {
  
  # 1. predict on the training data
  # continuous BART always computes the latent Gaussian regression function 
  # for both train and test sets during fitting, not just during prediction.
  # need to use predict for the multinomial case
  X_train <- X_train[, colnames(fit$varprob[[1]]), drop = FALSE] # drop covariates not used by BART
  pred <- predict(fit, newdata = X_train)   # no type="prob"
  K_model    <- pred$K
  
  # get matrix with most likely classes per holdout pixel (rows)
  # `pred` was generated using `fit`, 
  # so it's limited to the classes in `fit$cats`
  prob_ecol <- matrix(pred$prob.test.mean, ncol = K_model, byrow = TRUE)
  colnames(prob_ecol) <- present[fit$cats]
  
  # predictions
  pred_class <- present[max.col(prob_ecol, ties.method = "first")]

  # accuracy
  accuracy <- mean(pred_class == y, na.rm = TRUE)
  
  # entropy
  eps <- 1e-12
  entropy <- -rowSums(prob_ecol * log(pmax(prob_ecol, eps)))
  
  # return metrics only (not predictions)
  data.frame(
    covariate    = covariate,
    subbasin     = subbasin,
    year         = year,
    n_obs    = length(y),
    accuracy     = accuracy,
    mean_entropy = mean(entropy),
    top_var      = top_var,
    split        = "train"
  )
}
