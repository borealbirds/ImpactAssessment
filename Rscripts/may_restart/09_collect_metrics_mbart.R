collect_metrics_mbart <- function(fit, X_train, y, 
                                  covariate, subbasin, year, top_var, present) {
  
  # 1. predict on the training data
  # continuous BART always computes the latent Gaussian regression function 
  # for both train and test sets during fitting, not just during prediction.
  # need to use predict for the multinomial case
  pred <- predict(fit, newdata = X_train)   # no type="prob"
  K_model    <- pred$K
  cats_model <- fit$cats
  
  # pixel-major -> reshape into matrix (np × K_model)
  np <- length(pred$prob.test.mean) / K_model
  prob_mean_model <- matrix(
    pred$prob.test.mean,
    ncol = K_model,
    byrow = TRUE
  )
  
  # remap MBART class order → ecological class order
  K_present <- length(present)
  cols <- match(present, cats_model)
  
  prob_ecol <- matrix(0, nrow = np, ncol = K_present)
  for (i in seq_along(present)) {
    ci <- cols[i]
    if (!is.na(ci)) {
      prob_ecol[, i] <- prob_mean_model[, ci]
    } else {
      prob_ecol[, i] <- 0
    }
  }
  
  # predicted classes in ecological order
  pred_idx <- max.col(prob_ecol, ties.method = "first")
  pred_class <- present[pred_idx]
  
  # true ecological classes
  true_class <- present[y]
  
  # accuracy
  accuracy <- mean(pred_class == true_class, na.rm = TRUE)
  
  # entropy
  eps <- 1e-12
  entropy <- -rowSums(prob_ecol * log(pmax(prob_ecol, eps)))
  
  # return metrics only (not predictions)
  data.frame(
    covariate    = covariate,
    subbasin     = subbasin,
    year         = year,
    n_holdout    = length(y),
    accuracy     = accuracy,
    mean_entropy = mean(entropy),
    top_var      = top_var,
    split        = "train"
  )
}
