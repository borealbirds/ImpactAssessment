collect_holdout_metrics_mbart <- function(fit, X_holdout, y_holdout,
                                          covariate, subbasin, year, top_var) {
  
  # 1. Predict on holdout rows
  pred <- predict(fit, newdata = X_holdout)
  
  K      <- pred$K
  prob   <- pred$prob.test      # ndpost × (K*n)
  ndpost <- nrow(prob)
  np     <- ncol(prob) / K
  
  # 2. Reshape to ndpost × n × K
  prob_arr <- array(NA_real_, dim = c(ndpost, np, K))
  for (j in seq_len(K)) {
    cols <- seq(from = j, to = K * np, by = K)
    prob_arr[, , j] <- prob[, cols, drop = FALSE]
  }
  
  # 3. Posterior mean probabilities (n × K)
  prob_mean <- apply(prob_arr, c(2, 3), mean)
  
  # 4. Predicted class (MAP)
  pred_idx <- max.col(prob_mean, ties.method = "first")
  
  # True labels (must be integer-coded 1..K)
  y_mapped <- as.integer(y_holdout)
  
  accuracy <- mean(pred_idx == y_mapped, na.rm = TRUE)
  
  # 5. Entropy
  eps <- 1e-12
  entropy <- -rowSums(prob_mean * log(pmax(prob_mean, eps)))
  
  # 6. Return metrics
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
