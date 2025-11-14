collect_metrics_mbart <- function(fit, X_train, y, 
                                  covariate, subbasin, year) {
  
  # 1. Predict on the training data
  # continuous BART always computes the latent Gaussian regression function 
  # for both train and test sets during fitting, not just during prediction.
  # need to use predict for the multinomial case
  pred <- predict(fit, newdata = X_train)   # no type="prob"
  
  # 2. Extract K, n, ndpost
  K      <- pred$K
  prob   <- pred$prob.test        # ndpost × (K*n)
  ndpost <- nrow(prob)
  np     <- ncol(prob) / K        # number of training rows
  
  # 3. Reshape to 3D: ndpost × n × K
  #    BART flattens column-major by pixel then class
  prob_arr <- array(NA_real_, dim = c(ndpost, np, K))
  
  for (j in seq_len(K)) {
    cols <- seq(from = j, to = K * np, by = K)
    prob_arr[, , j] <- prob[, cols, drop = FALSE]
  }
  
  # 4. Posterior mean probs (n × K)
  prob_mean <- apply(prob_arr, c(2, 3), mean)
  
  # 5. Predicted class (MAP in mapped 1..K space)
  pred_idx <- max.col(prob_mean, ties.method = "first")
  
  # y is in original coding (e.g., 12,23,44)
  # pred_idx is in mapped coding (1..K)
  
  # SO: the calling code must map y to 1..K before calling this function
  #     (this matches usage in the main script)
  
  accuracy <- mean(pred_idx == y)
  
  # 6. Entropy
  eps <- 1e-12
  entropy <- -rowSums(prob_mean * log(pmax(prob_mean, eps)))
  
  # 7. Return metrics row
  data.frame(
    covariate    = covariate,
    subbasin     = subbasin,
    year         = year,
    ntrain       = length(y),
    accuracy     = accuracy,
    mean_entropy = mean(entropy)
  )
}