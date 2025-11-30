collect_metrics_mbart <- function(fit, X_train, y, 
                                  covariate, subbasin, year, top_var) {
  
  # 1. Predict on the training data
  # continuous BART always computes the latent Gaussian regression function 
  # for both train and test sets during fitting, not just during prediction.
  # need to use predict for the multinomial case
  pred <- predict(fit, newdata = X_train)   # no type="prob"
  
  # get most influential variable
  # collect metrics for the categorical branch
  extract_varprob_mbart <- function(fit) {
    # fit$varprob is list(ndpost) of numeric(p)
    vp <- fit$varprob
    ndpost <- length(vp)
    p      <- length(vp[[1]])
    
    # convert to ndpost × p matrix
    mat <- do.call(rbind, vp)  # ndpost × p
    
    # mean over posterior draws
    varprob_mean <- colMeans(mat)
    
    varprob_mean
  }
  
  varprob_mean <- extract_varprob_mbart(fit)
  top_var <- names(sort(varprob_mean, decreasing = TRUE))[1]
  
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
  
  # Ensure y is integer-coded 1..K
  y_mapped <- as.integer(y)
  accuracy <- mean(pred_idx == y_mapped, na.rm=TRUE)
  
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
    mean_entropy = mean(entropy),
    top_var      = top_var
  )
}
