deploy_gbart <- function(
    b,
    df_train,
    df_train_bart,
    df_backfill,
    df_backfill_bart,
    idx,
    biotic_cols,
    subbasin_index,
    year,
    metrics,
    out_layers,
    logp,
    cc,
    root
) {
  
  # source BART metrics summary functions
  if(!cc){
    source(file.path(getwd(), "Rscripts", "may_restart", "09_collect_metrics_gbart.R"))
    source(file.path(getwd(), "Rscripts", "may_restart", "09_collect_holdout_metrics_gbart.R"))
  }

  if(cc){
    source(file.path(root, "Rscripts",  "09_collect_metrics_gbart.R"))
    source(file.path(root, "Rscripts",  "09_collect_holdout_metrics_gbart.R"))
  }
  
  # define continuous response, and subset pixels to non-NAs for a given `b`
  # we use `df_train` instead of `df_train_bart` because the latter is purposely missing `b`
  y <- as.numeric(df_train[[b]][idx]) # idx ensures length(y) == length(df_train_bart)

  # reproducible per (year, subbasin, covariate)
  set.seed(abs(as.integer(sprintf("%d%03d", subbasin_index, which(biotic_cols==b)))))
      
  # split y into 90% training and 10% holdout
  holdout_idx <- sample(seq_len(length(y)), size = round(0.10 * length(y)))
  df_holdout    <- df_train_bart[holdout_idx, ] # get holdout rows before modifying df_train_bart
  df_train_bart <- df_train_bart[-holdout_idx, ]
  
  # check that y is not binary
  # this happened in subbasin 8 with SCANFIBalsamFir_5x5
  if (length(unique(log1p(y[-holdout_idx]))) <= 2) {
    
    df_backfill[[b]] <- mean(y[-holdout_idx])
    out_layers[[paste0(b, "_mean")]] <- rep(mean(y[-holdout_idx]), nrow(df_backfill))
    out_layers[[paste0(b, "_sd")]]   <- rep(0, nrow(df_backfill))
    
    logp("[%s] degenerate (≤2 unique values after log1p)", b)
    return(list(
    df_backfill = df_backfill,
    metrics     = metrics,
    out_layers  = out_layers
  ))
  }
  
  # check that there is >1 high HF pixel to backfill (subbasins 423,304,247,246,221,160)
  # if nrow(df_backfill_bart) == 1 then take the mean from the training pixels
  if (nrow(df_backfill_bart) == 1) {
    
    mu <- mean(log1p(y[-holdout_idx]))
    sd_y <- sd(log1p(y[-holdout_idx]))
    
    b_mean <- expm1(mu)
    b_sd   <- expm1(sd_y)
    
    df_backfill[[b]] <- b_mean
    out_layers[[paste0(b, "_mean")]] <- b_mean
    out_layers[[paste0(b, "_sd")]]   <- b_sd
    
    logp("[%s] single backfill pixel: skipped gbart()", b)
    return(list(
      df_backfill = df_backfill,
      metrics     = metrics,
      out_layers  = out_layers
    ))
  }
  
  # train and predict with BART
  fit <- BART::gbart(x.train = as.matrix(df_train_bart), 
                     y.train = log1p(y[-holdout_idx]), # log transform to avoid negative predictions
                     x.test = as.matrix(df_backfill_bart), # has the same columns as df_train_bart
                     type = "wbart",
                     k = 3, #shrinkage
                     ntree = 50L, 
                     ndpost = 700L, 
                     nskip = 300L, 
                     sparse = TRUE, # sampler focuses on informative predictors (not all predictors treated as informative)
                     sigest  = sd(log1p(y[-holdout_idx])), # the rough error standard deviation used in the prior
                     sigdf = 3) 
  
  logp("gbart fit to [%s]", b)
  
  # --------------------------------------
  # post-modelling PART 1
  # estimate posterior mean and sd
  # yhat.test are the predicted values of b at high HF locations (rows = posterior draws, columns = pixels)
  mu_test <- fit$yhat.test.mean # posterior mean (log1p scale)
  sigma2_test <- apply(fit$yhat.test, 2, var)  # posterior variance (log1p scale)
  
  b_mean <- expm1(mu_test) # E[Y] approximation for b at high HF locations
  b_sd   <- sqrt(expm1(sigma2_test) * exp(2*mu_test + sigma2_test))  # log-normal SD
  
  # replace backfilled b in backfilling dataset for using in next iteration (following `neworder`), and for outputs
  df_backfill[[b]] <- b_mean 
  out_layers[[paste0(b, "_mean")]] <- b_mean # add values of b to a list
  out_layers[[paste0(b, "_sd")]]   <- b_sd
  
  # --------------------------------------
  # post-modelling PART 2
  # posterior predictive check: get posterior predictions for training data
  # yhat.train are the predicted values of b at low HF locations (values are from the fitted model, not the training dataset)
  mu_train <- fit$yhat.train.mean
  sigma2_train <- apply(fit$yhat.train, 2, var)
  
  ppc_mean <- expm1(mu_train)
  ppc_sd   <- sqrt(mean(expm1(sigma2_train) * exp(2*mu_train + sigma2_train)))
  
  # p is the probability the test statistic in a replicated data set exceeds that in the original data
  # calculate Bayesian p-value (proportion of times observed statistic > predicted)
  # Bayesian p-values for mean and SD
  p_value_mean <- mean(apply(expm1(fit$yhat.train), 1, function(draw) mean(draw) > mean(y[-holdout_idx])))
  p_value_sd   <- mean(apply(expm1(fit$yhat.train), 1, function(draw) sd(draw) > sd(y[-holdout_idx])))
  
  
  # --------------------------------------
  # post-modelling PART 3
  # extract variable with highest importance 
  top_var <- NA_character_
  if ("varprob" %in% names(fit)) {
    top_var <- names(sort(colMeans(fit$varprob), decreasing=TRUE))[1]
  } else if ("varcount" %in% names(fit)) {
    vc <- colMeans(fit$varcount)
    top_var <- names(sort(vc, decreasing=TRUE))[1]
  }
  
  # collect metrics (one row per year x subbasin x covariate)
  metrics[[length(metrics) + 1L]] <- collect_metrics_gbart(
    fit = fit, 
    y = y[-holdout_idx],
    covariate = b,
    subbasin  = subbasin_index,
    year      = year,
    top_var   = top_var)
  
  
  # out-of-sample prediction
  # first, filter for covariates that actually contributed to `fit`
  df_holdout <- df_holdout[, colnames(fit$varprob), drop = FALSE]
  pred_holdout <- predict(fit, newdata = as.matrix(df_holdout)) # log1p scale
  yhat_holdout_mean <- expm1(colMeans(pred_holdout)) # get mean estimate per pixel from the posterior
  
  metrics[[length(metrics) + 1L]] <- collect_holdout_metrics_gbart(
    y_obs    = y[holdout_idx],
    yhat_mean = yhat_holdout_mean,
    covariate = b,
    subbasin  = subbasin_index,
    year      = year,
    top_var   = top_var)
  
  # send outputs to "train_and_backfill_subbasin_s"
  list(
    df_backfill = df_backfill,
    metrics     = metrics,
    out_layers  = out_layers
  )

} # close function 
