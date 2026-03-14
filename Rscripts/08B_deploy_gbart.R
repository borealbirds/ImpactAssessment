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
    ia_dir
) {

  # source BART metrics summary functions
  source(file.path(ia_dir, "Rscripts", "09_collect_metrics_gbart.R"))
  source(file.path(ia_dir, "Rscripts", "09_collect_holdout_metrics_gbart.R"))
  
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
    
    b_mean_deg <- mean(y[-holdout_idx])
    df_backfill[[b]] <- b_mean_deg
    out_layers[[paste0(b, "_mean")]] <- rep(b_mean_deg, nrow(df_backfill))
    draw_val <- log1p(b_mean_deg)
    for (d in seq_len(100L)) {
      out_layers[[paste0(b, "_draw_", sprintf("%03d", d))]] <- rep(draw_val, nrow(df_backfill))
    }

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
    
    y_tr   <- y[-holdout_idx]
    b_mean <- mean(y_tr)
    b_logmean <- mean(log1p(y_tr))

    df_backfill[[b]] <- b_mean
    out_layers[[paste0(b, "_mean")]] <- b_mean
    for (d in seq_len(100L)) {
      out_layers[[paste0(b, "_draw_", sprintf("%03d", d))]] <- b_logmean
    }
    
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
  # yhat.test: posterior draws [ndpost x n_px]; rows = draws, columns = pixels (log1p scale).
  # Sample 100 random draws and store them as {cov}_draw_001…100 layers.
  # 12B resamples from these draws directly, avoiding the Gaussian normality assumption.
  # Original-scale mean is kept as {cov}_mean for the biotic hierarchy cascade (08A).
  b_mean   <- colMeans(expm1(fit$yhat.test))  # original scale; used by cascade
  n_draws  <- 100L
  draw_idx <- sample(nrow(fit$yhat.test), n_draws)
  draw_mat <- fit$yhat.test[draw_idx, , drop = FALSE]  # [100 x n_px] in log1p space

  # replace backfilled b in backfilling dataset for using in next iteration (following `neworder`), and for outputs
  df_backfill[[b]] <- b_mean
  out_layers[[paste0(b, "_mean")]] <- b_mean
  for (d in seq_len(n_draws)) {
    out_layers[[paste0(b, "_draw_", sprintf("%03d", d))]] <- draw_mat[d, ]
  }
  
  # --------------------------------------
  # post-modelling PART 2
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
