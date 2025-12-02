train_and_backfill_subbasin_s <- function(
    subbasin_index, 
    year, 
    stack_y,     # SpatRaster 
    lowhf_mask,  # SpatRaster 
    highhf_mask, # SpatRaster 
    all_subbasins_subset, # SpatVector
    abiotic_vars, # tibble with column `predictor`
    biotic_vars,  # tibble with column `predictor`
    categorical_responses, # character vector of layer names
    ia_dir,
    neworder, 
    quiet = FALSE
) {
  
  # source BART metrics summary functions
  source(file.path(getwd(), "Rscripts", "may_restart", "09_collect_metrics_gbart.R"))
  source(file.path(getwd(), "Rscripts", "may_restart", "09_collect_metrics_mbart.R"))
  source(file.path(getwd(), "Rscripts", "may_restart", "09_collect_holdout_metrics_gbart.R"))
  source(file.path(getwd(), "Rscripts", "may_restart", "09_collect_holdout_metrics_mbart.R"))
  
  # for logging progress
  logfile <- file.path(ia_dir, "logs", sprintf("Y%d_S%d.log", year, subbasin_index))
  logp <- make_logger(logfile)
  on.exit(logp("done"), add = TRUE)
  
  logp("start")
  
  # isolate subbasin s
  subbasin_s <- all_subbasins_subset[subbasin_index]
  
  # crop covariate stack to subbasin
  cov_s <- 
    stack_y |>
    terra::crop(x = _, y = subbasin_s) |>
    terra::mask(x = _, mask = subbasin_s)
  
  logp("cropped/masked; ncell=%d", terra::ncell(cov_s))
  
 
  # for training: identify low HF pixels in subbasin_s 
  lowhf_mask_s  <- terra::resample(lowhf_mask,  cov_s, method="near")
  lowhf_mask_s  <- terra::mask(lowhf_mask_s, subbasin_s)
  
  # for backfilling: identify high HF pixels in subbasin_s 
  highhf_mask_s <- terra::resample(highhf_mask, cov_s, method="near")
  highhf_mask_s <- terra::mask(highhf_mask_s, subbasin_s)
  
  # for training and backfilling: select pixels marked with 1s 
  # convert covariate stacks to a dataframe for training/backfilling (rows are pixels, columns are covariates)
  # df_train will have NAs filtered per covariate using `idx` in `for (b in biotic_cols) {`
  cov_train_s <- terra::mask(x = cov_s, mask = lowhf_mask_s)
  df_train    <- terra::as.data.frame(cov_train_s, xy = TRUE, na.rm = FALSE)
  
  # put marked pixels (1s) from highhf_mask_s into a dataframe
  # keep track of non-NA pixels in `backfill_idx`
  # these will be used for re-populating an empty raster with backfill values 
  df_full <- terra::as.data.frame(cov_s, xy = TRUE, na.rm = FALSE, cells = TRUE)
  backfill_idx <- which(values(highhf_mask_s) == 1)
  df_backfill <- df_full[backfill_idx, , drop = FALSE]
  
  # coordinate scaling: define mean and variance of lat/long across subbasin_s
  # compute center/scale from rows with valid x,y (should be many unless something odd)
  xy_ok <- which(!is.na(df_full[["x"]]) & !is.na(df_full[["y"]]))
  xy_center <- colMeans(df_full[xy_ok, c("x", "y"), drop = FALSE], na.rm = TRUE)
  xy_scale  <- apply(df_full[xy_ok, c("x", "y"), drop = FALSE], 2, sd, na.rm = TRUE)
  
  df_train[, c("x","y")]    <- sweep(df_train[c("x","y")], 2, xy_center, "-") #subtract center from every coordinate
  df_train[, c("x","y")]    <- sweep(df_train[c("x","y")], 2, xy_scale, "/") # divide every coordinate by variance
  df_backfill[, c("x","y")] <- sweep(df_backfill[, c("x","y")], 2, xy_center, "-")
  df_backfill[, c("x","y")] <- sweep(df_backfill[, c("x","y")], 2, xy_scale,  "/")
  
  logp("backfill rows=%d", nrow(df_backfill))
  
  # define predictors and responses
  abiotic_cols <- intersect(names(df_train), abiotic_vars$predictor) # abiotic features present in subbasin_s
  biotic_cols  <- intersect(names(df_train), biotic_vars$predictor)  # biotic features present in subbasin_s
  biotic_cols  <- na.omit(biotic_cols[match(neworder, biotic_cols)]) # enforce hierarchy
  biotic_cols_cont <- setdiff(biotic_cols, categorical_responses) # identify continuous biotic features in subbasin_s
  
  # collect backfilled layers in `out_layers` list 
  # each element will be a numeric of the mean posterior values of a backfilled covariate 
  out_layers <- list()

  # store BART metrics here, each element will be a single row dataframe with model metrics as columns
  metrics <- list()
  
  # store confusion matrix information for this subbasin
  confusion <- list()
  
  # train a model for each vegetation feature
  # backfill each feature where human footprint is high
  # for (b in biotic_cols) {
  for (b in c("SCANFI_1km")) {
    
    t0 <- proc.time()[3]
    logp("[%s] in process: (%s)", b, if (b %in% categorical_responses) "categorical" else "continuous")
    
    # training data: keep only the rows where biotic feature `b` is observed
    idx <- which(!is.na(df_train[[b]]))
    
    # check that sample size of b > 1, and for b invariance
    if (length(idx) == 0 || length(unique(df_train[[b]][idx])) < 2) {
      
      const_val <- unique(df_train[[b]][idx])
      # if idx has no valid rows (empty), define const_val as NA (as it would be in V5 rasters)
      if (length(const_val) == 0) const_val <- NA
      
      # fill outputs
      df_backfill[[b]] <- const_val
      out_layers[[paste0(b, "_mean")]] <- rep(const_val, nrow(df_backfill))
      out_layers[[paste0(b, "_sd")]]   <- rep(0,          nrow(df_backfill))  # no uncertainty
      
      logp("[%s] constant in subbasin: backfill raster populated with %s", b, const_val)
      next
    }
  
    # biotic predictors that precede `b` in `neworder`
    if (!(b %in% categorical_responses)){
      b_before <- biotic_cols_cont[seq_len(match(b, biotic_cols_cont) - 1)]
    } else {
      b_before <- biotic_cols_cont 
    }
    
    # the predictors for b are abiotic_cols, b_before, and lat/long
    # note: categorical biotic features excluded from predicting continuous biotic features
    predictors <- c(abiotic_cols, b_before, "x", "y")
    
    # create training data frame (i.e. predictors) for b
    # columns are abiotic_cols, biotic_cols (except b), lat/long
    # rows are pixel locations with low HF (idx are non-NA cells for b)
    df_train_bart <- df_train[idx, predictors, drop = FALSE]
    
    # drop predictors that are all NAs
    df_train_bart <- df_train_bart[, colSums(!is.na(df_train_bart)) > 0]
    logp("training pixels = %d", nrow(df_train_bart))
    logp("training predictors = %d", ncol(df_train_bart))
    
    # drop predictors with zero variance
    # note: in subbasin 57 we had the unlikely case of only a single non-NA value for CAfire, which gives sd()=NA
    col_sd <- sapply(df_train_bart, function(x) sd(as.numeric(x), na.rm = TRUE))
    df_train_bart <- df_train_bart[, !is.na(col_sd) & col_sd > 0 , drop = FALSE] 
    
    # subset backfill dataframe to the same abiotic_cols, biotic_cols, lat/long
    # because we can't backfill covariates that we didn't train models for
    df_backfill_bart  <- df_backfill[, colnames(df_train_bart), drop = FALSE]
    
    # check that the backfill dataframe doesn't have any all-NA predictors (different pixels than training locations, so it's possible..)
    # if it does, remove them from both training and backfill dataframes
    df_backfill_bart <- df_backfill_bart[, colSums(!is.na(df_backfill_bart )) > 0]
    df_train_bart <- df_train_bart[, colnames(df_backfill_bart), drop = FALSE]
    
    # train: check if continuous or categorical
    if (!(b %in% categorical_responses)){
      
      # define continuous response, and subset pixels to non-NAs for a given `b`
      # we use `df_train` instead of `df_train_bart` because the latter is purposely missing `b`
      y <- as.numeric(df_train[[b]][idx]) # idx ensures length(y) == length(df_train_bart)

      # reproducible per (year, subbasin, covariate)
      set.seed(abs(as.integer(sprintf("%d%03d", subbasin_index, which(biotic_cols==b)))))
      
      # split y into 90% training and 10% holdout
      holdout_idx <- sample(seq_len(length(y)), size = round(0.10 * length(y)))
      df_holdout    <- df_train_bart[holdout_idx, ] # get holdout rows before modifying df_train_bart
      df_train_bart <- df_train_bart[-holdout_idx, ]
      
      # train and predict with BART
      fit <- BART::gbart(x.train = as.matrix(df_train_bart), 
                         y.train = y[-holdout_idx], 
                         x.test = as.matrix(df_backfill_bart), # has the same columns as df_train_bart
                         type = "wbart",
                         k = 3, #shrinkage
                         ntree = 50L, 
                         ndpost = 700L, 
                         nskip = 300L, 
                         sparse = TRUE, # sampler focuses on informative predictors (not all predictors treated as informative)
                         sigest  = sd(y), # the rough error standard deviation used in the prior
                         sigdf = 3) 
      
      logp("gbart fit to [%s]", b)
      
      # estimate posterior mean and sd
      # yhat.test are the predicted values of b at high HF locations (rows = posterior draws, columns = pixels)
      b_mean <- as.numeric(fit$yhat.test.mean) # mean prediction for b at high HF locations
      b_sd   <- apply(fit$yhat.test, 2, sd) # uncertainty in model-predicted values of b per pixel
      
      # replace backfilled b in backfilling dataset for using in next iteration (following `neworder`), and for outputs
      df_backfill[[b]] <- b_mean
      out_layers[[paste0(b, "_mean")]] <- b_mean # add values of b to a list
      out_layers[[paste0(b, "_sd")]]   <- b_sd
      
      # posterior predictive check: get posterior predictions for training data
      # yhat.train are the predicted values of b at low HF locations (values are from the fitted model, not the training dataset)
      ppc_mean <- mean(colMeans(fit$yhat.train))
      ppc_sd <- mean(apply(fit$yhat.train, 2, sd))
      
      # p is the probability the test statistic in a replicated data set exceeds that in the original data
      # calculate Bayesian p-value (proportion of times observed statistic > predicted)
      p_value_mean <- mean(apply(fit$yhat.train, 1, function(x) mean(x) > mean(y)))
      p_value_sd   <- mean(apply(fit$yhat.train, 1, function(x) sd(x) > sd(y)))
      
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
        fit, y[-holdout_idx],
        covariate = b,
        subbasin  = subbasin_index,
        year      = year,
        top_var   = top_var)
      
      
      # out-of-sample prediction
      pred_holdout <- predict(fit, newdata = as.matrix(df_holdout))
      yhat_holdout_mean <- colMeans(pred_holdout) # get mean estimate per pixel from the posterior
      
      metrics[[length(metrics) + 1L]] <- collect_holdout_metrics_gbart(
        y_obs    = y[holdout_idx],
        yhat_mean = yhat_holdout_mean,
        covariate = b,
        subbasin  = subbasin_index,
        year      = year,
        top_var   = top_var
      )
      
    } else { # if b is categorical, use the following protocol
      
      # BART::mbart needs the landcover classes as a consecutive sequence of integers 1..K' 
      # but actual classes could be any arbitrary non-consecutive sequence of integers 1..K
      # so, we convert 1..K -> 1..K' -> train -> backfill as 1..K' -> convert back to 1..K
      
      # the vector of original landcover codes for the selected training pixels. `idx` ensures no NAs present
      # don't use df_train_bart because it's missing `b`. 
      y_codes <- as.character(df_train[[b]][idx]) 

      # make sure y isn’t empty 
      if (!length(y_codes)) { logp("[%s] skip: empty y_codes", b); next } 
      
      # land cover classes that are *actually* present in the training subset
      present <- sort(unique(y_codes))
      K <- length(present) 
      if (K == 0) { logp("[%s] skip: no classes present", b); next }
      
      # pre-allocate outputs
      b_mode    <- rep(NA_real_, nrow(df_backfill_bart))
      b_maxprob <- rep(NA_real_, nrow(df_backfill_bart))
      b_entropy <- rep(NA_real_, nrow(df_backfill_bart))

      # build a forward and inverse map only over present classes (size K')
      fwd_map  <- setNames(seq_along(present), present)   # original code -> 1..K'
      inv_map  <- setNames(present, seq_along(present))   # 1..K' -> original code
      
      # determine complexity of landcover classes
      if (K < 2) {
        
        # trivial: only one class in training data
        b_mode    <- rep(present[1], nrow(df_backfill_bart)) # the one and only landcover class present
        b_maxprob <- rep(1, nrow(df_backfill_bart)) 
        b_entropy <- rep(0, nrow(df_backfill_bart))

      } else {
        
        # training labels mapped to 1..K'
        y_int <- unname(fwd_map[y_codes])  
        
        # reproducible per (year, subbasin, covariate)
        set.seed(abs(as.integer(sprintf("%d%03d", subbasin_index, which(biotic_cols==b)))))
        
        # split y into 90% training and 10% holdout
        holdout_idx <- sample(seq_len(length(y_int)), size = round(0.10 * length(y_int)))
        df_holdout    <- df_train_bart[holdout_idx, ] # get holdout rows before modifying df_train_bart
        df_train_bart <- df_train_bart[-holdout_idx, ]
        
        # multinomial bart keeps every 10th posterior sample (gbart keeps every sample)
        # so mbart will take at least 10 times longer
        # gbart assumes continuous residuals -> conjugate priors -> less autocorrelation -> no thinning
        # mbart uses latent variables -> non-conjugate updates -> more autocorrelation -> needs thinning
        # mbart runs K times (e.g. 11 classes will take 11 times longer than gbart)
        fit <- BART::mbart(x.train = as.matrix(df_train_bart),
                           y.train = y_int[-holdout_idx],
                           x.test  = as.matrix(df_backfill_bart),
                           type = "pbart",
                           k = 3,
                           ntree = 40L,
                           ndpost = 400L,
                           nskip = 200L,
                           keepevery = 10L,
                           printevery = 350,
                           sparse = TRUE)
        
        fit <- BART::mbart(x.train=as.matrix(df_train_bart),
                           y.train = y_int[-holdout_idx],
                           x.test=as.matrix(df_backfill_bart),
                           ndpost = 750L)
        
        logp("mbart fit to [%s]", b)
      
        # reshape `prob.test` into draws x pixels x model-classes
        pred <- predict(fit, newdata = as.matrix(df_backfill_bart))
        
        prob_test <- pred$prob.test   # ndpost × (npixels*K_model)
        K_model   <- fit$K
        cats_model <- fit$cats       # MBART internal class order
        ndpost    <- nrow(prob_test)
        npixels   <- ncol(prob_test) / K_model
        
        # Reshape from class-fastest → (draw × pixel × class)
        prob_array <- array(NA_real_, dim = c(ndpost, npixels, K_model))
        
        for (j in seq_len(K_model)) {
          cols_j <- ((j - 1) * npixels + 1):(j * npixels)
          prob_array[, , j] <- prob_test[, cols_j, drop = FALSE]
        }
        
        # Posterior mean across draws → (pixel × class_model)
        class_probs_model <- apply(prob_array, c(2, 3), mean)  # npixels × K_model
        
        # Match ecological classes present in df_train to MBART’s internal class order
        cols <- match(present, cats_model)
        
        # Build probability matrix in ECOLOGICAL order (pixel × length(present))
        class_probs <- matrix(0, nrow = npixels, ncol = length(present))
        
        for (i in seq_along(present)) {
          ci <- cols[i]
          if (!is.na(ci)) {
            # Probability for this ecological class = model column ci
            class_probs[, i] <- class_probs_model[, ci]
          } else {
            # If class not present in training → prob = 0
            class_probs[, i] <- 0
          }
        }
        
        # Most likely class index (in ecological order)
        pred_idx <- max.col(class_probs, ties.method = "first")
        b_mode   <- present[pred_idx]
        
        # Uncertainty: entropy
        eps <- 1e-12
        b_entropy <- -rowSums(class_probs * log(pmax(class_probs, eps)))
        
        # Probability of predicted class
        b_maxprob <- apply(class_probs, 1, max)
        
        # Variable importance
        varprob_mean <- colMeans(do.call(rbind, fit$varprob))
        top_var <- names(sort(varprob_mean, decreasing = TRUE))[1]
        
      } # close if single or multi-class
      
      # outputs
      df_backfill[[b]] <- b_mode
      out_layers[[b]]  <- b_mode
      out_layers[[paste0(b, "_maxprob")]]  <- b_maxprob
      out_layers[[paste0(b, "_entropy")]]  <- b_entropy
      
      if (K >= 2){
      metrics[[length(metrics) + 1L]] <- collect_metrics_mbart(
        fit,
        y         = y_int[-holdout_idx],
        X_train   = df_train_bart,
        covariate = b,
        subbasin  = subbasin_index,
        year      = year,
        top_var   = top_var)
      } 
      
      # out-of-sample prediction
      pred_holdout <- predict(fit, newdata = as.matrix(df_holdout))
      
      # MBART internal dimensions
      K_model <- pred_holdout$K
      cats_model <- fit$cats   # MBART internal class order
      
      prob   <- pred_holdout$prob.test
      ndpost <- nrow(prob)
      np     <- ncol(prob) / K_model   # number of holdout rows
      
      # reshape ndpost × np × K_model
      prob_arr <- array(NA_real_, dim = c(ndpost, np, K_model))
      for (j in seq_len(K_model)) {
        cols <- seq(from = j, to = K_model * np, by = K_model)
        prob_arr[, , j] <- prob[, cols, drop = FALSE]
      }
      
      # posterior mean probs (np × K_model)
      prob_holdout_mean <- apply(prob_arr, c(2, 3), mean)
      
      # Re-map MBART class order → ecological order
      K_present <- length(present)
      cols <- match(present, cats_model)
      
      prob_holdout_ecol <- matrix(0, nrow = np, ncol = K_present)
      
      for (i in seq_along(present)) {
        ci <- cols[i]
        if (!is.na(ci)) {
          prob_holdout_ecol[, i] <- prob_holdout_mean[, ci]
        } else {
          prob_holdout_ecol[, i] <- 0
        }
      }
      
      # predicted class (MAP)
      yhat_holdout_class <- max.col(prob_holdout_ecol, ties.method = "first")
      
      # entropy
      eps <- 1e-12
      b_entropy_holdout <- -rowSums(prob_holdout_ecol * log(pmax(prob_holdout_ecol, eps)))
      
      # collect holdout metrics
      metrics[[length(metrics) + 1L]] <- collect_holdout_metrics_mbart(
        fit       = fit,
        X_holdout = df_holdout,
        y_holdout = y_int[holdout_idx], 
        covariate = b,
        subbasin  = subbasin_index,
        year      = year,
        top_var   = top_var
      )
      
      # store data for confusion matrix
      # actual ecological codes for the holdout rows
      actual_holdout_class <- present[ y_int[holdout_idx] ]
      
      # store holdout predictions for later confusion matrix
      confusion[[length(confusion) + 1L]] <- data.frame(
        covariate  = b,
        subbasin   = subbasin_index,
        actual     = actual_holdout_class,
        predicted  = yhat_holdout_class,
        stringsAsFactors = FALSE
      )
    } # close if continuous or categorical
    
    logp("[%s] fit done in %.1fs", b, proc.time()[3] - t0)

  } # close for loop over biotic_cols
  
  # after backfilling all covariates, build and save new stack for subbasin_s
  # keep only the layers that were actually created (remove NULL layers)
  # build result raster with predicted layers only
  created <- names(out_layers)[!vapply(out_layers, is.null, logical(1))]
  if (!length(created)) return(invisible(NULL))
  
  template <- cov_s[[1]]
  result_raster <- terra::rast(template, nlyr = length(created))
  
  # fill only the high-HF cells
  for (j in seq_along(created)) {
    
    v <- rep(NA_real_, terra::ncell(template)) # create NAs for every cell
    vals <- as.numeric(out_layers[[ created[j] ]]) # fetch backfilled values from out_layers list
    
    # write correctly-aligned values
    v[backfill_idx] <- vals
    
    result_raster[[j]] <- v
  }
  
  # assign the actual layer names
  names(result_raster) <- created
  
  # write
  out_dir <- file.path(ia_dir, "bart_models", year, sprintf("subbasin_%s", subbasin_index))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  terra::writeRaster(
    result_raster,
    file.path(out_dir, sprintf("subbasin_%d_backfill.tif", subbasin_index)),
    overwrite = TRUE
  )
  
  logp("writing %d layers to %s", length(created), out_dir)
  
  
  # save metrics for post-processing later on
  saveRDS(metrics, file.path(out_dir, sprintf("subbasin_%d_metrics.rds", subbasin_index)))
  
  # save confusion matrix information
  saveRDS(confusion, file.path(out_dir, sprintf("subbasin_%d_confusion.rds", subbasin_index)))
  
  list(subbasin = subbasin_index, ok = TRUE)

} # close train_and_backfill_subbasin_s()


