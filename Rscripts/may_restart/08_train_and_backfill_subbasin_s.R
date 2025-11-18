train_and_backfill_subbasin_s <- function(
    subbasin_index, 
    year, 
    stack_y,     # SpatRaster (already loaded)
    lowhf_mask,  # SpatRaster (already loaded)
    highhf_mask, # SpatRaster (already loaded)
    all_subbasins_subset, # SpatVector(already loaded)
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
  
  
  # for logging progress
  logfile <- file.path(ia_dir, "logs", sprintf("Y%d_S%03d.log", year, subbasin_index))
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
  
  # for training: align low-HF to covariate stack s
  lowhf_mask_s <- 
    terra::crop(x = lowhf_mask, y = cov_s) |> 
    terra::resample(x = _, y = cov_s, method = "near")
  
  # for backfilling: find high HF pixels in subbasin_s
  highhf_mask_s <- 
    terra::crop(x = highhf_mask, y = subbasin_s) |> 
    terra::resample(x = _, y = cov_s, method = "near") |> 
    terra::mask(x = _, mask = cov_s[[1]])
  
  # for training: select pixels in covariate stack with low hf
  cov_train_s <- terra::mask(x = cov_s, mask = lowhf_mask_s)
  
  # convert covariate stacks to a dataframe for training/backfilling
  # rows are pixels, columns are covariates
  # df_train will have NAs filtered per covariate using `idx` in `for (b in biotic_cols) {`
  # df_backfill will have NAs filtered globally   using `backfill_idx` below
  df_train    <- terra::as.data.frame(cov_train_s, xy=TRUE, na.rm = FALSE, cells = TRUE)
  df_backfill <- terra::as.data.frame(cov_s, xy = TRUE, na.rm = FALSE, cells = TRUE) 
  
  # pixel indices in highhf_mask_s to backfill
  backfill_idx <- which(!is.na(terra::values(highhf_mask_s)))
  if (!length(backfill_idx)) {
    if (!quiet) message("No industry pixels in subbasin ", subbasin_index, " for year ", year, "; skipping.")
    return(invisible(NULL))
  }
  
  # subset rows to high-HF pixels
  df_backfill <- df_backfill[backfill_idx, , drop = FALSE]
  
  # coordinate scaling: learn on training, apply to backfill
  xy_train_scaled <- scale(df_train[, c("x","y")])
  cx <- attr(xy_train_scaled, "scaled:center")
  sx <- attr(xy_train_scaled, "scaled:scale")
  
  df_train[, c("x","y")] <- xy_train_scaled
  df_backfill[, c("x","y")] <- sweep(df_backfill[, c("x","y")], 2, cx, "-")
  df_backfill[, c("x","y")] <- sweep(df_backfill[, c("x","y")], 2, sx, "/")
  
  logp("backfill rows=%d", nrow(df_backfill))
  
  # define predictors and responses
  abiotic_cols <- intersect(names(df_train), abiotic_vars$predictor)
  biotic_cols  <- intersect(names(df_train), biotic_vars$predictor)
  biotic_cols  <- na.omit(biotic_cols[match(neworder, biotic_cols)])    # enforce hierarchy
  biotic_cols_cont <- setdiff(biotic_cols, categorical_responses)
  
  # collect backfilled layers in `out_layers` list 
  # each element will be a numeric of the mean posterior values of a backfilled covariate 
  out_layers <- list()

  # store BART metrics here
  metrics <- list()
  
  # train a model for each vegetation feature
  # backfill each feature where human footprint is high
  for (b in biotic_cols) {
    
    t0 <- proc.time()[3]
    logp("[%s] fit (%s)", b, if (b %in% categorical_responses) "categorical" else "continuous")
    
    # training data: keep only the rows where biotic feature `b` is observed
    idx <- which(!is.na(df_train[[b]]))
    if (length(idx) == 0) next  # skip `b` if absent from this subbasin
    
    # biotic predictors that precede `b` in `neworder`
    if (!(b %in% categorical_responses)){
      b_before <- biotic_cols_cont[seq_len(match(b, biotic_cols_cont) - 1)]
    } else {
      b_before <- biotic_cols_cont 
    }
    
    # the predictors for biotic_cols[b] are abiotic_cols, 
    # biotic_cols (except b and those before it in `neworder`), and lat/long
    # note: categorical biotic features excluded from predicting continuous biotic features
    predictors <- c(abiotic_cols, b_before, "x", "y")
    
    # subset the global predictors to those in the current subbasin
    predictors <- intersect(predictors, names(df_train)) 
    
    # training data frame (i.e. predictors) for b
    # columns are abiotic_cols, biotic_cols (except b), lat/long
    # rows are pixel locations with low HF
    df_train_bart <- df_train[idx, predictors, drop = FALSE]
    
    # drop predictors that are all NAs
    df_train_bart <- df_train_bart[, colSums(!is.na(df_train_bart)) > 0]
    logp("train rows=%d", nrow(df_train_bart))
    logp("train columns=%d", ncol(df_train_bart))
    
    # drop predictors with zero variance
    col_sd <- sapply(df_train_bart, function(x) sd(as.numeric(x), na.rm = TRUE))
    df_train_bart <- df_train_bart[, col_sd > 0, drop = FALSE] 
    
    # subset backfill dataframe to the same abiotic_cols, biotic_cols, lat/long
    # because we can't backfill covariates that we didn't train models for
    df_backfill_bart  <- df_backfill[, colnames(df_train_bart), drop = FALSE]
    
    # check that the backfill dataframe doesn't have any all-NA predictors
    # if it does, remove them from both training and backfill dataframes
    df_backfill_bart <- df_backfill_bart[, colSums(!is.na(df_backfill_bart )) > 0]
    df_train_bart <- df_train_bart[, colnames(df_backfill_bart), drop = FALSE]
    
    # train: check if continuous or categorical
    if (!(b %in% categorical_responses)){
      
      # define continuous response, and subset pixels to non-NAs for a given `b`
      # we use `df_train` instead of `df_train_bart` because the latter is purposely missing `b`
      y <- as.numeric(df_train[[b]][idx])

      # reproducible per (year, subbasin, covariate)
      set.seed(abs(as.integer(sprintf("%d%03d", subbasin_index, which(biotic_cols==b)))))
      
      # train and predict with BART
      fit <- BART::gbart(x.train = as.matrix(df_train_bart), 
                         y.train = y, 
                         x.test = as.matrix(df_backfill_bart),
                         type = "wbart",
                         k = 3, #shrinkage
                         ntree = 50L, 
                         ndpost = 700L, 
                         nskip = 300L, 
                         sparse = TRUE, # sampler focuses on informative predictors (not all predictors treated as informative)
                         sigest  = sd(y), # the rough error standard deviation used in the prior
                         sigdf = 3) 
      
      # estimate posterior mean and sd
      # `draws` is a matrix: rows = posterior draws, columns = pixels
      draws <- fit$yhat.test # yhat.test are the predicted values of b at high HF locations
      b_mean <- as.numeric(fit$yhat.test.mean) 
      b_sd   <- apply(draws, 2, sd) # uncertainty in model-predicted values of b per pixel
      
      # write for using in next iteration (following `neworder`), and for outputs
      df_backfill[[b]] <- b_mean
      out_layers[[paste0(b, "_mean")]] <- b_mean
      out_layers[[paste0(b, "_sd")]]   <- b_sd
      
      # posterior predictive check: get posterior predictions for training data
      post_pred_train <- fit$yhat.train # yhat.train are the predicted values of b at low HF locations (values are from the fitted model, not the training dataset)
      ppc_mean <- mean(colMeans(post_pred_train))
      ppc_sd <- mean(apply(post_pred_train, 2, sd))
      
      # p is the probability the test statistic in a replicated data set exceeds that in the original data
      # Calculate Bayesian p-value (proportion of times observed statistic > predicted)
      obs_mean <- mean(y)
      obs_sd <- sd(y)
      
      # compare each posterior draw with observed statistics
      p_value_mean <- mean(apply(post_pred_train, 1, function(x) mean(x) > obs_mean))
      p_value_sd <- mean(apply(post_pred_train, 1, function(x) sd(x) > obs_sd))
      
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
        fit, y,
        covariate = b,
        subbasin  = subbasin_index,
        year      = year,
        top_var   = top_var)
      
    } else { # if b is categorical, use the following protocol
      
      # BART::mbart needs the landcover classes as a consecutive sequence of integers 1..K' 
      # but actual classes could be any arbitrary non-consecutive sequence of integers 1..K
      # so, we convert 1..K -> 1..K' -> train -> backfill as 1..K' -> convert back to 1..K
      lut <- levels(cov_s[[b]])[[1]]$value  # LookUpTable: *possible* (not actual) original cell values
      y_codes <- lut[ as.integer(df_train[[b]][idx]) ] # the vector of original landcover codes for the selected training pixels. don't use df_train_bart because it's missing `b`
      
      # make sure y isn’t empty 
      if (!length(y_codes)) { logp("[%s] skip: empty y_codes", b); next } 
      
      # land cover classes that are *actually* present in the training subset
      # this is built from `df_train` so it may have some classes dropped during NA filtering
      # therefore `fit$K` may have fewer than K classes
      present <- sort(unique(y_codes))
      K <- length(present) 
      if (K == 0) { logp("[%s] skip: no classes present", b); next }
      
      # pre-allocate outputs
      b_mode    <- rep(NA_real_, nrow(df_backfill))
      b_maxprob <- rep(NA_real_, nrow(df_backfill))
      b_entropy <- rep(NA_real_, nrow(df_backfill))
      per_class <- list()  # will fill in if K >= 2
      
      # build a forward and inverse map only over present classes (size K')
      fwd_map  <- setNames(seq_along(present), present)   # original code -> 1..K'
      inv_map  <- setNames(present, seq_along(present))   # 1..K' -> original code
      
      # determine complexity of landcover classes
      if (K < 2) {
        
        # trivial: only one class in training data
        b_mode    <- rep(present[1], nrow(df_backfill)) # the one and only landcover class present
        b_maxprob <- rep(1, nrow(df_backfill)) 
        b_entropy <- rep(0, nrow(df_backfill))

        nm <- paste0(b, "_prob_", present[1])
        per_class[[nm]] <- rep(1, nrow(df_backfill))
        
      } else {
        
        # training labels mapped to 1..K'
        y_int <- unname(fwd_map[y_codes])  
        
        # reproducible per (year, subbasin, covariate)
        set.seed(abs(as.integer(sprintf("%d%03d", subbasin_index, which(biotic_cols==b)))))
        
        # multinomial bart keeps every 10th posterior sample (gbart keeps every sample)
        # so mbart will take at least 10 times longer
        # gbart assumes continuous residuals → conjugate priors → less autocorrelation → no thinning
        # mbart uses latent variables → non-conjugate updates → more autocorrelation → needs thinning
        fit <- BART::mbart(x.train = as.matrix(df_train_bart), 
                           y.train = y_int, 
                           x.test  = as.matrix(df_backfill_bart),
                           type = "pbart",
                           k = 3,
                           ntree = 50L,
                           ndpost = 500L,
                           nskip = 200L,
                           keepevery = 5L,
                           printevery = 450,
                           sparse = TRUE)
      
        # fit$prob.test has dimensions: posterior density (ndpost) x (# high HF pixels x # classes in training set)
        # therefore, for every high HF pixel we have a probability estimate for every class learned in from the training set
        # and this is multiplied by the density of the posterior distribution (ndpost)
        # however, we want to split from 2-D prob.test to 3-D `prob_array` (see below)
        ndpost <- nrow(fit$prob.test) # density of posterior distribution
        K_model <- fit$K # number of classes in df_train_bart
        npixels <- as.integer(ncol(fit$prob.test) / K_model) # number of high HF pixels
        
        # reshape `prob.test` into draws x pixels x model-classes
        prob_array <- array(fit$prob.test, dim = c(ndpost, npixels, K_model))
        
        # get the mean probability per pixel x class
        # each row is a pixel, each column is a class with probability of said class
        # class_probs_model is indexed in MBART’s internal class order which we'll convert to `class_probs` below
        class_probs_model <- apply(prob_array, c(2,3), mean) #
        
        # which model column corresponds to each desired code?
        cols <- match(present, fit$cats)   
        
        # build final class_probs in correct ecological order
        class_probs <- matrix(0, nrow = npixels, ncol = length(present))
        
        # match each ecological class to MBART’s output column
        # insert probability 0 for any class not present in training
        # produce a probability matrix in the correct ecological order
        for (i in seq_along(present)) {
          ci <- cols[i]
          if (!is.na(ci)) {
            class_probs[, i] <- class_probs_model[, ci]
          } else {
            # probability = 0 if class not seen in the training subset 
            class_probs[, i] <- 0
          }
        }
        
        # most likely ecological class
        pred_idx  <- max.col(class_probs, ties.method = "first")
        b_mode    <- present[pred_idx]   # back to original ecological codes
        
        # uncertainty
        eps <- 1e-12
        b_entropy <- -rowSums(class_probs * log(pmax(class_probs, eps)))
        
        # probability of the chosen class
        b_maxprob <- apply(class_probs, 1, max)
        
        # per-class probabilities 
        for (k in seq_len(K)) { 
          nm <- paste0(b, "_prob_", present[k]) 
          per_class[[nm]] <- class_probs[, k] 
          }
        
      } # close if single or multi-class
      
      # outputs
      df_backfill[[b]] <- b_mode
      out_layers[[b]]  <- b_mode
      out_layers[[paste0(b, "_maxprob")]]  <- b_maxprob
      out_layers[[paste0(b, "_entropy")]]  <- b_entropy
      
      # add per-class probability layers
      for (nm in names(per_class)) out_layers[[nm]] <- per_class[[nm]]
      
      if (K >= 2){
      metrics[[length(metrics) + 1L]] <- collect_metrics_mbart(
        fit,
        y         = y_int,
        X_train   = df_train_bart,
        covariate = b,
        subbasin  = subbasin_index,
        year      = year)
      } else {
        metrics[[length(metrics) + 1L]] <- NULL
      }
      
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
  names(result_raster) <- created
  
  # fill only the high-HF cells
  for (j in seq_along(created)) {
    v <- rep(NA_real_, terra::ncell(template))
    vals <- out_layers[[ created[j] ]] # already original codes (if landcover classes)
    v[backfill_idx] <- vals
    result_raster[[j]] <- v
  }
  
  # write
  out_dir <- file.path(ia_dir, "bart_models", year, sprintf("subbasin_%s", subbasin_index))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  terra::writeRaster(
    result_raster,
    file.path(out_dir, sprintf("subbasin_%03d_backfill.tif", subbasin_index)),
    overwrite = TRUE
  )
  
  logp("writing %d layers to %s", length(created), out_dir)
  
  invisible(TRUE)

} # close train_and_backfill_subbasin_s()