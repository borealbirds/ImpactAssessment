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
  
  # source BART metrics summary function
  source(file.path(getwd(), "Rscripts", "may_restart", "09_collect_metrics.R"))
  
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
  
  logp("train rows=%d; backfill rows=%d", nrow(df_train), nrow(df_backfill))
  
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
    b_before <- biotic_cols_cont[seq_len(match(b, biotic_cols_cont) - 1)]
    
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
    
    # subset backfill dataframe to the same abiotic_cols, biotic_cols, lat/long
    df_backfill_bart  <- df_backfill[, colnames(df_train_bart), drop = FALSE]
    
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
                         k=2, #shrinkage
                         ntree = 75L, 
                         ndpost = 700L, 
                         nskip = 300L, 
                         sparse = TRUE, # sampler focuses on informative predictors (not all predictors treated as informative)
                         sigest  = sd(y), # the rough error standard deviation used in the prior
                         sigdf=1.5) 
      
      # estimate posterior mean and sd
      # `draws` is a matrix: rows = posterior draws, columns = pixels
      draws <- fit$yhat.test
      pred_mean <- as.numeric(fit$yhat.test.mean) 
      pred_sd   <- apply(draws, 2, sd)
      
      # allocate full-length vectors (including NAs where inputs were NA)
      b_mean <- rep(NA_real_, nrow(df_backfill))
      b_sd   <- rep(NA_real_, nrow(df_backfill))
      b_mean <- pred_mean
      b_sd   <- pred_sd
      
      # write for hierarchical use + outputs
      df_backfill[[b]] <- b_mean
      out_layers[[paste0(b, "_mean")]] <- b_mean
      out_layers[[paste0(b, "_sd")]]   <- b_sd
      
      # posterior predictive check: get posterior predictions for training data
      post_pred_train <- fit$yhat.train
      ppc_mean <- mean(colMeans(post_pred_train))
      ppc_sd <- mean(apply(post_pred_train, 2, sd))
      
      # p = the probability the test statistic in a replicated data 
      # set exceeds that in the original data
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
      metrics[[length(metrics) + 1L]] <- collect_metrics(
        fit, y,
        covariate = b,
        subbasin  = subbasin_index,
        year      = year,
        top_var   = top_var)
      
    } else { # if b is categorical, use the following protocol
      
      # BART::mbart expects class integers 1..K
      # get original codes, even if the raster is a factor
      if (terra::is.factor(cov_s[[b]])) {
        lut <- levels(cov_s[[b]])[[1]]$value  # original cell values
        y_codes_all <- lut[ as.integer(df_train[[b]][idx]) ]
      } else {
        y_codes_all <- as.integer(df_train[[b]][idx])
      }
      
      # align y with X_train 
      y_codes <- y_codes_all
      if (!length(y_codes)) { logp("[%s] skip: empty y_codes after keep1", b); next } # make sure y isn’t empty 
      
      # land cover classes actually present in the training subset
      present <- sort(unique(y_codes))
      
      # pre-allocate outputs
      b_mode    <- rep(NA_real_, nrow(df_backfill))
      b_maxprob <- rep(NA_real_, nrow(df_backfill))
      b_entropy <- rep(NA_real_, nrow(df_backfill))
      per_class <- list()  # will fill below when K>=2
      
      # build a forward and inverse map only over present classes (size K')
      fwd_map  <- setNames(seq_along(present), present)   # original code -> 1..K'
      inv_map  <- setNames(present, seq_along(present))   # 1..K' -> original code
      K <- length(present)
      if (K == 0) { logp("[%s] skip: no classes present", b); next }
      
      if (K < 2) {
        # trivial: only one class in training data
        b_mode[keep2]    <- present[1]
        b_maxprob[keep2] <- 1
        b_entropy[keep2] <- 0

        nm <- paste0(b, "_prob_", present[1])
        v <- rep(NA_real_, nrow(df_backfill)); v[keep2] <- 1
        per_class[[nm]] <- v
        
      } else {
        y_int <- unname(fwd_map[as.character(y_codes)])  # 1..K'
        
        fit <- BART::mbart(x.train = as.matrix(X_train), y.train = y_int, x.test  = as.matrix(X_backfill))
      
        # prob.test: draws x ntest x K; average over draws -> ntest x K'
        class_probs <- apply(fit$prob.test, c(2, 3), mean)
        
        # map (1..K') and map back to original codes
        pred_idx <- max.col(class_probs, ties.method = "first")  # 1..K'
        mode_codes <- as.integer(inv_map[as.character(pred_idx)]) # back to original codes  
        
        # place into full-length vectors
        b_mode[keep2]    <- mode_codes
        b_maxprob[keep2] <- apply(class_probs, 1, max)
        # entropy with small epsilon to avoid log(0)
        eps <- 1e-12
        b_entropy[keep2] <- -rowSums(class_probs * log(pmax(class_probs, eps)))
        
        # per-class prob layers (named by original code)
        for (k in seq_len(K)) {
          nm <- paste0(b, "_prob_", present[k])   # name by original code
          v  <- rep(NA_real_, nrow(df_backfill))
          v[keep2] <- class_probs[, k]
          per_class[[nm]] <- v
        }
      } # close if single or multi-class
      
      # outputs
      df_backfill[[b]] <- b_mode
      out_layers[[b]]  <- b_mode
      out_layers[[paste0(b, "_maxprob")]]  <- b_maxprob
      out_layers[[paste0(b, "_entropy")]]  <- b_entropy
      
      # add per-class probability layers
      for (nm in names(per_class)) out_layers[[nm]] <- per_class[[nm]]
      
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