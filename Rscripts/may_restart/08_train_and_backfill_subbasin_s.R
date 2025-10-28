train_and_backfill_subbasin_s <- function(
    subbasin_index, 
    year, 
    stack_y, 
    lowhf_mask, 
    highhf_mask,
    all_subbasins_subset_path,
    abiotic_vars, 
    biotic_vars,
    categorical_responses, 
    ia_dir,
    neworder, 
    quiet = FALSE
) {
  
  # read in data from paths 
  stack_y      <- terra::rast(stack_y_path)
  lowhf_mask   <- terra::rast(lowhf_mask_path)
  highhf_mask  <- terra::rast(highhf_mask_path)
  all_subbasins_subset <- terra::vect(all_subbasins_subset_path)

  # isolate subbasin s
  subbasin_s <- all_subbasins_subset[subbasin_index]
  
  # crop covariate stack to subbasin
  cov_s <- 
    stack_y |>
    terra::crop(x = _, y = subbasin_s) |>
    terra::mask(x = _, mask = subbasin_s)
  
  # align low-HF to covariate stack s
  lowhf_mask_s <- terra::crop(lowhf_mask, cov_s, snap = "near")
  
  # for training: select pixels in covariate stack with low hf
  cov_train_s <- terra::mask(x = cov_s, mask = lowhf_mask_s)
  
  # convert covariate stack to a dataframe for modelling
  df_train <- terra::as.data.frame(cov_train_s, xy=TRUE) 
  df_train$x  <- df_train$x  / 1000 # longtidue
  df_train$y <- df_train$y / 1000 # latitude
  
  # define predictors and responses
  abiotic_cols <- intersect(names(df_train), abiotic_vars$predictor)
  biotic_cols <- intersect(names(df_train), biotic_vars$predictor)
  biotic_cols <- na.omit(biotic_cols[match(neworder, biotic_cols)])
  biotic_cols_cont <- setdiff(biotic_cols, categorical_responses)
  
  # for backfilling: find high HF pixels in subbasin_s
  highhf_s <- terra::crop(highhf_mask, subbasin_s, snap = "near")

  # pixel indices in highhf_s to backfill
  backfill_idx <- which(!is.na(terra::values(highhf_s)))
  if (!length(backfill_idx)) {
    if (!quiet) message("No industry pixels in subbasin ", subbasin_index, " for year ", year, "; skipping.")
    return(invisible(NULL))
  }
  
  # create dataframe where rows are pixels, columns are abiotic and biotic covariates
  # later we'll update this dataframe with backfilled features following the order of `neworder`
  df_backfill <- terra::as.data.frame(cov_s, xy = TRUE, na.rm = FALSE)
  df_backfill$x <- df_backfill$x  / 1000
  df_backfill$y <- df_backfill$y / 1000
  df_backfill <- df_backfill[backfill_idx, , drop = FALSE] # subset dataframe to pixels with high HF
  
  # set output directory for current year
  out_dir <- file.path(ia_dir, "bart_models", year, sprintf("subbasin_%s", subbasin_index))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # collect backfilled layers in `out_layers`
  out_layers <- vector("list", length(biotic_cols))  
  names(out_layers) <- biotic_cols
  
  # train a model for each vegetation feature, and backfill
  for (j in seq_along(biotic_cols)) {
    
    b <- biotic_cols[j] # for testing: b<-biotic_cols[1]
    
    # training data:
    # keep only the rows where the response is observed
    idx <- which(!is.na(df_train[[b]]))
    if (!length(idx)) next  # skip this covariate
    
    # predictors for biotic_cols[b]
    predictors_train <- c(abiotic_cols, intersect(setdiff(biotic_cols_cont, b), names(df_train)), "x","y")
    X <- df_train[idx, intersect(predictors_train, names(df_train)), drop = FALSE]
    X <- X[, colSums(!is.na(X)) > 0, drop = FALSE] # drop columns with all NAs
    if (!ncol(X)) next # skip if no predictor columns after drops
    train_cols <- colnames(X) # keep track of variables that were available (prediction needs to be limited to these)
    
    # set up backfill data
    X_backfill <- df_backfill[, intersect(train_cols, names(df_backfill)), drop = FALSE]
    miss <- setdiff(train_cols, colnames(X_backfill))
    if (length(miss)) X_backfill[miss] <- NA_real_
    X_backfill <- X_backfill[, train_cols, drop = FALSE]
    
    # train: check if continuous or categorical
    if (!(b %in% categorical_responses)){
      
      # log transform the response to prevent negative values in predictions
      y <- log(as.numeric(df_train[[b]][idx]) + 1)
      
      # train and predict with BART
      fit <- BART::gbart(x.train = as.matrix(X), y.train = y, x.test  = as.matrix(X_backfill))
      
      # back-transform predictions
      pred <- fit$yhat.test.mean
      vals <- pmax(exp(pred_mean) - 1, 0)
      
    } else { # if b is categorical, use the following protocol
      
      y_factor <- droplevels(as.factor(train_data[[b]]))
      n_classes <- nlevels(y_factor)
      
      # single class modelling
      if (n_classes < 2) {
        vals <- rep(as.integer(levels(y_factor)[1]), nrow(df_backfill))
      } else {
        X_train <- data.matrix(train_data[, predictors])
        X_test <- data.matrix(df_backfill[, predictors])
      
        fit <- BART::mbart(x.train = X_train, y.train = y_factor, x.test = X_test)
        
        # get most probable class
        class_probs <- apply(fit$prob.test, c(2, 3), mean)
        vals <- apply(class_probs, 1, which.max)
      }
    }
      
    # store results
    out_layers[[b]] <- vals
    
  } # close for loop
  
  # after backfilling all covariates, build and save new stack for subbasin_s
  # keep only the layers that were actually created (remove NULL layers)
  result_raster <- cov_s[[1]]
  for (b in names(out_layers)) {
    result_raster <- c(result_raster, result_raster * NA)
    terra::values(result_raster[[nlyr(result_raster)]])[backfill_idx] <- out_layers[[b]]
  }
  names(result_raster) <- names(out_layers)
  
  terra::writeRaster( result_raster, file.path(out_dir, sprintf("subbasin_%03d_backfill.tif", subbasin_index)), overwrite = TRUE)
  
  return(invisible(TRUE))
  
} # close train_and_backfill_subbasin_s()