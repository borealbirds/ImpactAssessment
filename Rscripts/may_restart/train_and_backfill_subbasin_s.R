train_and_backfill_subbasin_s <- function(
    subbasin_index, year, stack_y_path, lowhf_mask_y_path, abiotic_vars, biotic_vars,
    categorical_responses, combined_poly_path, subbasins_path, ia_dir,
    neworder, quiet = FALSE
) {
  
  # read in data from paths for furrr worker
  stack_y      <- terra::rast(stack_y_path)
  lowhf_mask_y <- terra::rast(lowhf_mask_y_path)
  all_subbasins_subset <- terra::vect(subbasins_path)
  combined_poly        <- terra::vect(combined_poly_path)
  
  # isolate subbasin s
  subbasin_s <- all_subbasins_subset[subbasin_index]
  
  # crop covariate stack to subbasin
  cov_s <- 
    stack_y |>
    terra::crop(x = _, y = subbasin_s) |>
    terra::mask(x = _, mask = subbasin_s)
  
  # align low-HF to covariate stack s
  lowhf_mask_s <- terra::crop(lowhf_mask_y, cov_s, snap = "near")
  
  # for training: select pixels in covariate stack with low hf
  cov_train_s <- terra::mask(x = cov_s, mask = lowhf_mask_s)
  
  # convert covariate stack to a dataframe for modelling
  df_train <- 
    terra::as.data.frame(cov_train_s, xy=TRUE) |> 
    dplyr::as_tibble() |> 
    dplyr::rename(easting=x, northing=y) |> 
    dplyr::mutate(
      easting  = easting/1000,
      northing = northing/1000)
  
  # define predictors and responses
  abiotic_cols <- base::intersect(names(df_train), abiotic_vars$predictor)
  biotic_cols <- base::intersect(names(df_train), biotic_vars$predictor)
  biotic_cols <- stats::na.omit(biotic_cols[match(neworder, biotic_cols)])
  biotic_cols_cont <- setdiff(biotic_cols, categorical_responses)
  
  # for backfilling: find industry pixels in subbasin_s and rasterize onto cov_s grid
  industry_s <- terra::crop(combined_poly, subbasin_s)
  industry_mask_s <- terra::rasterize(industry_s, cov_s[[1]], field = 1, background = NA, touches = TRUE)     
  
  # pixel indices to backfill
  backfill_idx <- which(!is.na(terra::values(industry_mask_s)))
  if (!length(backfill_idx)) {
    if (!quiet) message("No industry pixels in subbasin ", subbasin_index, " for year ", year, "; skipping.")
    return(invisible(NULL))  # or return an empty stack if you prefer
  }
  
  # later we'll update this dataframe with backfilled features 
  # following the order of `neworder`
  df_backfill <- 
    terra::as.data.frame(cov_s, xy = TRUE, na.rm = FALSE) |>
    dplyr::as_tibble() |>
    dplyr::rename(easting = x, northing = y) |>
    dplyr::mutate(
      easting  = easting/1000,
      northing = northing/1000) |> 
    dplyr::slice(backfill_idx)
  
  # set output directory for current year
  out_dir <- file.path(ia_dir, "xgboost_models",
                       sprintf("year=%d", year),
                       sprintf("subbasin=%s", subbasin_index))
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # collect backfilled layers here
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
    predictors_train <- c(abiotic_cols, intersect(setdiff(biotic_cols_cont, b), names(df_train)), "easting","northing")
    X <- dplyr::select(df_train[idx, , drop = FALSE], all_of(predictors_train))
    X <- X[, colSums(!is.na(X)) > 0, drop = FALSE] # drop columns with all NAs
    if (!ncol(X)) next # skip if no predictor columns after drops
    train_cols <- colnames(X) # keep track of variables that were available (prediction needs to be limited to these)
    
    
    # train: check if continuous or categorical
    if (!(b %in% categorical_responses)){
      
      # log transform the response to prevent negative values in predictions
      y <- log(as.numeric(df_train[[b]][idx]) + 1)
      
      # train
      dtrain_b <- xgboost::xgb.DMatrix(data = as.matrix(X), label = y, missing = NA)
      params_b <- xgboost::xgb.params(objective = "reg:squarederror", eval_metric = "rmse", max_depth = 3, eta = 0.2, nthread = 1, seed = 123)
      cv_b <-xgboost::xgb.cv(params = params_b, data = dtrain_b, nrounds=5000, nfold=5, early_stopping_rounds = 50, verbose=FALSE)
      model_b <- xgboost::xgb.train(params = params_b, data = dtrain_b, nrounds = cv_b$early_stop$best_iteration, verbose = 0)
      
      # set up backfill data
      X_backfill <- dplyr::select(df_backfill, dplyr::any_of(train_cols))
      
      # add any missing training columns as NA
      miss <- setdiff(train_cols, colnames(X_backfill))
      if (length(miss)) {
        X_backfill[miss] <- NA_real_
      }
      
      # order columns exactly like training 
      X_backfill <- X_backfill[, train_cols, drop = FALSE]
      
      # backfill
      pred <- predict(model_b, newdata = as.matrix(X_backfill))
      vals <- pmax(exp(pred) - 1, 0)
      df_backfill[[b]] <- vals  # make backfilled feature available for subsequent`predict()`
      
    } else { # if b is categorical, use the following protocol
      
      # multiclass xgboost labels need to start at 0
      # response as factor with compact levels, then 0..K-1 labels
      y_fac <- droplevels(df_train[[b]][idx])
      K <- nlevels(y_fac)
      y <- as.integer(y_fac) - 1
      
      # if there is only one land cover class, backfill with that class
      if (K < 2) {
        single_code <- as.numeric(levels(y_fac))[1]     # true code (e.g., 5 for "conifer" in SCANFI)
        vals <- rep(single_code, length(backfill_idx))  # singular prediction
        
        r <- cov_s[[1]] * NA_real_
        r[backfill_idx] <- vals
        out_layers[[b]] <- r
        
        if (!quiet) message("subbasin ", subbasin_index, " year ", year,
                            " — categorical '", b, "' has one class (", single_code,
                            "); filled without training.")
        next
      }
      
      # train (don't use as.matrix on X as we did with continuous features)
      dtrain_b <- xgboost::xgb.DMatrix(data = as.matrix(X), label = y, missing = NA)
      params_b <- xgboost::xgb.params(objective = "multi:softmax", eval_metric = "mlogloss", max_depth = 3, eta = 0.2, num_class = K, nthread = 1, seed = 123)
      cv_b <-xgboost::xgb.cv(params = params_b, data = dtrain_b, nrounds=5000, nfold=5, early_stopping_rounds = 50, verbose=FALSE)
      model_b <- xgboost::xgb.train(params = params_b, data = dtrain_b, nrounds = cv_b$early_stop$best_iteration, verbose = 0)
      
      # backfill (categoricals do not serve as predictors in the next iteration..too complicated)
      # set up backfill data
      X_backfill <- dplyr::select(df_backfill, dplyr::any_of(train_cols))
      
      # add any missing training columns as NA
      miss <- setdiff(train_cols, colnames(X_backfill))
      if (length(miss)) {
        X_backfill[miss] <- NA_real_
      }
      
      # order columns exactly like training 
      X_backfill <- X_backfill[, train_cols, drop = FALSE]
      pred_idx_1k <- as.integer(predict(model_b, newdata = as.matrix(X_backfill))) + 1  # 1..K
      
      # map predicted indexes to the true class codes (using training factor levels)
      vals <- as.numeric(levels(y_fac))[pred_idx_1k]  
      
    } # close if/else (continuous vs categorical training and backfilling)
    
    # make a layer for biotic_cols[b] for industry footprint pixels
    r <- cov_s[[1]] * NA_real_
    r[backfill_idx] <- vals
    out_layers[[b]] <- r
    
  } # close for loop
  
  # after backfilling all covariates, build and save new stack for subbasin_s
  # keep only the layers that were actually created (remove NULL layers)
  keep <- vapply(out_layers, inherits, logical(1), what = "SpatRaster")
  out_stack <- terra::rast(out_layers[keep]) 
  names(out_stack) <- names(out_layers[keep])
  
  terra::writeRaster(out_stack, file.path(out_dir, sprintf("backfilled_stack_subbasin-%03d.tif", subbasin_index)), overwrite = TRUE)
  
  if (!quiet) {
    message("finished subbasin ", subbasin_index, " for year ", year,
            " (", nlyr(out_stack), " layers written)")}
  
} # close train_and_backfill_subbasin_s()