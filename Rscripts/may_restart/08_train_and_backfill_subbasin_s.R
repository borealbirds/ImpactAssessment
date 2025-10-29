train_and_backfill_subbasin_s <- function(
    subbasin_index, 
    year, 
    stack_y,     # SpatRaster (already loaded)
    lowhf_mask,  # SpatRaster (already loaded)
    highhf_mask, # SpatRaster (already loaded)
    all_subbasins_subset, # SpatVector(already loaded)
    abiotic_vars, # tibble with $predictor
    biotic_vars,  # tibble with $predictor
    categorical_responses, # character vector of layer names
    ia_dir,
    neworder, 
    quiet = FALSE
) {
  
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
  highhf_s <- 
    terra::crop(x = highhf_mask, y = subbasin_s) |> 
    terra::resample(x = _, y = cov_s, method = "near")
  
  # for training: select pixels in covariate stack with low hf
  cov_train_s <- terra::mask(x = cov_s, mask = lowhf_mask_s)
  
  # convert covariate stacks to dataframes for training/backfilling
  df_train    <- terra::as.data.frame(cov_train_s, xy=TRUE)
  df_backfill <- terra::as.data.frame(cov_s, xy = TRUE) # will subset to high HF using backfill_idx below
  
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
  
  # pixel indices in highhf_s to backfill
  backfill_idx <- which(!is.na(terra::values(highhf_s)))
  if (!length(backfill_idx)) {
    if (!quiet) message("No industry pixels in subbasin ", subbasin_index, " for year ", year, "; skipping.")
    return(invisible(NULL))
  }
  
  # subset rows to high-HF pixels
  df_backfill <- df_backfill[backfill_idx, , drop = FALSE]
  
  # collect backfilled layers in `out_layers`
  out_layers <- vector("list", length(biotic_cols))  
  names(out_layers) <- biotic_cols
  
  # train a model for each vegetation feature, and backfill
  for (b in biotic_cols) {
    
    t0 <- proc.time()[3]
    logp("[%s] fit (%s)", b, if (b %in% categorical_responses) "categorical" else "continuous")
    
    # training data: keep only the rows where the response is observed
    idx <- which(!is.na(df_train[[b]]))
    if (length(idx) == 0) next  # skip this covariate
    
    # predictors for biotic_cols[b]: abiotic covariates, biotic covariates (except b), lat/long
    # categorical biotic features excluded
    predictors <- c(abiotic_cols, setdiff(biotic_cols_cont, b), "x", "y")
    predictors <- intersect(predictors, names(df_train))
    
    X_train    <- df_train[idx, predictors, drop = FALSE]
    keep       <- complete.cases(X_train)
    X_train    <- X_train[keep, , drop = FALSE]
    X_backfill <- df_backfill[, predictors, drop = FALSE]
    
    if (!ncol(X_train)) next
  
    # train: check if continuous or categorical
    if (!(b %in% categorical_responses)){
      
      # define continuous response
      y <- as.numeric(df_train[[b]][idx])[keep]
      
      # train and predict with BART
      fit <- BART::gbart(x.train = as.matrix(X_train), y.train = y, x.test  = as.matrix(X_backfill))
      
      # posterior mean 
      vals <- as.numeric(fit$yhat.test.mean)  
      
      # write into the working backfill rows so later targets can use it
      df_backfill[[b]] <- vals
      out_layers[[b]] <- vals
      
    } else { # if b is categorical, use the following protocol
      
      # BART::mbart expects class integers 1..K
      # get original codes, even if the raster is a factor
      if (terra::is.factor(cov_s[[b]])) {
        lut <- levels(cov_s[[b]])[[1]]$value  # original cell values
        y_codes <- lut[ as.integer(df_train[[b]][idx]) ]
      } else {
        y_codes <- as.integer(df_train[[b]][idx])
      }
      
      # land cover classes actually present in the training subset
      present <- sort(unique(y_codes))
      
      # build a forward and inverse map only over present classes (size K')
      fwd_map  <- setNames(seq_along(present), present)   # original code -> 1..K'
      inv_map  <- setNames(present, seq_along(present))   # 1..K' -> original code
      
      K <- length(present)
      if (K < 2) {
        # trivial: only one class in training data
        vals_codes <- rep.int(present[1], nrow(df_backfill))  # keep original code
        out_layers[[b]]  <- vals_codes
        df_backfill[[b]] <- vals_codes
      } else {
        y_int <- unname(fwd_map[as.character(y_codes)])  # 1..K'
        
        fit <- BART::mbart(x.train = as.matrix(X_train), y.train = y_int, x.test  = as.matrix(X_backfill))
      
        # prob.test: draws x ntest x K; average over draws -> ntest x K'
        class_probs <- apply(fit$prob.test, c(2, 3), mean)
        pred_idx <- max.col(class_probs, ties.method = "first")  # 1..K'
        vals_codes <- as.integer(inv_map[as.character(pred_idx)]) # back to original codes  
        
        out_layers[[b]]  <- vals_codes
        df_backfill[[b]] <- vals_codes
      } # close if single or multi-class
      
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