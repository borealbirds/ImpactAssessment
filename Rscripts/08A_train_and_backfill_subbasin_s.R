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
  quiet = FALSE,
  cc
) {
  
  # source BART metrics summary functions
  source(file.path(ia_dir, "Rscripts", "08B_deploy_gbart.R"))
  source(file.path(ia_dir, "Rscripts", "08B_deploy_mbart.R"))
  
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
  for (b in biotic_cols) {
    
    t0 <- proc.time()[3]
    logp("[%s] in process: (%s)", b, if (b %in% categorical_responses) "categorical" else "continuous")
    
    # training data: keep only the rows where biotic feature `b` is observed
    idx <- which(!is.na(df_train[[b]]))
    
    # check that sample size of b > 1, and for b invariance
    if (length(idx) == 0 || length(unique(df_train[[b]][idx])) < 2) {
      
      const_val <- unique(df_train[[b]][idx])
      # if idx has no valid rows (empty), define const_val as NA (as it would be in V5 rasters)
      if (length(const_val) == 0) const_val <- NA

      # fill outputs — naming must match deploy_gbart/mbart layer name conventions so that
      # 11_premosaic and 12C can locate the layers by their expected names
      df_backfill[[b]] <- const_val
      if (b %in% categorical_responses) {
        out_layers[[b]]                     <- rep(const_val, nrow(df_backfill))
        out_layers[[paste0(b, "_maxprob")]] <- rep(1,          nrow(df_backfill))
        out_layers[[paste0(b, "_entropy")]] <- rep(0,          nrow(df_backfill))
      } else {
        out_layers[[paste0(b, "_mean")]] <- rep(const_val, nrow(df_backfill))
        out_layers[[paste0(b, "_sd")]]   <- rep(0,          nrow(df_backfill))
      }
      
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
    
    # select BART method for a continuous or categorical feature
    if (!(b %in% categorical_responses)){
      
      res <- deploy_gbart(
        b               = b,
        df_train        = df_train,
        df_train_bart   = df_train_bart,
        df_backfill     = df_backfill,
        df_backfill_bart= df_backfill_bart,
        idx             = idx,
        biotic_cols     = biotic_cols,
        subbasin_index  = subbasin_index,
        year            = year,
        metrics         = metrics,
        out_layers      = out_layers,
        logp            = logp,
        cc              = cc,
        ia_dir          = ia_dir
      )
      
      df_backfill <- res$df_backfill
      metrics     <- res$metrics
      out_layers  <- res$out_layers
      
    } else { 
      
      res <- deploy_mbart(
        b               = b,
        df_train        = df_train,
        df_train_bart   = df_train_bart,
        df_backfill     = df_backfill,
        df_backfill_bart= df_backfill_bart,
        idx             = idx,
        biotic_cols     = biotic_cols,
        subbasin_index  = subbasin_index,
        year            = year,
        metrics         = metrics,
        confusion       = confusion,
        out_layers      = out_layers,
        logp            = logp,
        cc              = cc,
        ia_dir          = ia_dir
      )
      
      df_backfill <- res$df_backfill
      metrics     <- res$metrics
      confusion   <- res$confusion
      out_layers  <- res$out_layers
    
    } # close backfilling protocol
    
    logp("[%s] fit done in %.1fs", b, proc.time()[3] - t0)

  } # close loop over biotic_cols
  
  
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
  out_dir <- file.path(ia_dir, "data", "derived_data", "bart_models", year, sprintf("subbasin_%s", subbasin_index))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  terra::writeRaster(
    result_raster,
    file.path(out_dir, sprintf("subbasin_%d_backfill.tif", subbasin_index)),
    overwrite  = TRUE,
    datatype   = "FLT4S",
    gdal       = c("COMPRESS=DEFLATE", "PREDICTOR=3", "ZLEVEL=6", "BIGTIFF=YES")
  )
  
  logp("writing %d layers to %s", length(created), out_dir)
  
  
  # save metrics for post-processing later on
  saveRDS(metrics, file.path(out_dir, sprintf("subbasin_%d_metrics.rds", subbasin_index)))
  
  # save confusion matrix information
  saveRDS(confusion, file.path(out_dir, sprintf("subbasin_%d_confusion.rds", subbasin_index)))
  
  list(subbasin = subbasin_index, ok = TRUE)

} # close train_and_backfill_subbasin_s()


