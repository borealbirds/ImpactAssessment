# ---
# title: Impact Assessment: train and backfill one subbasin (sourced by 07)
# author: Mannfred Boehm
# ---

# a logger that appends time-stamped lines to `logfile`
make_logger <- function(logfile) {
  dir.create(dirname(logfile), recursive = TRUE, showWarnings = FALSE)
  function(fmt, ...) {
    cat(sprintf("[%s pid=%d host=%s] %s\n", format(Sys.time(), "%F %T"), Sys.getpid(),
                Sys.info()[["nodename"]], sprintf(fmt, ...)),
        file = logfile, append = TRUE)
  }
}


# BART design matrices for one response: the training rows' predictors and the backfill's,
# with the same columns in the same order. used by the amount models (gbart/mbart) and the
# blank model (pbart), each with its own training rows
prep_bart_design <- function(df_train_bart, df_backfill, logp) {

  # drop predictors that are all NA, or constant (sd is NA with one non-NA value: CAfire, subbasin 57)
  df_train_bart <- df_train_bart[, colSums(!is.na(df_train_bart)) > 0, drop = FALSE]
  logp("training pixels = %d", nrow(df_train_bart))
  logp("training predictors = %d", ncol(df_train_bart))
  col_sd <- sapply(df_train_bart, function(x) sd(as.numeric(x), na.rm = TRUE))
  df_train_bart <- df_train_bart[, !is.na(col_sd) & col_sd > 0, drop = FALSE]

  # drop predictors that are all NA on the backfilled pixels, unless they are partly NA on the
  # training pixels too: there an all-NA backfill is information (no trees, so no
  # SCANFIheightcv_1km), which the `_isNA` flag below carries
  df_backfill_bart <- df_backfill[, colnames(df_train_bart), drop = FALSE]
  flaggable <- colSums(is.na(df_train_bart)) > 0L & colSums(!is.na(df_train_bart)) > 0L
  df_backfill_bart <- df_backfill_bart[, colSums(!is.na(df_backfill_bart)) > 0 | flaggable, drop = FALSE]
  df_train_bart <- df_train_bart[, colnames(df_backfill_bart), drop = FALSE]

  # BART can't take NA (a pixel with one gets no posterior draw), so a predictor with NAs gets
  # the training median in their place and, where it varies in training, a `<col>_isNA` flag
  # BART can split on (as gbm splits on NA, e.g. greenup on water). a backfilled biotic
  # predictor is NA only where pbart modelled it blank (08A below)
  na_cols <- names(df_train_bart)[colSums(is.na(df_train_bart)) > 0L | colSums(is.na(df_backfill_bart)) > 0L]
  na_cols <- setdiff(na_cols, c("x", "y"))
  for (cna in na_cols) {
    fill <- stats::median(df_train_bart[[cna]], na.rm = TRUE)
    if (!is.finite(fill)) fill <- 0
    tr_flag <- as.integer(is.na(df_train_bart[[cna]]))
    if (length(unique(tr_flag)) > 1L) {
      df_train_bart[[paste0(cna, "_isNA")]]    <- tr_flag
      df_backfill_bart[[paste0(cna, "_isNA")]] <- as.integer(is.na(df_backfill_bart[[cna]]))
    }
    df_train_bart[[cna]][is.na(df_train_bart[[cna]])]       <- fill
    df_backfill_bart[[cna]][is.na(df_backfill_bart[[cna]])] <- fill
  }
  df_backfill_bart <- df_backfill_bart[, colnames(df_train_bart), drop = FALSE]
  if (length(na_cols) > 0L)
    logp("imputed+flagged %d partial-NA predictor(s): %s", length(na_cols), paste(na_cols, collapse = ", "))

  list(train = df_train_bart, backfill = df_backfill_bart)
}


# trains BART models for every biotic covariate of one subbasin and backfills its high HF
# pixels -> bart_models/{year}/subbasin_{i}/subbasin_{i}_{backfill.tif, metrics.rds, confusion.rds}
train_and_backfill_subbasin_s <- function(subbasin_index, year, stack_y, lowhf_mask, highhf_mask,
                                          all_subbasins_subset, abiotic_vars, footprint_vars,
                                          biotic_vars, categorical_responses, water_classes, ia_dir) {

  # BART_MODELS_DIR (default bart_models) sends a test run's outputs and logs elsewhere
  bart_dir <- Sys.getenv("BART_MODELS_DIR", "bart_models")
  logp <- make_logger(file.path(ia_dir, "logs", if (bart_dir == "bart_models") "" else bart_dir,
                                sprintf("Y%d_S%d.log", year, subbasin_index)))
  on.exit(logp("done"), add = TRUE)
  logp("start")

  subbasin_s <- all_subbasins_subset[subbasin_index]
  cov_s <-
    terra::crop(stack_y, subbasin_s) |>
    terra::mask(subbasin_s)
  logp("cropped/masked; ncell=%d", terra::ncell(cov_s))

  lowhf_mask_s  <- terra::mask(terra::resample(lowhf_mask,  cov_s, method = "near"), subbasin_s)
  highhf_mask_s <- terra::mask(terra::resample(highhf_mask, cov_s, method = "near"), subbasin_s)

  # rows are pixels, columns are covariates. df_train is NA outside the low HF pixels;
  # train_rows are its low HF rows, blank covariates included (for pbart)
  df_train    <- terra::as.data.frame(terra::mask(cov_s, lowhf_mask_s), xy = TRUE, na.rm = FALSE)
  train_rows  <- which(terra::values(lowhf_mask_s, mat = FALSE) == 1)
  df_full     <- terra::as.data.frame(cov_s, xy = TRUE, na.rm = FALSE, cells = TRUE)
  backfill_idx <- which(terra::values(highhf_mask_s) == 1)
  df_backfill <- df_full[backfill_idx, , drop = FALSE]

  # centre and scale the coordinates over the subbasin's bounding box
  xy_center <- colMeans(df_full[c("x", "y")])
  xy_scale  <- apply(df_full[c("x", "y")], 2, sd)
  df_train[, c("x", "y")]    <- sweep(sweep(df_train[c("x", "y")], 2, xy_center, "-"), 2, xy_scale, "/")
  df_backfill[, c("x", "y")] <- sweep(sweep(df_backfill[, c("x", "y")], 2, xy_center, "-"), 2, xy_scale, "/")
  logp("backfill rows=%d", nrow(df_backfill))

  # the counterfactual has no industry, so the backfill is predicted with the footprint
  # covariates at 0 (OL #9b; 12F does the same for the bird model). left at their observed
  # values, beyond the training range, trees would predict vegetation as if the road were
  # still there. the models still train on the observed values, which include 0
  footprint_cols <- intersect(names(df_backfill), footprint_vars)
  df_backfill[, footprint_cols] <- 0
  logp("footprint covariates set to 0 at backfilled pixels: %s",
       if (length(footprint_cols)) paste(footprint_cols, collapse = ", ") else "none")

  abiotic_cols     <- intersect(names(df_train), abiotic_vars)
  biotic_cols      <- intersect(biotic_vars, names(df_train))   # in hierarchy order (00)
  biotic_cols_cont <- setdiff(biotic_cols, categorical_responses)

  out_layers <- list()   # backfilled layers, one numeric per layer, in output order
  metrics    <- list()   # one data frame per model and split
  confusion  <- list()   # holdout classes of the land-cover models

  for (b in biotic_cols) {

    t0 <- proc.time()[3]
    is_cat <- b %in% categorical_responses
    logp("[%s] in process: (%s)", b, if (is_cat) "categorical" else "continuous")

    # training rows: where b is observed, minus open water for a land-cover layer, so a
    # footprint pixel's counterfactual is land (production 2026-09-28 made 94% of can13
    # subbasin 72's backfill SCANFI water)
    idx <- which(!is.na(df_train[[b]]))
    if (b %in% names(water_classes)) {
      is_water <- df_train[[b]][idx] %in% water_classes[[b]]
      if (any(is_water))
        logp("[%s] %d water training pixels (class %s) left out", b, sum(is_water),
             paste(water_classes[[b]], collapse = "/"))
      idx <- idx[!is_water]
    }

    # predictors: abiotic, the continuous biotic covariates above b in the hierarchy (all of
    # them for a land-cover layer), and the coordinates
    b_before   <- if (is_cat) biotic_cols_cont else biotic_cols_cont[seq_len(match(b, biotic_cols_cont) - 1)]
    predictors <- c(abiotic_cols, b_before, "x", "y")

    # reproducible per (subbasin, covariate)
    seed <- abs(as.integer(sprintf("%d%03d", subbasin_index, which(biotic_cols == b))))

    if (length(idx) == 0 || length(unique(df_train[[b]][idx])) < 2) {

      # b is constant on the training pixels (or absent: NA, as in the V5 rasters). 11 and 12F
      # look the layers up by these names
      const_val <- unique(df_train[[b]][idx])
      if (length(const_val) == 0) const_val <- NA
      n_bf <- nrow(df_backfill)
      df_backfill[[b]] <- const_val
      if (is_cat) {
        out_layers[[b]]                     <- rep(const_val, n_bf)
        out_layers[[paste0(b, "_maxprob")]] <- rep(1, n_bf)
        out_layers[[paste0(b, "_entropy")]] <- rep(0, n_bf)
        for (d in seq_len(100L))
          out_layers[[paste0(b, "_class_", sprintf("%03d", d))]] <- rep(const_val, n_bf)
      } else {
        out_layers[[paste0(b, "_mean")]] <- rep(const_val, n_bf)
        out_layers[[paste0(b, "_sd")]]   <- rep(0, n_bf)
      }
      logp("[%s] constant in subbasin: backfill raster populated with %s", b, const_val)

    } else {

      design <- prep_bart_design(df_train[idx, predictors, drop = FALSE], df_backfill, logp)
      res <- if (is_cat) deploy_mbart(b, df_train[[b]][idx], design, seed, subbasin_index, year, logp) else
                         deploy_gbart(b, df_train[[b]][idx], design, seed, subbasin_index, year, logp)
      df_backfill[[b]] <- res$backfill
      out_layers <- c(out_layers, res$layers)
      metrics    <- c(metrics, res$metrics)
      confusion  <- c(confusion, res$confusion)
    }

    # where the training land leaves b blank, model where the backfill is blank too (08B_deploy_pbart.R).
    # the amount layers keep BART's values; <b>_blank = 1 marks the pixels 12F passes to the bird
    # model as NA, and the cascade sees them as NA
    p_blank <- deploy_pbart_blank(b, df_train[train_rows, , drop = FALSE], predictors, df_backfill, seed, logp)
    if (!is.null(p_blank)) {
      df_backfill[[b]][p_blank > 0.5] <- NA
      out_layers[[paste0(b, "_pblank")]] <- p_blank
      out_layers[[paste0(b, "_blank")]]  <- as.integer(p_blank > 0.5)
    }

    logp("[%s] fit done in %.1fs", b, proc.time()[3] - t0)
  }

  if (!length(out_layers)) return(NULL)

  # fill one ncell x nlyr matrix and set the raster's values once. assigning layer by layer
  # (r[[j]] <- v) copies the whole stack each time and OOM'd 6 subbasins (2026-09-16).
  # each layer is dropped from out_layers as it goes in, to keep the peak down
  template  <- cov_s[[1]]
  lyr_names <- names(out_layers)
  M <- matrix(NA_real_, nrow = terra::ncell(template), ncol = length(lyr_names))
  for (j in seq_along(lyr_names)) {
    M[backfill_idx, j] <- as.numeric(out_layers[[lyr_names[j]]])
    out_layers[[lyr_names[j]]] <- NULL
  }
  result_raster <- terra::rast(template, nlyr = length(lyr_names))
  terra::values(result_raster) <- M
  rm(M)
  names(result_raster) <- lyr_names

  out_dir <- file.path(ia_dir, "data", "derived_data", bart_dir, year, sprintf("subbasin_%s", subbasin_index))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  terra::writeRaster(result_raster, file.path(out_dir, sprintf("subbasin_%d_backfill.tif", subbasin_index)),
                     overwrite = TRUE, datatype = "FLT4S",
                     gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3", "ZLEVEL=6", "BIGTIFF=YES"))
  logp("writing %d layers to %s", length(lyr_names), out_dir)

  saveRDS(metrics,   file.path(out_dir, sprintf("subbasin_%d_metrics.rds", subbasin_index)))
  saveRDS(confusion, file.path(out_dir, sprintf("subbasin_%d_confusion.rds", subbasin_index)))

  list(subbasin = subbasin_index, ok = TRUE)
}
