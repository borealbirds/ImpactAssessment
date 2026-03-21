# define density prediction function ------------------------------------------------------
predict_species_bcr <- function(species, year, all_subbasins_subset, sector_name, hirsh_dir) {
  
  # communicate
  message(Sys.time(), " | preparing for species=", species, " year=", year)
  
  # find all models (BCRs) for the current `species`
  rdata_files <- list.files(file.path(nm_root,"output/06_bootstraps", species), pattern = "can.*\\.Rdata$", full.names = TRUE)
  message(Sys.time(), " | ", species, " | found ", length(rdata_files), " BCR models")
  
  # this sub-function will work through all model-lists (one per BCR) for a given species
  results <- lapply(rdata_files, function(rdata_path) {
    
    message(Sys.time(), " | loading ", basename(rdata_path))
    e <- new.env(parent = emptyenv())
    
    # creates a `b.list` object in the environment
    load(rdata_path, envir = e)
    
    if (!exists("b.list", envir = e)) {
      stop("b.list not found in ", basename(rdata_path))
    }
    
    # loads b.list (32 boots) for some spp x bcr
    b.list <- e$b.list
    
    # get current bcr
    bcr_code <- attr(b.list[[1]], "bcr")
    message(Sys.time(), " | working on species=", species, " BCR=", bcr_code)
    
    # find subbasins in current BCR
    sub_ids <- 
      bcr_subbasins_ref |>
      dplyr::filter(bcr_code == !!bcr_code) |>
      dplyr::pull(sub_index) |>
      unique()
    
    message(Sys.time(), " | ", species, " ", bcr_code, " | subbasins=", length(sub_ids))
    
    if (length(sub_ids) == 0) {
      message(Sys.time(), " | ", species, " ", bcr_code, " | NO SUBBASINS — skipping")
      return(NULL)
    }
    
    # get observed environmental stack for current BCR
    stack_obs <- terra::rast(file.path(nm_root,"gis/stacks", paste0(bcr_code, "_", year, ".tif")))
    message(Sys.time(), " | done loading observed stack")

    # build sector mask projected to BCR prediction grid
    # "all_hf": pre-built "goodsectors" raster (CanHF >= 1 AND any *target* sector > 0).
    # "goodsectors" raster: pixels that are only forestry_harvest/night_lights/population_density/nav_water
    # are excluded and keep their observed landscape bird estimates.
    # individual sectors: sector > 0 AND CanHF >= 1 (original full HF mask)
    if (sector_name == "all_hf") {
      goodsectors_r <- terra::project(
                         terra::rast(file.path(hirsh_dir, "CanHF_1km_morethan1_goodsectors.tif")),
                         stack_obs, method = "near")
      sector_mask <- terra::ifel(!is.na(goodsectors_r), 1, NA)
    } else {
      canHF_r     <- terra::project(
                       terra::rast(file.path(hirsh_dir, "CanHF_1km_morethan1.tif")),
                       stack_obs, method = "near")
      sector_r    <- terra::project(
                       terra::rast(file.path(hirsh_dir, paste0(sector_name, ".tif"))),
                       stack_obs, method = "near")
      sector_mask <- terra::ifel((sector_r > 0) & (canHF_r >= 1), 1, NA)
    }
    message(Sys.time(), " | sector mask built for ", sector_name,
            " (", global(sector_mask, "notNA")[[1]], " active pixels)")
    
    # get mosaiced watershed stacks for current BCR
    stack_bf <- mosaic_backfilled_stacks(sub_ids, year)
    if (!is.null(stack_bf)) {
      stack_bf <- terra::resample(stack_bf, stack_obs, method = "near")
      bam_bcr_codes <- gsub("_", "", paste(bam_boundary$country, bam_boundary$subUnit, sep = "_"))
      bcr_poly <- bam_boundary[bam_bcr_codes == bcr_code, ]
      stack_bf <- terra::mask(stack_bf, bcr_poly)
    } else {
      message(Sys.time(), " | no backfilled stack — skipping BCR")
      return(NULL)
    }
    
    message(Sys.time(), " | done loading backfilled stack")
    
    # identify continuous biotic covariates with stored posterior draw layers
    draw_layer_names <- names(stack_bf)[grep("_draw_[0-9]{3}$", names(stack_bf))]
    draw_covs        <- unique(sub("_draw_[0-9]{3}$", "", draw_layer_names))
    n_draws          <- 100L
    n_scen           <- 100L

    # pre-compute bootstrap-invariant quantities
    # GBM bootstrap models share the same var.names; compute once rather than
    # repeating identically for every bootstrap iteration.
    model_vars_shared  <- b.list[[1]]$var.names
    cat_vars_shared    <- intersect(model_vars_shared, categorical_responses)
    dist_shared        <- intersect(disturbance_vars$predictor, model_vars_shared)
    biotic_cont_shared <- intersect(setdiff(model_vars_shared, categorical_responses),
                                    biotic_continuous_vars)

    missing_bf <- biotic_cont_shared[
      !paste0(biotic_cont_shared, "_draw_001") %in% names(stack_bf)]
    if (length(missing_bf) > 0) {
      message(Sys.time(), " | ", species, " ", bcr_code,
              " | missing backfilled cont vars: ", paste(missing_bf, collapse = ", "))
    }

    # rasterize subbasin polygons onto prediction grid (needed before sector extraction)
    subbasin_zone_r <- terra::rasterize(
      all_subbasins_subset[sub_ids, ], stack_obs[[1]],
      field = "first_HYBAS_ID"
    )

    # sector pixel cell indices and subbasin zone membership (computed once)
    sector_cell_idx <- which(!is.na(terra::values(sector_mask, mat = FALSE)))
    sector_zones    <- terra::values(subbasin_zone_r, mat = FALSE)[sector_cell_idx]
    message(Sys.time(), " | ", species, " ", bcr_code,
            " | sector pixels: ", length(sector_cell_idx),
            " of ", terra::ncell(stack_obs[[1]]), " BCR cells")

    # pre-extract observed covariate values at sector pixels
    obs_all_vals <- terra::values(stack_obs)                             # n_cells x n_layers
    X_obs_sector <- as.data.frame(obs_all_vals[sector_cell_idx, , drop = FALSE])
    rm(obs_all_vals); gc()

    # pre-extract BART posterior draw values at sector pixels (expm1-transformed, clamped >= 0)
    draw_vals_sector <- setNames(
      lapply(draw_covs, function(v) {
        sapply(seq_len(n_draws), function(d) {
          lyr_vals <- terra::values(stack_bf[[paste0(v, "_draw_", sprintf("%03d", d))]], mat = FALSE)
          pmax(expm1(lyr_vals[sector_cell_idx]), 0)
        })
      }),
      draw_covs
    )

    # pre-extract categorical backfilled values at sector pixels (constant across scenarios)
    cat_vals_sector <- setNames(
      lapply(cat_vars_shared, function(v)
        terra::values(stack_bf[[v]], mat = FALSE)[sector_cell_idx]),
      cat_vars_shared
    )

    # build K sector-pixel scenario data frames (replaces full-BCR X_bf_sets rasters)
    # each df has the same columns as stack_obs, with biotic/disturbance vars overwritten
    set.seed(sum(utf8ToInt(paste0(species, bcr_code))) %% .Machine$integer.max)
    X_sector_sets <- lapply(seq_len(n_scen), function(k) {
      chosen <- sample(n_draws, 1)
      df <- X_obs_sector
      for (v in draw_covs)       { if (v %in% names(df)) df[[v]] <- draw_vals_sector[[v]][, chosen] }
      for (v in cat_vars_shared) { if (v %in% names(df)) df[[v]] <- cat_vals_sector[[v]] }
      for (v in dist_shared)     { if (v %in% names(df)) df[[v]] <- 0 }
      df
    })
    message(Sys.time(), " | ", species, " ", bcr_code,
            " | zeroed out disturbance vars: ",
            paste(intersect(names(X_obs_sector), dist_shared), collapse = ", "))

    # observed covariate stack for terra::predict on full BCR (observed landscape)
    X_obs <- stack_obs[[intersect(model_vars_shared, names(stack_obs))]]

    # create containers: obs_preds are full-BCR rasters; bf_preds are sector-pixel vectors
    obs_preds <- vector("list", length(b.list))
    bf_preds  <- vector("list", n_scen)
    for (k in seq_len(n_scen)) bf_preds[[k]] <- vector("list", length(b.list))

    # --- bootstrap loop -------------------------------------------------------
    for (i in seq_along(b.list)) {

      model <- b.list[[i]]

      # predict bird density on observed landscape (full BCR raster — needed for obs totals)
      obs_preds[[i]] <- terra::predict(X_obs, model, type = "response", n.trees = model$n.trees)
      message(Sys.time(), " | ", species, " ", bcr_code, " bootstrap=", i,
              " | finished predict() on observed landscape")

      # predict bird density at sector pixels only (avoids full-BCR predict for 100 scenarios)
      for (k in seq_len(n_scen)) {
        bf_preds[[k]][[i]] <- gbm::predict.gbm(
          model, X_sector_sets[[k]], n.trees = model$n.trees, type = "response")
      }
      message(Sys.time(), " | ", species, " ", bcr_code, " bootstrap=", i,
              " | finished predict() on backfilled landscape")

      terra::tmpFiles(remove = TRUE)

    } # close loop over `b.list` (a single spp x BCR)
    
    # save mean and SD of prediction surfaces ------------------------------------------------------
    
    # observed rasters are sector-independent; write once to the shared predictions/ path
    obs_dir <- file.path(ia_dir, "data", "derived_data", "predictions", species, bcr_code, year)
    # sector-specific backfilled rasters go to predictions_sectors/
    bf_dir  <- file.path(ia_dir, "data", "derived_data", "predictions_sectors",
                         sector_name, species, bcr_code, year)
    dir.create(obs_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(bf_dir,  recursive = TRUE, showWarnings = FALSE)

    # compute mean + sd incrementally without stacking (too memory intense)
    # helper: compute raster mean and sd incrementally (no stacking)
    summarize_preds <- function(pred_list) {

      n <- length(pred_list)

      mean_r <- pred_list[[1]]
      m2_r   <- mean_r * 0

      for (i in seq_len(n)) {
        x <- pred_list[[i]]
        delta <- x - mean_r
        mean_r <- mean_r + delta / i
        m2_r <- m2_r + delta * (x - mean_r)
      }

      list(mean = mean_r, sd   = sqrt(m2_r / (n - 1)))
    }

    # observed
    obs_summary <- summarize_preds(obs_preds)
    obs_mean <- obs_summary$mean
    obs_sd   <- obs_summary$sd

    # backfilled: bf_preds are sector-pixel vectors; compute mean/SD then reconstruct
    # sparse rasters (NA at non-sector pixels) for 15B footprint-scale computations
    {
      all_bf_vecs <- unlist(bf_preds, recursive = FALSE)  # list of K*n_boot vectors
      n_total     <- length(all_bf_vecs)
      bf_sum_vec  <- Reduce("+", all_bf_vecs)
      bf_sq_vec   <- Reduce("+", lapply(all_bf_vecs, function(v) v * v))
      bf_m_vec    <- bf_sum_vec / n_total
      bf_sd_vec   <- sqrt(pmax((bf_sq_vec - n_total * bf_m_vec^2) / (n_total - 1), 0))
      rm(all_bf_vecs, bf_sum_vec, bf_sq_vec); gc()

      bf_mean_r <- terra::rast(stack_obs[[1]]); terra::values(bf_mean_r) <- NA_real_
      bf_sd_r   <- terra::rast(stack_obs[[1]]); terra::values(bf_sd_r)   <- NA_real_
      bf_mean_r[sector_cell_idx] <- bf_m_vec
      bf_sd_r[sector_cell_idx]   <- bf_sd_vec
      rm(bf_m_vec, bf_sd_vec); gc()
    }

    # observed rasters are sector-independent; skip if already written by a previous sector run
    if (!file.exists(file.path(obs_dir, "observed_mean.tif"))) {
      terra::writeRaster(obs_mean, file.path(obs_dir, "observed_mean.tif"), overwrite = TRUE)
      terra::writeRaster(obs_sd,   file.path(obs_dir, "observed_sd.tif"),   overwrite = TRUE)
    }

    # write pooled backfilled rasters
    terra::writeRaster(bf_mean_r, file.path(bf_dir, "backfilled_mean_mean.tif"), overwrite = TRUE)
    terra::writeRaster(bf_sd_r,   file.path(bf_dir, "backfilled_mean_sd.tif"),   overwrite = TRUE)
    
    
    
    # count birds from prediction surfaces intersecting our watersheds ------------------------------------------------------

    message(Sys.time(), " | ", species, " ", bcr_code, " | estimating population over ", length(sub_ids), " subbasins")

    agg <- function(obs_rasters, bf_vecs_list, sub_ids, sector_mask, subbasin_zone_r, sector_zones) {

      n_sub  <- length(sub_ids)
      n_boot <- length(obs_rasters)
      n_scen <- length(bf_vecs_list)

      # subbasin-level arrays: dimensions = subbasin x bootstrap x scenario
      pop_obs_total_arr <- array(NA, dim = c(n_sub, n_boot, n_scen))
      pop_obs_on_bf_arr <- array(NA, dim = c(n_sub, n_boot, n_scen))
      pop_bf_only_arr   <- array(NA, dim = c(n_sub, n_boot, n_scen))

      # BCR-level arrays: observed quantities don’t depend on scenario,
      # backfilled quantities do (mean/low/high BART posterior scenarios)
      pop_obs_bcr_arr          <- numeric(n_boot)
      pop_obs_on_sector_bcr    <- numeric(n_boot)
      pop_bf_on_sector_bcr_mat <- matrix(NA_real_, n_boot, n_scen)

      # precompute zone alignment once: zonal() returns zones sorted by value
      zone_ids  <- sort(unique(terra::values(subbasin_zone_r, na.rm = TRUE)))
      hybas_ids <- all_subbasins_subset$first_HYBAS_ID[sub_ids]
      idx       <- match(hybas_ids, zone_ids)

      # pre-filter sector pixels that fall within a subbasin (for bf vector aggregation)
      valid_px         <- !is.na(sector_zones)
      sector_zones_flt <- sector_zones[valid_px]

      # outer loop: bootstrap — obs quantities are scenario-independent, compute here
      for (i in seq_len(n_boot)) {

        obs_r     <- obs_rasters[[i]] * 100
        obs_on_fp <- terra::mask(obs_r, sector_mask)

        pop_obs_bcr_arr[i]       <- as.numeric(terra::global(obs_r,     "sum", na.rm = TRUE))
        pop_obs_on_sector_bcr[i] <- as.numeric(terra::global(obs_on_fp, "sum", na.rm = TRUE))

        sub_obs_total <- terra::zonal(obs_r,     subbasin_zone_r, "sum", na.rm = TRUE)
        sub_obs_on_fp <- terra::zonal(obs_on_fp, subbasin_zone_r, "sum", na.rm = TRUE)

        # inner loop: scenario — only backfilled quantities vary here
        for (s in seq_len(n_scen)) {

          # bf_vecs_list[[s]][[i]] is a sector-pixel vector (not a full BCR raster)
          bf_vals <- bf_vecs_list[[s]][[i]] * 100

          pop_bf_on_sector_bcr_mat[i, s] <- sum(bf_vals, na.rm = TRUE)

          # aggregate bf values by subbasin zone using tapply
          bf_by_zone  <- tapply(bf_vals[valid_px], sector_zones_flt, sum, na.rm = TRUE)
          sub_bf_vals <- as.numeric(bf_by_zone[match(hybas_ids, names(bf_by_zone))])
          sub_bf_vals[is.na(sub_bf_vals)] <- 0  # subbasins with no sector pixels

          pop_obs_total_arr[, i, s] <- sub_obs_total[[2]][idx]
          pop_obs_on_bf_arr[, i, s] <- sub_obs_on_fp[[2]][idx]
          pop_bf_only_arr[, i, s]   <- sub_bf_vals

        } # close for s in seq

      } # close for i in seq

      # BCR-level CF and sector impact across all bootstrap x scenario combinations
      obs_bcr_mat           <- matrix(pop_obs_bcr_arr,        n_boot, n_scen)
      obs_on_sector_bcr_mat <- matrix(pop_obs_on_sector_bcr,  n_boot, n_scen)
      cf_total_bcr_mat      <- obs_bcr_mat - obs_on_sector_bcr_mat + pop_bf_on_sector_bcr_mat
      sector_impact_bcr_mat <- pop_bf_on_sector_bcr_mat - obs_on_sector_bcr_mat

      # subbasin-level summary: collapse boot x scenario into mean and SD per subbasin
      combine_stats <- function(x) {
        mat <- matrix(x, nrow = n_sub, ncol = n_boot * n_scen)
        list(mean = rowMeans(mat, na.rm = TRUE),
             sd   = apply(mat, 1, sd, na.rm = TRUE))
      }

      # BCR-level summary: scalar mean and SD directly from the boot x scenario distribution
      bcr_scalar <- function(mat) {
        list(mean = mean(mat, na.rm = TRUE), sd = sd(as.vector(mat), na.rm = TRUE))
      }

      list(
        # subbasin-level
        pop_obs_total     = combine_stats(pop_obs_total_arr),
        pop_obs_on_bf     = combine_stats(pop_obs_on_bf_arr),
        pop_bf_only       = combine_stats(pop_bf_only_arr),
        # BCR-level (obs SD = bootstrap only; cf/sector SD = bootstrap + BART posterior)
        obs_total_bcr     = list(mean = mean(pop_obs_bcr_arr),       sd = sd(pop_obs_bcr_arr)),
        obs_on_sector_bcr = list(mean = mean(pop_obs_on_sector_bcr), sd = sd(pop_obs_on_sector_bcr)),
        bf_on_sector_bcr  = bcr_scalar(pop_bf_on_sector_bcr_mat),
        cf_total_bcr      = bcr_scalar(cf_total_bcr_mat),
        sector_impact_bcr = bcr_scalar(sector_impact_bcr_mat)
      )

    } # close agg()
    
    
    # run aggregation including posterior
    pop_lists <- agg(obs_preds, bf_preds, sub_ids, sector_mask, subbasin_zone_r, sector_zones)
    
    # free up some memory
    rm(obs_preds, bf_preds, stack_obs, stack_bf)
    gc()

    # summarize into tidy tibble
    out <- tibble(
      species  = species,
      subbasin = sub_ids,
      bcr      = bcr_code,
      sector   = sector_name,
      # subbasin-level: population within each subbasin polygon
      obs_total_mean         = pop_lists$pop_obs_total$mean,
      obs_total_sd           = pop_lists$pop_obs_total$sd,
      obs_on_bf_mean         = pop_lists$pop_obs_on_bf$mean,
      obs_on_bf_sd           = pop_lists$pop_obs_on_bf$sd,
      bf_only_mean           = pop_lists$pop_bf_only$mean,
      bf_only_sd             = pop_lists$pop_bf_only$sd,
      # BCR-level: scalar summaries, same value repeated for every subbasin row
      # obs SD = BRT bootstrap only; cf/sector/bf SDs = bootstrap + BART posterior
      obs_total_bcr_mean     = pop_lists$obs_total_bcr$mean,
      obs_total_bcr_sd       = pop_lists$obs_total_bcr$sd,
      obs_on_sector_bcr_mean = pop_lists$obs_on_sector_bcr$mean,
      obs_on_sector_bcr_sd   = pop_lists$obs_on_sector_bcr$sd,
      bf_on_sector_bcr_mean  = pop_lists$bf_on_sector_bcr$mean,
      bf_on_sector_bcr_sd    = pop_lists$bf_on_sector_bcr$sd,
      cf_total_bcr_mean      = pop_lists$cf_total_bcr$mean,
      cf_total_bcr_sd        = pop_lists$cf_total_bcr$sd,
      sector_impact_bcr_mean = pop_lists$sector_impact_bcr$mean,
      sector_impact_bcr_sd   = pop_lists$sector_impact_bcr$sd
    )
    
    message(Sys.time(), " | ", species, " ", bcr_code, " | returning ", nrow(out), " rows")
    out
    
  })
  dplyr::bind_rows(results)
} # close predict_species_bcr()