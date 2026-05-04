# ---
# title: Impact Assessment: predict bird density for one species x BCR under a given sector coalition
# author: Mannfred Boehm
# ---
# Sourced by 12B_repredict_birds.R. Defines predict_species_bcr(), which performs the core
# counterfactual prediction for one species across all of its BCRs under a specified coalition.
#
# For each BCR the function:
#   1. Reads the canonical observed bootstraps produced by 12A_observed.R
#   2. Builds a coalition mask: pixels where any sector in the coalition has footprint
#      AND CanHF >= 1 receive backfilled covariates; all other pixels use observed covariates
#   3. Runs joint BRT x BART sampling: for each BRT bootstrap iteration a fresh BART posterior
#      draw is made, so BART uncertainty is naturally nested inside BRT uncertainty
#   4. Returns subbasin-level density tables (bootstrap x BART-scenario arrays) that 12B
#      aggregates up to BCR-wide population totals
# ---

# define density prediction function ------------------------------------------------------
predict_species_bcr <- function(species, year, all_subbasins_subset, coalition, coalition_id, hirsh_dir, save_arrays = FALSE) {

  # communicate
  coalition_label <- if (length(coalition) == 0) "empty" else paste(coalition, collapse = "+")
  message(Sys.time(), " | preparing for species=", species, " year=", year,
          " coalition=", coalition_label)

  # find all models (BCRs) for the current `species`
  rdata_files <- list.files(file.path(nm_root,"output/06_bootstraps", species), pattern = "can.*\\.Rdata$", full.names = TRUE)
  message(Sys.time(), " | ", species, " | found ", length(rdata_files), " BCR models")

  # this sub-function will work through all model-lists (one per BCR) for a given species
  results <- lapply(rdata_files, function(rdata_path) {
    on.exit(gc())

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

    # ---- Build coalition mask (union of individual sector masks) ----
    # Empty coalition: no pixels are masked (trivial case, v = 0)
    if (length(coalition) == 0) {
      sector_mask <- terra::rast(stack_obs[[1]]); terra::values(sector_mask) <- NA
    } else {
      canHF_r <- terra::project(
        terra::rast(file.path(hirsh_dir, "CanHF_1km_morethan1.tif")),
        stack_obs, method = "near")

      # build union mask: any pixel where any coalition sector > 0 AND CanHF >= 1
      union_mask <- terra::rast(stack_obs[[1]]); terra::values(union_mask) <- 0L
      for (sec in coalition) {
        sec_r <- terra::project(
          terra::rast(file.path(hirsh_dir, paste0(sec, ".tif"))),
          stack_obs, method = "near")
        union_mask <- terra::ifel(sec_r > 0, 1L, union_mask)
      }
      sector_mask <- terra::ifel((union_mask == 1L) & (canHF_r >= 1), 1, NA)
      rm(canHF_r, union_mask); gc()
    }

    n_active <- terra::global(sector_mask, "notNA")[[1]]
    message(Sys.time(), " | coalition mask built for ", coalition_label,
            " (", n_active, " active pixels)")

    # trivial case: empty coalition or no active pixels
    if (n_active == 0) {
      message(Sys.time(), " | ", species, " ", bcr_code,
              " | no active pixels — returning zeros")
      return(list(
        table = tibble(
          species  = species,
          subbasin = sub_ids,
          bcr      = bcr_code,
          coalition_id = coalition_id,
          obs_total_mean         = 0, obs_total_sd         = 0,
          obs_on_coalition_mean  = 0, obs_on_coalition_sd  = 0,
          bf_on_coalition_mean   = 0, bf_on_coalition_sd   = 0
        ),
        n_dropped   = 0L,
        n_sector_px = 0L,
        arrays = NULL
      ))
    }

    # read pre-mosaiced BCR-wide backfilled stack (written by 11_premosaic_backfilled_stacks.R)
    bf_mosaic_path <- file.path(ia_dir, "data", "derived_data", "bart_models_mosaics",
                                year, paste0(bcr_code, "_backfilled.tif"))
    if (!file.exists(bf_mosaic_path)) {
      message(Sys.time(), " | no backfilled mosaic found for ", bcr_code, " — skipping BCR")
      return(NULL)
    }
    stack_bf <- terra::rast(bf_mosaic_path)
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
            " | coalition pixels: ", length(sector_cell_idx),
            " of ", terra::ncell(stack_obs[[1]]), " BCR cells")

    # pre-extract observed covariate values at coalition pixels
    obs_all_vals <- terra::values(stack_obs)                             # n_cells x n_layers
    X_obs_sector <- as.data.frame(obs_all_vals[sector_cell_idx, , drop = FALSE], check.names = FALSE)
    rm(obs_all_vals); gc()

    # pre-extract BART posterior draw values at coalition pixels (expm1-transformed, clamped >= 0)
    # Non-finite values (Inf from expm1 overflow, NA from edge pixels) are replaced with 0
    # before they can reach gbm's C prediction code, which cannot handle non-finite inputs.
    draw_vals_sector <- setNames(
      lapply(draw_covs, function(v) {
        lyr_names <- paste0(v, "_draw_", sprintf("%03d", seq_len(n_draws)))
        lyr_names <- intersect(lyr_names, names(stack_bf))
        mat_raw   <- terra::values(stack_bf[[lyr_names]], mat = TRUE)
        mat_sec   <- mat_raw[sector_cell_idx, , drop = FALSE]
        rm(mat_raw)
        mat_sec   <- expm1(mat_sec)
        mat_sec[!is.finite(mat_sec)] <- 0
        pmax(mat_sec, 0)
      }),
      draw_covs
    )

    # pre-extract categorical backfilled values at coalition pixels (constant across scenarios)
    cat_vals_sector <- setNames(
      lapply(cat_vars_shared, function(v)
        terra::values(stack_bf[[v]], mat = FALSE)[sector_cell_idx]),
      cat_vars_shared
    )

    # --- Read canonical observed bootstraps (from 12_observed.R) ---
    obs_dir <- file.path(ia_dir, "data", "derived_data", "predictions", species, bcr_code, year)
    obs_boot_path <- file.path(obs_dir, "observed_bootstraps.tif")

    if (!file.exists(obs_boot_path)) {
      stop("Canonical observed bootstraps not found at ", obs_boot_path,
           "\nRun 12_observed.R first (sbatch 12_observed.sh)")
    }

    obs_boot_stack <- terra::rast(obs_boot_path)
    n_boot <- terra::nlyr(obs_boot_stack)
    obs_preds <- lapply(seq_len(n_boot), function(i) obs_boot_stack[[i]])
    rm(obs_boot_stack)
    message(Sys.time(), " | ", species, " ", bcr_code,
            " | loaded ", n_boot, " canonical observed bootstraps")

    # --- Joint sampling: nest BART posterior scenario inside BRT bootstrap ---
    # Each (bootstrap, scenario) pair draws a fresh BART posterior realization
    # and predicts bird density, capturing the covariance between BART and BRT
    # uncertainty naturally.

    # containers: bf_preds[[boot]][[scen]] = vector of sector-pixel predictions
    bf_preds <- vector("list", n_boot)
    for (i in seq_len(n_boot)) bf_preds[[i]] <- vector("list", n_scen)

    n_dropped_bcr <- NA_integer_  # captured on first (i=1, k=1) iteration

    for (i in seq_along(b.list)) {
      model <- b.list[[i]]

      for (k in seq_len(n_scen)) {
        # deterministic seed per (species, bcr, boot, scen) for reproducibility
        set.seed((sum(utf8ToInt(paste0(species, bcr_code))) + i * 1000L + k) %%
                   .Machine$integer.max)
        chosen <- sample(n_draws, 1)

        # build scenario data frame: start from observed, overwrite biotic/disturbance
        X_k <- X_obs_sector
        for (v in draw_covs)       { if (v %in% names(X_k)) X_k[[v]] <- draw_vals_sector[[v]][, chosen] }
        for (v in cat_vars_shared) { if (v %in% names(X_k)) X_k[[v]] <- cat_vals_sector[[v]] }
        for (v in dist_shared)     { if (v %in% names(X_k)) X_k[[v]] <- 0 }

        # Exclude NA rows before prediction: gbm's C code segfaults on NAs.
        # Consistent with terra::predict() (propagates NAs) and
        # LandbirdModelsV5/analysis/12.Summarize.R (sum(na.rm=TRUE) excludes
        # NA pixels from population totals).
        complete_rows <- complete.cases(X_k[, model$var.names, drop = FALSE])

        # capture drop count on the first iteration; NA pattern is stable across
        # (i, k) because NAs come from fixed sources (X_obs_sector non-overwritten
        # columns; BART draws and categorical values are already clamped/non-NA).
        if (i == 1L && k == 1L) {
          n_dropped_bcr <- sum(!complete_rows)
        }

        pred_vec <- rep(NA_real_, nrow(X_k))
        if (any(complete_rows)) {
          pred_vec[complete_rows] <- gbm::predict.gbm(
            model, X_k[complete_rows, , drop = FALSE],
            n.trees = model$n.trees, type = "response")
        }
        bf_preds[[i]][[k]] <- pred_vec
      }

      message(Sys.time(), " | ", species, " ", bcr_code, " bootstrap=", i,
              " | finished predict() on backfilled landscape")
      gc()
    }

    message(sprintf(
      "%s | %s %s | incomplete-case pixels dropped: %d / %d (%.1f%%)",
      Sys.time(), species, bcr_code,
      n_dropped_bcr, length(sector_cell_idx),
      100 * n_dropped_bcr / max(length(sector_cell_idx), 1L)
    ))

    # ---- Save backfilled prediction rasters (mean/SD) for optional inspection ----

    bf_dir <- file.path(ia_dir, "data", "derived_data", "predictions_coalitions",
                        as.character(coalition_id), species, bcr_code, year)
    dir.create(bf_dir, recursive = TRUE, showWarnings = FALSE)

    {
      all_bf_vecs <- unlist(bf_preds, recursive = FALSE)
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

      terra::writeRaster(bf_mean_r, file.path(bf_dir, "backfilled_mean.tif"), overwrite = TRUE)
      terra::writeRaster(bf_sd_r,   file.path(bf_dir, "backfilled_sd.tif"),   overwrite = TRUE)
      rm(bf_m_vec, bf_sd_vec, bf_mean_r, bf_sd_r); gc()
    }

    # ---- Aggregate populations: subbasin-level (bottom-up) ----

    message(Sys.time(), " | ", species, " ", bcr_code,
            " | estimating population over ", length(sub_ids), " subbasins")

    agg <- function(obs_rasters, bf_vecs_list, sub_ids, sector_mask,
                    subbasin_zone_r, sector_zones, save_arrays = FALSE) {

      n_sub  <- length(sub_ids)
      n_boot <- length(obs_rasters)
      n_scen <- length(bf_vecs_list[[1]])

      # subbasin-level arrays: dimensions = subbasin x bootstrap x scenario
      pop_obs_total_arr     <- array(NA, dim = c(n_sub, n_boot, n_scen))
      pop_obs_on_coal_arr   <- array(NA, dim = c(n_sub, n_boot, n_scen))
      pop_bf_on_coal_arr    <- array(NA, dim = c(n_sub, n_boot, n_scen))

      # precompute zone alignment: zonal() returns zones sorted by value
      zone_ids  <- sort(unique(terra::values(subbasin_zone_r, na.rm = TRUE)))
      hybas_ids <- all_subbasins_subset$first_HYBAS_ID[sub_ids]
      idx       <- match(hybas_ids, zone_ids)

      # pre-filter sector pixels that fall within a subbasin
      valid_px         <- !is.na(sector_zones)
      sector_zones_flt <- sector_zones[valid_px]

      # Cap per-pixel BRT predictions at the maximum of the packaged mean raster.
      # Follows LandbirdModelsV5/analysis/12.Summarize.R, which uses the max of
      # 10.Package.R's cleaned mean layer as a species- and BCR-specific upper bound.
      # Falls back to Inf (no-op clamp) if the packaged raster is unavailable.
      pkg_path <- file.path(nm_root, "output", "10_packaged", species, bcr_code,
                            paste0(species, "_", bcr_code, "_", year, ".tif"))
      q99 <- if (file.exists(pkg_path)) {
        terra::global(terra::rast(pkg_path)[["mean"]], max, na.rm = TRUE)[, 1]
      } else {
        message("  WARNING: packaged raster not found for ", species, " ", bcr_code,
                " ", year, " — skipping prediction cap")
        Inf
      }

      for (i in seq_len(n_boot)) {

        obs_r     <- terra::clamp(obs_rasters[[i]], upper = q99) * 100
        obs_on_fp <- terra::mask(obs_r, sector_mask)

        sub_obs_total <- terra::zonal(obs_r,     subbasin_zone_r, "sum", na.rm = TRUE)
        sub_obs_on_fp <- terra::zonal(obs_on_fp, subbasin_zone_r, "sum", na.rm = TRUE)

        for (s in seq_len(n_scen)) {

          bf_vals <- bf_vecs_list[[i]][[s]] * 100

          # aggregate bf values by subbasin zone using tapply
          bf_by_zone  <- tapply(bf_vals[valid_px], sector_zones_flt, sum, na.rm = TRUE)
          sub_bf_vals <- as.numeric(bf_by_zone[match(hybas_ids, names(bf_by_zone))])
          sub_bf_vals[is.na(sub_bf_vals)] <- 0

          pop_obs_total_arr[, i, s]   <- sub_obs_total[[2]][idx]
          pop_obs_on_coal_arr[, i, s] <- sub_obs_on_fp[[2]][idx]
          pop_bf_on_coal_arr[, i, s]  <- sub_bf_vals

        }
      }

      # collapse boot x scenario into mean and SD per subbasin
      combine_stats <- function(x) {
        mat <- matrix(x, nrow = n_sub, ncol = n_boot * n_scen)
        list(mean = rowMeans(mat, na.rm = TRUE),
             sd   = apply(mat, 1, sd, na.rm = TRUE))
      }

      out <- list(
        pop_obs_total   = combine_stats(pop_obs_total_arr),
        pop_obs_on_coal = combine_stats(pop_obs_on_coal_arr),
        pop_bf_on_coal  = combine_stats(pop_bf_on_coal_arr)
      )

      if (save_arrays) {
        out$bcr_obs_total_mat   <- apply(pop_obs_total_arr,   c(2L, 3L), sum, na.rm = TRUE)
        out$bcr_obs_on_coal_mat <- apply(pop_obs_on_coal_arr, c(2L, 3L), sum, na.rm = TRUE)
        out$bcr_bf_on_coal_mat  <- apply(pop_bf_on_coal_arr,  c(2L, 3L), sum, na.rm = TRUE)
      }

      out
    }

    pop_lists <- agg(obs_preds, bf_preds, sub_ids, sector_mask,
                     subbasin_zone_r, sector_zones, save_arrays = save_arrays)

    rm(obs_preds, bf_preds, stack_obs, stack_bf)
    gc()

    # return subbasin-level tibble (no BCR scalars — 15B aggregates upward)
    out <- tibble(
      species  = species,
      subbasin = sub_ids,
      bcr      = bcr_code,
      coalition_id          = coalition_id,
      obs_total_mean        = pop_lists$pop_obs_total$mean,
      obs_total_sd          = pop_lists$pop_obs_total$sd,
      obs_on_coalition_mean = pop_lists$pop_obs_on_coal$mean,
      obs_on_coalition_sd   = pop_lists$pop_obs_on_coal$sd,
      bf_on_coalition_mean  = pop_lists$pop_bf_on_coal$mean,
      bf_on_coalition_sd    = pop_lists$pop_bf_on_coal$sd
    )

    message(Sys.time(), " | ", species, " ", bcr_code, " | returning ", nrow(out), " rows")
    list(
      table       = out,
      n_dropped   = n_dropped_bcr,
      n_sector_px = length(sector_cell_idx),
      arrays = if (save_arrays) list(
        obs_total_mat   = pop_lists$bcr_obs_total_mat,
        obs_on_coal_mat = pop_lists$bcr_obs_on_coal_mat,
        bf_on_coal_mat  = pop_lists$bcr_bf_on_coal_mat
      ) else NULL
    )

  })

  results_nonnull <- Filter(Negate(is.null), results)

  # species-wide dropped pixel summary (printed at end of each species job)
  total_dropped   <- sum(vapply(results_nonnull, function(r) if (is.null(r$n_dropped))   0L else r$n_dropped,   integer(1L)))
  total_sector_px <- sum(vapply(results_nonnull, function(r) if (is.null(r$n_sector_px)) 0L else r$n_sector_px, integer(1L)))
  message(sprintf(
    "%s | DROPPED PIXEL SUMMARY | species=%s coalition_id=%d | BCRs=%d | total_coalition_px=%d | dropped=%d (%.1f%%)",
    Sys.time(), species, coalition_id, length(results_nonnull),
    total_sector_px, total_dropped,
    100 * total_dropped / max(total_sector_px, 1L)
  ))

  national_arrays <- if (save_arrays) {
    mats <- Filter(Negate(is.null), lapply(results_nonnull, `[[`, "arrays"))
    if (length(mats) > 0) {
      list(
        obs_total_mat   = Reduce("+", lapply(mats, `[[`, "obs_total_mat")),
        obs_on_coal_mat = Reduce("+", lapply(mats, `[[`, "obs_on_coal_mat")),
        bf_on_coal_mat  = Reduce("+", lapply(mats, `[[`, "bf_on_coal_mat"))
      )
    } else NULL
  } else NULL

  list(
    table           = dplyr::bind_rows(lapply(results_nonnull, `[[`, "table")),
    national_arrays = national_arrays
  )
} # close predict_species_bcr()
