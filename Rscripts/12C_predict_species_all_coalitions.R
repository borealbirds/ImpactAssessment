# ---
# title: Impact Assessment: predict bird density for one species across ALL 255 coalitions
# author: Mannfred Boehm
# ---
# Restructured replacement for predict_species_bcr() (12C_predict_species_bcr.R).
# See DESIGN_12C_restructure.md for the full rationale.
#
# KEY INSIGHT: the per-pixel backfilled prediction is coalition-INDEPENDENT. The
# joint-sampling seed (set.seed below) keys only on species+bcr+bootstrap(i)+
# scenario(k); the predicted design-matrix row at a pixel keys only on that pixel.
# The coalition merely selects WHICH pixels are summed. So we compute the bf field
# ONCE per species x BCR over the SUPERSET (all-8-sectors mask = any sector
# footprint AND CanHF>=1), then reduce each of the 255 coalitions as a masked
# grouped sum (rowsum) — no extra raster I/O, no extra gbm predicts.
#
# This was verified BIT-IDENTICAL to looping the old per-coalition predict_species_bcr()
# over all 255 coalitions (verify_restructure.R, on CAWA can14 and can10, 2026-06-10);
# both that oracle and the verifier have since been retired.
#
# Sourced by 12B_repredict_all_coalitions.R. Relies on the same globals that
# 12B sets up: nm_root, ia_dir, year, bcr_subbasins_ref, categorical_responses,
# disturbance_vars, biotic_continuous_vars, q.out, l.out, and the 12E helpers
# canonical_sectors() / coalition_id_to_sectors() / sectors_to_coalition_id().

predict_species_all_coalitions <- function(species, year, all_subbasins_subset,
                                            hirsh_dir, save_arrays_ids = integer(0)) {

  sectors    <- canonical_sectors()                 # 8 sectors, canonical order
  n_sectors  <- length(sectors)
  all_cids   <- 2:(2L^n_sectors)                    # 2..256 (cid 1 = empty, skipped)

  message(Sys.time(), " | preparing ALL coalitions for species=", species, " year=", year)

  qsp <- q.out[q.out$spp == species, ]$q
  q0  <- l.out[l.out$spp == species, ]$denshthresh

  rdata_files <- list.files(file.path(nm_root, "output/06_bootstraps", species),
                            pattern = "can.*\\.Rdata$", full.names = TRUE)
  message(Sys.time(), " | ", species, " | found ", length(rdata_files), " BCR models")

  # ---- per-BCR work: build the superset field once, reduce all coalitions ----
  # returns a list keyed by coalition_id; each element is that coalition's
  # per-BCR density-table tibble (n_sub rows).
  per_bcr <- lapply(rdata_files, function(rdata_path) {
    on.exit(gc())

    message(Sys.time(), " | loading ", basename(rdata_path))
    e <- new.env(parent = emptyenv())
    load(rdata_path, envir = e)
    if (!exists("b.list", envir = e)) stop("b.list not found in ", basename(rdata_path))
    b.list <- e$b.list; rm(e)

    test_n_boot <- as.integer(Sys.getenv("TEST_N_BOOT", "0"))
    if (test_n_boot > 0L) {
      b.list <- b.list[seq_len(min(test_n_boot, length(b.list)))]
      message(Sys.time(), " | TEST_N_BOOT=", test_n_boot,
              " — truncating to ", length(b.list), " bootstrap(s)")
    }

    bcr_code <- attr(b.list[[1]], "bcr")
    message(Sys.time(), " | working on species=", species, " BCR=", bcr_code)

    test_bcr_env <- Sys.getenv("TEST_BCR", "")
    if (nchar(test_bcr_env) > 0 && !bcr_code %in% strsplit(test_bcr_env, ",")[[1]]) {
      message(Sys.time(), " | TEST_BCR=", test_bcr_env, " — skipping ", bcr_code)
      return(NULL)
    }

    sub_ids <- bcr_subbasins_ref |>
      dplyr::filter(bcr_code == !!bcr_code) |>
      dplyr::pull(sub_index) |>
      unique()
    message(Sys.time(), " | ", species, " ", bcr_code, " | subbasins=", length(sub_ids))
    if (length(sub_ids) == 0) {
      message(Sys.time(), " | ", species, " ", bcr_code, " | NO SUBBASINS — skipping")
      return(NULL)
    }

    stack_obs <- terra::rast(file.path(nm_root, "gis/stacks", paste0(bcr_code, "_", year, ".tif")))
    message(Sys.time(), " | done loading observed stack")

    # ---- SUPERSET mask: any of the 8 sectors > 0 AND CanHF >= 1 (= cid 256) ----
    canHF_r <- terra::project(
      terra::rast(file.path(hirsh_dir, "CanHF_1km_morethan1.tif")),
      stack_obs, method = "near")

    # per-sector reprojected footprint rasters, kept so coalition membership over
    # the superset pixels can be derived cheaply (no per-coalition reprojection).
    sec_rasters <- setNames(vector("list", n_sectors), sectors)
    union_mask  <- terra::rast(stack_obs[[1]]); terra::values(union_mask) <- 0L
    for (sec in sectors) {
      sec_r <- terra::project(
        terra::rast(file.path(hirsh_dir, paste0(sec, ".tif"))),
        stack_obs, method = "near")
      sec_rasters[[sec]] <- sec_r
      union_mask <- terra::ifel(sec_r > 0, 1L, union_mask)
    }
    super_mask <- terra::ifel((union_mask == 1L) & (canHF_r >= 1), 1, NA)
    rm(union_mask); gc()

    n_super_px <- terra::global(super_mask, "notNA")[[1]]
    message(Sys.time(), " | superset mask built (", n_super_px, " pixels)")
    if (n_super_px == 0) {
      message(Sys.time(), " | ", species, " ", bcr_code, " | superset empty — skipping BCR")
      return(NULL)
    }

    bf_mosaic_path <- file.path(ia_dir, "data", "derived_data", "bart_models_mosaics",
                                year, paste0(bcr_code, "_backfilled.tif"))
    if (!file.exists(bf_mosaic_path)) {
      message(Sys.time(), " | no backfilled mosaic found for ", bcr_code, " — skipping BCR")
      return(NULL)
    }
    stack_bf <- terra::rast(bf_mosaic_path)
    message(Sys.time(), " | done loading backfilled stack")

    draw_layer_names <- names(stack_bf)[grep("_draw_[0-9]{3}$", names(stack_bf))]
    draw_covs        <- unique(sub("_draw_[0-9]{3}$", "", draw_layer_names))
    n_draws          <- 100L
    n_scen           <- 100L

    # bootstrap-invariant quantities (identical to 12C) -------------------------
    model_vars_shared  <- b.list[[1]]$var.names
    cat_vars_shared    <- intersect(model_vars_shared, categorical_responses)
    dist_shared        <- intersect(disturbance_vars$predictor, model_vars_shared)
    biotic_cont_shared <- intersect(setdiff(model_vars_shared, categorical_responses),
                                    biotic_continuous_vars)

    cat_levels_shared <- setNames(
      lapply(cat_vars_shared, function(v)
        as.character(b.list[[1]]$var.levels[[match(v, model_vars_shared)]])),
      cat_vars_shared)
    null_cat_vars <- names(Filter(function(x) is.null(x) || length(x) == 0L, cat_levels_shared))
    if (length(null_cat_vars) > 0L)
      message(Sys.time(), " | WARNING: ", species, " ", bcr_code,
              " | var.levels empty for categorical var(s) ",
              paste(null_cat_vars, collapse = ", "), " — skipping factor conversion")

    # rasterize subbasins; zone field = first_HYBAS_ID (zone ids ARE hybas ids) --
    subbasin_zone_r <- terra::rasterize(all_subbasins_subset[sub_ids, ], stack_obs[[1]],
                                        field = "first_HYBAS_ID")
    zone_ids  <- sort(unique(terra::values(subbasin_zone_r, na.rm = TRUE)))
    hybas_ids <- all_subbasins_subset$first_HYBAS_ID[sub_ids]
    idx       <- match(hybas_ids, zone_ids)         # sub -> zone slot; NA if absent

    # superset pixel indices, their subbasin zone, and per-sector membership -----
    super_idx   <- which(!is.na(terra::values(super_mask, mat = FALSE)))
    super_zones <- terra::values(subbasin_zone_r, mat = FALSE)[super_idx]
    sec_member  <- setNames(lapply(sectors, function(sec)
      terra::values(sec_rasters[[sec]], mat = FALSE)[super_idx] > 0), sectors)
    sec_member  <- lapply(sec_member, function(v) { v[is.na(v)] <- FALSE; v })
    rm(sec_rasters, canHF_r); gc()
    message(Sys.time(), " | ", species, " ", bcr_code,
            " | superset pixels: ", length(super_idx),
            " of ", terra::ncell(stack_obs[[1]]), " BCR cells")

    # ---- Prediction weight (range x not-water x in-data-limit; built by 12A2) ---
    # Applied multiplicatively to BOTH observed and backfilled density so obs/bf
    # symmetry holds: w*bf - w*obs = w*(bf - obs). Replicates V5 10.Truncate
    # range/water/extent masking. Missing weight => UNMASKED (w == 1) with a warning
    # so the pipeline still runs. Built on the stack grid, so it aligns with
    # stack_obs and the observed bootstraps; resample only as a grid-drift guard.
    weight_path <- file.path(ia_dir, "data", "derived_data", "predictions",
                             species, bcr_code, year, "weight.tif")
    if (file.exists(weight_path)) {
      weight_r <- terra::rast(weight_path)
      if (!isTRUE(terra::compareGeom(weight_r, stack_obs[[1]], stopOnError = FALSE)))
        weight_r <- terra::resample(weight_r, stack_obs[[1]], method = "near")
    } else {
      message(Sys.time(), " | ", species, " ", bcr_code,
              " | no weight.tif — proceeding UNMASKED (run 12A2_build_prediction_weights)")
      weight_r <- terra::rast(stack_obs[[1]]); terra::values(weight_r) <- 1
    }
    weight_super <- terra::values(weight_r, mat = FALSE)[super_idx]

    # observed covariates at superset pixels (identical extraction to 12C) -------
    obs_all_vals  <- terra::values(stack_obs)
    X_obs_super   <- as.data.frame(obs_all_vals[super_idx, , drop = FALSE], check.names = FALSE)
    rm(obs_all_vals); gc()

    # BART posterior draws at superset pixels (expm1, non-finite -> NA, pmax 0) ---
    draw_vals_super <- setNames(
      lapply(draw_covs, function(v) {
        lyr_names <- paste0(v, "_draw_", sprintf("%03d", seq_len(n_draws)))
        lyr_names <- intersect(lyr_names, names(stack_bf))
        mat_raw   <- terra::values(stack_bf[[lyr_names]], mat = TRUE)
        mat_sec   <- mat_raw[super_idx, , drop = FALSE]; rm(mat_raw)
        mat_sec   <- expm1(mat_sec)
        mat_sec[!is.finite(mat_sec)] <- NA_real_
        pmax(mat_sec, 0)
      }), draw_covs)

    # categorical backfill at superset pixels (constant across scenarios) --------
    cat_vals_super <- setNames(
      lapply(cat_vars_shared, function(v) {
        lyr_name <- if (v %in% names(stack_bf)) v else paste0(v, "_mean")
        if (!lyr_name %in% names(stack_bf)) {
          message(Sys.time(), " | WARNING: no backfilled layer for categorical var ", v,
                  " — filling NA")
          return(rep(NA_integer_, length(super_idx)))
        }
        terra::values(stack_bf[[lyr_name]], mat = FALSE)[super_idx]
      }), cat_vars_shared)

    # canonical observed bootstraps (from 12A; always present now that the colleague's
    # observed prediction rasters are staged on the cluster — absence is a hard error,
    # not a bf-only fall-through).
    obs_dir       <- file.path(ia_dir, "data", "derived_data", "predictions", species, bcr_code, year)
    obs_boot_path <- file.path(obs_dir, "observed_bootstraps.tif")
    if (!file.exists(obs_boot_path))
      stop(species, " ", bcr_code, " | observed_bootstraps.tif missing at ", obs_boot_path,
           " — run 12A_observed.R and Globus-transfer the outputs before 12B")
    obs_boot_stack <- terra::rast(obs_boot_path)
    n_obs <- terra::nlyr(obs_boot_stack)
    obs_preds <- lapply(seq_len(n_obs), function(i) obs_boot_stack[[i]]); rm(obs_boot_stack)
    if (test_n_boot > 0L) obs_preds <- obs_preds[seq_len(min(test_n_boot, length(obs_preds)))]
    # apply prediction weight to observed density (range/water/extent masking).
    # weights BOTH O_super (obs_on source) and obs_total_byboot (zonal) downstream,
    # mirroring OLD 12C where a single weighted obs_preds fed both quantities.
    obs_preds <- lapply(obs_preds, function(r) r * weight_r)
    message(Sys.time(), " | ", species, " ", bcr_code, " | loaded ",
            length(obs_preds), " canonical observed bootstraps (weighted)")
    n_boot <- length(obs_preds)
    if (length(b.list) != n_boot)
      stop(species, " ", bcr_code, " | bootstrap count mismatch: b.list=",
           length(b.list), " obs=", n_boot)

    # ---- complete-case mask over the SUPERSET design matrix (coalition-free) ---
    X_rep <- X_obs_super
    for (v in cat_vars_shared) if (v %in% names(X_rep)) X_rep[[v]] <- cat_vals_super[[v]]
    for (v in dist_shared)     if (v %in% names(X_rep)) X_rep[[v]] <- 0
    for (v in cat_vars_shared) if (v %in% names(X_rep)) {
      lvls <- cat_levels_shared[[v]]
      if (!is.null(lvls) && length(lvls) > 0L)
        X_rep[[v]] <- factor(as.character(X_rep[[v]]), levels = lvls)
    }
    for (v in draw_covs) if (v %in% names(X_rep)) X_rep[[v]] <- draw_vals_super[[v]][, 1L]
    complete_mask <- stats::complete.cases(X_rep[, model_vars_shared, drop = FALSE])
    rm(X_rep)
    message(Sys.time(), " | ", species, " ", bcr_code, " | complete superset pixels: ",
            sum(complete_mask), " / ", length(super_idx),
            " (", round(100 * mean(complete_mask), 1), "%)")

    # ---- Joint BRT x BART sampling over the SUPERSET (the expensive part, ONCE)
    # Worker i returns an [n_super x n_scen] matrix of capped birds/ha predictions
    # (NA at incomplete-case pixels). Seed is coalition-free — identical to 12C.
    n_cores <- max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")))
    message(Sys.time(), " | ", species, " ", bcr_code, " | sampling ", n_boot,
            " bootstraps x ", n_scen, " scenarios on ", n_cores, " core(s)")

    boot_mats <- parallel::mclapply(seq_along(b.list), function(i) {
      model <- b.list[[i]]
      X_k <- X_obs_super
      for (v in cat_vars_shared) if (v %in% names(X_k)) X_k[[v]] <- cat_vals_super[[v]]
      for (v in dist_shared)     if (v %in% names(X_k)) X_k[[v]] <- 0
      for (v in cat_vars_shared) if (v %in% names(X_k)) {
        lvls <- cat_levels_shared[[v]]
        if (is.null(lvls) || length(lvls) == 0L) next
        X_k[[v]] <- factor(as.character(X_k[[v]]), levels = lvls)
      }
      sc <- matrix(NA_real_, nrow = nrow(X_k), ncol = n_scen)
      for (k in seq_len(n_scen)) {
        set.seed((sum(utf8ToInt(paste0(species, bcr_code))) + i * 1000L + k) %%
                   .Machine$integer.max)
        chosen <- sample(n_draws, 1)
        for (v in draw_covs) if (v %in% names(X_k)) X_k[[v]] <- draw_vals_super[[v]][, chosen]
        pred_vec <- rep(NA_real_, nrow(X_k))
        if (any(complete_mask))
          pred_vec[complete_mask] <- gbm::predict.gbm(
            model, X_k[complete_mask, , drop = FALSE],
            n.trees = model$n.trees, type = "response")
        sc[, k] <- pmin(pred_vec, qsp)
      }
      sc
    }, mc.cores = n_cores)

    failed <- vapply(boot_mats, inherits, logical(1L), "try-error")
    if (any(failed)) stop(sprintf("%s %s | %d/%d bootstrap workers failed",
                                  species, bcr_code, sum(failed), length(boot_mats)))

    # M: [n_super x (n_boot*n_scen)] capped birds/ha; col block per bootstrap.
    M <- do.call(cbind, boot_mats); rm(boot_mats)
    # apply prediction weight AFTER the in-worker qsp clamp (V5 order: truncate then
    # range-multiply). weight_super has length n_super = nrow(M); column-major
    # recycling multiplies row i of every column by weight_super[i]. Symmetric with
    # the observed-side weighting above, so bf - obs stays w*(bf - obs). The bf-only
    # arr stashed below is reshaped from this M, so it inherits the weighting too.
    M <- M * weight_super
    # O_super: observed birds/ha at superset pixels, [n_super x n_boot]
    O_super <- vapply(obs_preds, function(r) terra::values(r, mat = FALSE)[super_idx],
                      numeric(length(super_idx)))

    # ---- coalition-INDEPENDENT obs_total (whole-subbasin), computed ONCE -------
    # mirrors 12C: per bootstrap, zonal sum of obs_r*100 over subbasin zones.
    obs_total_byboot <- vapply(seq_len(n_boot), function(i) {
      z <- terra::zonal(obs_preds[[i]] * 100, subbasin_zone_r, "sum", na.rm = TRUE)
      z[[2]][idx]                                     # NA where sub's zone absent
    }, numeric(length(sub_ids)))
    if (is.null(dim(obs_total_byboot)))
      obs_total_byboot <- matrix(obs_total_byboot, nrow = length(sub_ids))

    rm(stack_obs, stack_bf, obs_preds, X_obs_super, draw_vals_super, cat_vals_super,
       b.list, cat_levels_shared); gc()

    # combine_stats: identical math to 12C (column order is irrelevant to mean/sd)
    combine_stats <- function(mat) {
      list(mean = rowMeans(mat, na.rm = TRUE),
           sd   = apply(mat, 1L, sd, na.rm = TRUE))
    }
    n_sub <- length(sub_ids)
    # broadcast a [n_sub x n_boot] matrix to [n_sub x (n_boot*n_scen)]
    bcols <- rep(seq_len(n_boot), each = n_scen)

    # NOTE: per-coalition inspection rasters (predictions_coalitions/.../backfilled_
    # {mean,sd}.tif) are intentionally NOT written here — no downstream script reads
    # them (12D/12F/14B consume only the .rds tables). If needed, emit them only for
    # the save_arrays_ids coalitions from M (birds/ha) + the q99/q0 caps, as in 12C.

    # ---- reduce every coalition (cheap: masked grouped sums) -------------------
    coalition_tables <- vector("list", length(all_cids))
    names(coalition_tables) <- as.character(all_cids)

    for (ci in seq_along(all_cids)) {
      cid       <- all_cids[ci]
      coalition <- coalition_id_to_sectors(cid, sectors)

      # coalition pixel membership over superset = OR of its sectors' footprints
      coal_mem <- Reduce(`|`, sec_member[coalition])

      # Reproduce 12C's n_active==0 short-circuit EXACTLY: a coalition with no
      # footprint pixels in this BCR returns an all-zeros table (incl. obs_total=0,
      # which is a coalition-invariant quantity 12C nonetheless zeroes here). Match
      # it for bit-identity; sum(coal_mem) == 12C's n_active (superset already ANDs
      # CanHF>=1, so coal_mem over super == coalition-mask notNA count).
      if (sum(coal_mem) == 0L) {
        # no footprint pixels for this coalition in this BCR: obs placeholders are 0
        # (bit-identical to OLD 12C's n_active==0 short-circuit), and bf is a genuine
        # 0 since no footprint here means no counterfactual change.
        coalition_tables[[ci]] <- tibble::tibble(
          species = species, subbasin = sub_ids, bcr = bcr_code, coalition_id = cid,
          obs_total_mean = 0, obs_total_sd = 0,
          obs_on_coalition_mean = 0, obs_on_coalition_sd = 0,
          bf_on_coalition_mean = 0, bf_on_coalition_sd = 0)
        next
      }

      keep   <- coal_mem & complete_mask & !is.na(super_zones)
      n_keep <- sum(keep)

      # bf_on_coalition: rowsum of M over kept pixels by subbasin zone, ->[n_sub x 3200]
      bf_mat <- matrix(0, nrow = n_sub, ncol = n_boot * n_scen)
      obs_on_mat <- matrix(NA_real_, n_sub, n_boot * n_scen)
      if (n_keep > 0) {
        zk     <- super_zones[keep]
        # x100: birds/ha -> birds/km2 (matches 12C `bf_vals * 100`)
        rs_bf  <- rowsum(M[keep, , drop = FALSE], group = zk, reorder = TRUE) * 100
        map_bf <- match(hybas_ids, rownames(rs_bf))           # bf: absent -> 0
        present <- !is.na(map_bf)
        if (any(present)) bf_mat[present, ] <- rs_bf[map_bf[present], , drop = FALSE]

        Ok <- O_super[keep, , drop = FALSE] * 100
        Ok[is.na(Ok)] <- 0                                  # na.rm-equivalent sum
        rs_obs <- rowsum(Ok, group = zk, reorder = TRUE)
        # obs_on is filled ONLY for subbasins that have kept (coalition AND
        # complete-case) pixels. Subbasins with no kept pixels — whether their
        # zone is absent from the raster (idx NA) or present-but-footprint-free —
        # stay NA, so combine_stats yields NaN mean / NA sd. This mirrors 12C
        # exactly: its terra::zonal over an all-NA zone returns NaN. Do NOT
        # pre-fill 0 here; 12C never assigns a real 0 to an empty obs_on cell.
        obs_on_byboot <- matrix(NA_real_, n_sub, n_boot)
        map_o <- match(hybas_ids, rownames(rs_obs))
        presO <- !is.na(map_o)
        if (any(presO)) obs_on_byboot[presO, ] <- rs_obs[map_o[presO], , drop = FALSE]
        obs_on_mat <- obs_on_byboot[, bcols, drop = FALSE]
      } else {
        # coalition has footprint in this BCR but no complete-case pixels at all:
        # no kept pixels anywhere -> every subbasin's obs_on is empty -> NA
        # (NaN mean / NA sd), matching 12C's all-NA zonal result.
        obs_on_mat <- matrix(NA_real_, n_sub, n_boot)[, bcols, drop = FALSE]
      }

      bf_stats      <- combine_stats(bf_mat)
      obs_total_mat <- obs_total_byboot[, bcols, drop = FALSE]
      ot_stats      <- combine_stats(obs_total_mat)
      ooc_stats     <- combine_stats(obs_on_mat)
      coalition_tables[[ci]] <- tibble::tibble(
        species = species, subbasin = sub_ids, bcr = bcr_code, coalition_id = cid,
        obs_total_mean = ot_stats$mean,  obs_total_sd = ot_stats$sd,
        obs_on_coalition_mean = ooc_stats$mean, obs_on_coalition_sd = ooc_stats$sd,
        bf_on_coalition_mean = bf_stats$mean,   bf_on_coalition_sd = bf_stats$sd)
    }

    list(coalition_tables = coalition_tables,
         bcr_code         = bcr_code)
  })

  per_bcr <- Filter(Negate(is.null), per_bcr)

  # ---- bind each coalition's per-BCR tables across BCRs -> 255 final tables ----
  cid_keys <- as.character(all_cids)
  tables_by_cid <- setNames(lapply(cid_keys, function(k) {
    parts <- lapply(per_bcr, function(x) x$coalition_tables[[k]])
    dplyr::bind_rows(Filter(Negate(is.null), parts))
  }), cid_keys)

  list(tables_by_cid = tables_by_cid)             # named list: "2".."256" -> tibble
}
