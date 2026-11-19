# ---
# title: Impact Assessment: predict bird density for one species x BCR across all 255 coalitions
# author: Mannfred Boehm
# ---

# sourced by 12D (and 12H, for combine_bcr_results() and coalition_task_table()). uses 12D's
# globals: nm_root, ia_dir, q.out, bcr_subbasins_ref, categorical_responses, footprint_vars,
# biotic_continuous_vars, and 12E and 12G's functions.
#
# a pixel's backfilled prediction doesn't depend on the coalition: its joint-sampling seed
# keys on species, BCR, bootstrap i and scenario k only, and its design row on the pixel only.
# so the backfilled field is predicted once, over the all-sector footprint (any sector > 0 and
# CanHF >= 1, = coalition 256), and each coalition is a masked sum of it (CLAUDE.md, "12F
# invariants").


# categorical columns as factors in THIS model's levels. bootstraps are fitted to different
# resamples, so their var.levels differ (CAWA can14: 28 of 32 know VLCE classes bootstrap 1
# never saw); bootstrap 1's levels sent those classes down the missing branch in the others
as_model_factors <- function(X, model, vars) {
  for (v in vars) if (v %in% names(X)) {
    lvls <- model$var.levels[[match(v, model$var.names)]]
    if (is.null(lvls) || length(lvls) == 0L) next
    X[[v]] <- factor(as.character(X[[v]]), levels = as.character(lvls))
  }
  X
}


# one species x BCR (rdata_path: its 06_bootstraps file). returns the coalition tables, the
# national-array pieces, the per-sample subbasin Shapley values and the BCR code; NULL when the
# BCR has nothing to predict (no subbasins, no footprint, or no backfill mosaic)
predict_species_all_coalitions <- function(species, year, rdata_path, all_subbasins_subset,
                                           hirsh_dir, save_arrays_ids = integer(0)) {

  sectors   <- canonical_sectors()
  n_sectors <- length(sectors)
  all_cids  <- 2:(2L^n_sectors)            # coalition 1 (no sector) has v = 0
  n_draws   <- 100L                        # BART draws per backfilled covariate (08B)
  n_scen    <- 100L                        # scenarios per bootstrap
  bcr_code  <- sub(paste0("^", species, "_(.*)\\.Rdata$"), "\\1", basename(rdata_path))
  test_n_boot <- as.integer(Sys.getenv("TEST_N_BOOT", "0"))
  message(Sys.time(), " | ", species, " ", bcr_code, " ", year)

  # ---- V5's two caps (12B), and the prediction weight (12C) ----
  # densmax: a species constant. q99: the 99.9th percentile of the observed landscape, frozen by
  # 12A: a counterfactual's own q99 would sit higher (no industry, more birds) and put part of
  # the cap into the obs - cf contrast. q.out has had no $q column since V5's 4e7fc83, and
  # pmin(x, NULL) is numeric(0) without a warning, hence the schema check
  if (!"densmax" %in% names(q.out))
    stop("q.out has no densmax column: restage SpeciesPredictionTruncationValues.Rdata (2026-06-04 or later)")
  qsp <- q.out[q.out$spp == species, ]$densmax
  if (length(qsp) != 1L || !is.finite(qsp)) stop("no usable densmax for ", species)

  params_path <- file.path(ia_dir, "data", "derived_data", "predictions", species, "truncation_params.rds")
  tp <- if (file.exists(params_path)) readRDS(params_path)[[paste(bcr_code, year, sep = "_")]]
  if (is.null(tp) || !is.finite(tp$q99)) stop(species, " ", bcr_code, ": no frozen q99 in ", params_path, " - run 12A")
  if (!isTRUE(all.equal(tp$densmax, qsp)))
    stop(species, " ", bcr_code, ": densmax is ", tp$densmax, " in ", params_path, " but ", qsp, " in q.out")
  q99 <- tp$q99
  message(Sys.time(), " | caps: densmax=", signif(qsp, 6), " q99=", signif(q99, 6))

  # the weight (range x not water x data limit x in BCR) is the ONLY masking: since 08A imputes
  # partial-NA predictors, water and out-of-range pixels no longer drop out by themselves, and
  # an unweighted run over-counted CAWA can10 7.8x. its band name is 12C's version stamp
  pred_dir    <- file.path(ia_dir, "data", "derived_data", "predictions", species, bcr_code, year)
  weight_path <- file.path(pred_dir, "weight.tif")
  if (!file.exists(weight_path)) stop(species, " ", bcr_code, ": no ", weight_path, " - run 12C")
  weight_r <- terra::rast(weight_path)
  if (!isTRUE(names(weight_r)[1] == "weight_v3_touches"))
    stop(species, " ", bcr_code, ": weight.tif is version '", names(weight_r)[1], "', not weight_v3_touches - rerun 12C")
  if (!isTRUE(terra::global(weight_r, "max", na.rm = TRUE)[[1]] > 0))
    stop(species, " ", bcr_code, ": weight.tif is 0 or NA over the whole BCR - rebuild it with 12C")

  obs_boot_path <- file.path(pred_dir, "observed_bootstraps.tif")
  if (!file.exists(obs_boot_path)) stop(species, " ", bcr_code, ": no ", obs_boot_path, " - run 12A")

  # ---- the bird models ----
  e <- new.env(parent = emptyenv())
  load(rdata_path, envir = e)
  b.list <- e$b.list
  rm(e)
  if (test_n_boot > 0L) b.list <- b.list[seq_len(min(test_n_boot, length(b.list)))]
  # 12A labels its outputs by the file name, so the two must agree
  if (!identical(attr(b.list[[1]], "bcr"), bcr_code))
    stop(basename(rdata_path), " holds a model for ", attr(b.list[[1]], "bcr"), ", not ", bcr_code)

  sub_ids <- unique(bcr_subbasins_ref$sub_index[bcr_subbasins_ref$bcr_code == bcr_code])
  n_sub   <- length(sub_ids)
  message(Sys.time(), " | subbasins=", n_sub)
  if (n_sub == 0) return(NULL)

  stack_obs <- terra::rast(file.path(nm_root, "gis/stacks", paste0(bcr_code, "_", year, ".tif")))
  if (!isTRUE(terra::compareGeom(weight_r, stack_obs[[1]], stopOnError = FALSE)))
    weight_r <- terra::resample(weight_r, stack_obs[[1]], method = "near")

  # ---- the all-sector footprint (the superset), and each sector's footprint and direct part ----
  on_grid <- function(f) terra::project(terra::rast(file.path(hirsh_dir, f)), stack_obs, method = "near")
  canHF_r    <- on_grid("CanHF_1km_morethan1.tif")
  sec_r      <- setNames(lapply(sectors, function(s) on_grid(paste0(s, ".tif"))), sectors)
  sec_dir_r  <- setNames(lapply(sectors, function(s) on_grid(paste0(s, "_direct.tif"))), sectors)
  union_mask <- terra::rast(stack_obs[[1]])
  terra::values(union_mask) <- 0L
  for (s in sectors) union_mask <- terra::ifel(sec_r[[s]] > 0, 1L, union_mask)
  super_mask <- terra::ifel((union_mask == 1L) & (canHF_r >= 1), 1, NA)
  super_idx  <- which(!is.na(terra::values(super_mask, mat = FALSE)))
  message(Sys.time(), " | superset: ", length(super_idx), " pixels")
  if (length(super_idx) == 0) return(NULL)

  bf_mosaic_path <- file.path(ia_dir, "data", "derived_data", "bart_models_mosaics", year, paste0(bcr_code, "_backfilled.tif"))
  if (!file.exists(bf_mosaic_path)) {
    message(Sys.time(), " | no backfill mosaic for ", bcr_code, " - skipping the BCR")
    return(NULL)
  }
  stack_bf <- terra::rast(bf_mosaic_path)

  # each superset pixel's subbasin, and whether each sector's footprint (> 0) and direct part (== 1) hold it
  subbasin_zone_r <- terra::rasterize(all_subbasins_subset[sub_ids, ], stack_obs[[1]], field = "first_HYBAS_ID")
  hybas_ids   <- all_subbasins_subset$first_HYBAS_ID[sub_ids]
  if (anyDuplicated(hybas_ids)) stop(species, " ", bcr_code, ": duplicate first_HYBAS_ID")
  super_zones <- terra::values(subbasin_zone_r, mat = FALSE)[super_idx]
  sec_member  <- lapply(sec_r, function(r) {
    v <- terra::values(r, mat = FALSE)[super_idx] > 0
    !is.na(v) & v
  })
  sec_direct  <- setNames(lapply(sectors, function(s) {
    v <- terra::values(sec_dir_r[[s]], mat = FALSE)[super_idx] == 1
    !is.na(v) & v & sec_member[[s]]
  }), sectors)
  rm(union_mask, super_mask, canHF_r, sec_r, sec_dir_r)
  gc()
  message(Sys.time(), " | direct footprint of any sector: ", sum(Reduce(`|`, sec_direct)), " pixels")

  # a weight 0 over the whole footprint is a range that misses it (BAWW, BBWA can41): impact 0
  weight_super <- terra::values(weight_r, mat = FALSE)[super_idx]
  if (sum(weight_super, na.rm = TRUE) == 0) message(Sys.time(), " | weight is 0 on the whole footprint: impact 0")

  # ---- covariates the bootstraps share ----
  model_vars <- b.list[[1]]$var.names
  if (!all(vapply(b.list, function(m) identical(m$var.names, model_vars), logical(1L))))
    stop(species, " ", bcr_code, ": the bootstraps differ in var.names")
  cat_vars    <- intersect(model_vars, categorical_responses)
  dist_vars   <- intersect(footprint_vars, model_vars)
  biotic_cont <- intersect(setdiff(model_vars, categorical_responses), biotic_continuous_vars)

  # ---- observed covariates at the superset pixels ----
  miss_vars <- setdiff(model_vars, names(stack_obs))
  if (length(miss_vars) > 0L) stop(species, " ", bcr_code, ": not in the observed stack: ", paste(miss_vars, collapse = ", "))
  obs_all_vals <- terra::values(stack_obs[[model_vars]])[super_idx, , drop = FALSE]
  # terra returns missing cells as NaN, but gbm sends only NA_real_ down a split's missing
  # branch (a NaN goes right). V5's observed predictions took the missing branch, so a NaN
  # here predicted the bf side by another rule (CAWA can11: ~+10% of the footprint's birds)
  obs_all_vals[is.nan(obs_all_vals)] <- NA_real_
  X_obs_super <- as.data.frame(obs_all_vals, check.names = FALSE)
  rm(obs_all_vals)
  gc()

  # ---- backfilled covariates at the superset pixels ----
  # continuous: the BART draws (08B; log1p scale, so expm1; non-finite -> NA, never 0), only
  # for covariates the models read. where 08A found b constant on a subbasin's training pixels
  # it wrote <b>_mean (raw scale) instead of draws: that constant fills every draw there (CAWA
  # can60: SCANFIBalsamFir_5x5 in 24 subbasins, 43% of weight > 0)
  draw_covs <- intersect(unique(sub("_draw_[0-9]{3}$", "", grep("_draw_[0-9]{3}$", names(stack_bf), value = TRUE))),
                         model_vars)
  draw_vals_super <- setNames(lapply(draw_covs, function(v) {
    lyr_names <- intersect(paste0(v, "_draw_", sprintf("%03d", seq_len(n_draws))), names(stack_bf))
    m <- expm1(terra::values(stack_bf[[lyr_names]], mat = TRUE)[super_idx, , drop = FALSE])
    m[!is.finite(m)] <- NA_real_
    mean_lyr <- paste0(v, "_mean")
    if (mean_lyr %in% names(stack_bf)) {
      const <- terra::values(stack_bf[[mean_lyr]], mat = FALSE)[super_idx]
      fill  <- rowSums(!is.na(m)) == 0L & is.finite(const)
      if (any(fill)) {
        m[fill, ] <- const[fill]
        message(Sys.time(), " | ", v, ": ", sum(fill), " superset pixels filled from ", mean_lyr)
      }
    }
    pmax(m, 0)
  }), draw_covs)

  # a covariate constant in EVERY subbasin of the BCR has only <v>_mean: carried as one draw
  mean_only_covs <- setdiff(intersect(sub("_mean$", "", grep("_mean$", names(stack_bf), value = TRUE)), biotic_cont),
                            draw_covs)
  for (v in mean_only_covs) {
    const <- terra::values(stack_bf[[paste0(v, "_mean")]], mat = FALSE)[super_idx]
    const[!is.finite(const)] <- NA_real_
    draw_vals_super[[v]] <- matrix(pmax(const, 0), ncol = 1L)
    message(Sys.time(), " | ", v, ": no draws in the BCR, using ", v, "_mean")
  }
  draw_covs <- c(draw_covs, mean_only_covs)

  # modelled blanks (08B_deploy_pbart.R): <v>_blank == 1 reaches gbm as NA_real_, its missing
  # branch, as on such land in V5's observed predictions, and counts as complete below
  blank_super <- setNames(lapply(union(draw_covs, cat_vars), function(v) {
    bl <- paste0(v, "_blank")
    if (!bl %in% names(stack_bf)) return(logical(length(super_idx)))
    terra::values(stack_bf[[bl]], mat = FALSE)[super_idx] %in% 1
  }), union(draw_covs, cat_vars))
  for (v in draw_covs) draw_vals_super[[v]][blank_super[[v]], ] <- NA_real_

  # categorical: the drawn classes <v>_class_001..100 (08B), never the modal class, which erases
  # minority classes. a pixel with no drawn class takes the modal one, as <v>_mean fills draws
  cat_vals_super <- setNames(lapply(cat_vars, function(v) {
    cls      <- paste0(v, "_class_", sprintf("%03d", seq_len(n_draws)))
    mode_lyr <- if (v %in% names(stack_bf)) v else paste0(v, "_mean")
    has_mode <- mode_lyr %in% names(stack_bf)
    mode_v   <- if (has_mode) terra::values(stack_bf[[mode_lyr]], mat = FALSE)[super_idx]
    if (all(cls %in% names(stack_bf))) {
      m <- terra::values(stack_bf[[cls]], mat = TRUE)[super_idx, , drop = FALSE]
      if (has_mode) {
        fill <- rowSums(!is.na(m)) == 0L & is.finite(mode_v)
        if (any(fill)) {
          m[fill, ] <- mode_v[fill]
          message(Sys.time(), " | ", v, ": ", sum(fill), " superset pixels filled from the modal class")
        }
      }
    } else if (has_mode) {
      message(Sys.time(), " | ", v, ": no drawn classes in the mosaic, using the modal class for every draw")
      m <- matrix(mode_v, ncol = 1L)
    } else {
      message(Sys.time(), " | WARNING: no backfilled layer for ", v, " - filling NA")
      return(matrix(NA_real_, length(super_idx), 1L))
    }
    m[blank_super[[v]], ] <- NA_real_
    m
  }), cat_vars)

  # ---- the observed side: 12A's capped V5 predictions, weighted ----
  obs_boot_stack <- terra::rast(obs_boot_path)
  obs_preds <- lapply(seq_len(terra::nlyr(obs_boot_stack)), function(i) obs_boot_stack[[i]] * weight_r)
  rm(obs_boot_stack)
  if (test_n_boot > 0L) obs_preds <- obs_preds[seq_len(min(test_n_boot, length(obs_preds)))]
  n_boot <- length(obs_preds)
  if (length(b.list) != n_boot) stop(species, " ", bcr_code, ": ", length(b.list), " models but ", n_boot, " observed layers")

  # ---- complete-case mask: ONE mask for the obs and bf sides ----
  # it gates only the backfilled covariates (draw 1 standing for all draws): an observed-only
  # covariate's NA goes down gbm's missing branch, as in V5's predictions, and a categorical is
  # tested raw (a class some bootstrap never saw takes that model's missing branch, as in V5)
  X_rep <- X_obs_super
  for (v in cat_vars)  X_rep[[v]] <- cat_vals_super[[v]][, 1L]
  for (v in dist_vars) X_rep[[v]] <- 0
  for (v in intersect(draw_covs, names(X_rep))) X_rep[[v]] <- draw_vals_super[[v]][, 1L]
  backfilled_vars <- intersect(c(draw_covs, cat_vars), model_vars)
  complete_mask   <- Reduce(`&`, lapply(backfilled_vars, function(v) !is.na(X_rep[[v]]) | blank_super[[v]]),
                            rep(TRUE, nrow(X_rep)))

  # coverage that matters is among weight > 0 pixels (the grid's 100 km buffer has weight 0)
  in_w <- is.finite(weight_super) & weight_super > 0
  message(Sys.time(), " | complete weight>0 superset pixels: ", sum(complete_mask & in_w), " / ", sum(in_w),
          " (", round(100 * sum(complete_mask & in_w) / max(1L, sum(in_w)), 1), "%)")
  lost <- in_w & !complete_mask
  if (any(lost)) {
    na_by_var <- vapply(backfilled_vars, function(v) sum(is.na(X_rep[[v]][lost])), integer(1))
    na_by_var <- sort(na_by_var[na_by_var > 0], decreasing = TRUE)
    message(Sys.time(), " | NA among lost weight>0 pixels, by covariate: ", paste0(names(na_by_var), "=", na_by_var, collapse = ", "))
  }
  rm(X_rep, in_w, lost)

  # a weight-0 pixel is 0 whatever gbm predicts, so only weight != 0 pixels are predicted (a
  # non-finite weight is still predicted, so it NAs the product and trips the audit below)
  pred_mask <- complete_mask & !(is.finite(weight_super) & weight_super == 0)
  pred_rows <- which(pred_mask)
  message(Sys.time(), " | predicting ", length(pred_rows), " complete pixels")

  # ---- identity gate: the workers' path must reproduce V5's own prediction ----
  # the OBSERVED design through the workers' path, against observed_bootstraps.tif: any gap is
  # a difference in HOW the two sides are predicted, which would land in bf - obs as a fake
  # impact (how the NaN/NA mismatch and the per-bootstrap levels were found). pixels with a
  # missing covariate or a class some bootstrap never saw are over-sampled, as that is where
  # such rules part. tolerance: the tif's float32 storage
  spread  <- function(p, k) p[unique(round(seq(1, length(p), length.out = min(k, length(p)))))]
  na_rows <- pred_rows[!stats::complete.cases(X_obs_super[pred_rows, , drop = FALSE])]
  odd_cls <- Reduce(`|`, lapply(cat_vars, function(v) {
    common <- Reduce(intersect, lapply(b.list, function(m) as.character(m$var.levels[[match(v, m$var.names)]])))
    x <- X_obs_super[[v]][pred_rows]
    !is.na(x) & !(as.character(x) %in% common)
  }), rep(FALSE, length(pred_rows)))
  id_rows <- sort(unique(c(spread(pred_rows, 1500L), spread(na_rows, 1500L), spread(pred_rows[odd_cls], 1500L))))
  if (length(id_rows) > 0L) {
    X_id    <- X_obs_super[id_rows, , drop = FALSE]
    obs_chk <- terra::rast(obs_boot_path)
    for (i in seq_len(n_boot)) {
      X_i  <- as_model_factors(X_id, b.list[[i]], cat_vars)
      p_id <- pmin(pmin(gbm_predict_design(b.list[[i]], gbm_design(b.list[[i]], X_i)), qsp), q99)
      o_id <- terra::values(obs_chk[[i]], mat = FALSE)[super_idx[id_rows]]
      bad  <- is.finite(o_id) & !(abs(p_id - o_id) <= 1e-5 * abs(o_id) + 1e-30)
      if (any(bad))
        stop(species, " ", bcr_code, ": identity gate failed for bootstrap ", i, ": ", sum(bad), " of ",
             length(id_rows), " pixels predicted from the observed covariates differ from observed_bootstraps.tif",
             " (max rel diff ", signif(max(abs(p_id - o_id)[bad] / abs(o_id[bad])), 3), ")")
    }
    message(Sys.time(), " | identity gate: all ", n_boot, " bootstraps reproduce V5 at ", length(id_rows),
            " pixels (", length(na_rows), " with a missing covariate, ", sum(odd_cls), " with a class some bootstrap never saw)")
    rm(X_id, X_i, obs_chk, p_id, o_id, bad)
  }

  # ---- coalition membership, once ----
  # the workers carry only the complete rows (c_rows), and keep_rows[[ci]] indexes into them.
  # two masks per coalition (OL #10): keep_rows is its whole footprint (any of its sectors > 0),
  # where the counterfactual has no footprint covariates (d0); direct_rows, a subset, is where
  # one of its sectors converted the land, the only place vegetation is backfilled (bf). so the
  # counterfactual is bf on direct_rows + d0 on the rest of keep_rows, and the split is exact:
  # direct part = d0 - obs on keep_rows, vegetation part = bf - d0 on direct_rows.
  # coal_state: 0 = no footprint in the BCR, 1 = footprint but no complete pixel, 2 = kept pixels
  c_rows      <- which(complete_mask)
  pred_in_c   <- pred_mask[c_rows]
  w_c         <- weight_super[c_rows]
  zones_c     <- super_zones[c_rows]
  zi_c        <- match(zones_c, hybas_ids)         # subbasin row of each complete row
  coal_state  <- integer(length(all_cids))
  keep_rows   <- vector("list", length(all_cids))
  direct_rows <- vector("list", length(all_cids))
  for (ci in seq_along(all_cids)) {
    coal_secs <- coalition_id_to_sectors(all_cids[ci], sectors)
    coal_mem  <- Reduce(`|`, sec_member[coal_secs])
    if (!any(coal_mem)) next
    kr <- which((coal_mem & !is.na(super_zones))[c_rows])
    coal_state[ci] <- if (length(kr) > 0L) 2L else 1L
    if (length(kr) > 0L) {
      keep_rows[[ci]]   <- kr
      direct_rows[[ci]] <- which((Reduce(`|`, sec_direct[coal_secs]) & !is.na(super_zones))[c_rows])
    }
  }

  # from here only the predicted rows are read; the BART draws dominate memory (CAWA can11:
  # 8.9 GB, 35% of it never predicted), so keep just those rows
  # rows are kept one matrix at a time, so the full matrices and their subsets never coexist
  X_pred <- X_obs_super[pred_rows, , drop = FALSE]
  rm(X_obs_super)
  cat_vals_pred <- lapply(cat_vals_super, function(v) v[pred_rows, , drop = FALSE])
  rm(cat_vals_super)
  draw_vals_pred <- draw_vals_super
  rm(draw_vals_super)
  for (v in names(draw_vals_pred)) draw_vals_pred[[v]] <- draw_vals_pred[[v]][pred_rows, , drop = FALSE]
  gc()

  # ---- joint BRT x BART sampling (the expensive part), one forked worker per bootstrap ----
  # worker i predicts its bootstrap's capped, weighted birds/ha at the complete rows for each
  # scenario, and reduces it to every coalition's [n_sub x n_scen] subbasin sums (coal_rowsum,
  # bit-identical to rowsum over the full matrix, which is never built)
  n_cores <- max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")))
  message(Sys.time(), " | sampling ", n_boot, " bootstraps x ", n_scen, " scenarios on ", n_cores, " core(s)")
  multi_draw <- any(vapply(c(draw_vals_pred[draw_covs], cat_vals_pred), ncol, integer(1L)) > 1L)
  ci_full    <- match(2L^n_sectors, all_cids)      # the all-sector coalition (256)

  boot_res <- parallel::mclapply(seq_along(b.list), function(i) {
    model <- b.list[[i]]
    n_c   <- length(c_rows)

    # each scenario re-seeds and draws one BART draw with replacement (only ~64 of 100 are
    # distinct), so each distinct draw is predicted once and copied to the scenarios that drew it
    chosen_k <- vapply(seq_len(n_scen), function(k) {
      set.seed((sum(utf8ToInt(paste0(species, bcr_code))) + i * 1000L + k) %% .Machine$integer.max)
      sample(n_draws, 1)
    }, integer(1L))
    key <- if (multi_draw) chosen_k else rep(1L, n_scen)

    # d0: this bootstrap with the observed vegetation and the footprint covariates at 0. d0 - obs
    # is the bird model's direct response to the footprint; it needs no BART draw
    d0_c <- numeric(n_c)                           # complete rows with weight 0 stay 0
    if (length(pred_rows) > 0L) {
      X_d <- X_pred
      for (v in dist_vars) X_d[[v]] <- 0
      d0_c[pred_in_c] <- gbm_predict_design(model, gbm_design(model, as_model_factors(X_d, model, cat_vars)))
      rm(X_d)
    }
    d0_c   <- pmin(pmin(d0_c, qsp), q99) * w_c
    n_bad0 <- sum(is.na(d0_c))
    d0 <- lapply(seq_along(all_cids), function(ci) {
      if (coal_state[ci] != 2L) return(NULL)
      coal_rowsum(matrix(d0_c, ncol = 1L), keep_rows[[ci]], zi_c, n_sub) * 100
    })
    # d0 on the footprint that is NOT direct, which keeps its observed vegetation
    d0_infl <- lapply(seq_along(all_cids), function(ci) {
      if (coal_state[ci] != 2L) return(NULL)
      d0[[ci]] - coal_rowsum(matrix(d0_c, ncol = 1L), direct_rows[[ci]], zi_c, n_sub) * 100
    })
    # the all-sector coalition's per-pixel counterfactual for the maps: d0, and on its direct
    # pixels the mean over this bootstrap's scenarios (filled in below)
    cf_full <- d0_c
    rm(d0_c)
    invisible(gc(full = FALSE))

    # the design matrix, built once per bootstrap; per draw only the draw columns change
    if (length(pred_rows) > 0L) {
      X_k <- X_pred
      for (v in cat_vars) {
        cv <- cat_vals_pred[[v]]
        X_k[[v]] <- cv[, if (ncol(cv) == 1L) 1L else chosen_k[1L]]
      }
      for (v in dist_vars) X_k[[v]] <- 0
      X_k  <- as_model_factors(X_k, model, cat_vars)
      dcol <- match(draw_covs, model$var.names)    # NA: not read by this model
      # drawn classes change per draw too: their 0-based codes in THIS model's levels, as
      # as_model_factors() + gbm_design() build them (an unseen class -> NA, the missing branch)
      ccol <- match(cat_vars, model$var.names)
      ccol[vapply(cat_vars, function(v) ncol(cat_vals_pred[[v]]) == 1L, logical(1L))] <- NA
      cat_code <- function(j, k) {
        lv <- model$var.levels[[ccol[j]]]
        x  <- cat_vals_pred[[cat_vars[j]]][, k]
        if (is.null(lv) || length(lv) == 0L) return(as.double(x))
        as.numeric(factor(as.character(x), levels = as.character(lv))) - 1
      }
      for (v in intersect(draw_covs, names(X_k))) {
        dv <- draw_vals_pred[[v]]
        X_k[[v]] <- dv[, if (ncol(dv) == 1L) 1L else chosen_k[1L]]
      }
      D <- gbm_design(model, X_k)
      gbm_check_fast(model, X_k, D, label = paste(species, bcr_code, "bootstrap", i))
      for (j in which(!is.na(ccol)))
        if (!identical(cat_code(j, chosen_k[1L]), D[, ccol[j]]))
          stop(species, " ", bcr_code, " bootstrap ", i, ": drawn-class codes for ", cat_vars[j], " differ from gbm_design()'s")
      rm(X_k)
    }

    # a forked worker inherits the parent's GC trigger (~10-15 GB after the draw reads on big
    # BCRs), so it could pile up that much garbage before collecting: 8 of them OOM'd at 64G.
    # hence only per-draw vectors here, and gc(full = FALSE), which scans only the worker's own objects
    sc    <- matrix(NA_real_, nrow = n_c, ncol = n_scen)
    n_bad <- 0
    for (d in unique(key)) {
      ks     <- which(key == d)
      pred_c <- numeric(n_c)                       # complete rows with weight 0 stay 0
      if (length(pred_rows) > 0L) {
        chosen <- chosen_k[ks[1L]]
        for (j in which(!is.na(dcol))) {
          dv <- draw_vals_pred[[draw_covs[j]]]
          D[, dcol[j]] <- dv[, if (ncol(dv) == 1L) 1L else chosen]
        }
        for (j in which(!is.na(ccol))) D[, ccol[j]] <- cat_code(j, chosen)
        pred_c[pred_in_c] <- gbm_predict_design(model, D)
      }
      # V5's order: densmax, then q99, then the weight; the same as the obs side, so bf - obs = w * (bf - obs)
      sc_d <- pmin(pmin(pred_c, qsp), q99) * w_c
      # coal_rowsum has no na.rm, so one NA would NA a subbasin's whole bf sum while the obs
      # side zeroes its NAs. the complete mask used draw 1, so another draw could still be NA
      n_bad <- n_bad + length(ks) * sum(is.na(sc_d))
      sc[, ks] <- sc_d
      rm(pred_c, sc_d)
      invisible(gc(full = FALSE))
    }

    # x 100: birds/ha -> birds per km2 pixel. a coalition's "bf" is its counterfactual on its
    # footprint: backfilled vegetation on its direct pixels, d0 (the same in every scenario) on the rest
    bf <- lapply(seq_along(all_cids), function(ci) {
      if (coal_state[ci] != 2L) return(NULL)
      coal_rowsum(sc, direct_rows[[ci]], zi_c, n_sub) * 100 + d0_infl[[ci]][, rep(1L, n_scen), drop = FALSE]
    })
    if (coal_state[ci_full] == 2L) {
      dr <- direct_rows[[ci_full]]
      cf_full[dr] <- rowMeans(sc[dr, , drop = FALSE])
    }
    list(bf = bf, d0 = d0, cf_full = cf_full, n_bad = as.numeric(n_bad + n_bad0))
  }, mc.cores = n_cores, mc.preschedule = FALSE)

  # a worker the kernel kills (out of memory) returns NULL
  died <- vapply(boot_res, is.null, logical(1L))
  if (any(died))
    stop(species, " ", bcr_code, ": ", sum(died), " of ", length(boot_res), " bootstrap workers died ",
         "(most likely out of memory: look for oom_kill in the job log and raise --mem)")
  failed <- vapply(boot_res, inherits, logical(1L), "try-error")
  if (any(failed))
    stop(species, " ", bcr_code, ": ", sum(failed), " bootstrap workers failed: ", as.character(boot_res[[which(failed)[1L]]]))
  n_bad <- sum(vapply(boot_res, function(r) r$n_bad, numeric(1L)))
  if (n_bad > 0)
    stop(species, " ", bcr_code, ": ", n_bad, " NA cells in the bf/d0 fields at complete pixels: some BART draw ",
         "other than draw 1 has thinner coverage than the complete mask assumes")

  # observed birds/ha at the complete superset pixels, [n_complete x n_boot]
  O_c <- vapply(obs_preds, function(r) terra::values(r, mat = FALSE)[super_idx[c_rows]], numeric(length(c_rows)))
  if (is.null(dim(O_c))) O_c <- matrix(O_c, ncol = n_boot)

  # each subbasin's whole observed population, per bootstrap (NA where its zone is absent)
  zone_slot <- match(hybas_ids, sort(unique(terra::values(subbasin_zone_r, na.rm = TRUE))))
  obs_total_byboot <- vapply(seq_len(n_boot), function(i) {
    z <- terra::zonal(obs_preds[[i]] * 100, subbasin_zone_r, "sum", na.rm = TRUE)
    z[[2]][zone_slot]
  }, numeric(n_sub))
  if (is.null(dim(obs_total_byboot))) obs_total_byboot <- matrix(obs_total_byboot, nrow = n_sub)

  rm(stack_obs, stack_bf, obs_preds, X_pred, draw_vals_pred, cat_vals_pred, b.list)
  gc()

  combine_stats <- function(mat) list(mean = rowMeans(mat, na.rm = TRUE), sd = apply(mat, 1L, sd, na.rm = TRUE))

  # [n_sub x n_boot] -> [n_sub x (n_boot * n_scen)], samples running scenario-within-bootstrap
  bcols         <- rep(seq_len(n_boot), each = n_scen)
  obs_total_mat <- obs_total_byboot[, bcols, drop = FALSE]
  ot_stats      <- combine_stats(obs_total_mat)

  # BCR totals per (bootstrap, scenario) for 15A's national arrays
  bcr_draws  <- function(mat) matrix(colSums(mat, na.rm = TRUE), nrow = n_boot, ncol = n_scen, byrow = TRUE)
  arr_cids   <- intersect(as.character(save_arrays_ids), as.character(all_cids))
  bcr_arrays <- setNames(vector("list", length(arr_cids)), arr_cids)

  # per-sample subbasin Shapley values for 14B: Shapley values are linear in v (phi = W v, 12E),
  # so they can be taken per (bootstrap, scenario) sample. summing those samples keeps the
  # covariances that propagated SDs lose: a BCR's subbasins share 32 bird models, and v(S) and
  # v(S + j) share samples. the direct part depends on the bootstrap only; the vegetation part
  # of sample (i, k) is phi_samples[, , (i, k)] - phi_direct[, , i]
  W_shap      <- shapley_weight_matrix(sectors)
  phi_samples <- array(0, c(n_sub, n_sectors, n_boot * n_scen), dimnames = list(NULL, sectors, NULL))
  phi_direct  <- array(0, c(n_sub, n_sectors, n_boot), dimnames = list(NULL, sectors, NULL))

  # ---- assemble every coalition's table ----
  coalition_tables <- setNames(vector("list", length(all_cids)), as.character(all_cids))

  for (ci in seq_along(all_cids)) {
    cid <- all_cids[ci]

    # no footprint in the BCR: all zeros, obs_total included (as the retired per-coalition code
    # wrote it; 14B takes obs_total from the other coalitions). the arrays keep the real obs_total
    if (coal_state[ci] == 0L) {
      coalition_tables[[ci]] <- tibble::tibble(
        species = species, subbasin = sub_ids, bcr = bcr_code, coalition_id = cid,
        obs_total_mean = 0, obs_total_sd = 0,
        obs_on_coalition_mean = 0, obs_on_coalition_sd = 0,
        bf_on_coalition_mean = 0, bf_on_coalition_sd = 0)
      if (as.character(cid) %in% arr_cids)
        bcr_arrays[[as.character(cid)]] <- list(obs_total_mat   = bcr_draws(obs_total_mat),
                                                obs_on_coal_mat = matrix(0, n_boot, n_scen),
                                                bf_on_coal_mat  = matrix(0, n_boot, n_scen))
      next
    }

    # a subbasin without kept pixels has obs_on NA (NaN mean, NA sd), never a 0
    bf_mat     <- matrix(0, nrow = n_sub, ncol = n_boot * n_scen)
    obs_on_mat <- matrix(NA_real_, n_sub, n_boot * n_scen)
    if (coal_state[ci] == 2L) {
      kr     <- keep_rows[[ci]]
      bf_mat <- do.call(cbind, lapply(boot_res, function(r) r$bf[[ci]]))   # bootstrap-major

      Ok <- O_c[kr, , drop = FALSE] * 100
      Ok[is.na(Ok)] <- 0
      rs_obs <- rowsum(Ok, group = zones_c[kr], reorder = TRUE)
      obs_on_byboot <- matrix(NA_real_, n_sub, n_boot)
      map_o <- match(hybas_ids, rownames(rs_obs))
      presO <- !is.na(map_o)
      if (any(presO)) obs_on_byboot[presO, ] <- rs_obs[map_o[presO], , drop = FALSE]
      obs_on_mat <- obs_on_byboot[, bcols, drop = FALSE]

      # direct part of v(S), per bootstrap
      d0_byboot <- do.call(cbind, lapply(boot_res, function(r) r$d0[[ci]]))   # [n_sub x n_boot]
      imp_d     <- d0_byboot - obs_on_byboot
      empty_d   <- is.na(obs_on_byboot)
      if (any(d0_byboot[empty_d] != 0)) stop(species, " ", bcr_code, " coalition ", cid, ": d0 != 0 in a subbasin without kept pixels")
      imp_d[empty_d] <- 0
      for (j in seq_len(n_sectors)) phi_direct[, j, ] <- phi_direct[, j, ] + W_shap[j, cid] * imp_d
      rm(d0_byboot, imp_d, empty_d)
    }

    # v(S) per sample; a subbasin without kept pixels has nothing backfilled, so v(S) = 0 there
    imp   <- bf_mat - obs_on_mat
    empty <- is.na(obs_on_mat)
    if (any(bf_mat[empty] != 0)) stop(species, " ", bcr_code, " coalition ", cid, ": bf != 0 in a subbasin without kept pixels")
    imp[empty] <- 0
    for (j in seq_len(n_sectors)) phi_samples[, j, ] <- phi_samples[, j, ] + W_shap[j, cid] * imp
    rm(imp, empty)

    bf_stats  <- combine_stats(bf_mat)
    ooc_stats <- combine_stats(obs_on_mat)
    if (as.character(cid) %in% arr_cids)
      bcr_arrays[[as.character(cid)]] <- list(obs_total_mat   = bcr_draws(obs_total_mat),
                                              obs_on_coal_mat = bcr_draws(obs_on_mat),
                                              bf_on_coal_mat  = bcr_draws(bf_mat))
    coalition_tables[[ci]] <- tibble::tibble(
      species = species, subbasin = sub_ids, bcr = bcr_code, coalition_id = cid,
      obs_total_mean = ot_stats$mean,  obs_total_sd = ot_stats$sd,
      obs_on_coalition_mean = ooc_stats$mean, obs_on_coalition_sd = ooc_stats$sd,
      bf_on_coalition_mean = bf_stats$mean,   bf_on_coalition_sd = bf_stats$sd)
  }

  # ---- per-pixel change of the all-sector coalition, for 18's maps ----
  # the mean observed and counterfactual birds/km2 (weighted) at each pixel the tables sum; per
  # subbasin they must add up to the coalition-256 table, so the maps add up to the tables
  if (coal_state[ci_full] == 2L) {
    kr_full <- keep_rows[[ci_full]]
    cf_px   <- Reduce(`+`, lapply(boot_res, `[[`, "cf_full"))[kr_full] / n_boot * 100
    O0      <- O_c[kr_full, , drop = FALSE]
    O0[is.na(O0)] <- 0
    obs_px  <- rowMeans(O0) * 100
    zk      <- factor(zi_c[kr_full], levels = seq_len(n_sub))
    sub_sum <- function(v) { s <- as.numeric(tapply(v, zk, sum)); s[is.na(s)] <- 0; s }
    tb      <- coalition_tables[[ci_full]]
    has_px  <- !is.nan(tb$obs_on_coalition_mean)
    if (!isTRUE(all.equal(sub_sum(cf_px), tb$bf_on_coalition_mean, tolerance = 1e-8)) ||
        !isTRUE(all.equal(sub_sum(obs_px)[has_px], tb$obs_on_coalition_mean[has_px], tolerance = 1e-8)))
      stop(species, " ", bcr_code, ": per-pixel means do not sum to the coalition-256 table")
    px_dir <- file.path(ia_dir, "data", "derived_data", "density_tables", "pixel_change")
    dir.create(px_dir, recursive = TRUE, showWarnings = FALSE)
    px_r   <- terra::rast(subbasin_zone_r, nlyrs = 2L)
    px_val <- matrix(NA_real_, terra::ncell(px_r), 2L)
    px_val[super_idx[c_rows[kr_full]], ] <- cbind(obs_px, cf_px)
    terra::values(px_r) <- px_val
    names(px_r) <- c("obs", "cf")
    px_f <- file.path(px_dir, paste0(species, "_", year, "_", bcr_code, ".tif"))
    terra::writeRaster(px_r, sub("\\.tif$", "_partial.tif", px_f), overwrite = TRUE, datatype = "FLT4S", gdal = "COMPRESS=DEFLATE")
    if (!file.rename(sub("\\.tif$", "_partial.tif", px_f), px_f)) stop("could not move ", px_f, " into place")
    message(Sys.time(), " | wrote ", px_f, " (", length(kr_full), " footprint pixels)")
  }

  list(coalition_tables = coalition_tables,
       bcr_arrays       = bcr_arrays,
       shapley_samples  = list(phi        = phi_samples,       # [n_sub x sector x sample]
                               phi_direct = phi_direct,        # [n_sub x sector x n_boot]
                               obs_total  = obs_total_byboot,  # [n_sub x n_boot]
                               subbasin = sub_ids, n_boot = n_boot, n_scen = n_scen),
       bcr_code         = bcr_code)
}


# binds per-BCR results (in the order given) into the 255 coalition tables, the national arrays
# and the per-sample Shapley values; 12H passes them in each species' 06_bootstraps order,
# which fixes the floating-point order of the national sums
combine_bcr_results <- function(per_bcr, save_arrays_ids = integer(0)) {
  per_bcr <- Filter(Negate(is.null), per_bcr)

  cid_keys <- as.character(2:(2L^length(canonical_sectors())))
  tables_by_cid <- setNames(lapply(cid_keys, function(k)
    dplyr::bind_rows(Filter(Negate(is.null), lapply(per_bcr, function(x) x$coalition_tables[[k]])))), cid_keys)

  # national [n_boot x n_scen] arrays for 15A: sums of the per-BCR arrays
  arr_keys <- intersect(as.character(save_arrays_ids), cid_keys)
  arrays_by_cid <- setNames(lapply(arr_keys, function(k) {
    parts <- Filter(Negate(is.null), lapply(per_bcr, function(x) x$bcr_arrays[[k]]))
    if (length(parts) == 0L) return(NULL)
    lapply(setNames(nm = c("obs_total_mat", "obs_on_coal_mat", "bf_on_coal_mat")),
           function(f) Reduce(`+`, lapply(parts, `[[`, f)))
  }), arr_keys)

  if (length(per_bcr) == 0L) return(list(tables_by_cid = tables_by_cid, arrays_by_cid = arrays_by_cid, shapley_samples = NULL))

  # per-sample Shapley values stacked across BCRs in table row order. 14B sums them sample by
  # sample, so every BCR needs the same (bootstrap, scenario) grid
  ss <- lapply(per_bcr, `[[`, "shapley_samples")
  n_samp <- vapply(ss, function(s) dim(s$phi)[3L], integer(1L))
  if (length(unique(n_samp)) > 1L)
    stop("per-BCR Shapley samples differ in size (", paste0(vapply(per_bcr, `[[`, character(1L), "bcr_code"), "=", n_samp, collapse = ", "), ")")
  n_rows <- vapply(ss, function(s) dim(s$phi)[1L], integer(1L))
  stack_rows <- function(f) {
    a <- array(0, c(sum(n_rows), dim(ss[[1L]][[f]])[2:3]), dimnames = list(NULL, dimnames(ss[[1L]][[f]])[[2L]], NULL))
    ends <- cumsum(n_rows)
    for (b in seq_along(ss)) a[(ends[b] - n_rows[b] + 1L):ends[b], , ] <- ss[[b]][[f]]
    a
  }
  list(tables_by_cid   = tables_by_cid,                       # "2".."256" -> tibble
       arrays_by_cid   = arrays_by_cid,                       # save_arrays_ids -> 3 national matrices
       shapley_samples = list(
         bcr        = rep(vapply(per_bcr, `[[`, character(1L), "bcr_code"), n_rows),
         subbasin   = unlist(lapply(ss, `[[`, "subbasin")),
         phi        = stack_rows("phi"),                      # [row x sector x sample]
         phi_direct = stack_rows("phi_direct"),               # [row x sector x bootstrap]
         obs_total  = do.call(rbind, lapply(ss, `[[`, "obs_total")),
         n_boot     = ss[[1L]]$n_boot, n_scen = ss[[1L]]$n_scen))
}


# coalitions whose national arrays are kept for 15A: all sectors, and each sector alone
coalition_array_ids <- function(sectors = canonical_sectors()) {
  c(sectors_to_coalition_id(sectors, sectors),
    vapply(sectors, function(s) sectors_to_coalition_id(s, sectors), numeric(1L)))
}


# 12D's array tasks: one row per species x BCR, each species' 06_bootstraps files in
# list.files() order (bcr_order, for 12H). TEST_BCR filters it. 12D and 12H must agree on it
coalition_task_table <- function(species_vec, nm_root, test_bcr = Sys.getenv("TEST_BCR", "")) {
  tt <- dplyr::bind_rows(lapply(species_vec, function(sp) {
    f <- list.files(file.path(nm_root, "output/06_bootstraps", sp), pattern = "can.*\\.Rdata$", full.names = TRUE)
    if (length(f) == 0L) return(NULL)
    tibble::tibble(species = sp, rdata_path = f,
                   bcr = sub(paste0("^", sp, "_(.*)\\.Rdata$"), "\\1", basename(f)),
                   bcr_order = seq_along(f))
  }))
  if (nchar(test_bcr) > 0L) tt <- tt[tt$bcr %in% strsplit(test_bcr, ",")[[1]], ]
  tt$task <- seq_len(nrow(tt))
  tt
}
