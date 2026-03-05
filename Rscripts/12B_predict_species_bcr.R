# define density prediction function ------------------------------------------------------
predict_species_bcr <- function(species, year, all_subbasins_subset, sector_name, hirsh_dir) {
  
  # communicate
  message(Sys.time(), " | preparing for species=", species, " year=", year)
  
  # find all models (BCRs) for the current `species`
  rdata_files <- list.files(file.path(nm_root,"output/06_bootstraps", species), pattern = "can14.*\\.Rdata$", full.names = TRUE)
  message(Sys.time(), " | ", species, " | found ", length(rdata_files), " BCR models")
  
  # this sub-function will work through all model-lists (one per BCR) for a given species
  do.call(dplyr::bind_rows, lapply(rdata_files, function(rdata_path) {
    
    message(Sys.time(), " | loading ", basename(rdata_path))
    e <- new.env(parent = emptyenv())
    
    # creates a `b.list` object in the environment
    load(rdata_path, envir = e)
    
    if (!exists("b.list", envir = e)) {
      stop("b.list not found in ", basename(rdata_path))
    }
    
    # loads b.list for some spp x bcr
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

    # build sector mask: sector > 0 AND CanHF >= 1, projected to BCR prediction grid
    sector_r  <- terra::project(
                   terra::rast(file.path(hirsh_dir, paste0(sector_name, ".tif"))),
                   stack_obs, method = "near")
    canHF_r   <- terra::project(
                   terra::rast(file.path(hirsh_dir, "CanHF_1km_morethan1.tif")),
                   stack_obs, method = "near")
    sector_mask <- terra::ifel((sector_r > 0) & (canHF_r >= 1), 1, NA)
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
    footprint_layer <- sector_mask
    message(Sys.time(), " | done loading backfilled stack")
    
    # split BART outputs (mean + empirical posterior quantiles)
    mu_stack <- stack_bf[[grep("_mean$", names(stack_bf))]]
    names(mu_stack)   <- sub("_mean$", "", names(mu_stack))

    q025_stack <- stack_bf[[grep("_q025$", names(stack_bf))]]
    names(q025_stack) <- sub("_q025$", "", names(q025_stack))

    q975_stack <- stack_bf[[grep("_q975$", names(stack_bf))]]
    names(q975_stack) <- sub("_q975$", "", names(q975_stack))

    # construct counterfactual scenarios using empirical posterior quantiles
    X_bf_sets <- list(
      mean = make_counterfactual_stack(stack_obs, mu_stack,    sector_mask),
      low  = make_counterfactual_stack(stack_obs, q025_stack,  sector_mask),
      high = make_counterfactual_stack(stack_obs, q975_stack,  sector_mask)
    )
    
    # create containers to store predictions
    obs_preds <- vector("list", length(b.list))
    bf_preds  <- list(
      mean = vector("list", length(b.list)),
      low  = vector("list", length(b.list)),
      high = vector("list", length(b.list))
    )
    
    # run predictions for 32 bootstraps in `b.list` (a single spp x BCR)
    for (i in seq_along(b.list)) {
      
      model <- b.list[[i]]
      model_vars <- model$var.names
      
      # cont_vars has abiotic *and* biotic continuous variables 
      cont_vars <- setdiff(model_vars, categorical_responses)
      cat_vars  <- intersect(model_vars, categorical_responses)
      
      # check for continuous biotic covariates missing from backfilled stack
      biotic_cont_vars <- intersect(cont_vars, biotic_continuous_vars)
      
      missing_bf <- biotic_cont_vars[!(paste0(biotic_cont_vars, "_mean") %in% names(stack_bf) &
                                         paste0(biotic_cont_vars, "_q025") %in% names(stack_bf) &
                                         paste0(biotic_cont_vars, "_q975") %in% names(stack_bf))]
      
      if (length(missing_bf) > 0) {
        message(Sys.time(), " | ", species, " ", bcr_code, " | bootstrap=", i, " | missing backfilled cont vars: ", paste(missing_bf, collapse = ", "))
      }
      
      # observed landscape for b.list[[i]]
      X_obs_i <- stack_obs[[intersect(model_vars, names(stack_obs))]]
      
      # overwrite categorical biotic covariates at sector pixels only
      # (stack_bf has backfilled values at all high-HF pixels; sector_mask restricts to sector)
      for (v in cat_vars) {
        X_bf_sets$mean[[v]] <- terra::ifel(sector_mask, stack_bf[[v]], stack_obs[[v]])
        X_bf_sets$low [[v]] <- terra::ifel(sector_mask, stack_bf[[v]], stack_obs[[v]])
        X_bf_sets$high[[v]] <- terra::ifel(sector_mask, stack_bf[[v]], stack_obs[[v]])
      }
      
      # set disturbance vars to zero at sector pixels only
      dist_i <- intersect(disturbance_vars$predictor, model_vars)

      for (k in names(X_bf_sets)) {
        for (v in dist_i) {
          X_bf_sets[[k]][[v]] <- terra::ifel(sector_mask, 0, stack_obs[[v]])
        }
      }
      
      zeroed <- intersect(names(X_bf_sets$mean), dist_i)
      message(Sys.time(), " | ", species, " ", bcr_code, " bootstrap=", i, " | zeroed out disturbance vars: ", paste(zeroed, collapse = ", "))

      # predict bird density on observed landscape
      obs_preds[[i]] <- terra::predict(X_obs_i, model, type = "response", n.trees = model$n.trees)
      message(Sys.time(), " | ", species, " ", bcr_code, " bootstrap=", i, " | finished predict() on observed landscape")

      # predict bird density on backfilled landscape
      for (k in names(X_bf_sets)) {
        bf_preds[[k]][[i]] <- terra::predict(X_bf_sets[[k]], model, type = "response", n.trees = model$n.trees)
      }
      message(Sys.time(), " | ", species, " ", bcr_code, " bootstrap=", i, " | finished predict() on backfilled landscape")
      
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

    # backfilled
    bf_summary <- lapply(bf_preds, summarize_preds)

    # write observed rasters
    terra::writeRaster(obs_mean, file.path(obs_dir, "observed_mean.tif"), overwrite = TRUE)
    terra::writeRaster(obs_sd,   file.path(obs_dir, "observed_sd.tif"),   overwrite = TRUE)

    # write backfilled rasters
    for (k in names(bf_summary)) {
      terra::writeRaster(bf_summary[[k]]$mean, file.path(bf_dir, paste0("backfilled_", k, "_mean.tif")), overwrite = TRUE)
      terra::writeRaster(bf_summary[[k]]$sd,   file.path(bf_dir, paste0("backfilled_", k, "_sd.tif")),   overwrite = TRUE)
    }
    
    
    
    # count birds from prediction surfaces intersecting our watersheds ------------------------------------------------------
    
    # force agg() to always return a matrix, even when there’s only one subbasin
    message(Sys.time(), " | ", species, " ", bcr_code, " | estimating population over ", length(sub_ids), " subbasins")
    
    agg <- function(obs_rasters, bf_rasters_list, sub_ids, footprint_layer) {
      
      n_sub  <- length(sub_ids)
      n_boot <- length(obs_rasters)
      n_scen <- length(bf_rasters_list) # e.g., mean/low/high
      
      # initialize arrays: sub x boot x scenario
      pop_obs_total_arr   <- array(NA, dim = c(n_sub, n_boot, n_scen))
      pop_obs_on_bf_arr   <- array(NA, dim = c(n_sub, n_boot, n_scen))
      pop_bf_only_arr     <- array(NA, dim = c(n_sub, n_boot, n_scen))
      
      scenario_names <- names(bf_rasters_list)
      
      for (s in seq_len(n_scen)) {
        
        bf_rasters <- bf_rasters_list[[s]]
        
        for (i in seq_len(n_boot)) {
          
          obs_r <- obs_rasters[[i]] * 100
          bf_r  <- bf_rasters[[i]] * 100
          
          for (j in seq_along(sub_ids)) {
            
            sub_idx <- sub_ids[j]
            poly <- all_subbasins_subset[sub_idx, ]
            
            # mask using the polygon directly
            pop_obs_total_arr[j, i, s] <- terra::extract(obs_r, poly, fun = sum, na.rm = TRUE)[,2]
            
            # restrict footprint to this subbasin
            fp_sub <- terra::mask(footprint_layer, poly)
            
            obs_fp <- terra::mask(obs_r, fp_sub)
            bf_fp  <- terra::mask(bf_r,  fp_sub)
            
            pop_obs_on_bf_arr[j, i, s] <- terra::extract(obs_fp, poly, fun = sum, na.rm = TRUE)[,2]
            pop_bf_only_arr[j, i, s] <- terra::extract(bf_fp, poly, fun = sum, na.rm = TRUE)[,2]
            
          } # close for j in seq
          
        } # close for i in seq
        
      } # close for s in seq
      
      # combine over bootstraps and scenarios
      combine_stats <- function(x) {
        
        # collapse to sub x (boot*scenario) matrix
        mat <- matrix(x, nrow = n_sub, ncol = n_boot * n_scen)
        list(mean = rowMeans(mat, na.rm = TRUE),
             sd   = apply(mat, 1, sd, na.rm = TRUE))
        
      } # close combine_stats()
      
      list(
        pop_obs_total = combine_stats(pop_obs_total_arr),
        pop_obs_on_bf = combine_stats(pop_obs_on_bf_arr),
        pop_bf_only   = combine_stats(pop_bf_only_arr)
      )
      
    } # close agg()
    
    
    # run aggregation including posterior
    pop_lists <- agg(obs_preds, bf_preds, sub_ids, footprint_layer)
    
    # for each subbasin, compute counterfactual total population:
    counterfactual_total <-
      pop_lists$pop_obs_total$mean -
      pop_lists$pop_obs_on_bf$mean +
      pop_lists$pop_bf_only$mean
    
    
    # free up some memory
    rm(obs_preds, bf_preds, stack_obs, stack_bf)
    gc()
    
    # summarize into tidy tibble including SD for obs_on_bf and bf_only
    out <- tibble(
      species = species,
      subbasin = sub_ids,
      bcr      = bcr_code,
      sector   = sector_name,
      obs_total_mean    = pop_lists$pop_obs_total$mean,
      obs_total_sd      = pop_lists$pop_obs_total$sd,
      obs_on_bf_mean    = pop_lists$pop_obs_on_bf$mean,
      obs_on_bf_sd      = pop_lists$pop_obs_on_bf$sd,
      bf_only_mean      = pop_lists$pop_bf_only$mean,
      bf_only_sd        = pop_lists$pop_bf_only$sd,
      counterfactual_mean = counterfactual_total,
      counterfactual_sd   = sqrt(pop_lists$pop_obs_total$sd^2 +
                                   pop_lists$pop_obs_on_bf$sd^2 +
                                   pop_lists$pop_bf_only$sd^2)  # approximate error propagation
    )    
    
    message(Sys.time(), " | ", species, " ", bcr_code, " | returning ", nrow(out), " rows")
    out
    
  }))
} # close predict_species_bcr()