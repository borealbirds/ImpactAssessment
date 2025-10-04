# ---
# title: Impact Assessment: re-predict bird densities after backfilling
# author: Mannfred Boehm
# created: September 28, 2025
# ---

library(furrr)
library(terra)
library(tidyverse)


# set root path
root <- "G:/Shared drives/BAM_NationalModels5"
ia_dir <- file.path(root, "data", "Extras", "sandbox_data", "impactassessment_sandbox")



# import data ------------------------------------------------------

# import BCR boundaries
bam_boundary <- terra::vect(file.path(root, "Regions", "BAM_BCR_NationalModel_Unbuffered.shp")) 

# import merged + subsetted subbasins
all_subbasins_subset <- vect(file.path(ia_dir, "hydrobasins_masked_merged_subset.gpkg"))

# import industry footprint pixels
combined_poly_path <- file.path(ia_dir, "combined_industry_footprint.gpkg")

# create a reference table for which subbasins are in which BCRs 
# (some subbasins will be in multiple BCRs, and that's OK)
bcr_subbasins_ref <-
  {
    # logical matrix: rows=subbasins, cols=BCRs
    hits <- terra::relate(centroids(all_subbasins_subset), bam_boundary, relation = "intersects")
    
    # row/col indexes of TRUE
    ij <- which(hits, arr.ind = TRUE)
    
    tibble(
      HYBAS_ID = all_subbasins_subset$first_HYBAS_ID[ij[, 1]],
      bcr_label = paste(bam_boundary$country[ij[, 2]], bam_boundary$subUnit[ij[, 2]], sep = "_"),
      bcr_code = gsub("_", "", bcr_label) # e.g., "can_14" -> "can14" for filenames like PIGR_can14.Rdata
    )
  }


# helper: mosaic (backfilled) subbasin stacks into a single BCR-wide stack
.mosaic_backfilled_stacks <- function(sub_ids, year, template = NULL) {
  
  sub_sv <- terra::vect(subbasins_path)
  
  ids_chr <- as.character(sub_ids)
  
  # lookup index used during training: position in the reference object
  idx_chr <- match(sub_ids, sub_sv$first_HYBAS_ID)
 
  ptab <- tibble(
    HYBAS_ID = ids_chr,
    subbasin_index = idx_chr,
    path = file.path(
      ia_dir, "xgboost_models", paste0("year=", as.character(year)),
      paste0("subbasin=", as.character(idx_chr)),
      paste0("backfilled_stack_subbasin-", sprintf("%03d", idx_chr), ".tif"))) |> 
    dplyr::filter(!is.na(subbasin_index), file.exists(path))
  
  if (nrow(ptab) == 0) return(NULL)
  
  stacks <- setNames(lapply(ptab$path, terra::rast), ptab$HYBAS_ID)
  
  # if a template is provided, resample each stack to it (nearest neighbor)
  if (!is.null(template)) {
    stacks <- lapply(stacks, \(r) terra::resample(r, template, method = "near"))
  } else {
    # else, use the first stack as the internal template
    ref <- stacks[[1]]
    stacks <- lapply(stacks, \(r) if (!terra::compareGeom(r, ref, stopOnError = FALSE)) terra::resample(r, ref, method = "near") else r)
  }
  
  
  # union of layer names across stacks
  vars <- sort(unique(unlist(lapply(stacks, names))))
  if (length(vars) == 0) return(NULL)
  
  # mosaic per-variable via cover()
  mosaic_one <- function(v) {
    rlist <- lapply(stacks, \(r) if (v %in% names(r)) r[[v]] else NULL) |> purrr::compact()
    if (length(rlist) == 1) rlist[[1]] else Reduce(terra::cover, rlist)
  }
  
  out <- terra::rast(lapply(vars, mosaic_one))
  names(out) <- vars
  out
}

# define disturbance variables, which we'll set to zero when re-predicting
disturbance_vars <-
  dplyr::tibble(BAMexploreR::predictor_metadata) |>
  dplyr::filter(version == "v5") |>
  dplyr::select(predictor, definition, predictor_class) |>
  dplyr::filter(predictor_class == "Disturbance")

# for running in parallel
bam_boundary_path <- file.path(root, "Regions", "BAM_BCR_NationalModel_Unbuffered.shp")
subbasins_path    <- file.path(ia_dir, "hydrobasins_masked_merged_subset.gpkg")
combined_poly_path <- file.path(ia_dir, "combined_industry_footprint.gpkg")


# define density prediction function -------------------------------
predict_bird_density_per_bcr_year <- function(species, year, bcr_subbasins_ref){
 
  # fetch the 22 BCRs that intersect our subbasins
  bcrs <- unique(bcr_subbasins_ref$bcr_label)
  
  # loop: for every BCR, get counterfactual species density estimates 
  do.call(rbind, lapply(bcrs, function(bcr_label) {
    
    message("working on ", bcr_label, " ", species, " ", year)
    
    tryCatch({
    
    # load locally for every worker
    bam_boundary         <- terra::vect(bam_boundary_path)
    all_subbasins_subset <- terra::vect(subbasins_path)
    
    # get BCR code
    bcr_code <- gsub("_", "", bcr_label)
  
    # subset bam_boundary to the current BCR
    bcr_subUnit <- as.integer(sub(".*_(\\d+)$", "\\1", bcr_label))
    if (is.na(bcr_subUnit)) return(NULL)  # bail out cleanly
    
    bcr_subUnit <- as.numeric(str_extract(bcr_code, "(\\d)+"))
    bcr_poly <- bam_boundary[bam_boundary$country == "can" & bam_boundary$subUnit == bcr_subUnit, ]
    bcr_poly <- if (length(bcr_poly) > 1) terra::aggregate(bcr_poly) else bcr_poly

    # subset subbasins to current BCR
    s <- which(terra::relate(all_subbasins_subset, bcr_poly, relation = "intersects"))
    if (length(s) == 0) return(NULL)
    subbasins_in_bcr <- all_subbasins_subset[s,]
    
      
    # observed stack --------------------------------
    obs_path <- file.path(root, "gis", "stacks", paste0(bcr_code, "_", year, ".tif"))
    if (!file.exists(obs_path)) return(NULL)  
    
    # subset industry footprints to current BCR (may be empty)
    industry_bcr <- tryCatch(
      terra::crop(terra::vect(combined_poly_path), bcr_poly) |> terra::mask(bcr_poly),
      error = function(e) NULL
      )
    
    roi <- if (!is.null(industry_bcr) && terra::nrow(industry_bcr) > 0) industry_bcr else bcr_poly
    
    # subset covariate stack to industry pixels
    stack_obs <- terra::rast(obs_path) |>
      (\(r) terra::mask(terra::crop(r, roi), roi))() 
    
    
    # backfilled stack --------------------------------
    
    stack_bf <-
      .mosaic_backfilled_stacks(subbasins_in_bcr$first_HYBAS_ID, year, template = stack_obs) |>
      (\(r) if (is.null(r)) return(NULL) else terra::mask(terra::crop(r, bcr_poly), bcr_poly))()
    
    if (is.null(stack_bf)) return(NULL)
    
    # fill out backfilled raster with abiotic and disturbance covariates
    abiotic_and_disturbance_vars <- setdiff(names(stack_obs), names(stack_bf))
    if (length(abiotic_and_disturbance_vars) > 0) {
      stack_bf <- c(stack_bf, stack_obs[[abiotic_and_disturbance_vars]])
    }
    
    # zero out disturbance layers in backfilled stack
    disturbance_layers <- intersect(names(stack_bf), disturbance_vars$predictor)
    if (length(disturbance_layers) > 0) {
      stack_bf[[disturbance_layers]] <- 0
      stack_bf[[disturbance_layers]] <- mask(stack_bf[[disturbance_layers]], stack_obs[[disturbance_layers]])
    }
                        
    
    # use V5 models for density estimation----------------
    
    # load bootstrap models for this species x BCR
    rdata_path <- file.path(root, "output", "06_bootstraps", species, paste0(species, "_", bcr_code, ".Rdata"))
    if (!file.exists(rdata_path)) return(NULL)
    
    loaded_names <- load(rdata_path)
    mdl_list <- get(loaded_names[1])  # assume the .Rdata contains a single object: list of 32 gbm models
    
    # variable set (order & match to model)
    # take the first model’s variable names; gbm stores them in attr(m, "var.names") or m$var.names
    varnames <- if (!is.null(attr(mdl_list[[1]], "var.names"))) attr(mdl_list[[1]], "var.names") else mdl_list[[1]]$var.names
    
    # restrict and order stacks to model vars (drop any extras)
    X_obs <- stack_obs[[intersect(varnames, names(stack_obs))]]
    X_bf  <- stack_bf [[intersect(varnames, names(stack_bf ))]]
    
    # align to the same var set across both stacks
    keep <- intersect(names(X_obs), names(X_bf))
    if (length(keep) == 0) return(NULL)
    X_obs <- X_obs[[keep]]
    X_bf  <- X_bf [[keep]]

    # create output folder for mean density rasters
    outdir <- file.path(ia_dir, "density_predictions", species, as.character(year))
    if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
    
    # per-bootstrap predictions and subbasin sums ----------------
    
    # (sum over pixels; if your model predicts density per unit area and pixel areas are equal, this is proportional to population)
    # replace type="response" if your gbm objects need different predict args
    pred_sum <- function(r_pred) {
      
      # get one-column result named 'sum' (no ID column)
      out <- terra::extract(r_pred * 100, subbasins_in_bcr, fun = sum, na.rm = TRUE, ID = FALSE)
      
      tibble::tibble(HYBAS_ID = subbasins_in_bcr$first_HYBAS_ID, sum = out$lyr1)
    }
    
    # define function that searches for n trees
    gbm_fun <- function(m, d) {
      nt <- if (!is.null(m$best.iter)) m$best.iter else
            if (!is.null(m$n.trees))  m$n.trees  else
            if (!is.null(m$gbm.call$ntrees)) m$gbm.call$ntrees else
            stop("n.trees (or best.iter) not found in GBM model")
      gbm::predict.gbm(m, d, n.trees = nt, type = "response")
    }
    
    # create density surfaces for 32 bootstraps
    obs_density <- lapply(mdl_list, \(m) terra::predict(X_obs, m, fun = gbm_fun))    
    bf_density  <- lapply(mdl_list, \(m) terra::predict(X_bf, m, fun = gbm_fun))
    
    # write mean prediction rasters
    obs_mean <- mean(terra::rast(obs_density)) * 100  # *100 to convert ha → km²
    bf_mean  <- mean(terra::rast(bf_density)) * 100
    
    obs_fn <- file.path(outdir, paste0(species, "_", bcr_code, "_", year, "_observed.tif"))
    bf_fn  <- file.path(outdir, paste0(species, "_", bcr_code, "_", year, "_backfilled.tif"))
    
    terra::writeRaster(obs_mean, obs_fn, overwrite = TRUE)
    terra::writeRaster(bf_mean , bf_fn , overwrite = TRUE)
    
    
    # get population estimates
    obs_sums <- lapply(obs_density, pred_sum)
    bf_sums <- lapply(bf_density, pred_sum)
    
    # bind by HYBAS_ID, compute mean/sd across bootstraps
    obs_mat <- reduce(obs_sums, left_join, by = "HYBAS_ID") |> arrange(HYBAS_ID)
    bf_mat  <- reduce(bf_sums , left_join, by = "HYBAS_ID") |> arrange(HYBAS_ID)
    
    # drop HYBAS_ID to get numeric matrices
    obs_vals <- as.matrix(obs_mat[,-1, drop = FALSE])
    bf_vals  <- as.matrix(bf_mat [,-1, drop = FALSE])
    
    tibble(
      species         = species,
      subbasin        = obs_mat$HYBAS_ID,
      bcr             = bcr_code,
      treatment       = "observed",
      population_mean = rowMeans(obs_vals, na.rm = TRUE),
      population_sd   = apply(obs_vals, 1, sd, na.rm = TRUE)) |>
      bind_rows(
        tibble(
          species         = species,
          subbasin        = bf_mat$HYBAS_ID,
          bcr             = bcr_code,
          treatment       = "backfilled",
          population_mean = rowMeans(bf_vals, na.rm = TRUE),
          population_sd   = apply(bf_vals, 1, sd, na.rm = TRUE)))
    },error = function(e) {
      message("Skipping ", bcr_label, " (", species, " ", year, "): ", conditionMessage(e))
      return(NULL)
    })
 })) 
}

# set up for parallel processing
#plan(multisession, workers = max(1, parallel::detectCores() - 1)) 
#plan(multisession, workers = 2)
plan(sequential)

# (optional) reproducible RNG across workers for any stochastic bits
opts <- furrr::furrr_options(seed = TRUE, scheduling = 1, packages = c("dplyr","purrr","terra","stringr","tibble","gbm"))

# run in parallel for 
res_cawa_osfl_2020 <-
  future_map_dfr(c("CAWA", "OSFL"),
    ~ predict_bird_density_per_bcr_year(.x, 2020, bcr_subbasins_ref),
    .options = opts,
    .progress = TRUE)
