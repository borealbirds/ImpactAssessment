# ---
# title: Impact Assessment: re-predict bird densities after backfilling
# author: Mannfred Boehm
# created: September 28, 2025
# ---

library(BAMexploreR)
library(parallel)
library(terra)
library(tidyverse)


# set paths ------------------------------------------------------
root <- "/home/mannfred/projects/def-ecknight/NationalModels"
ia_dir <- "/home/mannfred/scratch/impact_assessment"

bootstrap_dir     <- file.path(root, "output/06_bootstraps")
bart_backfill_dir <- file.path(ia_dir, "bart_models/2020")
importance_path   <- file.path(ia_dir, "bam_predictor_importance_v5.rds")
metrics_path      <- file.path(ia_dir, "continuous_holdout_metrics.rds")
industry_dir      <- file.path(ia_dir, "Rscripts/hirshpearson")


# import data ------------------------------------------------------

# import BCR boundaries
bam_boundary <- terra::vect(file.path(root, "Regions", "BAM_BCR_NationalModel_Unbuffered.shp")) 

# import merged + subsetted subbasins
all_subbasins_subset <- terra::vect(file.path(ia_dir, "hydrobasins_masked_merged_subset.gpkg"))

# import model performance metrics for continuous covariates
continuous_holdout_metrics <- readRDS("/home/mannfred/scratch/impact_assessment/Rscripts/continuous_holdout_metrics.rds")


# helper functions  ------------------------------------------------------

# create a reference table for which subbasins are in which BCRs 
# (some subbasins will be in multiple BCRs, and that's OK)
bcr_subbasins_ref <- {
  
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


# mosaic (backfilled) subbasin stacks into a single BCR-wide stack
mosaic_backfilled_stacks <- function(sub_ids) {
  
  paths <- file.path(
    bart_backfill_dir,
    paste0("subbasin_", sub_ids),
    paste0("subbasin_", sub_ids, "_backfill.tif")
  )
  
  paths <- paths[file.exists(paths)]
  if (length(paths) == 0) return(NULL)
  
  stacks <- lapply(paths, terra::rast)
  
  ref <- stacks[[1]]
  stacks <- lapply(stacks, \(r)
                   if (!terra::compareGeom(r, ref, stopOnError = FALSE))
                     terra::resample(r, ref, method = "near") else r
  )
  
  vars <- sort(unique(unlist(lapply(stacks, names))))
  out <- terra::rast(lapply(vars, function(v) {
    Reduce(terra::cover,
           lapply(stacks, \(r) if (v %in% names(r)) r[[v]] else NULL))
  }))
  
  names(out) <- vars
  out
}


# predictor importance bookkeeping  ------------------------------------------------------

# keep track of which covariates were 
# 1) requested by a bird model and 2)dropped when R^2 < 0.70
# later we’ll record ranked importance of dropped covariates
importance_tbl <- readRDS(importance_path)
metrics_tbl    <- readRDS(metrics_path)

good_backfill_vars <- metrics_tbl |>
  filter(r2 >= 0.7) |>
  pull(covariate) |>
  unique()

# import human footprint rasters ------------------------------------------------------

industry_rasters <- list.files(industry_dir, pattern = "\\.tif$", full.names = TRUE)

industry_names <- tools::file_path_sans_ext(basename(industry_rasters))
industry_stack <- setNames(lapply(industry_rasters, terra::rast),
                           industry_names)

# define disturbance variables, which we'll set to zero when re-predicting
disturbance_vars <-
  dplyr::tibble(BAMexploreR::predictor_metadata) |>
  dplyr::filter(version == "v5") |>
  dplyr::select(predictor, definition, predictor_class) |>
  dplyr::filter(predictor_class == "Disturbance")


# define density prediction function ------------------------------------------------------
predict_species_bcr <- function(species, year) {
  
  spp_dir <- file.path(bootstrap_dir, species)
  if (!dir.exists(spp_dir)) return(NULL)
  
  rdata_files <- list.files(spp_dir, pattern = "can.*\\.Rdata$", full.names = TRUE)
  
  do.call(bind_rows, lapply(rdata_files, function(rdata_path) {
    
    load(rdata_path)   # loads b.list
    mdl_list <- b.list
    bcr_code <- attr(mdl_list[[1]], "bcr")
    
    message(species, " ", bcr_code)
    
    # find subbasins in current BCR
    sub_ids <- bcr_subbasins_ref |>
      filter(bcr_code == !!bcr_code) |>
      pull(HYBAS_ID) |>
      unique()
    
    if (length(sub_ids) == 0) return(NULL)
    
    # get covariate stack for present BCR x year
    obs_path <- file.path(root, "gis/stacks", paste0(bcr_code, "_", year, ".tif"))
    if (!file.exists(obs_path)) return(NULL)
    stack_obs <- terra::rast(obs_path)
    
    # get mosaiced watershed stacks
    stack_bf <- mosaic_backfilled_stacks(sub_ids)
    if (is.null(stack_bf)) return(NULL)
    
    # get covariates from bootstrap models
    # NOTE: this is only considering boostrap 1. EXPAND TO ALL BOOTS
    varnames <- mdl_list[[1]]$var.names
    
    # get predictor names from observed landscape
    X_obs <- stack_obs[[intersect(varnames, names(stack_obs))]]
    
    # get backfilled predictor names with R^2 > 0.7
    cont_vars <- intersect(good_backfill_vars, varnames)
    cat_vars  <- setdiff(varnames, cont_vars)
    
    # populate bf_sets with mean and 95% CI of each continuous covariate
    bf_sets <- list(mean = X_obs, low  = X_obs, high = X_obs)
    
    for (v in cont_vars) {
      bf_sets$mean[[v]] <- stack_bf[[paste0(v, "_mean")]]
      bf_sets$low [[v]] <- stack_bf[[paste0(v, "_mean")]] -
        1.96 * stack_bf[[paste0(v, "_sd")]]
      bf_sets$high[[v]] <- stack_bf[[paste0(v, "_mean")]] +
        1.96 * stack_bf[[paste0(v, "_sd")]]
    }
    
    # IS THERE A PURPOSE OF SETting A LOW AND HIGH HERE???
    for (v in cat_vars) {
      if (v %in% names(stack_bf)) {
        bf_sets$mean[[v]] <- stack_bf[[v]]
        bf_sets$low [[v]] <- stack_bf[[v]]
        bf_sets$high[[v]] <- stack_bf[[v]]
      }
    }
    
    # define helper for running predict()
    pred_one <- function(m, X) {
      nt <- m$n.trees
      terra::predict(X, m, n.trees = nt, type = "response")
    }
    

    # run predictions for 32 bootstraps in `b.list`
    # observed landscape
    obs_preds <- lapply(mdl_list, pred_one, X = X_obs)
    
    # run predictions for 32 bootstraps x 3 landscape layers (mean, low, high)
    # backfilled landscape
    bf_preds <- lapply(bf_sets, function(X) lapply(mdl_list, pred_one, X = X))
    
  
    # count birds from prediction surfaces
    # WHY ALL_SUBBASINS_SUBSET???
    agg <- function(preds) {
      sapply(preds, function(r)
        terra::extract(r * 100, all_subbasins_subset,
                       fun = sum, na.rm = TRUE, ID = FALSE)[,1]
      )
    }
    
    obs_mat <- agg(obs_preds)
    bf_mat  <- lapply(bf_preds, agg)
    
    tibble(
      species = species,
      bcr     = bcr_code,
      treatment = "observed",
      mean = rowMeans(obs_mat, na.rm = TRUE),
      sd   = apply(obs_mat, 1, sd, na.rm = TRUE)
    ) |>
      bind_rows(
        lapply(names(bf_mat), function(k)
          tibble(
            species = species,
            bcr = bcr_code,
            treatment = paste0("backfilled_", k),
            mean = rowMeans(bf_mat[[k]], na.rm = TRUE),
            sd   = apply(bf_mat[[k]], 1, sd, na.rm = TRUE)
          )
        )
      )
  }))
}

# set up for parallel processing ------------------------------------------------------
cl <- makeCluster(max(1, detectCores() - 1))
clusterEvalQ(cl, {
  library(terra)
  library(tidyverse)
  library(BAMexploreR)
})

species_vec <- c("CAWA", "OSFL")


# run in parallel ------------------------------------------------------
res <- parLapply(cl, species_vec, predict_species_bcr, year = 2020)
stopCluster(cl)

res_all <- bind_rows(res)
