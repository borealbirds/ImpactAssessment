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


# import data ------------------------------------------------------

# import BCR boundaries
bam_boundary <- terra::vect(file.path(root, "Regions", "BAM_BCR_NationalModel_Unbuffered.shp")) 

# import merged + subsetted subbasins
all_subbasins_subset <- terra::vect(file.path(ia_dir, "hydrobasins_masked_merged_subset.gpkg"))


# ------------------------------------------------------
# create a reference table for which subbasins are in which BCRs 
# some subbasins will be in multiple BCRs, and that's OK

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

# ------------------------------------------------------
# define helper function:
# mosaic (backfilled) subbasin stacks into a single BCR-wide stack

mosaic_backfilled_stacks <- function(sub_ids) {
  
  paths <- file.path(
    file.path(ia_dir, "bart_models/2020"),
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
bam_predictor_importance_v5 <- readRDS(file.path(ia_dir, "bam_predictor_importance_v5.rds"))

# import model performance metrics for continuous covariates
continuous_holdout_metrics <- readRDS(file.path(ia_dir, "continuous_holdout_metrics.rds"))

# find covariates with R^2 >= 0.7
good_backfill_vars <- 
  continuous_holdout_metrics |>
  filter(r2 >= 0.7) |>
  pull(covariate) |>
  unique()

# define categorical covariates
categorical_responses = c("ABoVE_1km", "NLCD_1km","MODISLCC_1km", "MODISLCC_5x5","SCANFI_1km","VLCE_1km")

# import human footprint rasters ------------------------------------------------------

industry_rasters <- list.files(file.path(ia_dir, "Rscripts/hirshpearson"), pattern = "\\.tif$", full.names = TRUE)
industry_names <- sub(".tif", "", basename(industry_rasters))
industry_stack <- setNames(lapply(industry_rasters, terra::rast), industry_names)

# define disturbance variables, which we'll set to zero when re-predicting
disturbance_vars <-
  dplyr::tibble(BAMexploreR::predictor_metadata) |>
  dplyr::filter(version == "v5") |>
  dplyr::select(predictor, definition, predictor_class) |>
  dplyr::filter(predictor_class == "Disturbance")


# define density prediction function ------------------------------------------------------
predict_species_bcr <- function(species, year) {
  
  # find all models for all BCRs for spp_i 
  rdata_files <- list.files(file.path(root, "output/06_bootstraps/", species), 
                            pattern = "can.*\\.Rdata$", 
                            full.names = TRUE)
  
  # this sub-function will work through all BCRs for a given species x year
  do.call(dplyr::bind_rows, lapply(rdata_files, function(rdata_path) {
    
    load(rdata_path) # loads b.list for some spp x bcr
    bcr_code <- attr(b.list[[1]], "bcr") # get current bcr
    message(species, " ", bcr_code)
    
    # find subbasins in current BCR
    sub_ids <- 
      bcr_subbasins_ref |>
      dplyr::filter(bcr_code == !!bcr_code) |>
      dplyr::pull(HYBAS_ID) |>
      unique()
    
    # get observed environmental stack for current BCR x year
    stack_obs <- terra::rast(file.path(root, "gis/stacks", paste0(bcr_code, "_", year, ".tif")))
    
    # get mosaiced watershed stacks for current BCR x year
    stack_bf <- mosaic_backfilled_stacks(sub_ids)
    
    # get all covariates used by bootstrap models for this spp x BCR
    varnames <- unique(unlist(lapply(b.list, `[[`, "var.names")))
    
    # get predictor names from observed landscape
    X_obs <- stack_obs[[intersect(varnames, names(stack_obs))]]
    
    # get backfilled predictor names with R^2 > 0.7
    cont_vars <- intersect(good_backfill_vars, varnames)
    cat_vars  <- categorical_responses
    
    # create bf_sets to hold mean and 95% CI of each continuous covariate
    bf_sets <- list(mean = X_obs, low  = X_obs, high = X_obs)
    
    for (v in cont_vars) {
      bf_sets$mean[[v]] <- stack_bf[[paste0(v, "_mean")]]
      bf_sets$low [[v]] <- stack_bf[[paste0(v, "_mean")]] -
        1.96 * stack_bf[[paste0(v, "_sd")]]
      bf_sets$high[[v]] <- stack_bf[[paste0(v, "_mean")]] +
        1.96 * stack_bf[[paste0(v, "_sd")]]
    }
    
    # no uncertainty is propagated for cat_vars, here we are simply 
    # creating the same slots as cont_vars for downstream cohesion
    for (v in cat_vars) {
      if (v %in% names(stack_bf)) {
        bf_sets$mean[[v]] <- stack_bf[[v]]
        bf_sets$low [[v]] <- stack_bf[[v]]
        bf_sets$high[[v]] <- stack_bf[[v]]
      }
    }
    
    # run predictions for 32 bootstraps in `b.list`
    # observed landscape
    obs_preds <- lapply(b.list, predict(X, m, n.trees = m$n.trees, type = "response"), X = X_obs)
    
    # run predictions for 32 bootstraps x 3 landscape layers (mean, low, high)
    # backfilled landscape
    bf_preds <- lapply(bf_sets, function(X) lapply(b.list, predict(X, m, n.trees = m$n.trees, type = "response"), X = X))
    
    # count birds from prediction surfaces intersecting our watersheds
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
