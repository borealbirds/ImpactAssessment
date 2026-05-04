# ---
# title: Impact Assessment: canonical observed-landscape bird predictions
# author: Mannfred Boehm
# ---
# Run ONCE per species (SLURM array over species) BEFORE any coalition runs.
# For each BCR that the species has a BRT model for, compute observed bootstrap
# predictions and save:
#   predictions/{species}/{bcr_code}/{year}/observed_bootstraps.tif  (32 layers)
#   predictions/{species}/{bcr_code}/{year}/observed_mean.tif
#   predictions/{species}/{bcr_code}/{year}/observed_sd.tif
#
# Coalition runs (12A/12B) then read these files rather than recomputing them,
# eliminating floating-point drift across parallel SLURM jobs.
# ---

suppressPackageStartupMessages({
  library(BAMexploreR)
  library(gbm)
  library(terra)
  library(tidyverse)
})

# ---- Paths -------------------------------------------------------------------

nm_root <- "/home/mannfred/projects/def-ecknight/NationalModels"

cc    <- TRUE
local <- FALSE

if (cc)            { ia_dir <- "/home/mannfred/scratch/impact_assessment" }
if (!cc && local)  { ia_dir <- getwd() }
if (!cc && !local) { ia_dir <- file.path("G:/Shared drives/BAM_NationalModels5", "data", "Extras",
                                          "sandbox_data", "impactassessment_sandbox") }

# ---- Species from SLURM -----------------------------------------------------

species_vec <- c("CAWA")
# species_vec <- sort(c("BANS", "BARS", "BOBO", "CAWA", "EAWP", "EVGR", "GCTH", "GRSP", "GWWA", "LEYE", "OSFL"))
# species_vec <- sort(list.dirs(file.path(nm_root, "output/06_bootstraps"), full.names = FALSE, recursive = FALSE))
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
species <- species_vec[task_id]
year    <- 2020

message(Sys.time(), " | observed predictions for species=", species)

# ---- Find BCR models ---------------------------------------------------------

rdata_files <- list.files(file.path(nm_root, "output/06_bootstraps", species),
                          pattern = "can.*\\.Rdata$", full.names = TRUE)
message(Sys.time(), " | found ", length(rdata_files), " BCR models")

# ---- Helper: incremental mean/sd --------------------------------------------

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
  list(mean = mean_r, sd = sqrt(m2_r / (n - 1)))
}

# ---- Loop over BCRs ---------------------------------------------------------

for (rdata_path in rdata_files) {

  e <- new.env(parent = emptyenv())
  load(rdata_path, envir = e)
  if (!exists("b.list", envir = e)) {
    message("  b.list not found in ", basename(rdata_path), " — skipping")
    next
  }
  b.list   <- e$b.list
  bcr_code <- attr(b.list[[1]], "bcr")

  obs_dir <- file.path(ia_dir, "data/derived_data/predictions", species, bcr_code, year)
  obs_boot_path <- file.path(obs_dir, "observed_bootstraps.tif")

  # skip if already computed
  if (file.exists(obs_boot_path)) {
    message(Sys.time(), " | ", bcr_code, " | already exists — skipping")
    next
  }

  message(Sys.time(), " | ", bcr_code, " | predicting observed landscape")

  # load observed covariate stack
  stack_obs <- rast(file.path(nm_root, "gis/stacks", paste0(bcr_code, "_", year, ".tif")))
  model_vars <- b.list[[1]]$var.names
  X_obs <- stack_obs[[intersect(model_vars, names(stack_obs))]]

  # predict for each bootstrap
  obs_preds <- vector("list", length(b.list))
  for (i in seq_along(b.list)) {
    model <- b.list[[i]]
    obs_preds[[i]] <- terra::predict(X_obs, model, type = "response", n.trees = model$n.trees)
    message(Sys.time(), " | ", bcr_code, " | bootstrap ", i, "/", length(b.list))
    gc()
  }

  # Cap per-pixel BRT predictions at the maximum of the packaged mean raster.
  # Follows LandbirdModelsV5/analysis/12.Summarize.R, which uses the max of
  # 10.Package.R's cleaned mean layer as a species- and BCR-specific upper bound.
  pkg_path <- file.path(nm_root, "output", "10_packaged", species, bcr_code,
                        paste0(species, "_", bcr_code, "_", year, ".tif"))
  q99 <- if (file.exists(pkg_path)) {
    terra::global(terra::rast(pkg_path)[["mean"]], max, na.rm = TRUE)[, 1]
  } else {
    message("  WARNING: packaged raster not found for ", species, " ", bcr_code,
            " ", year, " — skipping prediction cap")
    Inf
  }
  obs_preds <- lapply(obs_preds, function(r) terra::clamp(r, upper = q99))

  # save full bootstrap stack + summary rasters
  dir.create(obs_dir, recursive = TRUE, showWarnings = FALSE)
  terra::writeRaster(rast(obs_preds), obs_boot_path, overwrite = TRUE)

  obs_summary <- summarize_preds(obs_preds)
  terra::writeRaster(obs_summary$mean, file.path(obs_dir, "observed_mean.tif"), overwrite = TRUE)
  terra::writeRaster(obs_summary$sd,   file.path(obs_dir, "observed_sd.tif"),   overwrite = TRUE)

  message(Sys.time(), " | ", bcr_code, " | done")
  rm(b.list, stack_obs, X_obs, obs_preds, obs_summary, e)
  gc()
}

message(Sys.time(), " | species=", species, " complete.")
