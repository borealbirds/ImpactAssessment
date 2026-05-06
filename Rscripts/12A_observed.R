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

# ---- Prediction thresholds --------------------------------------------------

load(file.path(ia_dir, "data", "raw_data", "SpeciesPredictionTruncationValues.Rdata"))

# ---- Species from SLURM -----------------------------------------------------

species_vec <- c("CAWA", "OVEN")
# species_vec <- sort(c("BANS", "BARS", "BOBO", "CAWA", "EAWP", "EVGR", "GCTH", "GRSP", "GWWA", "LEYE", "OSFL"))
# species_vec <- sort(list.dirs(file.path(nm_root, "output/06_bootstraps"), full.names = FALSE, recursive = FALSE))
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
species <- species_vec[task_id]
year    <- 2020
qsp <- q.out[q.out$spp == species, ]$q
q0  <- l.out[l.out$spp == species, ]$denshthresh

message(Sys.time(), " | observed predictions for species=", species)

# ---- Find BCR models ---------------------------------------------------------

rdata_files <- list.files(file.path(nm_root, "output/06_bootstraps", species),
                          pattern = "can.*\\.Rdata$", full.names = TRUE)
message(Sys.time(), " | found ", length(rdata_files), " BCR models")

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

  # read colleague's pre-computed unclamped bootstrap predictions (nm_root/output/07_predictions)
  raw_pred_path <- file.path(nm_root, "output/07_predictions", species,
                             paste0(species, "_", bcr_code, "_", year, ".tif"))
  if (!file.exists(raw_pred_path)) {
    message(Sys.time(), " | ", bcr_code, " | no raw predictions found — skipping")
    next
  }
  message(Sys.time(), " | ", bcr_code, " | reading pre-computed bootstraps")
  raw_stack <- terra::rast(raw_pred_path)
  obs_preds <- lapply(seq_len(terra::nlyr(raw_stack)), function(i) raw_stack[[i]])

  # Step 5 of 10.Package.R: clamp each bootstrap at the species-specific quantile
  obs_preds <- lapply(obs_preds, function(r) terra::clamp(r, upper = qsp))

  # save bootstrap stack (qsp-clamped; read by 12B/12C)
  dir.create(obs_dir, recursive = TRUE, showWarnings = FALSE)
  terra::writeRaster(rast(obs_preds), obs_boot_path, overwrite = TRUE)

  # Steps 6-7: mean, then secondary cap at 99.9th percentile of mean
  obs_stack <- terra::rast(obs_preds)
  mn_r  <- terra::app(obs_stack, mean, na.rm = TRUE)
  q99_r <- terra::global(mn_r, quantile, probs = 0.999, na.rm = TRUE)[1, 1]
  mn2_r <- terra::clamp(mn_r, upper = q99_r)

  # Step 8: SD from bootstraps also clamped at q99
  sd_r  <- terra::app(terra::clamp(obs_stack, upper = q99_r), sd, na.rm = TRUE)

  # Step 9: zero pixels below denshthresh
  mean_out <- terra::ifel(mn2_r < q0, 0, mn2_r)
  sd_out   <- terra::ifel(mn2_r < q0, 0, sd_r)

  terra::writeRaster(mean_out, file.path(obs_dir, "observed_mean.tif"), overwrite = TRUE)
  terra::writeRaster(sd_out,   file.path(obs_dir, "observed_sd.tif"),   overwrite = TRUE)

  message(Sys.time(), " | ", bcr_code, " | done")
  rm(b.list, raw_stack, obs_preds, obs_stack, mn_r, q99_r, mn2_r, sd_r, mean_out, sd_out, e)
  gc()
}

message(Sys.time(), " | species=", species, " complete.")
