# ---
# title: Impact Assessment: canonical observed-landscape bird predictions
# author: Mannfred Boehm
# created: May 15, 2026
# rewritten: September 11, 2026 to conform to V5 10.Truncate.R (commit f082866)
# ---
#
# Run once per species on the local machine. For each Canadian BCR the species has
# a BRT model for, read Elly's unclamped 32-bootstrap prediction surfaces from
#   G:/Shared drives/BAM_NationalModels5/output/07_predictions/{species}/{species}_{bcr}_{year}.tif
# apply V5's truncation transform, and write to
# ia_dir/data/derived_data/predictions/{species}/{bcr_code}/{year}/:
#   observed_bootstraps.tif  (32 layers, densmax- AND q99.9-clamped, UNWEIGHTED)
#   observed_mean.tif
#   observed_sd.tif
# plus one ia_dir/data/derived_data/predictions/{species}/truncation_params.rds
# per species.
#
# After running locally, Globus-transfer observed_bootstraps.tif AND
# truncation_params.rds to the same relative paths on the cluster. 12B/12C read
# these rather than recomputing, which also eliminates floating-point drift
# across parallel SLURM jobs.
#
# WHAT CHANGED (2026-09-11) and why it matters
# --------------------------------------------
# This script used to replicate steps 5-9 of V5's 10.Package.R. That script no
# longer exists: commit f082866 (2026-06-04) split it into 10.Truncate.R +
# 11.Package.R, and 4e7fc83 rewrote how the truncation values are derived.
#
# 1. q.out's schema changed from (spp, thresh, countmax, off, q) to
#    (spp, thresh, countmax, densmax). There is no $q column any more. The old
#    $q was a count quantile divided by a single null-QPAD offset; densmax is
#    the 99th percentile of per-observation density using the per-observation
#    corrections table. CAWA 1.1850 -> 0.8334, OVEN 0.8460 -> 1.3941.
#
# 2. Truncation has TWO upper stages and this script only ever applied the
#    first. Measured on can10 2020: densmax removes 1.29% (CAWA) / 0.32% (OVEN)
#    of total abundance, but the secondary 99.9th-percentile cap removes a
#    further 12.19% / 4.13%. q99.9 ~ 0.133 for CAWA against densmax 0.833 -- it
#    binds 6.3x lower and does nearly all the work. Both are now applied to the
#    saved stack, making observed_bootstraps.tif our EPSG:5072 analogue of V5's
#    output/10_truncated product.
#
# 3. The denshthresh / l.out low-density zeroing (old step 9) was deleted from
#    V5 entirely. It is gone here too.
#
# WHY q99 IS PERSISTED RATHER THAN RECOMPUTED
# -------------------------------------------
# q99 is the only data-dependent parameter in the transform. 12C must clamp the
# BACKFILLED predictions at the SAME value, frozen from the observed landscape.
# Removing industry raises density, so a counterfactual's own 99.9th percentile
# sits higher and would be clamped less than the observed landscape it is being
# differenced against -- putting a component of the cap itself into the obs/bf
# contrast. truncation_params.rds carries it to 12C for exactly this reason.
# (Known residual bias, to be quantified by the uncapped sensitivity pass: a
# frozen ceiling is hit more often by counterfactuals, which under-estimates
# impact in the highest-density pixels.)
#
# WHY THIS STAYS IN EPSG:5072 AND IS NOT MASKED HERE
# --------------------------------------------------
# V5 reprojects to EPSG:3978 before truncating, but 10.Truncate.R:19 documents
# that as a legacy artifact ("Future versions will not require this step"). It
# costs ~2.3% of abundance (bilinear resampling does not conserve sums) and,
# because clamp() does not commute with projection, it would break the exact
# superset -> masked-rowsum decomposition 12C relies on. The parameter itself is
# near CRS-invariant -- CAWA can10 q99.9 = 0.132457 in 5072 vs 0.132881 in 3978,
# 0.32% apart -- so the 5072-derived cap is effectively V5's cap.
#
# Range/water/data-limit masking is deliberately NOT applied here: it lives in
# weight.tif (12A2) and 12C:220 multiplies it into BOTH the observed and the
# backfilled side, which is what keeps w*bf - w*obs = w*(bf - obs) exact.
# Baking it in here would double-weight the observed side.
#
# Gate G1 (Rscripts/misc/verify_v5_truncate_port.R) confirms the shared transform
# in 12A0_v5_truncate.R reproduces V5's released output/10_truncated product.
# ---

suppressPackageStartupMessages({
  library(terra)
  library(tidyverse)
})

# terra::project() returns different VALUES under memory pressure (a 28% swing in
# a derived q99.9 was measured between memfrac 0.60 and 0.01). This script does
# not project, but pin it anyway so nothing here is memory-sensitive.
terraOptions(memfrac = 0.60)

# paths -------------------------------------------------------------------

nm_root <- "/home/mannfred/projects/def-ecknight/NationalModels"

cc    <- FALSE
local <- TRUE

if (cc)            { ia_dir <- "/home/mannfred/scratch/impact_assessment" }
if (!cc && local)  { ia_dir <- getwd() }
if (!cc && !local) { ia_dir <- file.path("G:/Shared drives/BAM_NationalModels5", "data", "Extras",
                                          "sandbox_data", "impactassessment_sandbox") }

# when running locally, BRT bootstrap models and raw prediction tifs are on G:
if (!cc) { nm_root <- "G:/Shared drives/BAM_NationalModels5" }

source(file.path(ia_dir, "Rscripts", "12A0_v5_truncate.R"))

# prediction thresholds --------------------------------------------------
# restaged from G:/Shared drives/BAM_NationalModels5/data/ on 2026-09-11
# (the 2026-06-04 rebuild; see the schema note above).

load(file.path(ia_dir, "data", "raw_data", "SpeciesPredictionTruncationValues.Rdata"))
if (!"densmax" %in% names(q.out))
  stop("q.out has no `densmax` column — SpeciesPredictionTruncationValues.Rdata is the ",
       "pre-2026-06-04 version. Restage it from G:/Shared drives/BAM_NationalModels5/data/.")

# species from SLURM -----------------------------------------------------

species_vec <- c("CAWA", "OVEN")
# species_vec <- sort(list.dirs(file.path(nm_root, "output/06_bootstraps"), full.names = FALSE, recursive = FALSE))
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
species <- species_vec[task_id]
year    <- 2020

densmax <- q.out[q.out$spp == species, ]$densmax
if (length(densmax) != 1 || !is.finite(densmax))
  stop("no usable densmax for species=", species)

message(Sys.time(), " | observed predictions for species=", species,
        " | densmax=", signif(densmax, 6))

# find models for relevant BCRs ------------------------------------------
# NOTE: this discovers every can* model, which for CAWA includes can40 — a model
# BAM withholds (review/ModelReleaseDecisions.xlsx, "remove" tab, AUC). It is
# kept here deliberately; the release filter belongs in 14B, not in the products.

rdata_files <- list.files(file.path(nm_root, "output/06_bootstraps", species),
                          pattern = "can.*\\.Rdata$", full.names = TRUE)

# same smoke-test filter 12B/12C honour: TEST_BCR=can10,can71
test_bcr <- Sys.getenv("TEST_BCR", "")
if (nchar(test_bcr) > 0) {
  want <- strsplit(test_bcr, ",")[[1]]
  rdata_files <- rdata_files[
    regmatches(basename(rdata_files), regexpr("can[0-9]+", basename(rdata_files))) %in% want]
  message(Sys.time(), " | TEST_BCR=", test_bcr, " -> ", length(rdata_files), " model(s)")
}
message(Sys.time(), " | found ", length(rdata_files), " BCR models")

params_path <- file.path(ia_dir, "data/derived_data/predictions", species,
                         "truncation_params.rds")
params <- if (file.exists(params_path)) readRDS(params_path) else list()

# make bird density predictions for applicable BCRs ----------------------

for (rdata_path in rdata_files) {

  e <- new.env(parent = emptyenv())
  load(rdata_path, envir = e)
  if (!exists("b.list", envir = e)) {
    message("  b.list not found in ", basename(rdata_path), " — skipping")
    next
  }
  bcr_code <- attr(e$b.list[[1]], "bcr")
  rm(e)

  obs_dir       <- file.path(ia_dir, "data/derived_data/predictions", species, bcr_code, year)
  obs_boot_path <- file.path(obs_dir, "observed_bootstraps.tif")

  # skip only if BOTH the stack and its frozen q99 are already on hand — a stack
  # without a recorded q99 is unusable by 12C and must be rebuilt.
  key <- paste(bcr_code, year, sep = "_")
  if (file.exists(obs_boot_path) && !is.null(params[[key]])) {
    message(Sys.time(), " | ", bcr_code, " | already exists — skipping")
    next
  }

  raw_pred_path <- file.path(nm_root, "output/07_predictions", species,
                             paste0(species, "_", bcr_code, "_", year, ".tif"))
  if (!file.exists(raw_pred_path)) {
    message(Sys.time(), " | ", bcr_code, " | no raw predictions found — skipping")
    next
  }
  message(Sys.time(), " | ", bcr_code, " | reading pre-computed bootstraps")

  # V5 10.Truncate.R steps 5-6, in native 5072, unweighted (see header).
  out <- v5_truncate(terra::rast(raw_pred_path),
                     spp = species, bcr = bcr_code, densmax = densmax,
                     q99 = NULL, project_to = NULL, apply_masks = FALSE)
  if (is.null(out)) {
    message(Sys.time(), " | ", bcr_code, " | all-NA mean, no q99 derivable — skipping")
    next
  }
  obs_stack <- out$stack

  # INTERLEAVE=BAND matters: 12C:215 pulls these 32 layers out one at a time, and
  # terra's default BIP layout makes each single-band read scan the whole file.
  dir.create(obs_dir, recursive = TRUE, showWarnings = FALSE)
  terra::writeRaster(obs_stack, obs_boot_path, overwrite = TRUE,
                     wopt = list(gdal = c("INTERLEAVE=BAND", "COMPRESS=DEFLATE",
                                          "BIGTIFF=YES", "TILED=NO")))

  # inspection products, matching 11.Package.R's mean/sd over the truncated stack
  # (V5 computes these on the MASKED stack; ours are unmasked, so they are for
  # inspection only — 12C reads observed_bootstraps.tif, never these).
  terra::writeRaster(terra::app(obs_stack, mean, na.rm = TRUE),
                     file.path(obs_dir, "observed_mean.tif"), overwrite = TRUE)
  terra::writeRaster(terra::app(obs_stack, sd, na.rm = TRUE),
                     file.path(obs_dir, "observed_sd.tif"), overwrite = TRUE)

  params[[key]] <- list(species = species, bcr = bcr_code, year = year,
                        densmax = out$densmax, q99 = out$q99,
                        crs = "EPSG:5072", masked = FALSE,
                        source = raw_pred_path, written = Sys.time())
  saveRDS(params, params_path)

  message(Sys.time(), " | ", bcr_code, " | done | q99=", signif(out$q99, 6),
          " | densmax/q99=", round(densmax / out$q99, 1), "x")
  rm(out, obs_stack); gc()
}

message(Sys.time(), " | species=", species, " complete | ", length(params),
        " BCR(s) in ", params_path)
