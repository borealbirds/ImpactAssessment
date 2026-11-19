# ---
# title: Impact Assessment: canonical observed-landscape bird predictions
# author: Mannfred Boehm
# created: May 15, 2026 (rewritten September 11, 2026 to follow V5's 10.Truncate.R, commit f082866)
# ---

# per species, each Canadian BCR it has a model for: Elly's unclamped 32-bootstrap predictions
# ({nm_root}/output/07_predictions/{species}/{species}_{bcr}_{year}.tif) through V5's two caps (12B)
#   -> data/derived_data/predictions/{species}/{bcr}/{year}/observed_{bootstraps,mean,sd}.tif
#   -> data/derived_data/predictions/{species}/truncation_params.rds (densmax, and q99 frozen for 12F)
# unmasked: 12C's weight masks both sides in 12F, which keeps w*bf - w*obs = w*(bf - obs). in
# EPSG:5072, not V5's legacy 3978: clamp() does not commute with projection (OL #8).
# models BAM withholds (CAWA can40) are run too; 14B drops them.
# on Fir: bash 12A_observed.sh OSFL GRSP   (one task per species runs all YEARS, which share
# one truncation_params.rds). locally, from G: SPECIES=OSFL Rscript Rscripts/12A_observed.R

library(terra)
terraOptions(memfrac = 0.60)   # project() returns different values under memory pressure

cc      <- nzchar(Sys.getenv("SLURM_JOB_ID"))
ia_dir  <- if (cc) "/home/mannfred/scratch/impact_assessment" else getwd()
nm_root <- if (cc) "/home/mannfred/projects/def-ecknight/NationalModels" else "G:/Shared drives/BAM_NationalModels5"

# terra sizes its blocks from the node's free RAM, not the job's allocation (see 11)
slurm_mem_mb <- suppressWarnings(as.numeric(Sys.getenv("SLURM_MEM_PER_NODE")))
if (is.finite(slurm_mem_mb)) terraOptions(memmax = 0.3 * slurm_mem_mb / 1024)

source(file.path(ia_dir, "Rscripts", "12B_v5_truncate.R"))

load(file.path(ia_dir, "data", "raw_data", "SpeciesPredictionTruncationValues.Rdata"))   # q.out
if (!"densmax" %in% names(q.out))
  stop("q.out has no densmax column: restage SpeciesPredictionTruncationValues.Rdata (2026-06-04 or later)")

# one species per array task (the .sh sets SPECIES); a local run takes exactly one
species_vec <- strsplit(Sys.getenv("SPECIES"), "[, ]+")[[1]]
if (!cc && length(species_vec) != 1L) stop("a local run takes one species: SPECIES=OSFL Rscript Rscripts/12A_observed.R")
species <- species_vec[as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))]
years   <- as.integer(strsplit(Sys.getenv("YEARS", "2020"), ",")[[1]])

densmax <- q.out[q.out$spp == species, ]$densmax
if (length(densmax) != 1 || !is.finite(densmax)) stop("no usable densmax for ", species)
message(Sys.time(), " | ", species, " | years ", paste(years, collapse = ","), " | densmax=", signif(densmax, 6))

rdata_files <- list.files(file.path(nm_root, "output/06_bootstraps", species), pattern = "can.*\\.Rdata$", full.names = TRUE)

# TEST_BCR=can10,can71: the smoke-test filter 12D and 12H honour too
test_bcr <- Sys.getenv("TEST_BCR", "")
if (nchar(test_bcr) > 0)
  rdata_files <- rdata_files[regmatches(basename(rdata_files), regexpr("can[0-9]+", basename(rdata_files))) %in%
                               strsplit(test_bcr, ",")[[1]]]
message(Sys.time(), " | ", length(rdata_files), " BCR models")

params_path <- file.path(ia_dir, "data/derived_data/predictions", species, "truncation_params.rds")
params <- if (file.exists(params_path)) readRDS(params_path) else list()



# ---------------------------------------------------
# truncate each BCR x year

for (rdata_path in rdata_files) for (year in years) {

  # the BCR from the file name, {species}_{bcr}.Rdata (loading the model to read it costs up to
  # 654 MB; 12F checks the model's own BCR against the name)
  bcr_code <- sub(paste0("^", species, "_(.*)\\.Rdata$"), "\\1", basename(rdata_path))
  if (!grepl("^can[0-9]+$", bcr_code)) stop("cannot read a BCR code from ", basename(rdata_path))
  key <- paste(bcr_code, year, sep = "_")

  obs_dir       <- file.path(ia_dir, "data/derived_data/predictions", species, bcr_code, year)
  obs_boot_path <- file.path(obs_dir, "observed_bootstraps.tif")

  # skip only if both the stack and its frozen q99 exist: 12F can't use one without the other
  if (file.exists(obs_boot_path) && !is.null(params[[key]])) {
    message(Sys.time(), " | ", key, " | done before - skipping")
    next
  }

  raw_pred_path <- file.path(nm_root, "output/07_predictions", species, paste0(species, "_", bcr_code, "_", year, ".tif"))
  if (!file.exists(raw_pred_path)) {
    message(Sys.time(), " | ", key, " | no raw predictions - skipping")
    next
  }

  # V5 10.Truncate.R steps 5-6: densmax, then q99 from this stack's bootstrap mean
  out <- v5_truncate(terra::rast(raw_pred_path), spp = species, bcr = bcr_code, densmax = densmax,
                     q99 = NULL, project_to = NULL, apply_masks = FALSE)
  if (is.null(out)) {
    message(Sys.time(), " | ", key, " | all-NA mean, no q99 - skipping")
    next
  }

  # INTERLEAVE=BAND: 12F reads the 32 layers one at a time, and terra's default (pixel
  # interleave) makes each single-band read scan the whole file
  dir.create(obs_dir, recursive = TRUE, showWarnings = FALSE)
  terra::writeRaster(out$stack, obs_boot_path, overwrite = TRUE,
                     wopt = list(gdal = c("INTERLEAVE=BAND", "COMPRESS=DEFLATE", "BIGTIFF=YES", "TILED=NO")))

  # for inspection only (12F reads observed_bootstraps.tif)
  terra::writeRaster(terra::app(out$stack, mean, na.rm = TRUE), file.path(obs_dir, "observed_mean.tif"), overwrite = TRUE)
  terra::writeRaster(terra::app(out$stack, sd, na.rm = TRUE), file.path(obs_dir, "observed_sd.tif"), overwrite = TRUE)

  params[[key]] <- list(species = species, bcr = bcr_code, year = year, densmax = out$densmax, q99 = out$q99,
                        crs = "EPSG:5072", masked = FALSE, source = raw_pred_path, written = Sys.time())
  saveRDS(params, params_path)

  message(Sys.time(), " | ", key, " | q99=", signif(out$q99, 6), " (", round(densmax / out$q99, 1), "x below densmax)")
  rm(out)
  gc()
}
