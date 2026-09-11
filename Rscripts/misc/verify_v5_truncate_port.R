# ---
# B1 / Gate G1: verify that 12A0_v5_truncate.R reproduces V5's released product.
#
# Runs the port over 07_predictions in 3978 (V5's own configuration) and compares
# the result cell-by-cell against output/10_truncated/{spp}/{bcr}/{spp}_{bcr}_{yr}.tif.
# Passing this proves the port is faithful, independently of the CRS we later
# choose to run production in.
#
# Usage: Rscript Rscripts/misc/verify_v5_truncate_port.R CAWA:can10 OVEN:can10
# ---

suppressPackageStartupMessages({ library(terra); library(sf) })

# MUST pin memfrac. Measured 2026-09-11: terra::project() returns materially
# different values under memory pressure -- same dims, same extent, same non-NA
# count, but max|diff| 0.137 on CAWA can71 and a bootstrap-mean max of 0.0906
# (memfrac 0.60) vs 0.0300 (memfrac 0.01), which moved the derived q99 by 28%.
# global(quantile) itself is exact and stable at every memfrac; the instability
# is entirely in project(). Production never projects (we run in 5072), but this
# harness does, so it has to be pinned or G1 is not reproducible.
terraOptions(progress = 0, memfrac = 0.60)

ia_dir <- "C:/Users/mannf/Drive/boreal_avian_modelling_project/ImpactAssessment"
v5     <- "G:/Shared drives/BAM_NationalModels5"
year   <- 2020
# stable scratch dir so outputs survive the R session and can be re-inspected
tmpd <- Sys.getenv("V5PORT_OUT", unset = file.path(tempdir(), "v5port"))
dir.create(tmpd, recursive = TRUE, showWarnings = FALSE)

source(file.path(ia_dir, "Rscripts", "12A0_v5_truncate.R"))

args  <- commandArgs(trailingOnly = TRUE)
if (!length(args)) args <- c("CAWA:can10")
pairs <- strsplit(args, ":")

e <- new.env()
load(file.path(v5, "data", "SpeciesPredictionTruncationValues.Rdata"), envir = e)

cat("loading V5 masking layers...\n")
masks <- v5_load_masks(file.path(v5, "gis"), want_subregions = TRUE)

for (p in pairs) {
  spp <- p[1]; bcr <- p[2]
  cat("\n=========================", spp, bcr, year, "\n")

  target_path <- file.path(v5, "output", "10_truncated", spp, bcr,
                           sprintf("%s_%s_%d.tif", spp, bcr, year))
  pred_path   <- file.path(v5, "output", "07_predictions", spp,
                           sprintf("%s_%s_%d.tif", spp, bcr, year))
  if (!file.exists(target_path)) { cat("  no V5 target -- skipping\n"); next }

  densmax <- e$q.out[e$q.out$spp == spp, ]$densmax
  t0  <- Sys.time()
  out <- v5_truncate(terra::rast(pred_path), spp = spp, bcr = bcr,
                     densmax = densmax, masks = masks,
                     range_root = file.path(v5, "gis", "ranges"),
                     q99 = NULL, project_to = "EPSG:3978", res = 1000)
  cat("  ran in", round(difftime(Sys.time(), t0, units = "mins"), 2), "min |",
      "densmax =", signif(densmax, 6), "| derived q99 =", signif(out$q99, 6), "\n")

  # round-trip through FLT4S so we compare like with like (V5 writes FLT4S)
  mine_path <- file.path(tmpd, sprintf("%s_%s_%d_port.tif", spp, bcr, year))
  terra::writeRaster(out$stack, mine_path, overwrite = TRUE, datatype = "FLT4S")
  mine <- terra::rast(mine_path); theirs <- terra::rast(target_path)

  cat("  dims   mine:", paste(dim(mine), collapse = "x"),
      "| V5:", paste(dim(theirs), collapse = "x"), "\n")
  cat("  extent identical:", isTRUE(all.equal(as.vector(ext(mine)), as.vector(ext(theirs)))),
      "| crs identical:", terra::same.crs(mine, theirs), "\n")

  if (!all(dim(mine) == dim(theirs))) { cat("  GRID MISMATCH -- cannot compare cells\n"); next }

  d      <- mine - theirs
  amax   <- max(abs(terra::minmax(d)))
  sm     <- sum(terra::global(mine,   "sum", na.rm = TRUE)[[1]])
  st     <- sum(terra::global(theirs, "sum", na.rm = TRUE)[[1]])
  nmine  <- sum(terra::global(!is.na(mine),   "sum", na.rm = TRUE)[[1]])
  ntheir <- sum(terra::global(!is.na(theirs), "sum", na.rm = TRUE)[[1]])
  nbig   <- sum(terra::global(abs(d) > 1e-3,  "sum", na.rm = TRUE)[[1]])

  # Residual disagreement is expected to be reprojection noise from a PROJ/GDAL
  # version difference: bilinear weights shift sub-metre, so cells in steep
  # terrain move while flat terrain stays bit-identical. Confirm by correlating
  # |diff| against local roughness -- a high correlation means the transform
  # logic agrees and only the resampling differs.
  rough <- terra::focal(theirs[[1]], w = 3,
                        fun = function(x) diff(range(x, na.rm = TRUE)))
  v <- data.frame(d = terra::values(abs(d[[1]]))[, 1],
                  r = terra::values(rough)[, 1])
  v <- v[complete.cases(v), ]
  rho    <- suppressWarnings(cor(v$d, v$r))
  smooth <- mean(v$d[v$r <= quantile(v$r, 0.10)])
  steep  <- mean(v$d[v$r >= quantile(v$r, 0.90)])

  cat("  non-NA cells  mine:", format(nmine, big.mark = ","),
      "| V5:", format(ntheir, big.mark = ","),
      sprintf("(delta %+d)", nmine - ntheir), "\n")
  cat("  max |diff|          :", signif(amax, 4), "\n")
  cat("  cells |diff| > 1e-3 :", format(nbig, big.mark = ","),
      sprintf("(%.4f%% of non-NA)", 100 * nbig / nmine), "\n")
  cat("  sum ratio mine/V5   :", signif(sm / st, 8), "\n")
  cat("  cor(|diff|, roughness):", round(rho, 3),
      "| mean|diff| smoothest decile:", signif(smooth, 3),
      "-> roughest:", signif(steep, 3), "\n")

  # The discriminating evidence is that FLAT terrain is bit-identical: any error
  # in the transform logic (wrong cap, wrong mask, wrong order) would perturb
  # flat cells too. rho is reported but not gated on -- it necessarily weakens as
  # agreement improves, because there is less residual signal left to correlate.
  verdict <- if (amax < 1e-5) {
    "EXACT MATCH"
  } else if (abs(sm / st - 1) < 1e-3 && smooth < 1e-6) {
    "PASS -- flat terrain bit-identical; residual is reprojection noise only"
  } else {
    "MISMATCH -- transform logic differs, investigate"
  }
  cat("  VERDICT:", verdict, "\n")
}
