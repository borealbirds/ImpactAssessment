# ---
# title: Gate G3/G4 -- does weight.tif reproduce V5's range/water/data-limit/BCR masking?
# author: Mannfred Boehm
# created: September 11, 2026
# ---
#
# Runs LOCALLY (needs G: for V5's own mask vectors). Rebuilds 12A2's weight on this
# machine and compares, per species x BCR:
#
#   A  our production path : sum(observed_bootstraps.tif * weight.tif * 100)
#   B  V5's own masking    : 10.Truncate.R steps 7-8 applied as vectors, same CRS
#   C  V5's own masking    : same, after V5's legacy projection to EPSG:3978
#
# A vs B isolates raster-weight-vs-vector masking; B vs C isolates the 5072-vs-3978
# crosswalk we deliberately skip (see 12A0_v5_truncate.R). C is the number that should
# match output/13_summary/BAMV5-results.xlsx.
#
# RESULT 2026-09-11 (CAWA can10 / CAWA can71 / OVEN can10):
#   A/B = 1.000011, 1.000002, 1.000001   -- float32 noise; masking parity is exact
#   B/C = 1.029,    1.029,    1.034      -- the accepted projection cost
#   C    = 19,389 / 161,922 / 54,062     -- matches published 0.019 / 0.162 / 0.054 M
#
# Two defects were found and fixed by this harness; both are guarded by the
# WEIGHT_VERSION stamp in 12A2, so re-running it after a change is the way to re-gate:
#   1. the missing crop to the BCR's own polygon (V5 grids are buffered far past the
#      subunit) -- CLAUDE.md Open Limitation #7;
#   2. rasterize(touches = FALSE), terra's default, vs the touches = TRUE behaviour of
#      V5's mask()/crop(). Worth 4.9-7.7% on the 29,545-polygon water mask.
# ---

suppressPackageStartupMessages({ library(terra) })
terraOptions(memfrac = 0.6)

ia_dir  <- getwd()
nm_root <- "G:/Shared drives/BAM_NationalModels5"
gis_dir <- file.path(ia_dir, "data", "raw_data", "v5_gis")
year    <- 2020

pairs <- list(c("CAWA", "can10"), c("CAWA", "can71"), c("OVEN", "can10"))

water_v <- vect(file.path(nm_root, "gis", "WaterMask_Canada.shp"))
limit_v <- vect(file.path(nm_root, "gis", "DataLimitationsMask.shp"))
sub_v   <- vect(file.path(nm_root, "gis", "Subregions_Mosaics_EPSG3978.shp"))

pop <- function(x) sapply(seq_len(nlyr(x)), function(i)
  global(x[[i]] * 100, "sum", na.rm = TRUE)[[1]])   # birds/ha -> birds/km2
fmt <- function(v) sprintf("%.0f [%.0f, %.0f]", median(v),
                           quantile(v, 0.05), quantile(v, 0.95))

# 12A2's weight, rebuilt here so the gate tests the definition rather than a stale file
build_weight <- function(spp, tmpl) {
  tc  <- crs(tmpl)
  ctg <- function(v) {
    tp <- as.polygons(ext(tmpl), crs = tc)
    project(crop(v, project(tp, crs(v))), tc)
  }
  rr <- rast(file.path(gis_dir, "ranges", paste0(spp, ".tif")))
  list(ctg = ctg,
       range = classify(project(rr, tmpl, method = "bilinear"), cbind(NA, 0)))
}

# V5 10.Truncate.R steps 7-8, in whatever CRS `r` is already in
v5_mask <- function(r, spp, bcr) {
  tc <- crs(r)
  rr <- rast(file.path(nm_root, "gis", "ranges", paste0(spp, ".tif")))
  rr <- if (crs(rr) == tc) resample(rr, r) else project(rr, r, method = "bilinear")
  m <- r * rr
  m[is.na(m)] <- 0
  m <- crop(m, project(sub_v[sub_v$bcr == bcr, ], tc), mask = TRUE)
  m <- crop(m, project(limit_v, tc), mask = TRUE)
  mask(m, project(water_v, tc), inverse = TRUE)
}

for (p in pairs) {
  spp <- p[1]; bcr <- p[2]
  obs  <- rast(file.path(ia_dir, "data/derived_data/predictions", spp, bcr, year,
                         "observed_bootstraps.tif"))
  tmpl <- rast(file.path(nm_root, "gis/stacks", paste0(bcr, "_", year, ".tif")))[[1]]
  if (!isTRUE(compareGeom(tmpl, obs[[1]], stopOnError = FALSE)))
    stop(spp, " ", bcr, ": stack template and observed_bootstraps.tif disagree")

  b  <- build_weight(spp, tmpl)
  tc <- crs(tmpl)
  # touches = TRUE on every vector term -- this is the V5 mask()/crop() behaviour and
  # NOT terra::rasterize()'s default. See 12A2 for why it matters so much on water.
  w <- b$range *
    (1 - rasterize(b$ctg(water_v), tmpl, field = 1, background = 0, touches = TRUE)) *
    rasterize(b$ctg(limit_v), tmpl, field = 1, background = 0, touches = TRUE) *
    rasterize(project(sub_v[sub_v$bcr == bcr, ], tc), tmpl, field = 1,
              background = 0, touches = TRUE)

  pA <- pop(obs * w)
  pB <- pop(v5_mask(obs, spp, bcr))
  pC <- pop(v5_mask(project(obs, "EPSG:3978", res = 1000), spp, bcr))

  cat(sprintf("\n=== %s %s ===\n", spp, bcr))
  cat("A  5072 + weight.tif  : ", fmt(pA), "\n")
  cat("B  5072 + V5 vectors  : ", fmt(pB), sprintf("  A/B = %.6f\n", median(pA) / median(pB)))
  cat("C  3978 + V5 vectors  : ", fmt(pC), sprintf("  B/C = %.4f\n", median(pB) / median(pC)))
}
