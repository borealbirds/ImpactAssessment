# ---
# title: Is the 5072-vs-3978 population gap a resampling loss or an area-units artifact?
# author: Mannfred Boehm
# created: September 14, 2026
# ---
#
# Runs LOCALLY. Needs only staged files (data/raw_data/v5_gis + Regions/), not G:.
#
# Background. Population is summed as `density * 100` (birds/ha -> birds/km2), which gives
# birds ONLY if every pixel is one square kilometre of GROUND. That holds in an equal-area
# projection and fails in a conformal one. EPSG:5072 is Albers Equal Area; EPSG:3978 (V5's
# delivery CRS) is Lambert Conformal Conic, whose scale factor dips below 1 between its
# standard parallels (49N, 77N) -- so its 1000 m cells cover MORE than 1 km2 of ground.
#
# This script computes, per species x BCR, both the naive sum (`* 100`, one km2 assumed per
# pixel) and the area-weighted sum (`* 100 * cellSize()`), in each CRS. If the gap is an
# area-units artifact the naive sums disagree and the area-weighted sums agree.
#
# RESULT 2026-09-14 (CAWA can10 / CAWA can71 / OVEN can10):
#   mean true cell area   5072 1.00000 km2         | 3978 1.02682 / 1.03085 / 1.02682 km2
#   naive sum 5072/3978   1.0288 / 1.0284 / 1.0346
#   area-wtd  5072/3978   1.0006 / 1.0005 / 1.0023
# So essentially the whole gap is area units; bilinear interpolation contributes ~0.1%.
# Our 5072 totals satisfy the `* 100` assumption and V5's 3978 totals do not.
# See CLAUDE.md Open Limitation #8.
# ---

suppressPackageStartupMessages({ library(terra) })
terraOptions(memfrac = 0.6)

ia_dir <- getwd()
gis    <- file.path(ia_dir, "data", "raw_data", "v5_gis")
year   <- 2020

water_v <- vect(file.path(gis, "WaterMask_Canada.shp"))
limit_v <- vect(file.path(gis, "DataLimitationsMask.shp"))
bcrp    <- vect(file.path(ia_dir, "data", "raw_data", "Regions",
                          "BAM_BCR_NationalModel_Unbuffered.shp"))
bcrp$code <- gsub("_", "", paste(bcrp$country, bcrp$subUnit, sep = "_"))

pairs <- list(c("CAWA", "can10"), c("CAWA", "can71"), c("OVEN", "can10"))

# V5 10.Truncate.R steps 7-8, in whatever CRS `r` is already in
v5_mask <- function(r, spp, bcr) {
  tc <- crs(r)
  rr <- rast(file.path(gis, "ranges", paste0(spp, ".tif")))
  rr <- if (crs(rr) == tc) resample(rr, r) else project(rr, r, method = "bilinear")
  m <- r * rr
  m[is.na(m)] <- 0
  m <- crop(m, project(bcrp[bcrp$code == bcr, ], tc), mask = TRUE)
  m <- crop(m, project(limit_v, tc), mask = TRUE)
  mask(m, project(water_v, tc), inverse = TRUE)
}

for (p in pairs) {
  spp <- p[1]; bcr <- p[2]
  obs <- rast(file.path(ia_dir, "data/derived_data/predictions", spp, bcr, year,
                        "observed_bootstraps.tif"))[[1]]
  B <- v5_mask(obs, spp, bcr)                                  # 5072, equal-area
  C <- v5_mask(project(obs, "EPSG:3978", res = 1000), spp, bcr) # 3978, conformal

  naiveB <- global(B * 100, "sum", na.rm = TRUE)[[1]]
  naiveC <- global(C * 100, "sum", na.rm = TRUE)[[1]]
  aB <- cellSize(B, unit = "km"); aC <- cellSize(C, unit = "km")
  trueB <- global(B * 100 * aB, "sum", na.rm = TRUE)[[1]]
  trueC <- global(C * 100 * aC, "sum", na.rm = TRUE)[[1]]
  mB <- global(ifel(is.na(B), NA, aB), "mean", na.rm = TRUE)[[1]]
  mC <- global(ifel(is.na(C), NA, aC), "mean", na.rm = TRUE)[[1]]
  nB <- global(!is.na(B), "sum", na.rm = TRUE)[[1]]
  nC <- global(!is.na(C), "sum", na.rm = TRUE)[[1]]

  cat(sprintf("\n=== %s %s ===\n", spp, bcr))
  cat(sprintf("  mean TRUE cell area  5072 %.5f km2 | 3978 %.5f km2  (3978/5072 = %.4f)\n",
              mB, mC, mC / mB))
  cat(sprintf("  n data cells         5072 %-8d | 3978 %-8d  (3978/5072 = %.4f)\n",
              nB, nC, nC / nB))
  cat(sprintf("  NAIVE  sum (x1 km2)  5072 %-9.0f| 3978 %-9.0f -> 5072/3978 = %.4f\n",
              naiveB, naiveC, naiveB / naiveC))
  cat(sprintf("  AREA-WEIGHTED sum    5072 %-9.0f| 3978 %-9.0f -> 5072/3978 = %.4f\n",
              trueB, trueC, trueB / trueC))
}
