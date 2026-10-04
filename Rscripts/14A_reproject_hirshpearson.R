# ---
# title: Build 1 km sector footprint layers from the raw 300 m Hirsh-Pearson rasters
# author: Mannfred Boehm
# ---

# Input: the 300 m sector rasters of Hirsh-Pearson et al. (2022), downloaded unchanged from
# Borealis (doi:10.5683/SP2/EVKAVL, V3) into data/raw_data/hirshpearson/raw_300m/.
# Output, per sector, on the grid of the CanHF masks (EPSG:5072, 1000 m, written by 03):
#   {sector}.tif         highest 300 m pressure score in the 1 km cell (0-10)
#   {sector}_direct.tif  1 if any 300 m cell in the 1 km cell is DIRECT footprint, else 0
#
# 12F uses both. A pixel is in sector j's footprint when {j}.tif > 0 and CanHF >= 1: there
# the footprint covariates are set to 0. Its vegetation is backfilled only where
# {j}_direct.tif is 1 as well (TODO C2f). A 300 m cell is direct footprint when:
#   built, crop, pasture   score > 0 (HP gives them no distance buffers)
#   roads, rail, dams      score >= 6: the 0-300 m band of secondary roads, rail and
#                          reservoirs (HP Tables 2-4). It also takes in the 300-600 m band of
#                          national/major highways and the 300-900 m band of the Trans-Canada,
#                          since a score does not identify the road type.
#   mines                  score >= 6: the 0-600 m band (Table 5), plus 600-1500 m of large mines.
#                          Most sites are large: median 79 300 m cells (~7 km2, a 1.5 km disc)
#   oil_gas                score 10: the core around the site, a disc of ~29 300 m cells (900 m
#                          radius, ~2.6 km2), not one cell (HP scores are integers, 10 at the
#                          site decaying to 1 near 5 km)
# HP mapped mines and oil/gas as points, so these discs, not the real pits or well fields, set
# the direct area; it is 4.8% (mines) and 7.7% (oil_gas) of their 1 km footprint. The "any"
# rule below enlarges each disc to whole 1 km cells (x2.2 mines, x2.9 oil_gas; 2026-10-01).
#
# Aggregation is "max", never bilinear: bilinear downsampling spreads each 300 m cell over
# 2x2 1 km cells, so the layers this script used to write (bilinear, until 2026-09-28) marked
# up to a 1 km ring of cells that hold none of the feature.

library(terra)
terraOptions(memfrac = 0.6)   # project() values can depend on memory (see memory notes)

ia_dir    <- getwd()
hirsh_dir <- file.path(ia_dir, "data/raw_data/hirshpearson")
raw_dir   <- file.path(hirsh_dir, "raw_300m")

template_r <- rast(file.path(hirsh_dir, "CanHF_1km_morethan1.tif"))

direct_rule <- list(
  built = function(x) x > 0, crop = function(x) x > 0, pasture = function(x) x > 0,
  roads = function(x) x >= 6, rail = function(x) x >= 6,
  dam_and_associated_reservoir = function(x) x >= 6,
  mines = function(x) x >= 6, oil_gas = function(x) x >= 10)

# optional: build only the sectors named on the command line (e.g. Rscript 14A... mines oil_gas)
only <- commandArgs(trailingOnly = TRUE)
if (length(only)) {
  if (!all(only %in% names(direct_rule)))
    stop("unknown sector(s): ", paste(setdiff(only, names(direct_rule)), collapse = ", "))
  direct_rule <- direct_rule[only]
}

for (sec in names(direct_rule)) {
  f <- file.path(raw_dir, paste0(sec, ".tif"))
  if (!file.exists(f)) stop("missing ", f, " - download it from doi:10.5683/SP2/EVKAVL")
  r <- rast(f)
  message(sec, ": ", paste(dim(r)[1:2], collapse = " x "), " cells at ",
          paste(round(res(r)), collapse = " x "), " m, ", crs(r, describe = TRUE)$name)

  score <- project(r, template_r, method = "max")
  names(score) <- sec
  writeRaster(score, file.path(hirsh_dir, paste0(sec, ".tif")), overwrite = TRUE,
              datatype = "FLT4S")

  direct_300 <- classify(direct_rule[[sec]](r), cbind(NA, 0))
  direct <- project(direct_300, template_r, method = "max")
  names(direct) <- paste0(sec, "_direct")
  writeRaster(direct, file.path(hirsh_dir, paste0(sec, "_direct.tif")), overwrite = TRUE,
              datatype = "INT1U")

  n_any <- global(score > 0, "sum", na.rm = TRUE)[[1]]
  n_dir <- global(direct == 1, "sum", na.rm = TRUE)[[1]]
  message(sprintf("  1 km cells with %s > 0: %d; direct: %d (%.1f%%)", sec, n_any, n_dir,
                  100 * n_dir / max(1, n_any)))
}
