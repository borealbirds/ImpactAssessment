# ---
# title: Build 1 km sector footprint layers from the raw 300 m Hirsh-Pearson rasters
# author: Mannfred Boehm
# ---

# input: the 300 m sector rasters of Hirsh-Pearson et al. (2022), unchanged from Borealis
# (doi:10.5683/SP2/EVKAVL, V3), in data/raw_data/hirshpearson/raw_300m/. output, per sector, on
# the grid of 03's CanHF masks:
#   {sector}.tif         highest 300 m pressure score in the 1 km cell (0-10)
#   {sector}_direct.tif  1 if any 300 m cell in the 1 km cell is direct footprint
# in 12F a pixel is in sector j's footprint where {j}.tif > 0 and CanHF >= 1 (footprint
# covariates set to 0), and its vegetation is backfilled only where {j}_direct.tif is 1 too
# (OL #10). a 300 m cell is direct footprint when:
#   built, crop, pasture   score > 0 (no distance buffers)
#   roads, rail, dams      score >= 6: the 0-300 m band of secondary roads, rail and reservoirs
#                          (HP Tables 2-4), plus 300-600 m of major highways and 300-900 m of the
#                          Trans-Canada, as a score does not say the road type
#   mines                  score >= 6: the 0-600 m band (Table 5), plus 600-1500 m of large mines
#   oil_gas                score 10: a ~900 m disc around the site (scores fall from 10 to 1 over ~5 km)
# mines and oil/gas are mapped as points, so their direct area is a disc, not the real pits or
# well fields: 4.8% and 7.7% of their 1 km footprint (x2.2 and x2.9 once enlarged to whole 1 km cells).
# aggregation by "max", never bilinear, which spreads each 300 m cell over 2 x 2 km cells
# optional: only the sectors named, e.g. Rscript Rscripts/14A_reproject_hirshpearson.R mines oil_gas

library(terra)
terraOptions(memfrac = 0.6)   # project() returns different values under memory pressure

hirsh_dir  <- file.path(getwd(), "data/raw_data/hirshpearson")
raw_dir    <- file.path(hirsh_dir, "raw_300m")
template_r <- rast(file.path(hirsh_dir, "CanHF_1km_morethan1.tif"))

direct_rule <- list(
  built = function(x) x > 0, crop = function(x) x > 0, pasture = function(x) x > 0,
  roads = function(x) x >= 6, rail = function(x) x >= 6,
  dam_and_associated_reservoir = function(x) x >= 6,
  mines = function(x) x >= 6, oil_gas = function(x) x >= 10)

only <- commandArgs(trailingOnly = TRUE)
if (length(only)) {
  if (!all(only %in% names(direct_rule))) stop("unknown sector(s): ", paste(setdiff(only, names(direct_rule)), collapse = ", "))
  direct_rule <- direct_rule[only]
}

for (sec in names(direct_rule)) {
  f <- file.path(raw_dir, paste0(sec, ".tif"))
  if (!file.exists(f)) stop("missing ", f, " - download it from doi:10.5683/SP2/EVKAVL")
  r <- rast(f)

  score <- project(r, template_r, method = "max")
  names(score) <- sec
  writeRaster(score, file.path(hirsh_dir, paste0(sec, ".tif")), overwrite = TRUE, datatype = "FLT4S")

  direct <- project(classify(direct_rule[[sec]](r), cbind(NA, 0)), template_r, method = "max")
  names(direct) <- paste0(sec, "_direct")
  writeRaster(direct, file.path(hirsh_dir, paste0(sec, "_direct.tif")), overwrite = TRUE, datatype = "INT1U")

  n_any <- global(score > 0, "sum", na.rm = TRUE)[[1]]
  n_dir <- global(direct == 1, "sum", na.rm = TRUE)[[1]]
  message(sprintf("%s: %d cells > 0, %d direct (%.1f%%)", sec, n_any, n_dir, 100 * n_dir / max(1, n_any)))
}
