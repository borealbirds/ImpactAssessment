# ---
# title: compare BCR backfill mosaics against a test build
# created: October 2, 2026
# ---
#
# Checks a change to 11_premosaic_backfilled_stacks.R: compares each
# bart_models_mosaics/{year}/{bcr}_backfilled.tif with the copy that a run with
# MOSAIC_DIR=<test dir> wrote, layer by layer (names, grid, identical() values).
# Prints one "COMPARE <bcr> MATCH|DIFFER|MISSING" line per BCR.
#
# Usage (cluster job, see compare_mosaics.sh): Rscript compare_mosaics.R can11 can42 ...
# YEARS (default 2020) and MOSAIC_DIR (default bart_models_mosaics_test) as for 11.

suppressPackageStartupMessages(library(terra))

ia_dir   <- "/home/mannfred/scratch/impact_assessment"
year     <- as.integer(Sys.getenv("YEARS", "2020"))
test_dir <- Sys.getenv("MOSAIC_DIR", "bart_models_mosaics_test")

for (bcr in commandArgs(trailingOnly = TRUE)) {
  fa <- file.path(ia_dir, "data", "derived_data", "bart_models_mosaics", year, paste0(bcr, "_backfilled.tif"))
  fb <- file.path(ia_dir, "data", "derived_data", test_dir, year, paste0(bcr, "_backfilled.tif"))
  if (!file.exists(fa) || !file.exists(fb)) {
    message("COMPARE ", bcr, " MISSING ", if (!file.exists(fa)) fa else fb)
    next
  }
  a <- rast(fa)
  b <- rast(fb)
  if (!identical(names(a), names(b)) || !compareGeom(a, b, stopOnError = FALSE)) {
    message("COMPARE ", bcr, " DIFFER names or grid (", nlyr(a), " vs ", nlyr(b), " layers)")
    next
  }
  bad <- character(0)
  for (i in seq_len(nlyr(a))) {
    if (!identical(values(a[[i]], mat = FALSE), values(b[[i]], mat = FALSE))) bad <- c(bad, names(a)[i])
  }
  message("COMPARE ", bcr, if (length(bad) == 0) " MATCH " else " DIFFER ", nlyr(a), " layers, ",
          length(bad), " differ", if (length(bad) > 0) paste0(" (first: ", paste(head(bad, 5), collapse = ", "), ")"))
}
