# Local de-risk for the CAfire NA-by-design fix (HANDOFF_cafire_backfill_fix.md §5).
# Recodes the EXISTING CAfire_2020_masked.tif (already the ysf_norm output of 02)
# NA -> 0 within the study area, instead of re-running 02's ~6h reprojection.
suppressMessages({library(terra)})

root   <- "G:/Shared drives/BAM_NationalModels5"
ia_dir <- file.path(root, "data", "Extras", "sandbox_data", "impactassessment_sandbox")
caf_path <- file.path(ia_dir, "CAfire", "CAfire_2020_masked.tif")

stopifnot(file.exists(caf_path))
caf <- terra::rast(caf_path)

ncell_tot <- terra::ncell(caf)
na_before <- terra::global(is.na(caf), "sum", na.rm = TRUE)[1, 1]
cat(sprintf("BEFORE: %d / %d cells NA (%.3f)\n", na_before, ncell_tot, na_before / ncell_tot))
cat("value summary (non-NA) BEFORE:\n"); print(summary(terra::values(caf)[, 1]))

# boundary used by 02 (buffered)
bam_boundary <- terra::vect(file.path(root, "Regions", "BAM_BCR_NationalModel_Buffered.shp"))
if (terra::crs(bam_boundary) != terra::crs(caf)) bam_boundary <- terra::project(bam_boundary, caf)

# in-study NA fraction (this is the cascade seed ~0.99 we expect)
caf_inb   <- terra::mask(caf, bam_boundary)
na_inb    <- terra::global(is.na(caf_inb), "sum", na.rm = TRUE)[1, 1]
ncell_inb <- terra::global(!is.na(terra::mask(caf * 0 + 1, bam_boundary)), "sum", na.rm = TRUE)[1, 1]
cat(sprintf("BEFORE (in-study): %d / %d in-boundary cells NA (%.3f)  <-- cascade seed\n",
            na_inb, ncell_inb, na_inb / ncell_inb))

# RECODE: NA -> 0, then re-mask so off-study-area returns to NA
caf_fix <- caf
caf_fix[is.na(caf_fix)] <- 0
caf_fix <- terra::mask(caf_fix, bam_boundary)

na_after_inb <- terra::global(is.na(caf_fix), "sum", na.rm = TRUE)[1, 1]
# in-study after: should be ~0 (only off-boundary remains NA, but we masked to boundary so global NA == off-boundary)
cat(sprintf("AFTER: %d cells NA globally (off-study only); in-study NA should be ~0\n", na_after_inb))
cat("value summary AFTER (full):\n"); print(summary(terra::values(caf_fix)[, 1]))

# back up original then overwrite in place (06 reads this exact filename)
bak <- sub("\\.tif$", "_PREFIX_orig.tif", caf_path)
if (!file.exists(bak)) {
  terra::writeRaster(caf, bak, overwrite = FALSE)
  cat("backed up original ->", bak, "\n")
} else {
  cat("backup already exists, not overwriting backup ->", bak, "\n")
}
terra::writeRaster(caf_fix, caf_path, overwrite = TRUE)
cat("wrote recoded CAfire ->", caf_path, "\n")
