# Decisive inspection of the subbasin-1 smoke output: get TRUE backfill cell count
# (high-HF cells) and per-draw finite/NaN/NA counts, independent of any heuristic.
suppressMessages({library(terra)})
terra::terraOptions(progress = 0)
ia_dir <- getwd()
year <- 2020; si <- 1L

stack_y <- terra::rast(file.path(ia_dir,"data","raw_data","covariates_mosaiced",
                                 sprintf("covariates_mosaiced_%d.tif", year)))
highhf <- terra::project(terra::rast(file.path(ia_dir,"data","raw_data","hirshpearson","CanHF_1km_morethan1.tif")), stack_y, method="near")
subs <- terra::project(terra::vect(file.path(ia_dir,"data","raw_data","hydrobasins_masked_merged_subset.gpkg")), stack_y)
sub_s <- subs[si]
cov_s <- terra::mask(terra::crop(stack_y, sub_s), sub_s)
hh <- terra::mask(terra::resample(highhf, cov_s, method="near"), sub_s)
backfill_idx <- which(terra::values(hh) == 1)
cat("TRUE high-HF backfill cells in subbasin 1 =", length(backfill_idx), "\n")

tif <- file.path(ia_dir,"data","derived_data","bart_models",year,
                 sprintf("subbasin_%s",si), sprintf("subbasin_%d_backfill.tif", si))
r <- terra::rast(tif)
cat("tif layers =", terra::nlyr(r), "  ncell =", terra::ncell(r), "\n")
dl <- grep("_draw_001$", names(r), value = TRUE)
cat("\nper continuous-biotic draw_001 layer over the", length(backfill_idx), "backfill cells:\n")
for (ln in head(dl, 6)) {
  v <- terra::values(r[[ln]], mat = FALSE)[backfill_idx]
  cat(sprintf("  %-26s finite=%d (%.1f%%)  NaN=%d  NA=%d\n",
              sub("_draw_001$","",ln), sum(is.finite(v)), 100*mean(is.finite(v)),
              sum(is.nan(v)), sum(is.na(v) & !is.nan(v))))
}
# also: a categorical backfilled layer (e.g. SCANFI_1km) coverage
catl <- intersect(c("SCANFI_1km","VLCE_1km","ABoVE_1km"), names(r))
for (ln in catl) {
  v <- terra::values(r[[ln]], mat = FALSE)[backfill_idx]
  cat(sprintf("  [cat] %-20s finite=%d (%.1f%%)  NA=%d\n",
              ln, sum(is.finite(v)), 100*mean(is.finite(v)), sum(is.na(v))))
}
