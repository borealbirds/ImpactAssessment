# diag_can71_v2.R — apples-to-apples coverage check.
#
# Refines diag_can71.R by separating two effects that both reduce mosaic
# pixel counts vs the per-subbasin union:
#   (A) BCR polygon clip: the mosaic is masked by bcr_poly; the original
#       union counted pixels within the BCR raster's bounding box. For
#       irregular / coastal BCRs this can be a large effect.
#   (B) A real pixel-drop bug in 11_premosaic_backfilled_stacks.R.
#
# This script reports:
#   poly_n             pixels inside the BCR polygon on the ref grid (cap)
#   mosaic_n           current mosaic _draw_001 footprint (polygon-clipped)
#   union_bbox_n       per-subbasin union over the BCR bbox (= old metric)
#   union_polygon_n    per-subbasin union, ALSO masked by bcr_poly
#
# Apples-to-apples ratio = mosaic_n / union_polygon_n.
#   ratio ~= 1.0  -> 11_premosaic is fine; the original 11.7x came from the
#                    polygon clip, not a bug.
#   ratio << 1.0  -> a real pixel-drop bug remains; localize next in the
#                    pre-resample / cover() chain.
#
# Run on a Fir login node:
#   Rscript /home/mannfred/scratch/impact_assessment/Rscripts/diag_can71_v2.R
suppressPackageStartupMessages({
  library(terra)
})
terraOptions(progress = 0)

ia_dir <- "/home/mannfred/scratch/impact_assessment"
base   <- file.path(ia_dir, "data", "derived_data")
sb_dir <- file.path(base, "bart_models", "2020")
mo_dir <- file.path(base, "bart_models_mosaics", "2020")

bam_boundary <- terra::vect(file.path(ia_dir, "data", "raw_data", "Regions",
                                      "BAM_BCR_NationalModel_Unbuffered.shp"))
bam_bcr_codes <- gsub("_", "", paste(bam_boundary$country,
                                     bam_boundary$subUnit, sep = "_"))

for (b in c("can71", "can61")) {
  f <- file.path(mo_dir, paste0(b, "_backfilled.tif"))
  if (!file.exists(f)) { cat(b, "mosaic MISSING\n"); next }

  mo <- rast(f)
  dl <- grep("_draw_001$", names(mo), value = TRUE)
  if (!length(dl)) { cat(b, "mosaic has no _draw_ layer\n"); next }
  ref      <- mo[[dl[1]]]
  mosaic_n <- terra::global(ref, "notNA")[1, 1]

  bcr_poly <- bam_boundary[bam_bcr_codes == b, ]

  # Polygon footprint on the ref grid (upper bound for any mosaic).
  poly_mask <- terra::rasterize(bcr_poly, ref, field = 1L, background = NA)
  poly_n    <- terra::global(!is.na(poly_mask), "sum", na.rm = TRUE)[1, 1]

  e2 <- ext(ref)

  counter <- ref            # empty 0/1 accumulator on the BCR grid
  counter[] <- 0L
  n_contrib <- 0L

  for (i in 1:674) {
    sf <- file.path(sb_dir, paste0("subbasin_", i),
                    paste0("subbasin_", i, "_backfill.tif"))
    if (!file.exists(sf)) next
    r   <- rast(sf)
    sdl <- grep("_draw_001$", names(r), value = TRUE)
    if (!length(sdl)) next
    e1 <- ext(r)
    if (!(e1$xmin < e2$xmax && e1$xmax > e2$xmin &&
          e1$ymin < e2$ymax && e1$ymax > e2$ymin)) next   # no overlap
    lc <- tryCatch(terra::crop(r[[sdl[1]]], e2), error = function(e) NULL)
    if (is.null(lc) || terra::ncell(lc) == 0) next
    lc <- terra::resample(lc, ref, method = "near")
    m  <- terra::ifel(is.na(lc), 0L, 1L)
    if (terra::global(m, "sum", na.rm = TRUE)[1, 1] > 0) {
      counter   <- counter + m
      n_contrib <- n_contrib + 1L
    }
  }

  union_bbox    <- counter > 0
  union_bbox_n  <- terra::global(union_bbox, "sum", na.rm = TRUE)[1, 1]

  union_poly    <- union_bbox & !is.na(poly_mask)
  union_poly_n  <- terra::global(union_poly, "sum", na.rm = TRUE)[1, 1]

  cat(sprintf(
    "%s : poly_n=%d | mosaic_n=%d | union_bbox=%d | union_polygon=%d | contrib=%d\n",
    b, poly_n, mosaic_n, union_bbox_n, union_poly_n, n_contrib))
  cat(sprintf(
    "       mosaic / union_polygon = %.3f  (apples-to-apples; ~1.0 = OK)\n",
    mosaic_n / union_poly_n))
  cat(sprintf(
    "       mosaic / poly          = %.3f  (frac of BCR polygon covered)\n",
    mosaic_n / poly_n))
  cat(sprintf(
    "       union_polygon / union_bbox = %.3f  (frac of bbox-union inside polygon)\n",
    union_poly_n / union_bbox_n))
}
cat("done.\n")
