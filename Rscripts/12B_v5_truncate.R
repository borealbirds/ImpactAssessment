# ---
# title: Impact Assessment: port of V5's analysis/10.Truncate.R (commit f082866)
# author: Mannfred Boehm
# created: September 11, 2026
# ---

# v5_truncate() applies V5's truncation and masking to a stack of bootstrap predictions, in
# 10.Truncate.R's order:
#   4. project to EPSG:3978 at 1 km     skipped when project_to = NULL (12A: V5 calls it legacy,
#                                       and clamp() does not commute with projection)
#   5. clamp at densmax                 a species constant, q.out$densmax
#   6. clamp at q99                     the 99.9th percentile of the bootstrap mean; pass q99 to
#                                       freeze it rather than derive it from this stack
#   7. x range, then NA -> 0            skipped when apply_masks = FALSE (12A: 12C's weight
#   8. crop to BCR and data limit,      masks both sides in 12F instead)
#      mask out water
# the masking path (v5_load_masks, steps 4, 7, 8) serves misc/harnesses/verify_v5_truncate_port.R,
# which checks the port against V5's released output/10_truncated (gate G1).
# returns list(stack, q99, densmax)

library(terra)   # sf is called as sf:: by v5_load_masks() only, so 12A runs where sf is absent

# V5's static masking layers. gis_root is G:/Shared drives/BAM_NationalModels5/gis, or
# data/raw_data/v5_gis (no Subregions_Mosaics shapefile there: want_subregions = FALSE)
v5_load_masks <- function(gis_root, want_subregions = TRUE) {
  list(water      = sf::read_sf(file.path(gis_root, "WaterMask_Canada.shp")),
       limit      = sf::st_transform(sf::read_sf(file.path(gis_root, "DataLimitationsMask.shp")), 3978),   # 10.Truncate.R:64-65
       subregions = if (want_subregions) sf::read_sf(file.path(gis_root, "Subregions_Mosaics_EPSG3978.shp")))
}

v5_truncate <- function(stack_in, spp, bcr, densmax, masks = NULL, range_root = NULL, q99 = NULL,
                        project_to = "EPSG:3978", res = 1000, apply_masks = TRUE) {

  stopifnot(length(densmax) == 1, is.finite(densmax))

  # step 4
  r <- stack_in
  if (!is.null(project_to)) r <- terra::project(r, project_to, res = res)

  # steps 5 and 6
  r <- terra::clamp(r, upper = densmax, values = TRUE)
  if (is.null(q99)) {
    q99 <- terra::global(terra::app(r, mean, na.rm = TRUE), quantile, probs = 0.999, na.rm = TRUE)[1, 1]
    if (is.na(q99)) return(NULL)   # 10.Truncate.R:131
  }
  r <- terra::clamp(r, upper = q99, values = TRUE)

  if (!apply_masks) return(list(stack = r, q99 = q99, densmax = densmax))
  if (is.null(masks) || is.null(range_root)) stop("apply_masks = TRUE needs masks and range_root")

  # step 7: the range as a continuous weight, then every NA to 0 (the crops below bring NA back
  # outside the kept shapes)
  r <- r * terra::resample(terra::rast(file.path(range_root, paste0(spp, ".tif"))), r)
  r[is.na(r)] <- 0

  # step 8
  if (!is.null(masks$subregions)) {
    sf_bcr <- masks$subregions[masks$subregions$bcr == bcr, ]
    if (nrow(sf_bcr) == 0) stop("no subregion polygon for bcr=", bcr)
    r <- terra::crop(r, terra::vect(sf_bcr), mask = TRUE)
  }
  r <- terra::crop(r, terra::vect(masks$limit), mask = TRUE)
  r <- terra::mask(r, terra::vect(masks$water), inverse = TRUE)

  list(stack = r, q99 = q99, densmax = densmax)
}
