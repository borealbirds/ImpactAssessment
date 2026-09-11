# ---
# title: Impact Assessment: faithful port of V5 `analysis/10.Truncate.R`
# author: Mannfred Boehm
# created: September 11, 2026
# ---
#
# Sourced helper. Defines v5_truncate(), a line-for-line port of the truncation /
# masking transform that BAM V5 applies to bird density predictions, as of repo
# commit f082866 ("split package script into two steps", 2026-06-04), which split
# the old 10.Package.R into 10.Truncate.R + 11.Package.R.
#
# V5 order of operations (10.Truncate.R steps 4-8):
#   4. project to EPSG:3978 @ 1 km            <- legacy delivery step, see `project_to`
#   5. clamp(stack, upper = densmax)          <- species constant, from q.out$densmax
#   6. q99 = 99.9th pctile of the MEAN of the densmax-clamped stack;
#      clamp(stack, upper = q99)              <- the dominant truncation
#   7. stack * range; then NA -> 0
#   8. crop to BCR polygon, crop to data-limit, mask out water
#
# Two deliberate generalisations, both required by the counterfactual use case:
#
#   `q99`  -- pass a value to FREEZE the secondary cap instead of deriving it from
#             the stack being transformed. This is mandatory when transforming a
#             counterfactual: q99 is the only data-dependent parameter in the
#             transform, and letting it adapt to each counterfactual landscape
#             would put a component of the cap itself into the obs/bf contrast.
#             Derive it once from the OBSERVED stack and pass it everywhere else.
#             NULL reproduces V5's own behaviour (derive from this stack).
#
#   `project_to` -- set NULL to skip step 4 and stay in the prediction's native
#             EPSG:5072. 10.Truncate.R:19 documents the 3978 reprojection as a
#             legacy artifact ("Future versions will not require this step").
#             It costs ~2.3% of total abundance (bilinear resampling does not
#             conserve sums) and, because clamp() does not commute with
#             projection, it would break the exact superset -> masked-rowsum
#             decomposition that 12C relies on. Production runs in 5072; 3978 is
#             used only to reproduce V5's released product as a verification gate.
#
#   `apply_masks` -- set FALSE to stop after step 6 (the two clamps) and skip
#             steps 7-8. Our pipeline carries V5's range/water/data-limit masking
#             as a separate precomputed weight.tif (12A2), which 12C multiplies
#             into BOTH the observed and the backfilled side so the obs/bf
#             symmetry w*bf - w*obs = w*(bf - obs) is exact. 12A must therefore
#             write a truncated-but-UNWEIGHTED stack, or 12C:220 would apply the
#             masking a second time. Note that q99 is derived BEFORE masking in
#             V5 too (step 6 precedes step 7), so the frozen parameter is
#             identical either way.
#
# Returns a list(stack, q99, densmax) so callers can persist the frozen params.
# ---

suppressPackageStartupMessages({
  library(terra)
  library(sf)
})

# Load the four static V5 masking layers once; reused across species and BCRs.
# `gis_root` is G:/Shared drives/BAM_NationalModels5/gis when validating against
# V5's own products, or ia_dir/data/raw_data/v5_gis for the staged cluster copy
# (which has no Subregions_Mosaics shapefile -- pass bcr_poly = NULL there).
v5_load_masks <- function(gis_root, want_subregions = TRUE) {
  water <- sf::read_sf(file.path(gis_root, "WaterMask_Canada.shp"))

  # 10.Truncate.R:64-65 -- the data-limit mask ships in 5072 and is transformed.
  limit <- sf::read_sf(file.path(gis_root, "DataLimitationsMask.shp")) |>
    sf::st_transform(3978)

  subregions <- NULL
  if (want_subregions) {
    subregions <- sf::read_sf(file.path(gis_root, "Subregions_Mosaics_EPSG3978.shp"))
  }

  list(water = water, limit = limit, subregions = subregions)
}

v5_truncate <- function(stack_in,
                        spp,
                        bcr,
                        densmax,
                        masks       = NULL,
                        range_root  = NULL,
                        q99         = NULL,
                        project_to  = "EPSG:3978",
                        res         = 1000,
                        apply_masks = TRUE) {

  stopifnot(length(densmax) == 1, is.finite(densmax))

  # -- step 4: project (legacy; skipped when project_to is NULL) ----------------
  r <- stack_in
  if (!is.null(project_to)) {
    r <- terra::project(r, project_to, res = res)
  }

  # -- step 5: species density cap ---------------------------------------------
  r <- terra::clamp(r, upper = densmax, values = TRUE)

  # -- step 6: secondary cap at the 99.9th pctile of the bootstrap mean --------
  if (is.null(q99)) {
    mn  <- terra::app(r, mean, na.rm = TRUE)
    q99 <- terra::global(mn, quantile, probs = 0.999, na.rm = TRUE)[1, 1]
    if (is.na(q99)) return(NULL)   # 10.Truncate.R:131
  }
  r <- terra::clamp(r, upper = q99, values = TRUE)

  if (!apply_masks) return(list(stack = r, q99 = q99, densmax = densmax))
  if (is.null(masks) || is.null(range_root))
    stop("apply_masks = TRUE requires both `masks` and `range_root`")

  # -- step 7: range as a multiplicative weight, then NA -> 0 -------------------
  # V5 applies the range continuously (not as a hard cutoff) and zeroes every NA,
  # including NAs originating in the prediction itself. The crops in step 8
  # reintroduce NA outside the retained shapes.
  range_r <- terra::rast(file.path(range_root, paste0(spp, ".tif"))) |>
    terra::resample(r)
  r <- r * range_r
  r[is.na(r)] <- 0

  # -- step 8: BCR extent, data limit, water -----------------------------------
  if (!is.null(masks$subregions)) {
    sf_bcr <- masks$subregions[masks$subregions$bcr == bcr, ]
    if (nrow(sf_bcr) == 0) stop("no subregion polygon for bcr=", bcr)
    r <- terra::crop(r, terra::vect(sf_bcr), mask = TRUE)
  }
  r <- terra::crop(r, terra::vect(masks$limit), mask = TRUE)
  r <- terra::mask(r, terra::vect(masks$water), inverse = TRUE)

  list(stack = r, q99 = q99, densmax = densmax)
}
