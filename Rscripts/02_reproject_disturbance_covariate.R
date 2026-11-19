# ---
# title: Impact Assessment: reproject time-since-disturbance to BAM crs
# author: Mannfred Boehm
# created: June 3, 2025
# ---

# year of last fire (Hermosilla et al. 2016, 30 m, Lambert conformal conic) -> aggregate to
# ~1 km -> reproject -> crop -> mask, then one "years since fire" layer per year
# downloaded from: https://opendata.nfis.org/mapserver/nfis-change_eng.html

library(terra)

root    <- "G:/Shared drives/BAM_NationalModels5"
raw_dir <- file.path(root, "data", "Extras", "sandbox_data", "impactassessment_sandbox", "CAfire")
out_dir <- file.path("data", "raw_data", "CAfire")

CAfire <- terra::rast(file.path(raw_dir, "CA_Forest_Fire_1985-2020.tif"))

# 33 x 33 cells of 30 m -> one ~1 km cell (~5 min; saves hours of reprojection)
fact <- round(1000 / terra::res(CAfire)[1])
CAfire_agg <- terra::aggregate(CAfire, fact = fact, fun = "modal", na.rm = TRUE)

# 1 km template in EPSG:5072 (projecting onto a template avoids GDAL errors)
bam_boundary <- terra::vect(file.path(root, "Regions", "BAM_BCR_NationalModel_Buffered.shp"))
bam_template <- terra::rast(terra::ext(bam_boundary), resolution = 1000, crs = terra::crs(bam_boundary))

CAfire_mask <-
  terra::project(CAfire_agg, bam_template, method = "near") |>   # years are categorical
  terra::crop(bam_boundary) |>
  terra::mask(bam_boundary)

# 0 = no fire on record
CAfire_mask[CAfire_mask == 0] <- NA

# years since fire in `current_year`, scaled to (0, 1]: 1 = burned this year, -> 0 for old fires
find_year_since_fire <- function(current_year) {
  ysf <- current_year - CAfire_mask
  ysf[ysf < 0] <- NA                  # burned after current_year
  ysf_norm <- 1 / (ysf + 1)

  # never burned is the ysf -> Inf limit, a real 0, not missing: left NA it dropped 99% of
  # pixels from the backfill (Open Limitation #5). re-mask so only off-study-area stays NA
  ysf_norm[is.na(ysf_norm)] <- 0
  ysf_norm <- terra::mask(ysf_norm, bam_boundary)
  names(ysf_norm) <- paste0("CAfire_", current_year)
  ysf_norm
}

for (y in seq(1990, 2020, by = 5)) {
  terra::writeRaster(find_year_since_fire(y), file.path(out_dir, paste0("CAfire_", y, "_masked.tif")),
                     overwrite = TRUE)
}
