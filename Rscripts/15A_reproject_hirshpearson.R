# ---
# title: Reproject and align Hirsh-Pearson sector rasters to pipeline CRS
# author: Mannfred Boehm
# ---

# The individual sector rasters in data/raw_data/hirshpearson/ (built, crop,
# forestry_harvest, etc.) may differ in CRS and extent from the rest of the
# pipeline, which uses EPSG:5072 (NAD83(NSRS2007) / Conus Albers) and the
# spatial extent of the Level 6 HydroBASINS study area.
#
# This script reprojects all non-CanHF sector rasters in-place to:
#   - CRS:        EPSG:5072 (taken from hydrobasins_masked_merged_subset.gpkg)
#   - Extent:     bounding box of the hydrobasins shapefile
#   - Resolution: 1000 m x 1000 m
#   - Method:     bilinear (scores are continuous)
#
# Run this script once before 15_sector_attribution.R.
# The CanHF* files are left untouched.
# ---

library(terra)

# ---- Execution context -------------------------------------------------------

cc    <- FALSE
local <- TRUE

if (cc)            { ia_dir <- "/home/mannfred/scratch/impact_assessment" }
if (!cc && local)  { ia_dir <- getwd() }
if (!cc && !local) { ia_dir <- file.path("G:/Shared drives/BAM_NationalModels5", "data", "Extras",
                                          "sandbox_data", "impactassessment_sandbox") }

# ---- Paths ------------------------------------------------------------------

hirsh_dir  <- file.path(ia_dir, "data/raw_data/hirshpearson")
basin_path <- file.path(ia_dir, "data/raw_data/hydrobasins_masked_merged_subset.gpkg")

# ---- Build spatial template from hydrobasins --------------------------------
# The hydrobasins shapefile is the authoritative spatial reference for the
# pipeline. Use its CRS and bounding box at 1 km resolution.

hydrobasins <- vect(basin_path)

message("Hydrobasins CRS: ", crs(hydrobasins, describe = TRUE)$code)
message("Hydrobasins extent: ", paste(as.vector(ext(hydrobasins)), collapse = ", "))

template_r <- rast(ext(hydrobasins), resolution = 1000, crs = crs(hydrobasins))

# ---- Identify sector files to reproject -------------------------------------

sector_files <- list.files(hirsh_dir, pattern = "\\.tif$", full.names = TRUE)
sector_files <- sector_files[!grepl("^CanHF", basename(sector_files))]

message("\nReprojecting ", length(sector_files), " sector rasters to EPSG:5072 (1 km):")

# ---- Reproject and overwrite ------------------------------------------------

for (f in sector_files) {
  message("  ", basename(f))

  r <- rast(f)
  message("    source CRS:    ", crs(r, describe = TRUE)$code,
          "  extent: ", paste(round(as.vector(ext(r)), 0), collapse = ", "))

  r_proj <- project(r, template_r, method = "bilinear")

  writeRaster(r_proj, filename = f, overwrite = TRUE)
  message("    written.")
}

message("\nDone. All sector rasters now conform to the hydrobasins grid.")
