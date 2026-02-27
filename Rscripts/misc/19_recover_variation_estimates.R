# this is a temporary script that recovers CV estimates 

library(terra)
library(tidyverse)

# stack_boot : SpatRaster with one layer per bootstrap (densities)
# subbasins  : SpatVector of subbasins (must contain id_col)
# industry   : SpatVector of industry polygons (any/multiple polygons ok)
# species    : character (e.g., "CAWA")
# treatment  : "observed" or "backfilled" (any label you prefer)
# id_col     : column name in `subbasins` that holds the subbasin id
# per_ha     : TRUE if densities are per hectare, FALSE if already per km^2
# nan_to_zero: TRUE to coerce NaN (all-NA cells in a subbasin) to 0, FALSE to NA
boot_stats_from_stack_industry <- function(stack_boot, subbasins, industry,
                                           species, treatment,
                                           id_col = "first_HYBAS_ID",
                                           per_ha = TRUE,
                                           nan_to_zero = FALSE) {
  # coerce inputs
  if (is.character(industry)) industry <- terra::vect(industry)
  stopifnot(inherits(stack_boot, "SpatRaster"), inherits(subbasins, "SpatVector"),
            inherits(industry, "SpatVector"))
  
  # align CRS to raster
  if (!terra::same.crs(industry, stack_boot))  industry  <- terra::project(industry,  stack_boot)
  if (!terra::same.crs(subbasins, stack_boot)) subbasins <- terra::project(subbasins, stack_boot)
  
  # mask to industry (ROI may be empty in some regions)
  r <- tryCatch(terra::mask(terra::crop(stack_boot, industry), industry),
                error = function(e) NULL)
  if (is.null(r) || terra::nlyr(r) == 0) {
    return(tibble(
      species, subbasin = subbasins[[id_col]], treatment,
      population_mean = NA_real_, population_sd = NA_real_, population_cv = NA_real_
    ))
  }
  
  # convert density → abundance per cell (km^2); multiply by 100 if input is per ha
  area_km2 <- terra::cellSize(r, unit = "km")
  r_abund  <- (if (per_ha) r * 100 else r) * area_km2
  
  # per-subbasin bootstrap totals (matrix: n_subbasins × n_bootstraps)
  mat <- terra::extract(r_abund, subbasins, fun = sum, na.rm = TRUE, ID = FALSE) |> as.matrix()
  
  means <- rowMeans(mat, na.rm = TRUE)
  sds   <- apply(mat, 1, sd, na.rm = TRUE)
  cvs   <- sds / means
  
  # NaN handling
  means[is.nan(means)] <- if (nan_to_zero) 0 else NA_real_
  sds[is.nan(sds)]     <- NA_real_
  cvs[is.nan(cvs)]     <- NA_real_
  
  tibble(
    species         = species,
    subbasin        = subbasins[[id_col]],
    treatment       = treatment,
    population_mean = means,
    population_sd   = sds,
    population_cv   = cvs
  )
}

# inputs you already have
test <- terra::rast(file.path(root, "output", "08_mosaics_can", "OSFL_2020.tif"))
subbasins <- vect(file.path(ia_dir, "hydrobasins_masked_merged_subset.gpkg"))
industry  <- file.path(ia_dir, "combined_industry_footprint.gpkg")  # all industry footprints nationwide

res_cawa_obs <- boot_stats_from_stack_industry(
  stack_boot = test,
  subbasins  = subbasins,
  industry   = industry,
  species    = "CAWA",
  treatment  = "observed",
  id_col     = "first_HYBAS_ID",
  per_ha     = TRUE,        # set FALSE if your stack is already per km^2
  nan_to_zero = FALSE       # or TRUE if you want empty subbasins to be 0
)

regional_cv <- res_cawa_obs %>%
  summarise(
    total_mean = sum(population_mean, na.rm = TRUE),
    total_var  = sum(population_sd^2, na.rm = TRUE),
    total_sd   = sqrt(total_var),
    total_cv   = total_sd / total_mean
  )

regional_cv
