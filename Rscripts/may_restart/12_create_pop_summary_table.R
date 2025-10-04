library(terra)
library(tidyverse)

rebuild_population_from_mean_rasters_industry <- function(species_vec, year) {
  out_root <- file.path(ia_dir, "density_predictions")
  subbasins <- terra::vect(file.path(ia_dir, "hydrobasins_masked_merged_subset.gpkg"))
  industry  <- terra::vect(file.path(ia_dir, "combined_industry_footprint.gpkg"))
  
  bind_rows(lapply(species_vec, function(sp) {
    dir_sp <- file.path(out_root, sp, as.character(year))
    if (!dir.exists(dir_sp)) return(NULL)
    
    fns <- list.files(
      dir_sp,
      pattern = paste0("^", sp, "_[A-Za-z]+\\d+_", year, "_(observed|backfilled)\\.tif$"),
      full.names = TRUE
    )
    if (length(fns) == 0) return(NULL)
    
    bind_rows(lapply(fns, function(f) {
      # parse metadata from filename
      base <- basename(f)
      bcr_code  <- sub(paste0("^", sp, "_([A-Za-z]+\\d+)_", year, "_(observed|backfilled)\\.tif$"), "\\1", base)
      treatment <- sub(paste0("^", sp, "_([A-Za-z]+\\d+)_", year, "_(observed|backfilled)\\.tif$"), "\\2", base)
      
      r_mean <- terra::rast(f)                         # units: individuals / km^2
      # Always mask to **industry only**; if no industry overlaps this raster, skip
      ind_r  <- tryCatch(terra::mask(terra::crop(r_mean, industry), industry),
                         error = function(e) NULL)
      if (is.null(ind_r) || terra::nlyr(ind_r) == 0) return(NULL)
      
      area_km2 <- terra::cellSize(ind_r, unit = "km")
      r_abund  <- ind_r * area_km2
      
      vals <- terra::extract(r_abund, subbasins, fun = sum, na.rm = TRUE, ID = FALSE) |>
        as_tibble()
      
      tibble(
        species         = sp,
        subbasin        = subbasins$first_HYBAS_ID,
        bcr             = bcr_code,
        treatment       = treatment,
        population_mean = as.numeric(vals[[1]]),
        population_sd   = NA_real_
      )
    }))
  })) |>
    # Merge duplicates across BCRs: sum per species × treatment × subbasin (industry-only)
    dplyr::group_by(species, treatment, subbasin) |>
    dplyr::summarise(population_mean = sum(population_mean, na.rm = TRUE),
                     population_sd   = NA_real_, .groups = "drop")
}
# set root path
root <- "G:/Shared drives/BAM_NationalModels5"
ia_dir <- file.path(root, "data", "Extras", "sandbox_data", "impactassessment_sandbox")

# Example: rebuild for CAWA & OSFL in 2020
res_cawa_osfl_2020_rebuilt <- rebuild_population_from_mean_rasters_industry(c("CAWA","OSFL"), 2020)

pop_summary <- res_cawa_osfl_2020_rebuilt %>%
  group_by(species, treatment) %>%
  summarise(
    total_population = sum(population_mean, na.rm = TRUE),
    .groups = "drop"
  )

pop_summary
