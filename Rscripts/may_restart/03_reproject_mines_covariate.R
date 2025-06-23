# ---
# title: Impact Assessment: spatialize the MinCan dataset
# author: Mannfred Boehm
# created: June 10, 2025
# ---

# Clara Dalliare-Fortier. A comprehensive historical and geolocalized database of 
# mining activities in Canada. 2025. https://doi.org/10.1038/s41597-024-03116-3


library(terra)
library(tidyverse)


root <- "G:/Shared drives/BAM_NationalModels5"

# define variables of interest
mine_vars <- c("namemine", "latitude", "longitude",  
               "open1", "close1", "open2", "close2", "open3", "close3", 
               "commodity1", "commodity2", "commodity3", "commodity4",
               "commodity5", "commodity6", "commodity7", "commodity8")
  
# here, we set the latest closing date as the only closing date
# this simplifies cases where mines open and close repeatedly, without sufficient time for the landscape to recover
# this approach may not be optimal if e.g. a mine closed just before our analysis begin date (1985), e.g. 1980, and 
# remained closed for *most* of our date range, but then re-opened close to our analysis end date 
# (2020), e.g. 2015. In this case, we'd be assigning a close date of 2015 from 1985-2010, when 1980 would have been
# more accurate for that time period. 
mines_df <- 
  readxl::read_excel(file.path(root, "gis", "other_landscape_covariates", "mincan_dataset.xlsx"), sheet="Data") |> 
  dplyr::select(all_of(mine_vars)) |> 
  dplyr::mutate(open = open1, 
                close = coalesce(close3, close2, close1)) |>  # returns the first non-NA value in a sequence of columns
  dplyr::mutate(close = ifelse(close == "open", 2020, close)) |> 
  dplyr::mutate(close = as.numeric(close)) |> 
  dplyr::select(-c(open1, open2, open3, close1, close2, close3)) # clean up df    

# define a mine impact score function:
# for a given analysis year (1985, 1990,...,2020) this function 
# A) assigns an impact score of 0 if the current year is before the mine opened,
# B) assigns an impact score of 1 if the mine is operational in the current year,
# C) assigns an impact score between 0 and 1 if the mine closed before the current year.
# the impact score is a time decay model exp(-k * (year - close))
# recovery rate is a best case scenario as reported by: https://doi.org/10.3390/su151411287
# r = 0.01 assumes full recovery at 100 years, which is generous
compute_impact <- function(open, close, year, r = 0.01) {
  if (year < open) {
    return(0)
  } else if (year <= close) {
    return(1)
  } else {
    recovery <- r * (year - close) # quantify the amount of recovery
    impact <- pmax(0, 1 - recovery) # estimate remaining impact
    return(impact)
  }
}


# define a year-specific impact raster function:
# for a given analysis year (1985, 1990,...,2020) this function 
# deploys `compute_impact` on the CanMin dataset to estimate the impact of all CanMin mines at year_i
# it then removes mines with zero impact (opened after year_i, or fully recovered)
# and transforms the dataframe to a SpatVect
generate_mine_rasters <- function(mines_df, year) {
  mines_df <- 
    mines_df |>
    dplyr::mutate(impact = mapply(compute_impact, year, open, close)) |>
    dplyr::filter(impact > 0)  
  
  mines_vec <- 
    terra::vect(mines_df, geom = c("longitude", "latitude"), crs = "epsg:4326") |>
    terra::project(x=_, y=template_raster) |> 
    terra::buffer(x=_, width = max(res(template_raster)))
  
  mines_rast <-
    terra::rasterize(mines_vec, template_raster, field = "impact", background = 0, fun = "max") |> 
    terra::crop(x=_, y=bam_boundary) |> 
    terra::mask(x=_, mask=bam_boundary)
    
  return(mines_rast)
}


# import necessary reference data
# set CRS to match BAM data, then crop and mask (some mines are far north of BAM data)
bam_boundary <- terra::vect(file.path(root, "Regions", "BAM_BCR_NationalModel_Buffered.shp"))
template_raster <- terra::rast(file.path(root, "PredictionRasters", "Biomass", "SCANFI", "1km", "SCANFIBalsamFir_1km_2020.tif"))

# generate mines layers for every analysis year
test <- lapply()
  
  

terra::writeVector(mines_vector_mask, file.path(root, "gis", "other_landscape_covariates", "mincan_dataset_masked.gpkg"), overwrite=TRUE)
