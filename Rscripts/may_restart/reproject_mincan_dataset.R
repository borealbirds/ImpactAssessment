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
ia_dir <- file.path(root, "data", "Extras", "sandbox_data", "impactassessment_sandbox")


# define variables of interest from mining dataset
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
  readxl::read_excel(file.path(ia_dir, "mincan_dataset.xlsx"), sheet="Data") |> 
  dplyr::select(all_of(mine_vars)) |> 
  dplyr::mutate(open = open1, 
                close = coalesce(close3, close2, close1)) |>  # returns the first non-NA value in a sequence of columns
  dplyr::mutate(close = ifelse(close == "open", 2020, close)) |> 
  dplyr::mutate(close = as.numeric(close)) |> 
  dplyr::select(-c(open1, open2, open3, close1, close2, close3)) # clean up df    

# define a year-specific impact raster function:
# for a given analysis year (1985, 1990,...,2020) this function 
# determines if a mine was open on or before the current year
# it then removes mines that opened after the current year (from the future)
# and transforms the dataframe to a SpatRaster
generate_mine_rasters <- function(mines_df, year) {
  
  # filter for mines opened on or before the current year
  mine_presence <- 
    dplyr::filter(mines_df, !is.na(open) & open <= year) |> 
    dplyr::mutate(presence = 1)
  
  mines_vec <- 
    terra::vect(mine_presence, geom = c("longitude", "latitude"), crs = "epsg:4326") |>
    terra::project(x=_, y=template_raster) |> 
    terra::buffer(x=_, width = 5*res(template_raster)[1]) # add a 5km buffer around the mine
  
  # exclude mines outside `bam_boundary`
  mines_vec <- mines_vec[bam_boundary, ]
  
  mines_rast <-
    terra::rasterize(mines_vec, template_raster, field = "presence", background = 0, fun = "max") |> 
    terra::crop(x=_, y=bam_boundary) |> 
    terra::mask(x=_, mask=bam_boundary)
  
  # remove buffers that overlap with urbanized areas
  # remove oil and gas buffer pixels that overlap with urban
  # set oil&gas = 0 wherever urban == 1 (handles NA in urban by leaving mines as-is)
  mines_no_urban <- mines_rast * ((urban_reproj != 1) | is.na(urban_reproj))
    
  return(mines_no_urban)
}


# import necessary reference data
# set CRS to match BAM data, then crop and mask (some mines are far north of BAM data)
bam_boundary <- terra::vect(file.path(root, "Regions", "BAM_BCR_NationalModel_UnBuffered.shp"))
template_raster <- terra::rast(file.path(root, "PredictionRasters", "Biomass", "SCANFI", "1km", "SCANFIBalsamFir_1km_2020.tif"))


# import Hirsh-Pearson population density layer
# Creston: 664, Saulte Ste. Marie: 324; Fort Mac 1304
popden <- terra::rast(file.path(ia_dir, "hirshpearson_population_density.tif"))

# define urban mask (TRUE for >=400 persons/km^2)
# Stats Can threshold for urbanized is 400: https://www150.statcan.gc.ca/n1/en/catalogue/92-164-X
# Hirsh-Pearson transformed pop density by: threshold == 3.333 * log10(popden + 1) 
density_threshold <- 3.333 * log10(400 + 1) 
urban <- popden >= density_threshold
urban[urban == 0] <- NA  # keep only TRUE as 1, rest NA

# now that the popden raster is simplified, project
urban_reproj <- 
  terra::project(x = urban, y = template_raster, method = "near") |> 
  terra::crop(x = _, y = bam_boundary) |> 
  terra::mask(x = _, mask = bam_boundary)

urban_reproj<- terra::ifel(urban_reproj, 1, 0)



# generate mines layers for every analysis year
years <- seq(1990, 2020, 5)
mine_rasters <- purrr::map(.x = years, .f = ~generate_mine_rasters(mines_df, .x)) # ~ begins a formula-style anonymous function
names(mine_rasters) <- paste0("mines_", years)
mine_rasters <- purrr::map(.x = mine_rasters, .f = function(x){ terra::varnames(x) <- "mine_presence"; x})

# save each element as a separate raster
# iwalk uses the names of the list elements as .y 
purrr::iwalk(mine_rasters, ~ {
  terra::writeRaster(.x,
                     filename = file.path(ia_dir, paste0("mincan_", .y, "_masked.tif")),
                     overwrite = TRUE)})

