# ---
# title: Impact Assessment: test if presence of mines is associated with SCANFI "rock" and/or reduced biomass
# author: Mannfred Boehm
# created: June 11, 2025
# ---

library(terra)
library(tidyverse)

root <- "G:/Shared drives/BAM_NationalModels5"

# define covariates of interest (biotic vars affected by mines and disturbance features)
vars <- 
  readr::read_csv(file.path(root, "data", "Extras", "sandbox_data", "impactassessment_sandbox", "nice_var_names_v5.csv")) |> 
  unique() |> 
  dplyr::filter(var_class %in% c("Greenup", "Biomass", "Wetland", "Landcover", "Disturbance", "Road")) |> 
  dplyr::filter(!(var %in% c("WetOccur_1km", "WetOccur_5x5", "WetRecur_1km", "WetSeason_1km", "SurfaceWater_1km", "LakeEdge_1km"))) |> # keep peatland but discard other water variables
  dplyr::filter(!(var %in% c("usroad_1km", "usroad_5x5"))) |> # exclude US covariates
  dplyr::pull(var)

# index covariate stacks
stacks <- 
  list.files(file.path(root, "gis", "stacks"), pattern = "\\.tif$", full.names = TRUE) |> 
  tibble(file_name = _) |> 
  dplyr::filter(grepl(pattern="can", file_name)) |> # keep only Canadian BCRs (no mining data for US) 
  dplyr::filter(!grepl(pattern="1985", file_name)) # exclude 1985 data

# index mining data
mine_years <-
  list.files(file.path(root, "gis", "other_landscape_covariates"), pattern = "nobuffer_masked", full.names = TRUE)


# FUNCTION:
# import raster stack for BCR i x Year j
# import mining data for Year j and crop to BCR i
# for every mine x covariate tuple, calculate pixel values in the 9x9 square surrounding the
# mine center and enter 9 rows into a dataframe with columns:
# year, bcr, var, mine, pixel_distance,  value 
# enter the dataframe as an element in a list of 133 (number of canadian stacks excluding 1985)

# for testing: covariate_stack_path <- stacks$file_name[1]
estimate_mine_impact <- function(covariate_stack_path){
  
  # extract bcr and year metadata from current covariate stack
  file_parts <- str_match(basename(covariate_stack_path), "(can\\d+)_(\\d{4})\\.tif")
  bcr <- file_parts[2]
  year <- as.numeric(file_parts[3])
  
  # load covariate stack_i
  stack_i <- terra::rast(covariate_stack_path)
  
  # only keep covariate names that are present in the stack
  valid_vars <- intersect(vars, names(stack_i))
  stack_i <- stack_i[[valid_vars]]

  # find the matching mining raster by year and load appropriate mine layer
  mine_raster_path <- mine_years[str_detect(mine_years, as.character(year))]
  mine_raster <- 
    terra::rast(mine_raster_path) |> 
    terra::crop(x=_, y=stack_i$SCANFIprcC_1km) |> # choosing a covariate that is likely in every cov stack
    terra::mask(x=_, mask=stack_i$SCANFIprcC_1km) 
  
  
  # get the row and column index of each cell in the raster with mines
  rc <- terra::rowColFromCell(mine_raster, cell = which(values(mine_raster) == 1))
  
}
stack_bcr10_2020 <- terra::rast(file.path(root, "gis", "stacks", "can10_2020.tif"))
stack_bcr10_1990 <- terra::rast(file.path(root, "gis", "stacks", "can10_1990.tif"))

# import 1990 mining data
mines_1990 <- 
  terra::rast(file.path(root, "gis", "other_landscape_covariates", "mincan_mines_1990_masked.tif")) |> 
  terra::crop(x=_, y=ext(stack_bcr10_1990)) |> 
  terra::mask(x=_, mask=stack_bcr10_1990$year)


