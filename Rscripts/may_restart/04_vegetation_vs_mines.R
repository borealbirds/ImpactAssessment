# ---
# title: Impact Assessment: test if presence of mines is associated with SCANFI "rock" and/or reduced biomass
# author: Mannfred Boehm
# created: June 11, 2025
# ---

library(progress)
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
# enter the dataframe as an element in a list of 133 (number of Canadian BCR stacks excluding 1985)

# for testing: covariate_stack_path <- stacks$file_name[1]
estimate_mine_impact <- function(covariate_stack_path){
  
  # extract bcr and year metadata from current covariate stack
  file_parts <- str_match(basename(covariate_stack_path), "(can\\d+)_(\\d{4})\\.tif")
  bcr <- file_parts[2]
  year <- as.numeric(file_parts[3])
  
  # load covariate stack_i
  stack_ij <- terra::rast(covariate_stack_path)
  
  # only keep covariate names that are present in the stack
  valid_vars <- intersect(vars, names(stack_ij))
  if (length(valid_vars) == 0) return(NULL)
  stack_ij <- stack_ij[[valid_vars]]

  # find the matching mining raster by year and load appropriate mine layer
  mine_raster_path <- mine_years[str_detect(mine_years, as.character(year))]
  mine_raster <- 
    terra::rast(mine_raster_path) |> 
    terra::crop(x=_, y=stack_ij$SCANFIprcC_1km) |> # choosing a covariate that is likely in every cov stack
    terra::mask(x=_, mask=stack_ij$SCANFIprcC_1km) 
  
  
  # get the row and column index (location) of each cell in the raster with mines
  mine_cells <- which(values(mine_raster) == 1)
  rc <- terra::rowColFromCell(mine_raster, cell = mine_cells)
  
  # define a function that searches an arbitrary grid around every mine
  # for the pixel values of every covariate in `var`
  # this function works on the current loaded mine raster (for a given BCR x year tuple)
  # it then calls an internal function that loops through all covariates for that mine
  # the output is a dataframe with columns `year`, `bcr`, `var`, `mine`, `pixel_distance`, `value`
  pb <- progress::progress_bar$new(total = length(mine_cells))
  bcr_year_df <- purrr::map_dfr(seq_along(mine_cells), function(mine_k) {
    
    # display progress
    pb$tick()
    
    # index the row and column (location) for the current mine i
    row <- rc[mine_k, 1]
    col <- rc[mine_k, 2]
    
    # define search window
    r_range <- (row - 10):(row + 10)
    c_range <- (col - 10):(col + 10)
    
    grid_coords_k <- expand.grid(row = r_range, col = c_range)
  
    # extract values for covariate m from the grid surrounding mine k
    purrr::map_dfr(vars, function(var_m) {
      
      # check that var_m exists in the current stack
      if (!(var_m %in% names(stack_ij))) return(NULL)
      
      # for testing: var_m <- vars[47]
      # isolate one covariate layer (m) from raster stack ij (BCR x year)
      cov_m <- stack_ij[[var_m]]
      
      # index the cells from the grid surrounding grid
      cells_k <- terra::cellFromRowCol(cov_m, row = grid_coords_k$row, col = grid_coords_k$col)
      
            # find the pixel values from the grid surrounding mine k
      values_m <- terra::extract(x = cov_m, y = cells_k)
      if (nrow(values_m) == 0) return(NULL)
      
      # define center pixel and convert to lat-long
      center_cell <- terra::cellFromRowCol(object = cov_m, row = row, col = col)
      center_xy <- terra::xyFromCell(cov_m, center_cell)
      
      # calculate distances from center pixel
      # note: Conus Albers is a projected coordinate system designed
      # to represent distances over North America on a flat (2D Euclidean) plane so 
      # we can use the Euclidian distance as the actual distance between pixels in meters
      dist_vals <- 
        terra::xyFromCell(object = cov_m, cell = cells_k) |> 
        dplyr::as_tibble() |> 
        dplyr::mutate(pixel_distance = sqrt((x - center_xy[1])^2 + (y - center_xy[2])^2))
      
      # index distance vs covariate values along with current BCR x year x mine info
      # Combine with covariate values
      tibble(
        year = year,
        bcr = bcr,
        var = var_m,
        mine_id_x = center_xy[1],
        mine_id_y = center_xy[2],
        pixel_distance = dist_vals$pixel_distance,
        value = values_m[[var_m]])
    }) # finish iterating over all covariates for the current mine
        
  }) # finish iterating over all mines in the current BCR x year

  return(bcr_year_df) # keep all mine data for the current BCR x year
  
} # close function


mine_impacts_list <- purrr::map(.x = stacks$file_name, .f = estimate_mine_impact)
combined_df <- bind_rows(results_list)
