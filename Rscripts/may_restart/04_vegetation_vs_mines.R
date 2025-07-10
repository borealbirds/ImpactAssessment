# ---
# title: Impact Assessment: test if presence of mines is associated with SCANFI "rock" and/or reduced biomass
# author: Mannfred Boehm
# created: June 11, 2025
# ---


#1. attach packages----
library(furrr)
library(progress)
library(terra)
library(tidyverse)

# set up parallel computing
future::plan(multisession, workers = parallel::detectCores() - 2)
options(future.globals.maxSize = 4 * 1024^3)



#2. import and index covariate data----

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



#3. define mine impacts function----
# FUNCTION:
# import raster stack for BCR i x Year j
# import mining data for Year j and crop to BCR i
# for every mine x covariate tuple, calculate pixel values in a  n x n grid surrounding the
# mine center and enter n^2 rows into a dataframe with columns:
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
    
    grid_coords_k <- distinct(expand.grid(row = r_range, col = c_range))
  
    # extract values for covariate m from the grid surrounding mine k
    purrr::map_dfr(vars, function(var_m) {
      
      # check that var_m exists in the current stack
      if (!(var_m %in% names(stack_ij))) return(NULL)
      
      # for testing: var_m <- vars[47]
      # isolate one covariate layer (m) from raster stack ij (BCR x year)
      cov_m <- stack_ij[[var_m]]
      
      # index the cells from the surrounding grid
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
      # combine with covariate values
      tibble(
        year = year,
        bcr = bcr,
        var = var_m,
        mine_id_x = center_xy[1],
        mine_id_y = center_xy[2],
        pixel_distance = dist_vals$pixel_distance,
        value = values_m[[var_m]]) |> 
      distinct()
      
    }) # finish iterating over all covariates for the current mine
        
  }) # finish iterating over all mines in the current BCR x year

  return(bcr_year_df) # keep all mine data for the current BCR x year
  
} # close function


# wrap `estimate_mine_impact` function to return a value instead of an error
safe_estimate_mine_impact <- purrr::possibly(estimate_mine_impact, otherwise = NULL)

# run mine impact estimation in parallel
mine_impacts_list <- furrr::future_map(.x = stacks$file_name, .f = safe_estimate_mine_impact, .progress = TRUE)
saveRDS(mine_impacts_list, file=file.path("C:/Users/mannf/Downloads/mine_impacts_list.rds"))




#4. convert mine impacts list (every element is a bcr x year) into a dataframe----
mine_impacts_df <- 
  mine_impacts_list |> 
  dplyr::bind_rows() |> 
  tidyr::drop_na() |> 
  dplyr::group_by(var) |> 
  dplyr::mutate(value_scaled = (value - min(value, na.rm = TRUE)) /
                  (max(value, na.rm = TRUE) - min(value, na.rm = TRUE))) |> 
  dplyr::ungroup() 


mine_impacts_df_cont <- 
  dplyr::filter(mine_impacts_df, !(var %in% c("SCANFI_1km", "VLCE_1km", "MODISLCC_1km", "MODISLCC_5x5"))) |> 
  tidyr::drop_na()
# saveRDS(mine_impacts_df_cont, file="C:/Users/mannf/Downloads/mine_impacts_df_cont.rds")
mine_impacts_df_cont <- readRDS(file="C:/Users/mannf/Downloads/mine_impacts_df_cont.rds")


mine_impacts_df_cat <- 
  dplyr::filter(mine_impacts_df, var %in% c("SCANFI_1km", "VLCE_1km", "MODISLCC_1km", "MODISLCC_5x5")) |> 
  tidyr::drop_na()
# saveRDS(mine_impacts_df_cat, file="C:/Users/mannf/Downloads/mine_impacts_df_cat.rds")
mine_impacts_df_cat <- readRDS(file="C:/Users/mannf/Downloads/mine_impacts_df_cat.rds")



#5. quantify mine impacts by regression---- 

# estimate slope of distance from mine vs covariate value
# do this for every mine and every year (do not average across years)
quantify_impact_by_mine <- function(impact_data) {
  
  # store the var name for the current covariate data unit 
  # (note: `group_split` is used to run this function in parallel over each covariate)
  var_name <- unique(impact_data$var)   
  
  impact_summary <-
    impact_data |> 
    dplyr::group_by(year, mine_id_x, mine_id_y) |> 
    dplyr::group_modify(~{  # define a function that can iterate on the grouped tibble
      
      df <- .x
      df <- filter(df, !is.na(value_scaled), !is.na(pixel_distance))
      
      if (nrow(df) < 3 || length(unique(df$pixel_distance)) == 1) {
        return(tibble(beta = NA, p_value = NA, adj_r2 = NA))
      }
      
      model <- lm(value_scaled ~ pixel_distance, data = df)
      model_summary <- summary(model)
      tibble(beta = coef(model_summary)["pixel_distance", "Estimate"],
             p_value = coef(model_summary)["pixel_distance", "Pr(>|t|)"],
             adj_r2  = model_summary$adj.r.squared)
      
      }) |> 
    dplyr::ungroup() |> 
    dplyr::group_by(mine_id_x, mine_id_y) |> 
    dplyr::summarise(
      beta_mean = mean(beta, na.rm = TRUE),
      beta_sd = sd(beta, na.rm = TRUE),
      p_mean = mean(p_value, na.rm = TRUE),
      p_sd = sd(p_value, na.rm = TRUE),
      adj_r2_mean = mean(adj_r2, na.rm = TRUE),
      adj_r2_sd = sd(adj_r2, na.rm =TRUE),
      .groups = "drop"
    ) |> 
    dplyr::mutate(var = var_name) 
  
  return(impact_summary)
}

# split dataset by number of covariates (e.g., 39 vars -> 39 units to compute through)
mine_impacts_split <- dplyr::group_split(mine_impacts_df_cont, var)

# fit lm to impact vs pixel distance 
# compute in parallel for every covariate 
mine_impact_summary <- 
  furrr::future_map(.x = mine_impacts_split, .f = quantify_impact_by_mine, .progress = TRUE) |> 
  dplyr::bind_rows() |> 
  dplyr::arrange(desc(abs(beta_mean)))

write_csv(x = mine_impact_summary, file = file.path(getwd(), "/data/derived_data/mine_impact_summary.csv"))


# inspect mines of interest by converting coordinates to EPSG:4326
xy_proj <- vect(matrix(c(-2182500,	3997500), ncol = 2), crs = "EPSG:5070")

# Reproject to decimal degrees (lat/lon, EPSG:4326)
xy_lonlat <- project(xy_proj, "EPSG:4326")

# See the coordinates
crds(xy_lonlat)



#6. plot impact vs distance from mine (for selected mines)----


# filter for some specific mine x year combination 
plot_data <- dplyr::filter(mine_impacts_df_cont, 
                           dplyr::near(mine_id_x, -1754500) & 
                           dplyr::near(mine_id_y, 3151500) & 
                           var == "SCANFILodgepolePine_1km")
  
  
# exported at 1000 x 650
ggplot(plot_data, aes(x = pixel_distance, y = value_scaled, color = factor(year), group = year)) +
  geom_smooth(se = FALSE, method = "loess", formula = y ~ x, span = 0.75, alpha = 0.1) +
  geom_point(alpha = 0.1, size = 2) + 
  #facet_wrap(~var) +
  labs(
    x = "distance from mine (m)",
    y = "scaled covariate value",
    color = "year",
    title = paste0("mine impacts on ", plot_data$var[1], " for mine ", 
                   plot_data$mine_id_x[1], "-", plot_data$mine_id_y[1])
  ) +
  scale_color_viridis_d(option = "plasma") + 
  theme_minimal()




# -----------------------------------------------------------
# estimate slope of distance from mine vs covariate value
# take the average slope across years
mine_impact_summary <- 
  mine_impacts_df_cont |> 
  dplyr::group_by(var, year) |> 
  dplyr::group_modify(~{  # define a function that can iterate on the grouped tibble
    
    df <- .x
    df <- filter(df, !is.na(value_scaled), !is.na(pixel_distance))
    model <- lm(value_scaled ~ pixel_distance, data = df)
    model_summary <- summary(model)
    tibble(beta = coef(model_summary)["pixel_distance", "Estimate"],
           p_value = coef(model_summary)["pixel_distance", "Pr(>|t|)"],
           adj_r2  = model_summary$adj.r.squared)
    
  }) |> 
  dplyr::ungroup() |> 
  dplyr::group_by(var) |> 
  dplyr::summarise(
    beta_mean = mean(beta, na.rm = TRUE),
    beta_sd = sd(beta, na.rm = TRUE),
    p_mean = mean(p_value, na.rm = TRUE),
    p_sd = sd(p_value, na.rm = TRUE),
    adj_r2_mean = mean(adj_r2, na.rm = TRUE),
    adj_r2_sd = sd(adj_r2, na.rm =TRUE),
    .groups = "drop"
  )

saveRDS(mine_impact_summary, file="C:/Users/mannf/Downloads/mine_impact_summary.rds")


