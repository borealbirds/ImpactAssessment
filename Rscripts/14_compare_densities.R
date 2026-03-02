library(dplyr)
library(tibble)
library(terra)

dir <- "C:/Users/mannf/Downloads"
files <- list.files(dir, pattern = "\\.rds$", full.names = TRUE)

test <- readRDS(files)
test2 <- rast(files)

# Function to summarize a single species
summarize_species <- function(file) {
  
  species_data <- readRDS(file) |> as_tibble()
  
  # compute total mean and SD for each treatment
  summary_tbl <- 
    species_data |> 
    group_by(treatment) |> 
    summarise(
      total_mean = sum(mean),
      total_sd = sqrt(sum(sd^2)),  # sum variances, then sqrt
      .groups = "drop"
    )
  
  # Extract observed
  obs <- summary_tbl |> filter(treatment == "observed")
  
  # Extract backfilled scenarios
  bf_mean  <- summary_tbl |> filter(treatment == "backfilled_mean")
  bf_low   <- summary_tbl |> filter(treatment == "backfilled_low")
  bf_high  <- summary_tbl |> filter(treatment == "backfilled_high")
  
  tibble(
    species = unique(species_data$species),
    
    obs_mean = obs$total_mean,
    obs_sd   = obs$total_sd,
    
    bf_mean = bf_mean$total_mean,
    bf_mean_sd = bf_mean$total_sd,
    population_change_mean = bf_mean$total_mean - obs$total_mean,
    population_change_mean_sd = sqrt(obs$total_sd^2 + bf_mean$total_sd^2),
    
    bf_low = bf_low$total_mean,
    bf_low_sd = bf_low$total_sd,
    population_change_low = bf_low$total_mean - obs$total_mean,
    population_change_low_sd = sqrt(obs$total_sd^2 + bf_low$total_sd^2),
    
    bf_high = bf_high$total_mean,
    bf_high_sd = bf_high$total_sd,
    population_change_high = bf_high$total_mean - obs$total_mean,
    population_change_high_sd = sqrt(obs$total_sd^2 + bf_high$total_sd^2)
  )
}

# Apply to all species
results <- map_dfr(files, summarize_species)

results
