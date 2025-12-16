# ---
# title: Impact Assessment: train models per subbasin and backfill industry footprints
# author: Mannfred Boehm
# created: December 15, 2025
# ---

library(terra)
library(tidyverse)

root <- "G:/Shared drives/BAM_NationalModels5"
ia_dir <- file.path(root, "data", "Extras", "sandbox_data", "impactassessment_sandbox")

# -------------------------------------------------------
# inspect and synthesize confusion matrices

import_matrix <- function(i, year){
  readRDS(file.path(ia_dir, "bart_models", year, sprintf("subbasin_%d/subbasin_%d_confusion.rds", i, i)))
}

subbasin_indices <- c(3,4,7,57)
cms <- lapply(X = subbasin_indices, FUN = import_matrix, year = 2020)
names(cms) <- paste("subbasin_", subbasin_indices)

# pool results and create a confusion matrix per covariate
VLCE_1km_all <- 
  do.call(rbind, unlist(cms, recursive = FALSE)) |> 
  dplyr::filter(covariate == "VLCE_1km")

all_levels <- sort(unique(c(VLCE_1km_all$actual, VLCE_1km_all$predicted)))

VLCE_1km_cm <- table(factor(VLCE_1km_all$actual,    levels = all_levels),
                    factor(VLCE_1km_all$predicted, levels = all_levels))
VLCE_1km_cm

# -------------------------------------------------------
# inspect and synthesize model metrics

# -------------------------------------------------------
# inspect individual backfill rasters

