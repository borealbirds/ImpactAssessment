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

# merge confusion tables across subbasins
subbasin_indices <- c(1:57)
confusion_tables <- lapply(X = subbasin_indices, FUN = import_matrix, year = 2020)
names(confusion_tables) <- paste("subbasin_", subbasin_indices)
confusion_tables_merged <- do.call(rbind, unlist(confusion_tables, recursive = FALSE))

# for every categorical covariate: 
# pool results and create a confusion matrix
categorical_responses = c("ABoVE_1km", "NLCD_1km","MODISLCC_1km", "MODISLCC_5x5","SCANFI_1km","VLCE_1km")

confusion_matrices <- list()

create_cm <- function(cat){
  
  merged_tables <- 
    confusion_tables_merged |> 
    dplyr::filter(covariate == cat)
  
  all_levels <- sort(unique(c(merged_tables$actual, merged_tables$predicted)))
  
  confusion_matrices[[cat]] <- 
    table(factor(merged_tables$actual, levels = all_levels),
          factor(merged_tables$predicted, levels = all_levels))
}

confusion_matrices <- lapply(X = categorical_responses, FUN = create_cm)
names(confusion_matrices) <- categorical_responses

saveRDS(confusion_matrices, file.path(getwd(), "data/derived_data/rds_files/confusion_matrices.rds"))

# -------------------------------------------------------
# inspect and synthesize model metrics

categorical_responses = c("ABoVE_1km", "NLCD_1km","MODISLCC_1km", "MODISLCC_5x5","SCANFI_1km","VLCE_1km")

import_metrics <- function(i, year){
  readRDS(file.path(ia_dir, "bart_models", year, sprintf("subbasin_%d/subbasin_%d_metrics.rds", i, i)))
}

# every RDS file is a list of dataframes containing subbasin metrics
# 
subbasin_indices <- c(1:57)
metrics <- lapply(X = subbasin_indices, FUN = import_metrics, year = 2020)
names(metrics) <- paste("subbasin_", subbasin_indices)

metrics_merged <- tibble(dplyr::bind_rows(metrics))

continuous_train_metrics <- 
  metrics_merged |> 
  dplyr::filter(!(covariate %in% categorical_responses) & split == "train") |> 
  dplyr::select(1:14)
#saveRDS(continuous_train_metrics, file.path(getwd(), "data/derived_data/rds_files/model_metrics/continuous_train_metrics.rds"))

continuous_holdout_metrics <-
  metrics_merged |> 
  dplyr::filter(!(covariate %in% categorical_responses) & split == "holdout") |> 
  dplyr::select(1:3,5:6,13:16)
#saveRDS(file.path(continuous_holdout_metrics, getwd(), "data/derived_data/rds_files/model_metrics/continuous_holdout_metrics.rds"))

categorical_train_metrics <-
  metrics_merged |> 
  dplyr::filter((covariate %in% categorical_responses) & split == "train") |> 
  dplyr::select(1:3,13:15,17:18)
saveRDS(categorical_train_metrics, file.path(getwd(), "data/derived_data/rds_files/model_metrics/categorical_train_metrics.rds"))

categorical_holdout_metrics <-
  metrics_merged |> 
  dplyr::filter((covariate %in% categorical_responses) & split == "holdout") |> 
  dplyr::select(1:3,13:15,17:18)
saveRDS(categorical_holdout_metrics, file.path(getwd(), "data/derived_data/rds_files/model_metrics/categorical_holdout_metrics.rds"))

# -------------------------------------------------------
# plot model performance by covariate
continuous_train_metrics <- readRDS(file.path(getwd(), "data/derived_data/rds_files/model_metrics/continuous_train_metrics.rds"))
continuous_holdout_metrics <- readRDS(file.path(getwd(), "data/derived_data/rds_files/model_metrics/continuous_holdout_metrics.rds"))

continuous_holdout_metrics |> 
  mutate(covariate = reorder(covariate, mae, FUN = median, na.rm = TRUE, decreasing = FALSE)) |> 
  ggplot(aes(x = covariate, y = mae, fill = covariate, colour = covariate)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.2) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1, shape = 16) +
  scale_fill_manual(values = rep_len(c("#999999", "#E69F00", "#56B4E9"), length(unique(continuous_holdout_metrics$covariate)))) +
  scale_colour_manual(values = rep_len(c("#999999", "#E69F00", "#56B4E9"), length(unique(continuous_holdout_metrics$covariate)))) +
  #coord_flip() +
  #coord_cartesian(ylim = c(-2, 1)) + 
  theme_classic() +
  labs(
    x = "landscape feature",
    #y = bquote("coefficient of determination"~(R^2))
    y = "mean absolute error"
  ) +
  theme(axis.text.x = element_text(angle = 70, hjust = 1)) +
  theme(legend.position = "none")
  



# -------------------------------------------------------
# inspect individual backfill rasters

