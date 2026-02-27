# ---
# title: Impact Assessment: process and synthesize model metrics for "inspect_backfill_metrics.R"
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
subbasin_indices <- c(1:674)
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

#saveRDS(confusion_matrices, file.path(getwd(), "data/derived_data/rds_files/confusion_matrices.rds"))
#confusion_matrices <- readRDS(file.path(getwd(), "data/derived_data/rds_files/confusion_matrices.rds"))

# write each confusion matrix as a CSV
for (nm in names(confusion_matrices)) {
  write.csv(confusion_matrices[[nm]], file.path(getwd(), "data/derived_data/rds_files/", paste0(nm, ".csv")), row.names = TRUE)
}

# create companion matrices that give row-normalized prediction accuracy
accuracy_matrices <- lapply(confusion_matrices, function(cm) {
  sweep(cm, 1, rowSums(cm), "/") *100
})

# write each accuracy matrix as a CSV
for (nm in names(accuracy_matrices)) {
  write.csv(accuracy_matrices[[nm]], file.path(getwd(), "data/derived_data/rds_files/", paste0(nm, ".csv")), row.names = TRUE)
}

# -------------------------------------------------------
# inspect and synthesize holdout metrics 

categorical_responses = c("ABoVE_1km", "NLCD_1km","MODISLCC_1km", "MODISLCC_5x5","SCANFI_1km","VLCE_1km")

import_metrics <- function(i, year){
  readRDS(file.path(ia_dir, "bart_models", year, sprintf("subbasin_%d/subbasin_%d_metrics.rds", i, i)))
}

# every RDS file is a list of dataframes containing subbasin metrics
subbasin_indices <- c(1:674)
metrics <- lapply(X = subbasin_indices, FUN = import_metrics, year = 2020)
names(metrics) <- paste("subbasin_", subbasin_indices)

metrics_merged <- tibble(dplyr::bind_rows(metrics))

continuous_train_metrics <- 
  metrics_merged |> 
  dplyr::filter(!(covariate %in% categorical_responses) & split == "train") |> 
  dplyr::select(1:14)
saveRDS(continuous_train_metrics, file.path(getwd(), "data/derived_data/rds_files/model_metrics/continuous_train_metrics.rds"))

continuous_holdout_metrics <-
  metrics_merged |> 
  dplyr::filter(!(covariate %in% categorical_responses) & split == "holdout") |> 
  dplyr::select(1:3,5:6,13:16)
saveRDS(continuous_holdout_metrics, file.path(getwd(), "data/derived_data/rds_files/model_metrics/continuous_holdout_metrics.rds"))

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
# inspect continuous holdout metrics

# how many models have R^2 above some threshold?
continuous_holdout_metrics |> 
  summarise(
    n_total = n(),
    n_above = sum(r2 >= 0.7),
    prop_above = mean(r2 >= 0.7)
  )

# how many models have a rmse/mae <= 1.5?
continuous_holdout_metrics |> 
  filter(r2 >= 0.7) |> 
  summarise(
    n_total = n(),
    mean_ratio = mean(rmse/mae),
    sd_ratio = sd(rmse/mae)
    )

# visualize distribution of rmse/mae
test <- 
  filter(continuous_holdout_metrics) |> 
  mutate(ratio = rmse/mae)

ggplot(test, aes(x = ratio)) + 
  geom_histogram(binwidth = 0.05, fill="#69b3a2", color="#69b3a2", alpha=1) +
  geom_vline(aes(xintercept = mean(ratio, na.rm = TRUE))) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  xlim(c(0.7,11)) +
  labs(x = "RMSE/MAE") + 
  theme_classic()


# -------------------------------------------------------
# inspect categorical holdout metrics


categorical_holdout_metrics |> 
  group_by(covariate) |> 
  summarise(
    mean_accuracy = mean(accuracy),
    sd_accuracy = sd(accuracy)
  )
