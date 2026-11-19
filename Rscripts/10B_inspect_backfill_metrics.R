# ---
# title: Impact Assessment: inspect model metrics from 10A
# author: Mannfred Boehm
# created: December 19, 2025
# ---

library(terra)
library(tidyverse)

ia_dir     <- getwd()
metric_dir <- file.path(ia_dir, "data", "derived_data", "rds_files", "model_metrics")

continuous_train_metrics    <- readRDS(file.path(metric_dir, "continuous_train_metrics.rds"))
continuous_holdout_metrics  <- readRDS(file.path(metric_dir, "continuous_holdout_metrics.rds"))
categorical_holdout_metrics <- readRDS(file.path(metric_dir, "categorical_holdout_metrics.rds"))

three_cols <- function(d) rep_len(c("#999999", "#E69F00", "#56B4E9"), n_distinct(d$covariate))



# -------------------------------------------------------
# continuous covariates: model performance by covariate

continuous_train_metrics |>
  mutate(covariate = reorder(covariate, r2_median, FUN = median, na.rm = TRUE, decreasing = TRUE)) |>
  ggplot(aes(x = covariate, y = r2_median, fill = covariate, colour = covariate)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.2) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1, shape = 16) +
  scale_fill_manual(values = three_cols(continuous_train_metrics)) +
  scale_colour_manual(values = three_cols(continuous_train_metrics)) +
  labs(x = "landscape feature", y = bquote("Bayesian"~R^2)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 70, hjust = 1), legend.position = "none")

continuous_holdout_metrics |>
  mutate(covariate = reorder(covariate, r2, FUN = median, na.rm = TRUE, decreasing = TRUE)) |>
  ggplot(aes(x = covariate, y = r2, fill = covariate, colour = covariate)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.2, lwd = 0.3) +
  geom_jitter(width = 0.2, alpha = 0.2, size = 1, shape = 16) +
  scale_fill_manual(values = three_cols(continuous_holdout_metrics)) +
  scale_colour_manual(values = three_cols(continuous_holdout_metrics)) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "landscape feature", y = bquote("coefficient of determination"~(R^2))) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 70, hjust = 1), legend.position = "none")



# -------------------------------------------------------
# categorical covariates: holdout accuracy

categorical_holdout_metrics |>
  group_by(covariate) |>
  summarise(mean_accuracy = mean(accuracy, na.rm = TRUE),
            sd_accuracy   = sd(accuracy, na.rm = TRUE),
            n             = n(),
            .groups = "drop")



# -------------------------------------------------------
# is accuracy correlated with the size of the training set?

stack_y <- terra::rast(file.path(ia_dir, "data", "raw_data", "covariates_mosaiced", "covariates_mosaiced_2020.tif"))

lowhf_mask <-
  terra::rast(file.path(ia_dir, "data", "raw_data", "hirshpearson", "CanHF_1km_lessthan1.tif")) |>
  terra::project(stack_y, method = "near")

all_subbasins_subset <-
  terra::vect(file.path(ia_dir, "data", "raw_data", "hydrobasins_masked_merged_subset.gpkg")) |>
  terra::project(stack_y)

for (cov in c("VLCE_1km", "SCANFI_1km", "MODISLCC_1km")) {

  # training pixels per subbasin: low HF, with the covariate observed
  train_ok <- terra::ifel(lowhf_mask == 1 & !is.na(stack_y[[cov]]), 1, NA)
  n_train  <- terra::extract(train_ok, all_subbasins_subset[1:57], fun = sum, na.rm = TRUE)
  names(n_train) <- c("subbasin", "n_train")

  d <- left_join(filter(categorical_holdout_metrics, covariate == cov), n_train, by = "subbasin")

  print(ggplot(d, aes(x = accuracy, y = n_train)) + geom_point() + labs(title = cov) + theme_classic())
  print(summary(lm(accuracy ~ n_train, data = d)))
}
