# ---
# title: Impact Assessment: process and synthesize model metrics for 10B
# author: Mannfred Boehm
# created: December 15, 2025
# ---

# pools the 674 subbasins' BART metrics and land-cover confusion tables (07):
#   output_tables/{covariate}_confusion.xlsx
#   data/derived_data/rds_files/model_metrics/{continuous,categorical}_{train,holdout}_metrics.rds

library(openxlsx)
library(tidyverse)

ia_dir     <- getwd()
bart_dir   <- file.path(ia_dir, "data", "derived_data", "bart_models", "2020")
metric_dir <- file.path(ia_dir, "data", "derived_data", "rds_files", "model_metrics")

categorical_responses <- c("ABoVE_1km", "NLCD_1km", "MODISLCC_1km", "MODISLCC_5x5", "SCANFI_1km", "VLCE_1km")

read_subbasins <- function(what) {
  lapply(1:674, function(i) readRDS(file.path(bart_dir, sprintf("subbasin_%d/subbasin_%d_%s.rds", i, i, what))))
}



# -------------------------------------------------------
# confusion matrices, pooled over subbasins

confusion <- bind_rows(unlist(read_subbasins("confusion"), recursive = FALSE))

# one sheet per land-cover covariate: cells read "n (row %)", header row and column and the
# diagonal in bold
bold <- createStyle(textDecoration = "bold")

for (cov in unique(confusion$covariate)) {
  d    <- filter(confusion, covariate == cov)
  lvls <- sort(unique(c(d$actual, d$predicted)))
  cm   <- table(factor(d$actual, levels = lvls), factor(d$predicted, levels = lvls))
  acc  <- sweep(cm, 1, rowSums(cm), "/") * 100
  n    <- length(lvls)

  cells <- matrix(ifelse(cm == 0, "0", sprintf("%d (%.1f%%)", cm, acc)), n, n)
  combined <- data.frame(lvls, cells)
  colnames(combined) <- c("actual \\ predicted", lvls)

  wb <- createWorkbook()
  addWorksheet(wb, cov)
  writeData(wb, cov, combined, rowNames = FALSE)
  addStyle(wb, cov, bold, rows = 1, cols = seq_len(n + 1), gridExpand = TRUE)
  addStyle(wb, cov, bold, rows = seq(2, n + 1), cols = 1, gridExpand = TRUE)
  for (k in seq_len(n)) addStyle(wb, cov, bold, rows = k + 1, cols = k + 1)
  saveWorkbook(wb, file.path(ia_dir, "output_tables", paste0(cov, "_confusion.xlsx")), overwrite = TRUE)
}



# -------------------------------------------------------
# training and holdout metrics

metrics <- tibble(bind_rows(read_subbasins("metrics")))

continuous_train_metrics <-
  metrics |>
  filter(!covariate %in% categorical_responses, split == "train") |>
  select(covariate, subbasin, year, ntrain, rmse, mae, r2_median, r2_lo, r2_hi,
         sigma_mean, sigma_sd, coverage95, top_var, split)

continuous_holdout_metrics <-
  metrics |>
  filter(!covariate %in% categorical_responses, split == "holdout") |>
  select(covariate, subbasin, year, rmse, mae, top_var, split, n_holdout, r2)

categorical_train_metrics <-
  metrics |>
  filter(covariate %in% categorical_responses, split == "train") |>
  select(covariate, subbasin, year, top_var, split, n_holdout, n_obs, accuracy)

categorical_holdout_metrics <-
  metrics |>
  filter(covariate %in% categorical_responses, split == "holdout") |>
  select(covariate, subbasin, year, top_var, split, n_holdout, n_obs, accuracy)

saveRDS(continuous_train_metrics,    file.path(metric_dir, "continuous_train_metrics.rds"))
saveRDS(continuous_holdout_metrics,  file.path(metric_dir, "continuous_holdout_metrics.rds"))
saveRDS(categorical_train_metrics,   file.path(metric_dir, "categorical_train_metrics.rds"))
saveRDS(categorical_holdout_metrics, file.path(metric_dir, "categorical_holdout_metrics.rds"))



# -------------------------------------------------------
# continuous covariates: holdout R^2 and RMSE/MAE

continuous_holdout_metrics |>
  summarise(n_total = n(), n_above = sum(r2 >= 0.7), prop_above = mean(r2 >= 0.7))

# RMSE/MAE among the models with R^2 >= 0.7
continuous_holdout_metrics |>
  filter(r2 >= 0.7) |>
  summarise(n_total = n(), mean_ratio = mean(rmse / mae), sd_ratio = sd(rmse / mae))

p_rmse_mae <-
  continuous_holdout_metrics |>
  mutate(ratio = rmse / mae) |>
  ggplot(aes(x = ratio)) +
  geom_histogram(binwidth = 0.05, fill = "#69b3a2", color = "#69b3a2", alpha = 1) +
  geom_vline(aes(xintercept = mean(ratio, na.rm = TRUE))) +
  geom_vline(aes(xintercept = mean(ratio, na.rm = TRUE) - sd(ratio, na.rm = TRUE)), linetype = "dotted") +
  geom_vline(aes(xintercept = mean(ratio, na.rm = TRUE) + sd(ratio, na.rm = TRUE)), linetype = "dotted") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0), breaks = 0:11, limits = c(0.7, 11)) +
  labs(x = "RMSE/MAE") +
  theme_classic()

ggsave(file.path(ia_dir, "output_figures", "rmse_mae_ratio.png"), p_rmse_mae, width = 8, height = 5)



# -------------------------------------------------------
# categorical covariates: holdout accuracy

categorical_holdout_metrics |>
  group_by(covariate) |>
  summarise(mean_accuracy = mean(accuracy), sd_accuracy = sd(accuracy))
