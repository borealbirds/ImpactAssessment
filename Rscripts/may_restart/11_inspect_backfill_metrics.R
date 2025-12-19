# ---
# title: Impact Assessment: inspect model metrics from "process_backfill_metrics.R"
# author: Mannfred Boehm
# created: December 19, 2025
# ---

root <- "G:/Shared drives/BAM_NationalModels5"
ia_dir <- file.path(root, "data", "Extras", "sandbox_data", "impactassessment_sandbox")

# -------------------------------------------------------
# continuous features: plot model performance by covariate

continuous_train_metrics <- readRDS(file.path(getwd(), "data/derived_data/rds_files/model_metrics/continuous_train_metrics.rds"))

continuous_train_metrics |> 
  mutate(covariate = reorder(covariate, r2_median, FUN = median, na.rm = TRUE, decreasing = TRUE)) |> 
  ggplot(aes(x = covariate, y = r2_median, fill = covariate, colour = covariate)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.2) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1, shape = 16) +
  scale_fill_manual(values = rep_len(c("#999999", "#E69F00", "#56B4E9"), length(unique(continuous_train_metrics$covariate)))) +
  scale_colour_manual(values = rep_len(c("#999999", "#E69F00", "#56B4E9"), length(unique(continuous_train_metrics$covariate)))) +
  theme_classic() +
  labs(
    x = "landscape feature",
    y = bquote("Bayesian"~R^2)
  ) +
  theme(axis.text.x = element_text(angle = 70, hjust = 1)) +
  theme(legend.position = "none")


continuous_holdout_metrics <- readRDS(file.path(getwd(), "data/derived_data/rds_files/model_metrics/continuous_holdout_metrics.rds"))

continuous_holdout_metrics |> 
  mutate(covariate = reorder(covariate, r2, FUN = median, na.rm = TRUE, decreasing = TRUE)) |> 
  ggplot(aes(x = covariate, y = r2, fill = covariate, colour = covariate)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.2) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1, shape = 16) +
  scale_fill_manual(values = rep_len(c("#999999", "#E69F00", "#56B4E9"), length(unique(continuous_holdout_metrics$covariate)))) +
  scale_colour_manual(values = rep_len(c("#999999", "#E69F00", "#56B4E9"), length(unique(continuous_holdout_metrics$covariate)))) +
  coord_cartesian(ylim = c(-3, 1)) + 
  theme_classic() +
  labs(
    x = "landscape feature",
    y = bquote("coefficient of determination"~(R^2))
  ) +
  theme(axis.text.x = element_text(angle = 70, hjust = 1)) +
  theme(legend.position = "none")


# -------------------------------------------------------
# categorical features: inspect confusion matrices and accuracy


categorical_holdout_metrics <- readRDS(file.path(getwd(), "data/derived_data/rds_files/model_metrics/categorical_holdout_metrics.rds"))

# get mean abd SD accuracy per categorical covariate
test1 <- categorical_holdout_metrics |>
  group_by(covariate) |>
  summarise(
    mean_accuracy = mean(accuracy, na.rm = TRUE),
    sd_accuracy   = sd(accuracy, na.rm = TRUE),
    n             = n(),
    .groups = "drop"
  )

# inspect accuracy per subbasin for a given covariate
test2 <- categorical_holdout_metrics |>
  filter(covariate == "VLCE_1km")
group_by(subbasin) |>
  summarise(
    mean_accuracy = mean(accuracy, na.rm = TRUE),
    sd_accuracy   = sd(accuracy, na.rm = TRUE),
    n             = n(),
    .groups = "drop"
  )

# -------------------------------------------------------
# is accuracy correlated with sample size of the training set?

# import covariate stacks
year <- 2020
stack_y <- terra::rast(file.path(ia_dir, sprintf("covariates_mosaiced_%d.tif", year)))

# import training set
lowhf_mask <- terra::rast(file.path(ia_dir, "hirshpearson", "CanHF_1km_lessthan1.tif"))
lowhf_mask <- terra::project(x=lowhf_mask, y=stack_y, method = "near")

# import subbasin boundaries and project to current stack
all_subbasins_subset <- terra::vect(file.path(ia_dir, "hydrobasins_masked_merged_subset.gpkg"))
all_subbasins_subset <- terra::project(x=all_subbasins_subset, y=stack_y)

# VLCE: get lowHF pixels per subbasin
train_ok1   <- terra::ifel(lowhf_mask == 1 & !is.na(stack_y$VLCE_1km), 1, NA)
samplesize1 <- terra::extract(train_ok,  all_subbasins_subset[1:57], fun = sum, na.rm = TRUE)
samplesize1 <- dplyr::rename(samplesize, subbasin = ID)

# VLCE: get accuracy per subbasin 
accuracy_vlce <- dplyr::filter(categorical_holdout_metrics, covariate == "VLCE_1km") 
vlce_df <- left_join(accuracy_vlce, samplesize1, by = "subbasin")

# VLCE: plot sample size vs accuracy 
ggplot(data = vlce_df, mapping = aes(x = accuracy, y = CanHF_1km)) +
  geom_point()+
  theme_classic()

# no effect of sample size on accuracy
summary(lm(accuracy ~ CanHF_1km, data = vlce_df))



# SCANFI: get lowHF pixels per subbasin
train_ok2   <- terra::ifel(lowhf_mask == 1 & !is.na(stack_y$SCANFI_1km), 1, NA)
samplesize2 <- terra::extract(train_ok2,  all_subbasins_subset[1:57], fun = sum, na.rm = TRUE)
samplesize2 <- dplyr::rename(samplesize2, subbasin = ID)

# SCANFI: get accuracy per subbasin 
accuracy_scanfi <- dplyr::filter(categorical_holdout_metrics, covariate == "SCANFI_1km") 
scanfi_df <- left_join(accuracy_scanfi, samplesize2, by = "subbasin")

# for a given covariate, plot sample size vs accuracy 
ggplot(data = scanfi_df, mapping = aes(x = accuracy, y = CanHF_1km)) +
  geom_point()+
  theme_classic()

# no effect of sample size on accuracy
summary(lm(accuracy ~ CanHF_1km, data = scanfi_df))



# MODIS: get lowHF pixels per subbasin
train_ok3   <- terra::ifel(lowhf_mask == 1 & !is.na(stack_y$MODISLCC_1km), 1, NA)
samplesize3 <- terra::extract(train_ok3,  all_subbasins_subset[1:57], fun = sum, na.rm = TRUE)
samplesize3 <- dplyr::rename(samplesize3, subbasin = ID)

# SCANFI: get accuracy per subbasin 
accuracy_modis <- dplyr::filter(categorical_holdout_metrics, covariate == "MODISLCC_1km") 
modis_df <- left_join(accuracy_modis, samplesize3, by = "subbasin")

# for a given covariate, plot sample size vs accuracy 
ggplot(data = modis_df, mapping = aes(x = accuracy, y = CanHF_1km)) +
  geom_point()+
  theme_classic()

# no effect of sample size on accuracy
summary(lm(accuracy ~ CanHF_1km, data = modis_df))
