# ---
# title: Impact Assessment: low HF is not a proxy for pre-industrialization
# author: Mannfred Boehm
# created: December 9, 2025
# ---

library(BAMexploreR)
library(ranger)
library(terra)
library(tidyterra)
library(tidyverse)

root <- "G:/Shared drives/BAM_NationalModels5"
ia_dir <- file.path(root, "data", "Extras", "sandbox_data", "impactassessment_sandbox")
year <- 2020

# ---------------------------------------------------
# import pre-mosaiced covariate stack for year_y
stack_y <- terra::rast(file.path(ia_dir, sprintf("covariates_mosaiced_%d.tif", year)))

# import subbasin boundaries and project to current stack
all_subbasins_subset <- terra::vect(file.path(ia_dir, "hydrobasins_masked_merged_subset.gpkg"))
all_subbasins_subset <- terra::project(x=all_subbasins_subset, y=stack_y)

# subset to a random sample of subbasins
set.seed(123)
subbasin_index <- sample(1:length(all_subbasins_subset), size = 20)
subbasin_s <- all_subbasins_subset[subbasin_index]

# crop covariate stack to subbasin
cov_s <- 
  stack_y |>
  terra::crop(x = _, y = subbasin_s) |>
  terra::mask(x = _, mask = subbasin_s)


# ---------------------------------------------------
# visualize randomly selected subbasins

# get BCR boundaries to orient our subbasins in space
bam_boundary <- terra::vect(file.path(root, "Regions", "BAM_BCR_NationalModel_UnBuffered.shp"))
bam_boundary <- bam_boundary[bam_boundary$subUnit != 23, ]

# plot (exported at 1000 x 751)
ggplot() +
  
  # BCR and basin outlines
  geom_spatvector(data = bam_boundary, fill = NA, colour = "grey") +
  geom_spatvector(data = subbasin_s, fill = NA, color = "#D55E00", linewidth = 0.3) +
  
  coord_sf(crs = crs(subbasin_s)) +
  theme_minimal()



# ---------------------------------------------------
# prepare dataframes from randomly selected subbasins

# convert stacks from random subbasins to dataframe
cov_df <- terra::as.data.frame(cov_s)
saveRDS(cov_df, file=file.path(getwd(), "data/derived_data/rds_files/covariate_dataframe_random_subbasins_2020.rds"))


# import low hf layer and project to current stack
lowhf_mask <- terra::rast(file.path(ia_dir, "hirshpearson", "CanHF_1km_lessthan1.tif"))
lowhf_mask <- terra::project(x=lowhf_mask, y=cov_s, method = "near")
cov_lowhf  <- terra::mask(cov_s, lowhf_mask) 
cov_lowhf_df <- terra::as.data.frame(cov_lowhf)
saveRDS(cov_lowhf_df, file=file.path(getwd(), "data/derived_data/rds_files/lowhf_covariate_df_random_subbasins_2020.rds"))
cov_lowhf_df <- readRDS(file=file.path(getwd(), "data/derived_data/rds_files/lowhf_covariate_df_random_subbasins_2020.rds"))


# import high hf layer and project to current stack
highhf_mask <- terra::rast(file.path(ia_dir, "hirshpearson", "CanHF_1km_morethan1.tif"))
highhf_mask <- terra::project(x=highhf_mask, y=cov_s, method = "near")
cov_highhf  <- terra::mask(cov_s, highhf_mask) 
cov_highhf_df <- terra::as.data.frame(cov_highhf)
saveRDS(cov_highhf_df, file=file.path(getwd(), "data/derived_data/rds_files/highhf_covariate_df_random_subbasins_2020.rds"))
cov_highhf_df <- readRDS(file=file.path(getwd(), "data/derived_data/rds_files/highhf_covariate_df_random_subbasins_2020.rds"))

# ---------------------------------------------------
# define abiotic environmental features 
# (don't want to compare biotic features between low and high HF pixels)
predictor_metadata <-
  dplyr::tibble(BAMexploreR::predictor_metadata) |>
  dplyr::filter(version == "v5") |>
  dplyr::select(predictor, definition, predictor_class) |>
  dplyr::mutate(dplyr::across('predictor', stringr::str_replace, 'Year', 'year')) |>
  dplyr::mutate(dplyr::across('predictor', stringr::str_replace, 'Method','method'))

soil_covs <- tibble::tibble(predictor = c("cec_0-5cm_mean_1000", "cec_100-200cm_mean_1000",
                                          "cec_15-30cm_mean_1000", "cec_30-60cm_mean_1000",  
                                          "cec_5-15cm_mean_1000", "cec_60-100cm_mean_1000", 
                                          "soc_0-5cm_mean_1000",  "soc_100-200cm_mean_1000",
                                          "soc_15-30cm_mean_1000", "soc_30-60cm_mean_1000",  
                                          "soc_5-15cm_mean_1000", "soc_60-100cm_mean_1000"),
                            predictor_class = rep("Soil Properties", 12))

# convert some abiotic variables to biotic variables
actually_biotic_what <- c("StandardDormancy_1km", "StandardGreenup_1km", "Peatland_5x5", "Peatland_1km")

# define abiotic variables (V5 abiotic + CAfire + soil properties)
abiotic_vars <-
  predictor_metadata |> 
  dplyr::filter(predictor_class %in% c("Annual Climate", "Climate Normals", "Topography", "Wetland", "Disturbance")) |> 
  tibble::add_row(predictor = "CAfire", predictor_class ="Time Since Disturbance") |> 
  dplyr::filter(!(predictor %in% actually_biotic_what)) |> 
  dplyr::bind_rows(soil_covs)


# ---------------------------------------------------
# test for differences in environmental features between datasets

# combine low and high HF data (start with 88 covariates) 
df <- rbind(cov_highhf_df, cov_lowhf_df)
group <- factor(c(rep("high", nrow(cov_highhf_df)),
                  rep("low",  nrow(cov_lowhf_df))))

# remove invariant features
inv_covs <- sapply(df, function(x) length(unique(x)) == 1)
df2 <- df[, !inv_covs] # 3 covariates removed

# remove variables that aren't in `abiotic_vars`
df3 <- dplyr::select(df2, any_of(abiotic_vars$predictor)) 

# remove disturbances because they will differ by definition
disturbances <- c("CCNL_1km", "CanHF_1km", "CanHF_5x5", "canroad_1km", "canroad_5x5" )
df4 <- dplyr::select(df3, -all_of(disturbances)) 

# add `group` column to identify low vs high HF rows
# ensure all remaining covariates are numerics
#df5 <- data.frame(group = group, df4)
#saveRDS(df5, file.path(getwd(), "data/derived_data/rds_files/combined_covariate_df_random_subbasins_2020.rds"))
df5 <- readRDS(file.path(getwd(), "data/derived_data/rds_files/combined_covariate_df_random_subbasins_2020.rds"))

# no complete cases, so we can't use PCA on raw values (will use RF embeddings)
sum(complete.cases(df5))

# are low and high HF pixels environmentally distinct?
m1 <- ranger::ranger(x = df5[, setdiff(names(df5), "group")], 
                    y = df5$group, 
                    probability = TRUE,
                    importance = "impurity_corrected",
                    num.trees = 300)
saveRDS(m1, file.path(getwd(), "data/derived_data/rds_files/ranger_model1_random_subbasins_2020.rds"))

# get environmental features driving variance
top_vars <- head(sort(m1$variable.importance, decreasing = TRUE), 3)

head(sort(m1$variable.importance, decreasing = TRUE))
# ERATavesm_1km  ERATavewt_1km     ERAMAT_1km ERATavesmt_1km 
# 4074.229       3363.591       3291.618       2934.801 
# ERAPPTwt_1km       DD18_1km 
# 2498.760       1808.353 

# now re-run without importance (`predict()` suggests this for some reason)
m2 <-  ranger::ranger(x = df5[, setdiff(names(df5), "group")], 
                      y = df5$group, 
                      probability = TRUE,
                      importance = "none",
                      num.trees = 300)

saveRDS(m2, file.path(getwd(), "data/derived_data/rds_files/ranger_model2_random_subbasins_2020.rds"))



# get terminal node something..embeddings?
set.seed(32)
sample_idx <- sample(nrow(df5), 2000)
terminal_nodes <- predict(m2, data = df5, type = "terminalNodes")$predictions
terminal_nodes_sample <- terminal_nodes[sample_idx,]

# --------------

library(Matrix)
leaf_mat <- Matrix::sparse.model.matrix(~.-1, data = as.data.frame(terminal_nodes_sample))

# PCA on this embedding
pca <- prcomp(leaf_mat, rank. = 2)

# create a dataframe for ggplot
pca_df <- data.frame(
  pc1 = pca$x[,1],
  pc2 = pca$x[,2],
  group = df5$group[sample_idx]
)

ggplot(pca_df, aes(x = pc1, y = pc2, color = group)) +
  geom_point(size = 2, alpha=0.5) +
  scale_color_manual(values = c("#D55E00", "#56B4E9")) +
  theme_classic() +
  labs(color = "group")

# correlations of each variable with PC1 and PC2
cors1 <- apply(df5[sample_idx, setdiff(names(df5), "group")], 2, function(v) cor(v, pca_df$pc1, use="complete.obs"))
cors2 <- apply(df5[sample_idx, setdiff(names(df5), "group")], 2, function(v) cor(v, pca_df$pc2, use="complete.obs"))


arrow_df <- data.frame(var = names(cors1), pc1 = cors1, pc2 = cors2)
arrow_df2 <- arrow_df[arrow_df$var %in% top_vars, ]

# scale arrows to look reasonable on the PCA plot
arrow_df2$pc1s <- arrow_df2$pc1 * 200000   # scaling factor
arrow_df2$pc2s <- arrow_df2$pc2 * 200000

ggplot(pca_df, aes(x = pc1, y = pc2, color = group)) +
  
  geom_point(size = 2, alpha=0.5) +
  scale_color_manual(values = c("#D55E00", "#56B4E9")) +
  labs(color = "group") +
  
  geom_segment(data = arrow_df2,
               aes(x = 0, y = 0, xend = pc1s, yend = pc2s),
               arrow = arrow(length = unit(0.3, "cm")),
               inherit.aes = FALSE,
               color = "black", linewidth = 1) +
  
  geom_text(data = arrow_df2,
            aes(x = pc1s, y = pc2s, label = var),
            color = "black",
            hjust = 0.5, vjust = -0.5) +
  
  theme_classic() 

