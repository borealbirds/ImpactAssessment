# ---
# title: Impact Assessment: low HF is not a proxy for pre-industrialization
# author: Mannfred Boehm
# created: December 9, 2025
# ---

library(BAMexploreR)
library(Matrix)
library(patchwork)
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

# subset to a haphazard sample of subbasins
# see " .png" for subbasin IDs
subbasin_index <- c(101, 670, 609, 113, 94, 587, 426, 10, 328, 264, 188, 274, 589, 154)
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

# create a raster of subbasin IDs to add to each row of dataframes (below)
subbasin_id_rast <- terra::rasterize(subbasin_s, cov_s, field = "first_HYBAS_ID")
cov_s <- c(cov_s, subbasin_id_rast)


# import low hf layer and project/mask to current stack
lowhf_mask <- terra::rast(file.path(ia_dir, "hirshpearson", "CanHF_1km_lessthan1.tif"))
lowhf_mask <- terra::project(x=lowhf_mask, y=cov_s, method = "near")
cov_lowhf  <- terra::mask(cov_s, lowhf_mask) 
cov_lowhf_df <- terra::as.data.frame(cov_lowhf)
saveRDS(cov_lowhf_df, file=file.path(getwd(), "data/derived_data/rds_files/lowhf_covariate_df_random_subbasins_2020.rds"))
cov_lowhf_df <- readRDS(file=file.path(getwd(), "data/derived_data/rds_files/lowhf_covariate_df_random_subbasins_2020.rds"))


# import high hf layer and project/mask to current stack
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
# prepare data for testing environmental features between low and high HF pixels

# combine low and high HF data (start with 88 covariates) 
df <- rbind(cov_highhf_df, cov_lowhf_df)
group <- factor(c(rep("high", nrow(cov_highhf_df)),
                  rep("low",  nrow(cov_lowhf_df))))

# remove invariant features
inv_covs <- sapply(df, function(x) length(unique(x)) == 1)
df2 <- df[, !inv_covs] # 3 covariates removed

# remove variables that aren't in `abiotic_vars`
df3 <- dplyr::select(df2, any_of(c(abiotic_vars$predictor, "first_HYBAS_ID")))

# remove disturbances because they will differ by definition
disturbances <- c("CCNL_1km", "CanHF_1km", "CanHF_5x5", "canroad_1km", "canroad_5x5" )
df4 <- dplyr::select(df3, -all_of(disturbances)) 

# add `group` column to identify low vs high HF rows
#df5 <- data.frame(group = group, df4)
#df5 <- df5[!is.na(df5$first_HYBAS_ID),]
#saveRDS(df5, file.path(getwd(), "data/derived_data/rds_files/combined_covariate_df_random_subbasins_2020.rds"))
df5 <- readRDS(file.path(getwd(), "data/derived_data/rds_files/combined_covariate_df_random_subbasins_2020.rds"))

# no complete cases, so we can't use PCA on raw values (will use RF embeddings)
sum(complete.cases(df5))


# ---------------------------------------------------
# fit models to each subbasin

results <- list()

for(s in unique(df5$first_HYBAS_ID)) {
  
  # subset to this basin
  df_s <- dplyr::filter(df5, first_HYBAS_ID == s)
  
  # get sample size
  n_high <- sum(df_s$group == "high")
  n_low  <- sum(df_s$group == "low")
  
  # are low and high HF pixels environmentally distinct?
  m1 <- ranger::ranger(x = df_s[, setdiff(names(df_s), "group")], 
                    y = df_s$group, 
                    probability = TRUE,
                    importance = "impurity_corrected",
                    num.trees = 100)
  
  # get environmental features driving variance
  top_vars <- head(sort(m1$variable.importance, decreasing = TRUE), 3)

  # now re-run without importance (`predict()` suggests this for some reason)
  m2 <-  ranger::ranger(x = df_s[, setdiff(names(df_s), "group")], 
                      y = df_s$group, 
                      probability = TRUE,
                      importance = "none",
                      num.trees = 100)

  # get terminal node embeddings
  terminal_nodes <- predict(m2, data = df_s, type = "terminalNodes")$predictions

  # PCA on embeddings
  leaf_mat <- Matrix::sparse.model.matrix(~.-1, data = as.data.frame(terminal_nodes))
  pca <- prcomp(leaf_mat, rank. = 2)
  
  # create a dataframe for ggplot
  pca_df <- data.frame(pc1 = pca$x[,1], pc2 = pca$x[,2], group = df_s$group)
  
  # correlations of each top variables with PC1 and PC2
  df_s_top <- df_s[, names(top_vars), drop = FALSE]
  cors1 <- sapply(df_s_top, function(v) cor(v, pca_df$pc1, use="pairwise.complete.obs"))
  cors2 <- sapply(df_s_top, function(v) cor(v, pca_df$pc2, use="pairwise.complete.obs"))
  
  arrow_df2 <- data.frame(var = names(top_vars), pc1 = cors1, pc2 = cors2)
  
  # scale arrows to look reasonable on the PCA plot
  arrow_df2$pc1s <- arrow_df2$pc1 * max(pca_df$pc1)   # scaling factor
  arrow_df2$pc2s <- arrow_df2$pc2 * max(pca_df$pc2)

  # store results
  results[[as.character(s)]] <- list(
    subbasin = s,
    n_high = n_high,
    n_low  = n_low,
    model1 = m1,
    model2 = m2,
    top_vars = top_vars,
    pca_df = pca_df,
    arrow_df2 = arrow_df2
  )
  
}


saveRDS(results, file.path(getwd(), "data/derived_data/rds_files/rf_models.rds"))
#results <-readRDS(file.path(getwd(), "data/derived_data/rds_files/rf_models.rds"))


# summarise results
summary_df <- map_df(results, function(res) {
  
  # extract variables
  m1  <- res$model1
  tvars <- res$top_vars
  
  # convert top_vars named vector into tidy form
  top_tbl <- tibble(
    top1     = names(tvars)[1],
    top1_def = predictor_metadata$definition[match(names(tvars)[1], predictor_metadata$predictor)],
    top1_imp = tvars[1],
    top2     = names(tvars)[2],
    top2_def = predictor_metadata$definition[match(names(tvars)[2], predictor_metadata$predictor)],
    top2_imp = tvars[2],
    top3     = names(tvars)[3],
    top3_def = predictor_metadata$definition[match(names(tvars)[3], predictor_metadata$predictor)],
    top3_imp = tvars[3]
  )
  
  # get the original data for this subbasin
  df_s <- df5[df5$first_HYBAS_ID == s, ]
  
  tibble(
    subbasin = res$subbasin,
    n_high   = res$n_high,  
    n_low    = res$n_low,
    oob_error = m1$prediction.error
  ) %>% 
    bind_cols(top_tbl)
})


saveRDS(summary_df, file.path(getwd(), "data/derived_data/rds_files/rf_summary.rds"))
write.csv(summary_df, file.path(getwd(), "data/derived_data/rf_summary.csv"))

# ---------------------------------------------------
# visualize embedding space per subbasin
plot_per_subbasin <- function(pca_df, arrow_df2, subbasin_id){
  
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
              color = "black") +
    
    ggtitle(subbasin_id) +
    
    theme_classic() +
    theme(axis.ticks = element_blank(),
          axis.text  = element_blank())
}

plots <- lapply(X = results, 
                FUN = function(result_i) plot_per_subbasin(result_i$pca_df, result_i$arrow_df2, result_i$subbasin))

patchwork::wrap_plots(plots, guides = "collect") 

