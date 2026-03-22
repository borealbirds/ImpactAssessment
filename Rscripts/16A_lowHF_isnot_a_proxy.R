# ---
# title: Impact Assessment: low HF is not a proxy for pre-industrialization
# author: Mannfred Boehm
# created: December 9, 2025
# ---

library(BAMexploreR)
library(ggrepel)
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
# define abiotic environmental features first so we can subset the raster
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
# import pre-mosaiced covariate stack for year_y (subset to abiotic layers only)
stack_y_full <- terra::rast(file.path(ia_dir, sprintf("covariates_mosaiced_%d.tif", year)))
needed_layers <- intersect(abiotic_vars$predictor, names(stack_y_full))
stack_y <- stack_y_full[[needed_layers]]
rm(stack_y_full)

# import subbasin boundaries and project to current stack
all_subbasins_subset <- terra::vect(file.path(ia_dir, "hydrobasins_masked_merged_subset.gpkg"))
all_subbasins_subset <- terra::project(x=all_subbasins_subset, y=stack_y)

# load HF masks early for pre-screening
lowhf_mask_raw  <- terra::rast(file.path(ia_dir, "hirshpearson", "CanHF_1km_lessthan1.tif"))
highhf_mask_raw <- terra::rast(file.path(ia_dir, "hirshpearson", "CanHF_1km_morethan1.tif"))
lowhf_mask_raw  <- terra::project(lowhf_mask_raw,  stack_y, method = "near")
highhf_mask_raw <- terra::project(highhf_mask_raw, stack_y, method = "near")

# sample subbasins until 16 have >= 100 high-HF and >= 100 low-HF pixels
min_n <- 100
n_target <- 64
set.seed(42)
candidate_order <- sample(seq_len(nrow(all_subbasins_subset)))
subbasin_index  <- integer(0)

for (idx in candidate_order) {
  sub_poly <- all_subbasins_subset[idx]
  lo_sub <- tryCatch(terra::crop(lowhf_mask_raw,  sub_poly) |> terra::mask(sub_poly), error = function(e) NULL)
  hi_sub <- tryCatch(terra::crop(highhf_mask_raw, sub_poly) |> terra::mask(sub_poly), error = function(e) NULL)
  n_lo <- if (!is.null(lo_sub)) terra::global(lo_sub, "notNA")[[1]] else 0
  n_hi <- if (!is.null(hi_sub)) terra::global(hi_sub, "notNA")[[1]] else 0
  if (n_lo >= min_n && n_hi >= min_n) subbasin_index <- c(subbasin_index, idx)
  if (length(subbasin_index) == n_target) break
}
cat(sprintf("Selected %d subbasins (all with >= %d low-HF and >= %d high-HF pixels)\n",
            length(subbasin_index), min_n, min_n))

subbasin_s      <- all_subbasins_subset[subbasin_index]
subbasin_labels <- seq_along(subbasin_index)          # simplified 1..n labels for plotting


# ---------------------------------------------------
# visualize randomly selected subbasins

# get BCR boundaries to orient our subbasins in space
bam_boundary <- terra::vect(file.path(root, "Regions", "BAM_BCR_NationalModel_UnBuffered.shp"))
bam_boundary <- bam_boundary[bam_boundary$subUnit != 23, ]

# get centroids for plotting labels; attach simplified labels as attribute
sub_centroids       <- terra::centroids(all_subbasins_subset[subbasin_index])
sub_centroids$label <- subbasin_labels
centroid_df <- data.frame(
  label = sub_centroids$label,
  x     = terra::crds(sub_centroids)[, 1],
  y     = terra::crds(sub_centroids)[, 2]
)

out_fig_dir <- file.path(getwd(), "output_figures", "lowHF_isnotaproxy")
dir.create(out_fig_dir, recursive = TRUE, showWarnings = FALSE)

p_map <- ggplot() +

  # BCR and basin outlines
  geom_spatvector(data = bam_boundary, fill = NA, colour = "grey") +
  geom_spatvector(data = subbasin_s, fill = NA, color = "#D55E00", linewidth = 0.5) +

  # repelled labels with leader lines
  ggrepel::geom_text_repel(
    data          = centroid_df,
    aes(x = x, y = y, label = label),
    size          = 4,
    max.overlaps  = Inf,
    segment.color = "grey40",
    segment.size  = 0.3,
    box.padding   = 0.4,
    force         = 5) +

  coord_sf(crs = crs(subbasin_s)) +
  theme_minimal()

ggsave(file.path(out_fig_dir, "sampled_subbasins_map.png"),
       p_map, width = 1500, height = 1000, units = "px", dpi = 150)



# ---------------------------------------------------
# fit models to each subbasin (covariates extracted per-subbasin to limit memory)

disturbances <- c("CCNL_1km", "CanHF_1km", "CanHF_5x5", "canroad_1km", "canroad_5x5")

results        <- list()
subbasin_counter <- 0

for (idx_i in seq_along(subbasin_index)) {
  subbasin_counter <- subbasin_counter + 1
  sub_poly <- all_subbasins_subset[subbasin_index[idx_i]]

  # extract covariate stack for this subbasin only
  cov_sub <- tryCatch(
    stack_y |> terra::crop(sub_poly) |> terra::mask(sub_poly),
    error = function(e) NULL
  )
  if (is.null(cov_sub)) { message("Skipping subbasin ", subbasin_counter, " (crop error)"); next }

  lo_mask <- tryCatch(terra::crop(lowhf_mask_raw,  cov_sub) |> terra::mask(sub_poly), error = function(e) NULL)
  hi_mask <- tryCatch(terra::crop(highhf_mask_raw, cov_sub) |> terra::mask(sub_poly), error = function(e) NULL)

  lo_df <- terra::as.data.frame(terra::mask(cov_sub, lo_mask))
  hi_df <- terra::as.data.frame(terra::mask(cov_sub, hi_mask))
  rm(cov_sub, lo_mask, hi_mask)
  gc()

  n_high <- nrow(hi_df)
  n_low  <- nrow(lo_df)

  # build combined df
  df_s <- rbind(data.frame(group = "high", hi_df),
                data.frame(group = "low",  lo_df))
  df_s$group <- factor(df_s$group)
  rm(lo_df, hi_df)

  # remove invariant columns, keep only abiotic, drop disturbances
  inv_cols <- sapply(df_s[, setdiff(names(df_s), "group")], function(x) length(unique(x)) == 1)
  df_s <- df_s[, c(TRUE, !inv_cols)]
  df_s <- dplyr::select(df_s, any_of(c("group", abiotic_vars$predictor)))
  df_s <- dplyr::select(df_s, -any_of(disturbances))

  # are low and high HF pixels environmentally distinct?
  m1 <- ranger::ranger(
    x         = df_s[, setdiff(names(df_s), "group")],
    y         = df_s$group,
    probability = TRUE,
    importance  = "impurity_corrected",
    num.trees   = 100
  )
  top_vars <- head(sort(m1$variable.importance, decreasing = TRUE), 3)

  # re-run without importance for terminal node prediction (ranger advises this)
  m2 <- ranger::ranger(
    x         = df_s[, setdiff(names(df_s), "group")],
    y         = df_s$group,
    probability = TRUE,
    importance  = "none",
    num.trees   = 100
  )

  # get terminal node embeddings
  terminal_nodes <- predict(m2, data = df_s[, setdiff(names(df_s), "group")],
                            type = "terminalNodes")$predictions
  leaf_mat <- Matrix::sparse.model.matrix(~.-1, data = as.data.frame(terminal_nodes))
  pca      <- prcomp(leaf_mat, rank. = 2)
  pca_df   <- data.frame(pc1 = pca$x[, 1], pc2 = pca$x[, 2], group = df_s$group)
  rm(terminal_nodes, leaf_mat, pca)

  # correlations of top variables with PC1 and PC2
  df_s_top <- df_s[, names(top_vars), drop = FALSE]
  cors1 <- sapply(df_s_top, function(v) cor(v, pca_df$pc1, use = "pairwise.complete.obs"))
  cors2 <- sapply(df_s_top, function(v) cor(v, pca_df$pc2, use = "pairwise.complete.obs"))
  arrow_df2      <- data.frame(var = names(top_vars), pc1 = cors1, pc2 = cors2)
  arrow_df2$pc1s <- arrow_df2$pc1 * max(pca_df$pc1)
  arrow_df2$pc2s <- arrow_df2$pc2 * max(pca_df$pc2)

  # overlap metric 1: mean RF P(high) for high-HF pixels
  rf_probs_high  <- predict(m2,
                            data = df_s[df_s$group == "high", setdiff(names(df_s), "group")],
                            type = "response")$predictions[, "high"]
  mean_prob_high <- mean(rf_probs_high, na.rm = TRUE)

  # overlap metric 2: fraction of high-HF pixels outside 95th-pctile KDE envelope of low-HF
  lo_pcs <- pca_df[pca_df$group == "low",  c("pc1", "pc2")]
  hi_pcs <- pca_df[pca_df$group == "high", c("pc1", "pc2")]
  kde_lo <- MASS::kde2d(lo_pcs$pc1, lo_pcs$pc2, n = 100,
                        lims = c(range(pca_df$pc1), range(pca_df$pc2)))
  interp_density <- function(kd, x, y) {
    ix <- findInterval(x, kd$x)
    iy <- findInterval(y, kd$y)
    ix <- pmax(1, pmin(ix, length(kd$x) - 1))
    iy <- pmax(1, pmin(iy, length(kd$y) - 1))
    kd$z[cbind(ix, iy)]
  }
  lo_density_vals <- interp_density(kde_lo, lo_pcs$pc1, lo_pcs$pc2)
  hi_density_vals <- interp_density(kde_lo, hi_pcs$pc1, hi_pcs$pc2)
  kde_threshold   <- quantile(lo_density_vals, 0.05)
  frac_outside    <- mean(hi_density_vals < kde_threshold, na.rm = TRUE)

  results[[as.character(subbasin_index[idx_i])]] <- list(
    subbasin       = subbasin_index[idx_i],
    label          = subbasin_counter,
    n_high         = n_high,
    n_low          = n_low,
    model1         = m1,
    model2         = m2,
    top_vars       = top_vars,
    pca_df         = pca_df,
    arrow_df2      = arrow_df2,
    mean_prob_high = mean_prob_high,
    frac_outside   = frac_outside
  )

  rm(df_s)
  gc()
}


saveRDS(results, file.path(getwd(), "data/derived_data/rds_files/rf_models.rds"))


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
  
  tibble(
    label          = res$label,
    subbasin       = res$subbasin,
    n_high         = res$n_high,
    n_low          = res$n_low,
    oob_error      = m1$prediction.error,
    mean_prob_high = res$mean_prob_high,   # mean RF P(high) for high-HF pixels; 0.5 = chance
    frac_outside   = res$frac_outside      # fraction of high-HF pixels outside 95th-pctile KDE envelope of low-HF
  ) %>%
    bind_cols(top_tbl)
})

saveRDS(summary_df, file.path(getwd(), "data/derived_data/rds_files/rf_summary.rds"))
write.csv(summary_df, file.path(getwd(), "data/derived_data/rf_summary.csv"), row.names = FALSE)

# ---------------------------------------------------
# visualize embedding space per subbasin
plot_per_subbasin <- function(pca_df, arrow_df2, label){

  hull_data <- do.call(rbind, lapply(split(pca_df, pca_df$group), function(g) {
    idx <- grDevices::chull(g$pc1, g$pc2)
    g[idx, ]
  }))

  ggplot(pca_df, aes(x = pc1, y = pc2, fill = group, color = group)) +

    geom_polygon(data = hull_data, alpha = 0.3) +
    scale_fill_manual(values  = c("high" = "#D55E00", "low" = "#56B4E9")) +
    scale_color_manual(values = c("high" = "#D55E00", "low" = "#56B4E9")) +

    geom_segment(data = head(arrow_df2, 1),
                 aes(x = 0, y = 0, xend = pc1s, yend = pc2s),
                 arrow = arrow(length = unit(0.3, "cm")),
                 inherit.aes = FALSE,
                 color = "black", linewidth = 1) +

    geom_text(data = head(arrow_df2, 1),
              aes(x = pc1s, y = pc2s, label = var),
              inherit.aes = FALSE,
              color = "black") +

    ggtitle(label) +

    theme_classic() +
    theme(axis.ticks       = element_blank(),
          axis.text        = element_blank(),
          legend.position  = "none")
}

plots <- lapply(X = results,
                FUN = function(result_i) plot_per_subbasin(result_i$pca_df, result_i$arrow_df2, result_i$label))

for (chunk_i in seq_len(4)) {
  idx_start <- (chunk_i - 1) * 16 + 1
  idx_end   <- chunk_i * 16
  chunk_plots <- plots[idx_start:idx_end]
  fig_chunk   <- patchwork::wrap_plots(chunk_plots, ncol = 4)
  ggsave(
    file.path(out_fig_dir, sprintf("environmental_embedding_space_%d.png", chunk_i)),
    fig_chunk, width = 2000, height = 2000, units = "px", dpi = 150
  )
}

