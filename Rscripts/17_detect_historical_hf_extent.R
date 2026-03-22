##############################################################################
# 16_detect_historical_hf_extent.R
#
# Detect whether the 2020 human footprint (HF) mask was smaller in earlier
# years by comparing observed vs backfilled biotic covariate values within
# the HF mask using a multivariate (PCA) approach.
#
# Approach:
#   1. Fit PCA on biotic feature space using low-HF pixels (the "natural" reference)
#   2. Build a null distribution of PCA-space distances by simulating BART
#      prediction noise on low-HF pixels (using per-covariate holdout RMSE)
#   3. For each year T with completed backfilling:
#      a. For each subbasin, extract observed and backfilled biotic vectors
#         at HF pixels
#      b. Project both into PCA space, compute Euclidean distance
#      c. Classify pixels as "disturbed" (distance > 95th pctile of null)
#         vs "not yet disturbed"
#   4. Summarise proportion of 2020 HF that existed at each time step
#
# Uses _1km biotic covariates only (drops _5x5 to avoid double-weighting).
#
# Prerequisites:
#   - Backfilling (07_train_and_backfill.R) completed for each year
#   - Covariate stacks (covariates_mosaiced_{year}.tif) available
#   - Per-subbasin metrics files available
##############################################################################

library(terra)
library(tidyverse)

# --- Configuration --------------------------------------------------------
cc    <- FALSE
local <- TRUE

if (cc)            { ia_dir <- "/home/mannfred/scratch/impact_assessment" }
if (!cc && local)  { ia_dir <- getwd() }
if (!cc && !local) { ia_dir <- file.path("G:/Shared drives/BAM_NationalModels5", "data", "Extras", "sandbox_data", "impactassessment_sandbox") }

# Use only _1km continuous biotic covariates (no _5x5 to avoid redundancy)
# Excluded: SCANFIbiomass_1km (91% NA in low-HF), SCANFIclosure_1km (60% NA),
#           SCANFIPonderosaPine_1km (near-zero variance in low-HF)
biotic_covs <- c(
  "Peatland_1km",
  "SCANFIheight_1km",
  "SCANFIheightcv_1km",
  "SCANFIprcC_1km",
  "SCANFIprcD_1km",
  "SCANFIBalsamFir_1km",
  "SCANFIBlackSpruce_1km",
  "SCANFIDouglasFir_1km",
  "SCANFIJackPine_1km",
  "SCANFILodgepolePine_1km",
  "SCANFITamarack_1km",
  "SCANFIWhiteRedPine_1km"
)
# Number of PCs to retain (will be set after PCA based on cumulative variance)
cumvar_threshold <- 0.95  # retain PCs explaining this fraction of variance

# Years to process
years <- seq(1990, 2020, by = 5)

# Number of subbasins
n_subbasins <- 674

# Sample size for PCA fitting and null distribution
n_pca_sample <- 50000
n_null_sims  <- 50000

# --- Paths ----------------------------------------------------------------
raw_dir     <- file.path(ia_dir, "data", "raw_data")
derived_dir <- file.path(ia_dir, "data", "derived_data")
out_dir     <- file.path(derived_dir, "expanding_HF")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

lo_hf_path  <- file.path(raw_dir, "hirshpearson", "CanHF_1km_lessthan1.tif")
hi_hf_path  <- file.path(raw_dir, "hirshpearson", "CanHF_1km_morethan1.tif")
basins_path <- file.path(raw_dir, "hydrobasins_masked_merged_subset.gpkg")

# --- Load spatial references ----------------------------------------------
cat("Loading spatial references...\n")
lo_mask   <- rast(lo_hf_path)
hi_mask   <- rast(hi_hf_path)
subbasins <- vect(basins_path)

# Use the 2020 covariate stack as the reference for PCA fitting
cov_2020_path <- file.path(raw_dir, "covariates_mosaiced", "covariates_mosaiced_2020.tif")
stopifnot(file.exists(cov_2020_path))
cov_2020 <- rast(cov_2020_path)

# =========================================================================
# STEP 1: Fit PCA on low-HF pixels
# =========================================================================
cat("\n=== STEP 1: Fitting PCA on low-HF biotic feature space ===\n")

set.seed(42)

# Sample low-HF pixel coordinates
lo_cells <- cells(lo_mask)
lo_idx   <- lo_cells[sample(length(lo_cells), min(n_pca_sample, length(lo_cells)))]
lo_xy    <- xyFromCell(lo_mask, lo_idx)

# Extract biotic covariates at low-HF locations
lo_biotic <- terra::extract(cov_2020[[biotic_covs]], lo_xy)
lo_biotic <- as.data.frame(lo_biotic)

# Impute remaining NAs with column medians (preserves sample size)
col_medians <- sapply(lo_biotic, median, na.rm = TRUE)
for (v in names(lo_biotic)) {
  lo_biotic[[v]][is.na(lo_biotic[[v]])] <- col_medians[v]
}

# Drop any columns with zero variance (can't scale)
col_sds <- sapply(lo_biotic, sd, na.rm = TRUE)
zero_var <- names(col_sds)[col_sds < 1e-8]
if (length(zero_var) > 0) {
  cat("  Dropping zero-variance columns:", paste(zero_var, collapse = ", "), "\n")
  lo_biotic <- lo_biotic[, !(names(lo_biotic) %in% zero_var), drop = FALSE]
  biotic_covs <- setdiff(biotic_covs, zero_var)
}
cat("  Samples for PCA:", nrow(lo_biotic), "across", ncol(lo_biotic), "covariates\n")

# Fit PCA (centered and scaled)
pca_fit <- prcomp(lo_biotic, center = TRUE, scale. = TRUE)

# Backfill layer names (must be after potential zero-var filtering)
biotic_bf_lyrs <- paste0(biotic_covs, "_mean")

# Determine number of PCs to retain
cumvar  <- cumsum(pca_fit$sdev^2) / sum(pca_fit$sdev^2)
n_pcs   <- which(cumvar >= cumvar_threshold)[1]
cat("  Retaining", n_pcs, "PCs explaining",
    round(cumvar[n_pcs] * 100, 1), "% of variance\n")

# Print PC importance
cat("  Variance explained by each PC:\n")
var_pct <- round(pca_fit$sdev^2 / sum(pca_fit$sdev^2) * 100, 1)
for (j in seq_len(min(n_pcs + 2, length(var_pct)))) {
  cat(sprintf("    PC%d: %.1f%% (cumulative: %.1f%%)\n",
              j, var_pct[j], cumvar[j] * 100))
}

# Print top loadings per retained PC
cat("\n  Top 3 loadings per retained PC:\n")
for (j in seq_len(n_pcs)) {
  loads <- sort(abs(pca_fit$rotation[, j]), decreasing = TRUE)
  top3  <- names(loads)[1:3]
  signs <- ifelse(pca_fit$rotation[top3, j] > 0, "+", "-")
  cat(sprintf("    PC%d: %s\n", j,
              paste(sprintf("%s%s (%.2f)", signs, top3,
                            pca_fit$rotation[top3, j]), collapse = ", ")))
}

# =========================================================================
# STEP 2: Build null distribution using BART holdout RMSE
# =========================================================================
cat("\n=== STEP 2: Building null distribution of PCA distances ===\n")

# Collect per-covariate median holdout RMSE across all subbasins (from 2020)
cat("  Collecting per-covariate holdout RMSE from 2020 backfilling...\n")

rmse_by_cov <- setNames(rep(NA_real_, length(biotic_covs)), biotic_covs)

for (cov_name in biotic_covs) {
  cov_rmses <- numeric(0)
  for (i in seq_len(n_subbasins)) {
    metrics_path <- file.path(derived_dir, "bart_models", 2020,
                              sprintf("subbasin_%d", i),
                              sprintf("subbasin_%d_metrics.rds", i))
    if (!file.exists(metrics_path)) next
    m <- readRDS(metrics_path)
    for (entry in m) {
      if (!is.null(entry$covariate) && entry$covariate == cov_name &&
          !is.null(entry$split) && entry$split == "holdout" &&
          !is.null(entry$rmse) && !is.na(entry$rmse)) {
        cov_rmses <- c(cov_rmses, entry$rmse)
        break
      }
    }
  }
  if (length(cov_rmses) > 0) {
    rmse_by_cov[cov_name] <- median(cov_rmses)
  }
}

cat("  Per-covariate median holdout RMSE:\n")
for (cov_name in biotic_covs) {
  cat(sprintf("    %s: %.4f\n", cov_name, rmse_by_cov[cov_name]))
}

# Simulate null: for low-HF pixels, backfilled ≈ observed + N(0, RMSE²)
# Project both into PCA space, compute distance → null distribution
cat("  Simulating null distances...\n")

# Use the same low-HF sample (already in lo_biotic)
null_sample <- lo_biotic[sample(nrow(lo_biotic), min(n_null_sims, nrow(lo_biotic))), ]

# Scale to PCA space: "observed" projection
obs_scaled  <- scale(null_sample, center = pca_fit$center, scale = pca_fit$scale)
obs_pc      <- obs_scaled %*% pca_fit$rotation[, 1:n_pcs]

# Add noise proportional to holdout RMSE for each covariate
noise <- matrix(0, nrow = nrow(null_sample), ncol = length(biotic_covs))
colnames(noise) <- biotic_covs
for (j in seq_along(biotic_covs)) {
  rmse_j <- rmse_by_cov[biotic_covs[j]]
  if (!is.na(rmse_j)) {
    noise[, j] <- rnorm(nrow(null_sample), mean = 0, sd = rmse_j)
  }
}

bf_simulated <- as.matrix(null_sample) + noise
bf_scaled    <- scale(bf_simulated, center = pca_fit$center, scale = pca_fit$scale)
bf_pc        <- bf_scaled %*% pca_fit$rotation[, 1:n_pcs]

# Euclidean distance in PCA space
null_distances <- sqrt(rowSums((bf_pc - obs_pc)^2))

# Threshold from null distribution
threshold_95 <- quantile(null_distances, 0.95)
threshold_99 <- quantile(null_distances, 0.99)

cat(sprintf("  Null distribution: median = %.3f, 95th = %.3f, 99th = %.3f\n",
            median(null_distances), threshold_95, threshold_99))

# Use 95th percentile as the threshold
pca_threshold <- threshold_95

# =========================================================================
# STEP 3: Process each year — subbasin by subbasin
# =========================================================================

# Project subbasins to match HF mask
subbasins_proj <- project(subbasins, hi_mask)

all_results    <- list()
all_pixel_data <- list()  # store per-pixel results for diagnostics

for (yr in years) {
  cat("\n====== Processing year", yr, "======\n")

  bf_dir <- file.path(derived_dir, "bart_models", yr)
  if (!dir.exists(bf_dir)) {
    cat("  Skipping: no backfill directory for", yr, "\n")
    next
  }

  # Check available subbasins
  available <- sum(dir.exists(
    file.path(bf_dir, sprintf("subbasin_%d", 1:n_subbasins))
  ))
  cat("  Available subbasins:", available, "of", n_subbasins, "\n")
  if (available < 10) {
    cat("  Skipping: too few subbasins\n")
    next
  }

  # Load observed covariate stack for this year
  cov_path <- file.path(raw_dir, "covariates_mosaiced",
                         sprintf("covariates_mosaiced_%d.tif", yr))
  if (!file.exists(cov_path)) {
    cat("  Skipping: covariate stack not found for", yr, "\n")
    next
  }
  obs_stack <- rast(cov_path)

  # Process subbasin by subbasin
  pixel_distances  <- numeric(0)
  pixel_subbasin_v <- integer(0)
  pixel_xy_list    <- list()
  n_processed <- 0
  n_skipped   <- 0

  for (i in seq_len(n_subbasins)) {
    bf_path <- file.path(bf_dir, sprintf("subbasin_%d", i),
                         sprintf("subbasin_%d_backfill.tif", i))
    if (!file.exists(bf_path)) {
      n_skipped <- n_skipped + 1
      next
    }

    # Load backfill raster
    bf_rast <- tryCatch(rast(bf_path), error = function(e) NULL)
    if (is.null(bf_rast)) { n_skipped <- n_skipped + 1; next }

    # Check that all biotic _mean layers are present
    available_lyrs <- intersect(biotic_bf_lyrs, names(bf_rast))
    if (length(available_lyrs) < length(biotic_bf_lyrs) * 0.5) {
      n_skipped <- n_skipped + 1
      next
    }

    # Get the subbasin polygon
    sub_poly <- subbasins_proj[i]

    # Crop HF mask to this subbasin
    hi_sub <- tryCatch({
      terra::crop(hi_mask, sub_poly) |> terra::mask(sub_poly)
    }, error = function(e) NULL)
    if (is.null(hi_sub)) { n_skipped <- n_skipped + 1; next }

    hf_cells <- cells(hi_sub)
    if (length(hf_cells) == 0) { n_skipped <- n_skipped + 1; next }

    hf_xy <- xyFromCell(hi_sub, hf_cells)

    # Extract backfilled biotic values at HF pixels
    bf_vals <- terra::extract(bf_rast[[available_lyrs]], hf_xy)
    # Strip the "_mean" suffix to match covariate names
    colnames(bf_vals) <- sub("_mean$", "", colnames(bf_vals))

    # Extract observed biotic values at same locations
    obs_covs_available <- intersect(sub("_mean$", "", available_lyrs), names(obs_stack))
    obs_vals <- terra::extract(obs_stack[[obs_covs_available]], hf_xy)

    # Align columns: use only covariates present in both
    shared_covs <- intersect(colnames(bf_vals), colnames(obs_vals))
    shared_covs <- intersect(shared_covs, biotic_covs)  # must be in PCA set
    if (length(shared_covs) < length(biotic_covs) * 0.5) {
      n_skipped <- n_skipped + 1
      next
    }

    bf_mat  <- as.matrix(bf_vals[, shared_covs, drop = FALSE])
    obs_mat <- as.matrix(obs_vals[, shared_covs, drop = FALSE])

    # For covariates not in this subbasin's backfill, fill with observed
    # (so PCA distance contribution from those dimensions is zero)
    missing_covs <- setdiff(biotic_covs, shared_covs)
    if (length(missing_covs) > 0) {
      fill_obs <- terra::extract(obs_stack[[missing_covs]], hf_xy)
      fill_mat <- as.matrix(fill_obs)
      bf_full  <- cbind(bf_mat,  fill_mat)[, biotic_covs, drop = FALSE]
      obs_full <- cbind(obs_mat, fill_mat)[, biotic_covs, drop = FALSE]
    } else {
      bf_full  <- bf_mat[, biotic_covs, drop = FALSE]
      obs_full <- obs_mat[, biotic_covs, drop = FALSE]
    }

    # Impute NAs with PCA-fitted column medians (col_medians from step 1)
    for (cv in biotic_covs) {
      bf_full[is.na(bf_full[, cv]), cv]  <- col_medians[cv]
      obs_full[is.na(obs_full[, cv]), cv] <- col_medians[cv]
    }

    # Drop rows where ALL covariates were NA (completely outside data extent)
    valid_rows <- rowSums(is.na(bf_full)) < ncol(bf_full)
    if (sum(valid_rows) == 0) { n_skipped <- n_skipped + 1; next }

    bf_full     <- bf_full[valid_rows, , drop = FALSE]
    obs_full    <- obs_full[valid_rows, , drop = FALSE]
    hf_xy_valid <- hf_xy[valid_rows, , drop = FALSE]

    # Project into PCA space
    obs_sc <- scale(obs_full, center = pca_fit$center, scale = pca_fit$scale)
    bf_sc  <- scale(bf_full,  center = pca_fit$center, scale = pca_fit$scale)

    obs_pcs <- obs_sc %*% pca_fit$rotation[, 1:n_pcs]
    bf_pcs  <- bf_sc  %*% pca_fit$rotation[, 1:n_pcs]

    # Euclidean distance
    dists <- sqrt(rowSums((bf_pcs - obs_pcs)^2))

    pixel_distances  <- c(pixel_distances, dists)
    pixel_subbasin_v <- c(pixel_subbasin_v, rep(i, length(dists)))
    pixel_xy_list[[length(pixel_xy_list) + 1]] <- hf_xy_valid

    n_processed <- n_processed + 1
    if (n_processed %% 50 == 0) {
      cat(sprintf("    Processed %d subbasins (%d pixels so far)...\n",
                  n_processed, length(pixel_distances)))
    }
  }

  cat(sprintf("  Done: %d subbasins processed, %d skipped, %d HF pixels\n",
              n_processed, n_skipped, length(pixel_distances)))

  if (length(pixel_distances) == 0) {
    cat("  Skipping: no valid pixels\n")
    next
  }

  # Step 4: Classify using the null-derived threshold
  already_disturbed <- pixel_distances > pca_threshold

  n_total     <- length(pixel_distances)
  n_disturbed <- sum(already_disturbed)
  pct_dist    <- n_disturbed / n_total * 100

  cat(sprintf("\n  Results for %d:\n", yr))
  cat(sprintf("    Total HF pixels:       %d\n", n_total))
  cat(sprintf("    Already disturbed:      %d (%.1f%%)\n", n_disturbed, pct_dist))
  cat(sprintf("    Not yet disturbed:      %d (%.1f%%)\n",
              n_total - n_disturbed, 100 - pct_dist))
  cat(sprintf("    Mean PCA distance (disturbed):   %.3f\n",
              mean(pixel_distances[already_disturbed], na.rm = TRUE)))
  cat(sprintf("    Mean PCA distance (undisturbed): %.3f\n",
              mean(pixel_distances[!already_disturbed], na.rm = TRUE)))
  cat(sprintf("    Threshold (95th null):           %.3f\n", pca_threshold))

  all_results[[as.character(yr)]] <- tibble(
    year              = yr,
    n_hf_pixels       = n_total,
    n_disturbed       = n_disturbed,
    n_not_disturbed   = n_total - n_disturbed,
    pct_disturbed     = round(pct_dist, 2),
    pct_not_disturbed = round(100 - pct_dist, 2),
    mean_dist         = round(mean(pixel_distances), 4),
    median_dist       = round(median(pixel_distances), 4),
    pca_threshold_95  = round(pca_threshold, 4),
    n_pcs_used        = n_pcs,
    cumvar_explained  = round(cumvar[n_pcs], 4)
  )

  # Save per-pixel distances for this year (for diagnostics / mapping)
  pixel_xy_all <- do.call(rbind, pixel_xy_list)
  pixel_df <- tibble(
    x = pixel_xy_all[, 1],
    y = pixel_xy_all[, 2],
    subbasin   = pixel_subbasin_v,
    pca_dist   = round(pixel_distances, 4),
    disturbed  = already_disturbed
  )
  pixel_out <- file.path(out_dir, sprintf("hf_pixel_distances_%d.csv", yr))
  write_csv(pixel_df, pixel_out)
  cat("  Per-pixel distances saved to:", pixel_out, "\n")
}

# =========================================================================
# STEP 4: Combine and save results
# =========================================================================
if (length(all_results) == 0) {
  cat("\nNo years processed. Ensure backfilling has been completed.\n")
} else {
  summary_df <- bind_rows(all_results)

  cat("\n\n====== HISTORICAL HF EXTENT SUMMARY (PCA) ======\n")
  cat("Biotic covariates used:", length(biotic_covs), "\n")
  cat("PCs retained:", n_pcs, "(", round(cumvar[n_pcs] * 100, 1), "% variance)\n")
  cat("Null threshold (95th pctile):", round(pca_threshold, 3), "\n\n")
  print(summary_df, n = 20)

  out_path <- file.path(out_dir, "historical_hf_extent_pca.csv")
  write_csv(summary_df, out_path)
  cat("\nSummary saved to:", out_path, "\n")

  # Save PCA object and null distribution for reproducibility
  saveRDS(
    list(pca_fit        = pca_fit,
         n_pcs          = n_pcs,
         cumvar         = cumvar,
         biotic_covs    = biotic_covs,
         rmse_by_cov    = rmse_by_cov,
         null_distances = null_distances,
         pca_threshold  = pca_threshold),
    file.path(out_dir, "hf_pca_model.rds")
  )
  cat("PCA model and null saved to:", file.path(out_dir, "hf_pca_model.rds"), "\n")

  # --- Plots ---------------------------------------------------------------
  # Distribution of PCA distances with threshold
  if (exists("pixel_distances") && length(pixel_distances) > 0) {
    dist_plot_df <- tibble(
      distance = c(null_distances, pixel_distances),
      group    = c(rep("Null (BART noise)", length(null_distances)),
                   rep("HF pixels", length(pixel_distances)))
    )

    p1 <- ggplot(dist_plot_df, aes(x = distance, fill = group)) +
      geom_density(alpha = 0.5) +
      geom_vline(xintercept = pca_threshold, linetype = "dashed", colour = "red") +
      annotate("text", x = pca_threshold, y = 0, label = "95th pctile threshold",
               hjust = -0.1, vjust = -1, colour = "red", size = 3.5) +
      scale_fill_manual(values = c("Null (BART noise)" = "grey60",
                                    "HF pixels" = "steelblue")) +
      labs(title = "PCA Distance: Null (BART noise) vs HF Pixels",
           x = "Euclidean distance in PCA space", y = "Density", fill = NULL) +
      theme_minimal(base_size = 12)

    ggsave(file.path(out_dir, "hf_pca_distance_distribution.png"), p1,
           width = 9, height = 5, dpi = 150)
  }

  # Temporal trend
  if (nrow(summary_df) > 1) {
    p2 <- ggplot(summary_df, aes(x = year, y = pct_disturbed)) +
      geom_line(linewidth = 1, colour = "firebrick") +
      geom_point(size = 3, colour = "firebrick") +
      geom_text(aes(label = paste0(pct_disturbed, "%")),
                vjust = -1.2, size = 3.5) +
      scale_x_continuous(breaks = summary_df$year) +
      scale_y_continuous(limits = c(0, 100)) +
      labs(
        title = "Estimated Human Footprint Extent Over Time (PCA)",
        subtitle = sprintf("%d biotic covariates, %d PCs (%.0f%% var), threshold=%.2f",
                           length(biotic_covs), n_pcs, cumvar[n_pcs] * 100, pca_threshold),
        x = "Year",
        y = "% of 2020 HF pixels already disturbed"
      ) +
      theme_minimal(base_size = 13)

    ggsave(file.path(out_dir, "historical_hf_extent_pca.png"), p2,
           width = 9, height = 5, dpi = 150)
  }

  cat("Plots saved.\n")
}
