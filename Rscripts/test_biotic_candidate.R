##############################################################################
# test_biotic_candidate.R
#
# Empirically identify which biotic covariate best discriminates high-HF
# from low-HF pixels — i.e., which covariate is most suppressed by human
# footprint and therefore best suited as a proxy for detecting historical
# HF extent via backfilling.
#
# Approach:
#   1. Sample pixels from low-HF and high-HF masks
#   2. For each continuous biotic covariate, compute:
#      - Mean in low-HF vs high-HF pixels
#      - Cohen's d (standardised effect size)
#      - Spearman correlation with continuous HF intensity (CanHF_1km)
#      - Proportion of high-HF pixels with near-zero values
#   3. Rank covariates and produce summary table + plots
##############################################################################

library(terra)
library(tidyverse)

# --- Paths ----------------------------------------------------------------
ia_dir <- getwd()
cov_path   <- file.path(ia_dir, "data/raw_data/covariates_mosaiced/covariates_mosaiced_2020.tif")
lo_hf_path <- file.path(ia_dir, "data/raw_data/hirshpearson/CanHF_1km_lessthan1.tif")
hi_hf_path <- file.path(ia_dir, "data/raw_data/hirshpearson/CanHF_1km_morethan1.tif")

stopifnot(file.exists(cov_path), file.exists(lo_hf_path), file.exists(hi_hf_path))

# --- Load rasters ---------------------------------------------------------
cov_stack <- rast(cov_path)
lo_mask   <- rast(lo_hf_path)
hi_mask   <- rast(hi_hf_path)

# --- Identify continuous biotic covariates in the stack --------------------
hierarchy <- readRDS(file.path(ia_dir, "data/raw_data/biotic_variable_hierarchy.rds"))

# Categorical covariates to exclude (effect size is meaningless)
categorical <- c("ABoVE_1km", "NLCD_1km", "MODISLCC_1km", "MODISLCC_5x5",
                 "SCANFI_1km", "VLCE_1km")

# Keep only continuous biotic covariates that are actually in the stack
biotic_continuous <- setdiff(hierarchy, categorical)
biotic_in_stack   <- intersect(biotic_continuous, names(cov_stack))

cat("Continuous biotic covariates in stack:", length(biotic_in_stack), "\n")
cat(paste(" ", biotic_in_stack), sep = "\n")

# --- Also grab CanHF_1km for correlation ----------------------------------
stopifnot("CanHF_1km" %in% names(cov_stack))

# --- Sample pixels --------------------------------------------------------
# Use the masks to identify cell indices, then sample for tractability
set.seed(42)
n_sample <- 50000

cat("\nIdentifying low-HF and high-HF cell indices...\n")

# Low-HF cells: non-NA cells in each mask
lo_cells <- cells(lo_mask)
hi_cells <- cells(hi_mask)

cat("  Low-HF pixels:", length(lo_cells), "\n")
cat("  High-HF pixels:", length(hi_cells), "\n")

# Sample
lo_idx <- lo_cells[sample(length(lo_cells), min(n_sample, length(lo_cells)))]
hi_idx <- hi_cells[sample(length(hi_cells), min(n_sample, length(hi_cells)))]

# --- Extract values -------------------------------------------------------
cat("Extracting covariate values for sampled pixels...\n")

layers_to_extract <- c(biotic_in_stack, "CanHF_1km")

# The masks and covariate stack may not share the same grid, so use xy coordinates
lo_xy <- xyFromCell(lo_mask, lo_idx)
hi_xy <- xyFromCell(hi_mask, hi_idx)

lo_vals <- terra::extract(cov_stack[[layers_to_extract]], lo_xy)
hi_vals <- terra::extract(cov_stack[[layers_to_extract]], hi_xy)

# Tag source
lo_vals$group <- "low_HF"
hi_vals$group <- "high_HF"

all_vals <- bind_rows(lo_vals, hi_vals)

cat("  Extracted", nrow(lo_vals), "low-HF and", nrow(hi_vals), "high-HF samples\n")

# --- Compute metrics per biotic covariate ---------------------------------
cat("\nComputing discrimination metrics...\n")

results <- map_dfr(biotic_in_stack, function(v) {
  lo <- lo_vals[[v]]
  hi <- hi_vals[[v]]

  # Drop NAs
  lo <- lo[!is.na(lo)]
  hi <- hi[!is.na(hi)]

  if (length(lo) < 100 || length(hi) < 100) {
    return(tibble(covariate = v, n_lo = length(lo), n_hi = length(hi),
                  mean_lo = NA, mean_hi = NA, cohens_d = NA,
                  spearman_r = NA, pct_hi_near_zero = NA))
  }

  # Means
  m_lo <- mean(lo)
  m_hi <- mean(hi)

  # Pooled SD for Cohen's d
  sd_lo <- sd(lo)
  sd_hi <- sd(hi)
  sd_pool <- sqrt(((length(lo) - 1) * sd_lo^2 + (length(hi) - 1) * sd_hi^2) /
                    (length(lo) + length(hi) - 2))

  d <- (m_lo - m_hi) / sd_pool  # positive = higher in low-HF (natural)

  # Spearman correlation with CanHF_1km (across ALL sampled pixels)
  combined_v   <- all_vals[[v]]
  combined_hf  <- all_vals[["CanHF_1km"]]
  valid        <- !is.na(combined_v) & !is.na(combined_hf)
  rho <- cor(combined_v[valid], combined_hf[valid], method = "spearman")

  # Proportion of high-HF pixels with near-zero values
  # "Near zero" = below 5th percentile of the low-HF distribution
  thresh <- quantile(lo, 0.05)
  pct_near_zero <- mean(hi <= thresh)

  tibble(
    covariate       = v,
    n_lo            = length(lo),
    n_hi            = length(hi),
    mean_lo         = round(m_lo, 3),
    mean_hi         = round(m_hi, 3),
    mean_diff       = round(m_lo - m_hi, 3),
    cohens_d        = round(d, 3),
    spearman_r      = round(rho, 3),
    pct_hi_near_zero = round(pct_near_zero, 3)
  )
})

# --- Rank -----------------------------------------------------------------
# Composite score: weight Cohen's d, |spearman_r|, and pct_hi_near_zero
results <- results %>%
  filter(!is.na(cohens_d)) %>%
  mutate(
    rank_d     = rank(-cohens_d),            # higher d = better → lower rank
    rank_rho   = rank(spearman_r),           # more negative = better → lower rank
    rank_zero  = rank(-pct_hi_near_zero),    # higher pct = better → lower rank
    composite  = rank_d + rank_rho + rank_zero
  ) %>%
  arrange(composite)

# --- Print results --------------------------------------------------------
cat("\n=== BIOTIC COVARIATE DISCRIMINATION RANKING ===\n\n")
cat("Positive Cohen's d → higher values in low-HF (natural) pixels\n")
cat("Negative Spearman r → decreases with human footprint intensity\n")
cat("pct_hi_near_zero → fraction of HF pixels below 5th pctile of natural\n\n")

print(results %>% select(covariate, mean_lo, mean_hi, mean_diff, cohens_d,
                          spearman_r, pct_hi_near_zero, composite),
      n = 40)

# --- Save results ---------------------------------------------------------
out_dir <- file.path(ia_dir, "data/derived_data/expanding_HF")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

write_csv(results, file.path(out_dir, "biotic_candidate_ranking.csv"))
cat("\nResults saved to", file.path(out_dir, "biotic_candidate_ranking.csv"), "\n")

# --- Plot top candidates --------------------------------------------------
top_n <- min(10, nrow(results))
top_covs <- results$covariate[1:top_n]

# Density plots: low-HF vs high-HF
plot_data <- all_vals %>%
  select(all_of(c(top_covs, "group"))) %>%
  pivot_longer(-group, names_to = "covariate", values_to = "value") %>%
  filter(!is.na(value)) %>%
  mutate(covariate = factor(covariate, levels = top_covs))

p1 <- ggplot(plot_data, aes(x = value, fill = group)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~covariate, scales = "free", ncol = 2) +
  scale_fill_manual(values = c("low_HF" = "forestgreen", "high_HF" = "firebrick")) +
  labs(title = "Distribution of Top Biotic Covariates: Low-HF vs High-HF Pixels",
       x = "Covariate Value", y = "Density", fill = "Group") +
  theme_minimal(base_size = 11)

ggsave(file.path(out_dir, "biotic_candidate_densities.png"), p1,
       width = 10, height = 12, dpi = 150)

# Bar plot of Cohen's d
p2 <- ggplot(results[1:top_n, ], aes(x = reorder(covariate, cohens_d), y = cohens_d)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(title = "Cohen's d: Low-HF vs High-HF (positive = higher in natural areas)",
       x = NULL, y = "Cohen's d") +
  theme_minimal(base_size = 12)

ggsave(file.path(out_dir, "biotic_candidate_cohens_d.png"), p2,
       width = 8, height = 6, dpi = 150)

cat("Plots saved.\n")
cat("\n=== RECOMMENDATION ===\n")
cat("Best candidate:", results$covariate[1], "\n")
cat("  Cohen's d =", results$cohens_d[1], "\n")
cat("  Spearman r =", results$spearman_r[1], "\n")
cat("  % HF near zero =", results$pct_hi_near_zero[1], "\n")
