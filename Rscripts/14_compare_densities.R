# ---
# title: Summarise observed vs counterfactual bird populations
# author: Mannfred Boehm
# ---
# Reads per-subbasin density tables (one .rds per species x year from 12A) and
# collapses them to a species x BCR summary with population totals and signed
# percent change at three spatial scales:
#
#   footprint  — high-HF pixels only; (bf_only_j - obs_on_bf_j) / obs_on_bf_j
#                per subbasin, averaged (each subbasin equal weight)
#   subbasin   — full subbasin population; (counterfactual_j - obs_total_j) / obs_total_j
#                per subbasin, averaged (each subbasin equal weight)
#   BCR        — full BCR population; ratio computed on BCR totals
#                (size-weighted across subbasins)
#
# All three share the same numerator structure (bf_only_j - obs_on_bf_j).
# Only the denominator grows as we zoom out: obs_on_bf_j → obs_total_j → sum(obs_total_j).
# This guarantees footprint >= subbasin per subbasin (obs_on_bf_j <= obs_total_j),
# while BCR vs subbasin can vary with subbasin size × impact correlation.
# ---

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
})

cc    <- FALSE
local <- TRUE

if (cc)            { ia_dir <- "/home/mannfred/scratch/impact_assessment" }
if (!cc && local)  { ia_dir <- getwd() }
if (!cc && !local) { ia_dir <- file.path("G:/Shared drives/BAM_NationalModels5", "data", "Extras", "sandbox_data", "impactassessment_sandbox") }

year <- 2020

# load and stack all per-subbasin density tables ----------------------------

density_files <- list.files(
  file.path(ia_dir, "data", "derived_data", "density_tables"),
  pattern = "\\.rds$",
  full.names = TRUE
)
density_long <- purrr::map_dfr(density_files, readRDS)

# collapse to species x BCR -------------------------------------------------
# All per-subbasin % changes and their propagated SDs are computed inline.
# Delta-method SD for a per-subbasin ratio r_j = (a_j - b_j) / b_j:
#   SD(r_j) = (1/b_j) * sqrt(sd_a_j^2 + (a_j/b_j)^2 * sd_b_j^2)
# SD of the unweighted mean across n subbasins (assuming independence):
#   SD(mean(r_j)) = sqrt(sum(SD(r_j)^2)) / n

summary_tbl <- density_long |>
  group_by(species, bcr) |>
  summarise(
    n_subbasins = n(),

    # BCR-level population totals
    obs_total               = sum(obs_total_mean),
    obs_total_sd            = sqrt(sum(obs_total_sd^2)),
    counterfactual_total    = sum(counterfactual_mean),
    counterfactual_total_sd = sqrt(sum(counterfactual_sd^2)),

    # footprint scale: mean of (bf_only_j - obs_on_bf_j) / obs_on_bf_j
    pct_change_footprint = mean(
      ifelse(obs_on_bf_mean > 0,
             (bf_only_mean - obs_on_bf_mean) / obs_on_bf_mean * 100,
             NA_real_),
      na.rm = TRUE),
    pct_change_footprint_sd = {
      valid <- obs_on_bf_mean > 0
      sd_j  <- ifelse(valid,
                      (100 / obs_on_bf_mean) * sqrt(
                        bf_only_sd^2 + (bf_only_mean / obs_on_bf_mean)^2 * obs_on_bf_sd^2),
                      NA_real_)
      sqrt(sum(sd_j^2, na.rm = TRUE)) / sum(valid)
    },

    # subbasin scale: mean of (counterfactual_j - obs_total_j) / obs_total_j
    pct_change_subbasin = mean(
      ifelse(obs_total_mean > 0,
             (counterfactual_mean - obs_total_mean) / obs_total_mean * 100,
             NA_real_),
      na.rm = TRUE),
    pct_change_subbasin_sd = {
      valid <- obs_total_mean > 0
      sd_j  <- ifelse(valid,
                      (100 / obs_total_mean) * sqrt(
                        counterfactual_sd^2 + (counterfactual_mean / obs_total_mean)^2 * obs_total_sd^2),
                      NA_real_)
      sqrt(sum(sd_j^2, na.rm = TRUE)) / sum(valid)
    },

    .groups = "drop"
  ) |>
  mutate(
    # BCR scale: ratio of BCR totals (size-weighted across subbasins)
    pct_change_bcr    = (counterfactual_total - obs_total) / obs_total * 100,
    pct_change_bcr_sd = 100 * sqrt(
      (counterfactual_total_sd / obs_total)^2 +
      (counterfactual_total * obs_total_sd / obs_total^2)^2)
  ) |>
  select(
    species, bcr, n_subbasins,
    obs_total, obs_total_sd,
    counterfactual_total, counterfactual_total_sd,
    pct_change_footprint, pct_change_footprint_sd,
    pct_change_subbasin,  pct_change_subbasin_sd,
    pct_change_bcr,       pct_change_bcr_sd
  )

# write outputs -------------------------------------------------------------

out_dir <- file.path(ia_dir, "data", "derived_data", "rds_files")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

saveRDS(summary_tbl, file.path(out_dir, paste0("population_summary_", year, ".rds")))
write.csv(summary_tbl,
          file.path(out_dir, paste0("population_summary_", year, ".csv")),
          row.names = FALSE)

message(Sys.time(), " | wrote population_summary_", year,
        " (", nrow(summary_tbl), " rows)")
summary_tbl
