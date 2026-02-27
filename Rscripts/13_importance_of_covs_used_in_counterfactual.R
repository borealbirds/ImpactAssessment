# ---
# title: Impact Assessment: keep track of percentile importance of covariates used in re-prediction
# author: Mannfred Boehm
# created: February 25, 2025
# ---

library(BAMexploreR)
library(gbm)
library(terra)
library(tidyverse)

# set paths ------------------------------------------------------

# colleague's NationalModels directory — paths here must not change
nm_root <- "/home/mannfred/projects/def-ecknight/NationalModels"

cc    <- TRUE    # TRUE = Compute Canada cluster
local <- FALSE   # TRUE = local RProject machine

if (cc)            { ia_dir <- "/home/mannfred/scratch/impact_assessment" }
if (!cc && local)  { ia_dir <- getwd() }
if (!cc && !local) { ia_dir <- file.path("G:/Shared drives/BAM_NationalModels5", "data", "Extras", "sandbox_data", "impactassessment_sandbox") }


# import data ------------------------------------------------------

# import BCR boundaries
bam_boundary <- terra::vect(file.path(ia_dir, "data", "raw_data", "Regions", "BAM_BCR_NationalModel_Unbuffered.shp"))

# import merged + subsetted subbasins
all_subbasins_subset <- terra::vect(file.path(ia_dir, "data", "raw_data", "hydrobasins_masked_merged_subset.gpkg"))

# import predictor importance from V5 models
bam_predictor_importance_v5 <- readRDS(file.path(ia_dir, "data", "derived_data", "rds_files", "bam_predictor_importance_v5.rds"))

# import model performance metrics for continuous covariates
continuous_holdout_metrics <- readRDS(file.path(ia_dir, "data", "derived_data", "rds_files", "continuous_holdout_metrics.rds"))

# ------------------------------------------------------
# create a reference table for which subbasins are in which BCRs 
# some subbasins will be in multiple BCRs, and that's OK

bcr_subbasins_ref <- {
  
  # logical matrix: rows=subbasins, cols=BCRs
  hits <- terra::relate(centroids(all_subbasins_subset), bam_boundary, relation = "intersects")
  
  # row/col indexes of TRUE
  ij <- which(hits, arr.ind = TRUE)
  
  tibble(
    sub_index = ij[, 1],  
    HYBAS_ID  = all_subbasins_subset$first_HYBAS_ID[ij[, 1]],
    bcr_label = paste(bam_boundary$country[ij[, 2]],
                      bam_boundary$subUnit[ij[, 2]], sep = "_"),
    bcr_code  = gsub("_", "", bcr_label))
}

# predictor importance bookkeeping  ------------------------------------------------------

# build backfilled covariate R^2 x percentile importance table
# 1) subbasin-BCR lookup
subbasin_bcr_lookup <-
  bcr_subbasins_ref |>
  dplyr::select(subbasin = sub_index, bcr = bcr_code) |>
  distinct()

# 2) add BCR info to continuous_holdout_metrics
covariate_backfill_perf <-
  continuous_holdout_metrics |>
  dplyr::select(subbasin, covariate, r2) |>
  left_join(subbasin_bcr_lookup, by = "subbasin")

# 3) get percentile importance per species x BCR
covariate_percentile_importance <-
  bam_predictor_importance_v5 |>
  dplyr::select(
    species = spp,
    bcr,
    covariate = predictor,
    mean_rel_inf
  ) |>
  group_by(species, bcr) |>
  mutate(
    importance_pct = 1 - percent_rank(mean_rel_inf)
  ) |>
  ungroup()

# 4) join percentile importance to backfill performance
covariate_backfill_importance <-
  covariate_percentile_importance |>
  left_join(
    covariate_backfill_perf,
    by = c("bcr", "covariate")
  )

saveRDS(covariate_backfill_importance, file.path(ia_dir, "data", "derived_data", "rds_files", "continuous_holdout_metrics_w_importance.rds"))
