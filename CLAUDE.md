# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a boreal bird impact assessment pipeline. The goal is to estimate the counterfactual effect of industrial human footprint (HF) on bird populations across Canada. The approach:

1. Use low-HF pixels within each hydrological subbasin to train BART models predicting biotic (vegetation) covariates from abiotic covariates
2. Apply those models to backfill (predict) biotic covariates in high-HF pixels — removing the imprint of industry
3. Re-predict bird densities using the backfilled covariates to estimate populations under a no-industry scenario

## Execution Contexts

Scripts toggle between local (Google Drive) and Compute Canada cluster using a `cc` flag:

```r
cc <- TRUE  # cluster
if (!cc) { root <- "G:/Shared drives/BAM_NationalModels5" }
if (cc)  { root <- "/home/mannfred/scratch/impact_assessment" }
ia_dir <- file.path(root, "data", "Extras", "sandbox_data", "impactassessment_sandbox")  # local
ia_dir <- root  # cluster
```

Most prep scripts (01–06) run locally. Compute-heavy scripts (07, 12A) run on the cluster via SLURM array jobs.

## Pipeline

Scripts are numbered in execution order:

| Script | Where | Purpose |
|--------|-------|---------|
| `00_set_backfilling_order.R` | local | Define hierarchy order for biotic covariate backfilling; saves `biotic_variable_hierarchy.rds` |
| `01_reproject_soil_covariates.R` | local | Crop and reproject ISRIC soil data to BAM boundary |
| `02_reproject_disturbance_covariate.R` | local | Prepare CAfire time-since-disturbance layer |
| `03_reproject_humanfootprint.R` | local | Create low-HF (`CanHF_1km_lessthan1.tif`) and high-HF (`CanHF_1km_morethan1.tif`) masks |
| `04_reproject_and_crop_hydrobasins.R` | local | Crop Level 6 HydroBASINS to BAM study area |
| `05_merge_low_density_subbasins.R` | local | Merge data-sparse subbasins so every unit has ≥Q25 low-HF pixels; subset to subbasins with any high-HF pixels → `hydrobasins_masked_merged_subset.gpkg` (674 subbasins) |
| `06_build_covariate_stacks.R` | local | Mosaic BCR covariate stacks by year and add soil + CAfire layers → `covariates_mosaiced_{year}.tif` |
| `07_train_and_backfill.R` + `.sh` | cluster | Entry point: for a given SLURM array index (subbasin), train BART models and backfill high-HF pixels |
| `08A_train_and_backfill_subbasin_s.R` | sourced | Core `train_and_backfill_subbasin_s()` function; loops over biotic covariates in hierarchy order |
| `08B_deploy_gbart.R` | sourced | `deploy_gbart()`: Gaussian BART for continuous biotic covariates (log1p-transformed, 90/10 train/holdout split) |
| `08B_deploy_mbart.R` | sourced | `deploy_mbart()`: Multinomial BART for categorical land-cover covariates |
| `09_collect_metrics_gbart/mbart.R` | sourced | Collect in-sample BART metrics |
| `09_collect_holdout_metrics_gbart/mbart.R` | sourced | Collect holdout BART metrics |
| `10_process_backfill_metrics.R` | local | Aggregate per-subbasin metrics and confusion matrices into CSVs |
| `11_inspect_backfill_metrics.R` | local | Inspect and visualize model accuracy |
| `12A_repredict_birds.R` + `.sh` | cluster | Re-predict bird densities using observed + backfilled covariate stacks |
| `12B_predict_species_bcr.R` | sourced | `predict_species_bcr()`: mosaics backfilled subbasin stacks per BCR and runs BRT density models |
| `12C_repredict_birds_check.R` | local | Sanity checks on re-prediction outputs |
| `13_importance_of_covs_used_in_counterfactual.R` | cluster | Assess percentile importance of backfilled covariates in V5 bird models |
| `misc/` | local | Downstream analysis: population summaries, visualization, vegetation vs. mines comparisons |

## Submitting Cluster Jobs

Backfilling (SLURM array, one job per subbasin index):
```bash
sbatch 07_train_and_backfill.sh
# --array=1-674 for full run; set specific indices to rerun failures
```

Re-predicting birds (SLURM array, one job per species):
```bash
sbatch 12_repredict_birds.sh
# --array=1-N where N = number of species
```

## Key Architecture Decisions

**Biotic variable hierarchy**: Continuous biotic covariates are backfilled in a fixed order (defined in `00_set_backfilling_order.R`). Each covariate's BART model uses all preceding biotic covariates as additional predictors, so backfilled values cascade down the hierarchy.

**Categorical covariates**: Six land-cover layers (`ABoVE_1km`, `NLCD_1km`, `MODISLCC_1km`, `MODISLCC_5x5`, `SCANFI_1km`, `VLCE_1km`) are treated as categorical throughout and use `mbart()` instead of `gbart()`. They are excluded from predicting continuous features.

**Coordinate scaling**: Within each subbasin, lat/lon are centered and scaled before being added as predictors to improve BART performance.

**Output per subbasin**: `bart_models/{year}/subbasin_{i}/subbasin_{i}_backfill.tif` (mean and SD layers for each covariate), `_metrics.rds`, `_confusion.rds`.

**Re-prediction**: `12A` mosaics backfilled subbasin rasters back to BCR-level stacks, then overlays them on the observed covariate stack to produce an "observed" prediction and a "backfilled" (counterfactual) prediction for each species × BCR.

## Key Packages

- `BART`: `gbart()` (Gaussian) and `mbart()` (multinomial) — core backfilling models
- `BAMexploreR`: installed from https://github.com/borealbirds/BAMexploreR. provides `predictor_metadata`, `bam_get_bcr()`
- `terra`: All raster/vector spatial operations (CRS: EPSG:5072, Canada Albers)
- `tidyverse`: Data manipulation throughout

## Data Notes

Large spatial files (`.tif`, `.rds`, `.gpkg`, `.shp`) are gitignored. Versioned outputs are the CSV accuracy/confusion matrices in `data/derived_data/rds_files/`.
