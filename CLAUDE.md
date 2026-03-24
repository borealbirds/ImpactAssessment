# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a boreal bird impact assessment pipeline. The goal is to estimate the counterfactual effect of industrial human footprint (HF) on bird populations across Canada. The approach:

1. Use low-HF pixels within each hydrological subbasin to train BART models predicting biotic (vegetation) covariates from abiotic covariates
2. Apply those models to backfill (predict) biotic covariates in high-HF pixels — removing the imprint of industry
3. Re-predict bird densities using the backfilled covariates to estimate populations under a no-industry scenario

## Execution Contexts

Scripts support three execution contexts controlled by `cc` and `local` flags:

```r
cc    <- TRUE   # TRUE = Compute Canada cluster
local <- FALSE  # TRUE = local RProject machine

if (cc)            { ia_dir <- "/home/mannfred/scratch/impact_assessment" }
if (!cc && local)  { ia_dir <- getwd() }
if (!cc && !local) { ia_dir <- file.path("G:/Shared drives/BAM_NationalModels5", "data", "Extras", "sandbox_data", "impactassessment_sandbox") }
```

The cluster and local RProject share the same subdirectory layout (see **Data Directory Structure** below). The Google Drive path is legacy and retains the old flat layout.

**Colleague's NationalModels directory** (`nm_root`): scripts 12A–12C and 13 also reference
`/home/mannfred/projects/def-ecknight/NationalModels` for BRT bootstrap models and BCR
covariate stacks. This path is cluster-only and must not be changed.

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
| `10C_abiotic_extrapolation_diagnostics.R` | cluster | Flag subbasins where BART may be extrapolating (KS + Mahalanobis diagnostics) → `extrapolation_flags.csv` |
| `11_inspect_backfill_metrics.R` | local | Inspect and visualize model accuracy |
| `11_premosaic_backfilled_stacks.R` + `.sh` | cluster | Mosaic per-subbasin BART backfill rasters into BCR-wide stacks (run before 12) |
| `12_observed.R` + `.sh` | cluster | Phase 1: canonical observed-landscape BRT predictions (one job per species, run BEFORE coalition jobs) |
| `12A_repredict_birds.R` + `.sh` | cluster | Phase 2: re-predict bird densities for a given coalition of sectors (SLURM array over species, one job per coalition) |
| `12B_predict_species_bcr.R` | sourced | `predict_species_bcr()`: reads canonical observed bootstraps, builds coalition mask, runs joint BRT×BART sampling, returns subbasin-level density tables |
| `12C_repredict_birds_check.R` | local | Sanity checks on re-prediction outputs (legacy) |
| `13_importance_of_covs_used_in_counterfactual.R` | cluster | Assess percentile importance of backfilled covariates in V5 bird models |
| `15B_sector_attribution.R` | local | Reads coalition density tables, computes exact Shapley values per sector, aggregates bottom-up (subbasin → BCR → national) → `sector_effects/shapley_*.csv` |
| `shapley_utils.R` | sourced | Coalition enumeration, Shapley value computation utilities |
| `misc/` | local | Downstream analysis: population summaries, visualization, vegetation vs. mines comparisons |

## Submitting Cluster Jobs

Backfilling (SLURM array, one job per subbasin index):
```bash
sbatch 07_train_and_backfill.sh
# --array=1-674 for full run; set specific indices to rerun failures
```

Re-predicting birds (two-phase submission):
```bash
# Phase 1: canonical observed predictions (one array job per species)
OBS_JOB_ID=$(sbatch --parsable 12_observed.sh)

# Phase 2: coalition counterfactual predictions (255 jobs, one per non-empty coalition)
# Each coalition job is a SLURM array over species, with dependency on Phase 1
bash 12_repredict_birds.sh $OBS_JOB_ID

# To run only single-sector coalitions (equivalent to old one-at-a-time):
# for CID in 2 3 5 9 17 33 65 129; do
#   sbatch --export=ALL,COALITION_ID=$CID 12_repredict_birds.sh
# done
```

## Key Architecture Decisions

**Biotic variable hierarchy**: Continuous biotic covariates are backfilled in a fixed order (defined in `00_set_backfilling_order.R`). Each covariate's BART model uses all preceding biotic covariates as additional predictors, so backfilled values cascade down the hierarchy.

**Categorical covariates**: Six land-cover layers (`ABoVE_1km`, `NLCD_1km`, `MODISLCC_1km`, `MODISLCC_5x5`, `SCANFI_1km`, `VLCE_1km`) are treated as categorical throughout and use `mbart()` instead of `gbart()`. They are excluded from predicting continuous features.

**Coordinate scaling**: Within each subbasin, lat/lon are centered and scaled before being added as predictors to improve BART performance.

**Output per subbasin**: `data/derived_data/bart_models/{year}/subbasin_{i}/subbasin_{i}_backfill.tif` (mean and SD layers for each covariate), `_metrics.rds`, `_confusion.rds`.

**Re-prediction**: `11_premosaic` mosaics backfilled subbasin rasters into BCR-wide stacks. `12_observed` produces canonical observed bootstrap predictions once per species. `12A/12B` then run coalition-based counterfactual predictions: for a given coalition S of sectors, pixels where any sector in S has footprint (AND CanHF ≥ 1) use backfilled covariates; all other pixels use observed. Joint BART×BRT sampling nests BART posterior draws inside BRT bootstrap iterations.

**Shapley attribution**: 8 sectors → 256 coalitions (2^8). Each coalition is a SLURM job. `15B` computes exact Shapley values from the 256 coalition density tables. Shapley values sum exactly to the total HF impact. `shapley_utils.R` provides coalition enumeration and the Shapley formula.

**Bottom-up aggregation**: The subbasin is the atomic spatial unit. BCR totals = sum of subbasin values within the BCR. National totals = sum of BCR values. Uncertainty propagates under subbasin independence within BCR and BCR independence nationally.

**Abiotic extrapolation diagnostics**: `10C` computes per-subbasin KS statistics and Mahalanobis exceedance rates comparing low-HF vs high-HF abiotic distributions. Flagged subbasins are annotated in `15B` output.

## Key Packages

- `BART`: `gbart()` (Gaussian) and `mbart()` (multinomial) — core backfilling models
- `BAMexploreR`: installed from https://github.com/borealbirds/BAMexploreR. provides `predictor_metadata`, `bam_get_bcr()`
- `terra`: All raster/vector spatial operations (CRS: EPSG:5072, Canada Albers)
- `tidyverse`: Data manipulation throughout

## Data Directory Structure

```
data/
├── raw_data/
│   ├── biotic_variable_hierarchy.rds
│   ├── covariates_mosaiced/
│   │   └── covariates_mosaiced_{year}.tif
│   ├── hirshpearson/
│   │   ├── CanHF_1km_lessthan1.tif
│   │   ├── CanHF_1km_morethan1.tif
│   │   └── {footprint_type}.tif      # built, crop, forestry_harvest, mines, etc.
│   ├── hydrobasins_masked_merged_subset.gpkg
│   └── Regions/
│       └── BAM_BCR_NationalModel_Unbuffered.shp
└── derived_data/
    ├── bart_models/
    │   └── {year}/
    │       └── subbasin_{i}/
    │           ├── subbasin_{i}_backfill.tif
    │           ├── subbasin_{i}_metrics.rds
    │           └── subbasin_{i}_confusion.rds
    ├── bart_models_mosaics/
    │   └── {year}/
    │       └── {bcr_code}_backfilled.tif        # BCR-wide mosaic of backfilled subbasin stacks
    ├── density_tables/
    │   └── {species}_{year}_coalition_{id}.rds   # subbasin-level population arrays per coalition
    ├── predictions/
    │   └── {species}/{bcr_code}/{year}/
    │       ├── observed_bootstraps.tif            # canonical 32-layer bootstrap stack
    │       ├── observed_mean.tif / observed_sd.tif
    ├── predictions_coalitions/
    │   └── {coalition_id}/{species}/{bcr_code}/{year}/
    │       └── backfilled_mean.tif / backfilled_sd.tif
    ├── sector_effects/
    │   ├── shapley_subbasin.csv
    │   ├── shapley_bcr.csv
    │   └── shapley_national.csv
    └── rds_files/
        ├── model_metrics/{year}/
        │   ├── continuous_train_metrics.rds
        │   ├── continuous_holdout_metrics.rds
        │   ├── categorical_train_metrics.rds
        │   └── categorical_holdout_metrics.rds
        ├── accuracy_matrices/{year}/{covariate}.csv
        ├── confusion_matrices/{year}/{covariate}.csv
        ├── bam_predictor_importance_v5.rds
        ├── continuous_holdout_metrics.rds
        └── continuous_holdout_metrics_w_importance.rds
```

Large spatial files (`.tif`, `.gpkg`, `.shp`) and most `.rds` files are gitignored. Versioned outputs are the CSV accuracy/confusion matrices in `data/derived_data/rds_files/`.

## Open Limitations (updated 2026-03-23)

### Conceptual

1. ~~**Abiotic overlap / extrapolation risk**~~: **Addressed** by `10C_abiotic_extrapolation_diagnostics.R`. Per-subbasin KS and Mahalanobis diagnostics flag subbasins where BART extrapolates. Flags are annotated in 15B output. Remaining gap: diagnostics are post-hoc; they don't correct the extrapolation, only flag it.

2. ~~**Sector impacts are not additive**~~: **Addressed** by Shapley value attribution. 12A/12B now run all 2^8 = 256 sector coalitions. 15B computes exact Shapley values that sum to total HF impact. Remaining gap: Shapley values assume the coalition value function v(S) is well-estimated for all S; coalitions with large combined footprints may have more extrapolation risk.

3. **No spatial spillover**: The formula `cf = obs_on_non_coalition + backfilled_on_coalition` assumes removing a coalition's footprint only affects birds on those pixels. Edge effects, area sensitivity, and functional connectivity mean impacts extend beyond the footprint boundary (especially important for linear features like roads and seismic lines).

### Computational

4. ~~**Uncertainty combination is ad hoc**~~: **Addressed** by joint sampling in 12B. Each (bootstrap, scenario) pair draws a fresh BART posterior realization inside the BRT bootstrap loop, capturing the covariance between BART and BRT uncertainty naturally. No post-hoc variance combination.

## Instructions from Masa
1.Always ignore the directory /Rscripts/misc when thinking. It's not immediately relevant to the project.
2.Note that population estimates made by boosted regression trees via `gbm` predict bird density at the hectare scale.
In order to match the scale of our covariate rasters (1km^2), we always multiply the pixel estimates by 100 to get the 1km^2 
density estimate. 
