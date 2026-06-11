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

**Colleague's NationalModels directory** (`nm_root`): scripts 12B–12C and 13 reference
`/home/mannfred/projects/def-ecknight/NationalModels` for BRT bootstrap models and BCR
covariate stacks. This path is cluster-only and must not be changed. **12A is an exception**:
when `cc=FALSE`, `nm_root` is overridden to `G:/Shared drives/BAM_NationalModels5` so the
script reads Elly's bootstrap models and raw prediction tifs from the local Google Drive.

Most prep scripts (01–06) run locally. Compute-heavy scripts (07, 12B) run on the cluster via SLURM array jobs.

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
| `12A_observed.R` | local | Run locally (not on cluster). Reads Elly's unclamped 32-bootstrap prediction tifs from `G:/Shared drives/BAM_NationalModels5/output/07_predictions/{species}/` and bootstrap model `.Rdata` files from `G:/Shared drives/BAM_NationalModels5/output/06_bootstraps/{species}/`. Applies species-specific clamping (Steps 5-9 of `10.Package.R`) and writes `observed_bootstraps.tif` (32-layer clamped stack), `observed_mean.tif`, and `observed_sd.tif` to `data/derived_data/predictions/{species}/{bcr_code}/{year}/`. After running, Globus-transfer the `observed_bootstraps.tif` files to the same relative path on the cluster before running 12B. Only Canadian BCRs (`can*`) are processed. |
| `12A2_build_prediction_weights.R` + `.sh` | cluster | Build per-species×BCR `weight.tif` (= range membership × not-water × inside-data-limit) replicating V5 `10.Truncate` range/water/extent masking. Reads source masks from `data/raw_data/v5_gis/` (no G: access); grid template is the BCR stack. Run ONCE before 12B. 12C multiplies BOTH observed and backfilled density by this weight, preserving obs/bf symmetry (`w·bf − w·obs = w·(bf − obs)`). |
| `12B_repredict_all_coalitions.R` + `.sh` | cluster | Entry point: re-predict bird densities for ALL 255 coalitions in one job (SLURM array, one task per species: 1=CAWA, 2=OVEN). Sources `12C_predict_species_all_coalitions.R`; writes `density_tables/{species}_{year}_coalition_{cid}.rds` (cid 2..256). Requires `observed_bootstraps.tif` present (hard error if missing — observed rasters are now always staged on the cluster). |
| `12C_predict_species_all_coalitions.R` | sourced | `predict_species_all_coalitions()`: builds the backfilled field ONCE per species×BCR over the all-8-sectors superset, then reduces all 255 coalitions as cheap masked `rowsum`s. Verified bit-identical to the old per-coalition path (since retired). Runs joint BRT×BART sampling; returns subbasin-level density tables. |
| `13_importance_of_covs_used_in_counterfactual.R` | cluster | Assess percentile importance of backfilled covariates in V5 bird models |
| `14B_sector_attribution.R` | local | Reads coalition density tables, computes exact Shapley values per sector, aggregates bottom-up (subbasin → BCR → national) → `sector_effects/shapley_*.csv` |
| `12E_shapley_utils.R` | sourced | Coalition enumeration, Shapley value computation utilities |
| `misc/` | local | Downstream analysis: population summaries, visualization, vegetation vs. mines comparisons |

## Transferring Files to the Cluster (Globus)

Use the Globus CLI at `C:\Users\mannf\AppData\Local\Python\pythoncore-3.14-64\Scripts\globus.exe`.

**IMPORTANT: never use `--batch` mode** — it always returns "permission denied". Submit one `globus transfer` call per file instead.

Endpoint IDs:
- Local: `a7878ccc-747b-11ef-b4b8-8fef73a45f39` (local paths must start with `/C/`)
- Cluster (Fir): `8dec4129-9ab4-451d-a45f-5b4b8471f7a3`

Example (one file at a time):
```powershell
$globus = "C:\Users\mannf\AppData\Local\Python\pythoncore-3.14-64\Scripts\globus.exe"
& $globus transfer `
  "a7878ccc-747b-11ef-b4b8-8fef73a45f39:/C/Users/mannf/Drive/.../Rscripts/12B_repredict_all_coalitions.R" `
  "8dec4129-9ab4-451d-a45f-5b4b8471f7a3:/home/mannfred/scratch/impact_assessment/Rscripts/12B_repredict_all_coalitions.R" `
  --label "my file" --sync-level checksum
```

## Submitting Cluster Jobs

Backfilling (SLURM array, one job per subbasin index):
```bash
sbatch 07_train_and_backfill.sh
# --array=1-674 for full run; set specific indices to rerun failures
```

Re-predicting birds:
```bash
# Phase 1: Run 12A_observed.R LOCALLY (Rscript Rscripts/12A_observed.R, array task 1 then 2).
# It reads from G: drive and writes observed_bootstraps.tif to data/derived_data/predictions/.
# Globus-transfer those tifs to the cluster BEFORE running 12B — they are now a hard
# dependency (12C stops with an error if observed_bootstraps.tif is missing).

# Phase 1b: Build prediction weights ONCE before 12B (range/water/extent masking).
# Reads data/raw_data/v5_gis (staged from G:); writes weight.tif per species x BCR.
sbatch 12A2_build_prediction_weights.sh   # --array=1-<n_species>
# If weight.tif is absent, 12C runs UNMASKED (with a warning) — so this is a
# correctness step, not a hard dependency for the pipeline to execute.

# Phase 2: ONE job per species computes ALL 255 coalitions in a single pass
# (superset restructure — see DESIGN_12C_restructure.md). No per-coalition fan-out,
# no resource tiers: 12B_repredict_all_coalitions.sh fixes --array=1-2 (1=CAWA, 2=OVEN)
# at 384G / 24:00:00 each.
#
# For a FRESH full re-run, delete stale density tables first (the ~510 pre-2026-06-05
# tables are UNMASKED and must be overwritten):
rm -f ../data/derived_data/density_tables/*.rds
sbatch 12B_repredict_all_coalitions.sh    # writes {species}_{year}_coalition_{cid}.rds, cid 2..256

# Smoke test (one species, one BCR, 2 bootstraps):
sbatch --array=1 --time=01:00:00 --mem=192G --export=ALL,TEST_BCR=can60,TEST_N_BOOT=2 12B_repredict_all_coalitions.sh
```

## Fir Cluster Specifications

**Cluster**: Fir, operated by Simon Fraser University (SFU), part of the Digital Research Alliance of Canada.
- Login: `fir.alliancecan.ca`
- Automation (non-interactive scripts): `robot.fir.alliancecan.ca`
- rsync/scp: use login node (dedicated DTN not yet available)
- Globus: collection `alliancecan#fir-globus` (endpoint ID already in Globus section above)

**CPU nodes** (what our jobs use):
- 864 nodes × 192 cores, 750 GB RAM per node
- 2× AMD EPYC 9655 (Zen 5) @ 2.7 GHz; 24 CCDs per node (8 cores/CCD, 8 NUMA nodes)
- Local NVMe: 7.84 TB per node
- Interconnect: InfiniBand NDR, 27:5 blocking over islands of 216 nodes

**SLURM limits and recommended directives (CPU jobs)**:
- Max walltime: `--time=7-00:00:00` (168 h); min 1 h (test jobs: 5 min)
- Max per node: 192 cores, 750 GB RAM → `--mem=750G` or use `--mem-per-cpu`
- Keep threads within one CCD for best cache locality: `--cpus-per-task=8`
- To fill a full node: `--ntasks-per-node=24 --cpus-per-task=8`
- Always specify `--account=def-ecknight` (or the relevant allocation account)
- No partition flag needed for standard CPU jobs (default partition is used)

**Storage** (all on 51 PB DDN Lustre):
- `$HOME` — small quota, cannot be increased; avoid writing job output here
- `$HOME/scratch` (`$SCRATCH`) — large quota, **no backup**, files purged automatically; use for all job I/O (`ia_dir` lives here)
- `$HOME/project/${def-project-id}` — large adjustable quota, daily backup; use for outputs to keep long-term
- All three mounts share the same Lustre filesystem; storage access is non-blocking

**Site-specific policies**:
- `crontab` is not supported — use `robot.fir.alliancecan.ca` for automation
- Compute nodes have full internet access
- SCRATCH files are purged automatically — copy important outputs to PROJECT promptly

## Key Architecture Decisions

**Biotic variable hierarchy**: Continuous biotic covariates are backfilled in a fixed order (defined in `00_set_backfilling_order.R`). Each covariate's BART model uses all preceding biotic covariates as additional predictors, so backfilled values cascade down the hierarchy.

**Categorical covariates**: Six land-cover layers (`ABoVE_1km`, `NLCD_1km`, `MODISLCC_1km`, `MODISLCC_5x5`, `SCANFI_1km`, `VLCE_1km`) are treated as categorical throughout and use `mbart()` instead of `gbart()`. They are excluded from predicting continuous features.

**Coordinate scaling**: Within each subbasin, lat/lon are centered and scaled before being added as predictors to improve BART performance.

**Output per subbasin**: `data/derived_data/bart_models/{year}/subbasin_{i}/subbasin_{i}_backfill.tif` (mean and SD layers for each covariate), `_metrics.rds`, `_confusion.rds`.

**Re-prediction**: `11_premosaic` mosaics backfilled subbasin rasters into BCR-wide stacks. `12A_observed.R` runs locally, reading Elly's unclamped 32-bootstrap prediction tifs from `G:/Shared drives/BAM_NationalModels5/output/07_predictions/` and bootstrap model files from `G:/Shared drives/BAM_NationalModels5/output/06_bootstraps/`, applying species-specific clamping (Steps 5-9 of `10.Package.R`), and writing `observed_bootstraps.tif` (32-layer clamped stack) to `data/derived_data/predictions/{species}/{bcr_code}/{year}/`. These are then Globus-transferred to the cluster, where `observed_bootstraps.tif` is now a hard dependency for `12B/12C` (12C stops with an error if it is missing). `12B_repredict_all_coalitions.R` runs ONE job per species and computes all 255 coalitions in a single pass: it sources `12C_predict_species_all_coalitions.R`, which builds the backfilled field ONCE per species×BCR over the all-8-sectors superset and reduces every coalition as a cheap masked `rowsum`. For a given coalition S of sectors, pixels where any sector in S has footprint (AND CanHF ≥ 1) use backfilled covariates; all other pixels use observed. Joint BART×BRT sampling nests BART posterior draws inside BRT bootstrap iterations. (This superset restructure was verified bit-identical to the older per-coalition `predict_species_bcr()` + `12D_combine` flow, both since retired.)

**Shapley attribution**: 8 sectors → 256 coalitions (2^8; cid 1 = empty is skipped → 255 computed). All 255 are produced by a single SLURM job per species (12B). `14B_sector_attribution.R` computes exact Shapley values from the coalition density tables. Shapley values sum exactly to the total HF impact. `12E_shapley_utils.R` provides coalition enumeration and the Shapley formula.

**Bottom-up aggregation**: The subbasin is the atomic spatial unit. BCR totals = sum of subbasin values within the BCR. National totals = sum of BCR values. Uncertainty propagates under subbasin independence within BCR and BCR independence nationally.

**Abiotic extrapolation diagnostics**: `10C` computes per-subbasin KS statistics and Mahalanobis exceedance rates comparing low-HF vs high-HF abiotic distributions. Flagged subbasins are annotated in `14B` output.

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
│   ├── v5_gis/                          # V5 masking layers (staged from G:, Canada-only)
│   │   ├── ranges/{species}.tif         # continuous range-membership rasters (EPSG:3978)
│   │   ├── WaterMask_Canada.{shp,shx,dbf,prj}
│   │   └── DataLimitationsMask.{shp,shx,dbf,prj,cpg}
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
    │       ├── observed_bootstraps.tif            # canonical 32-layer bootstrap stack (UNweighted)
    │       ├── observed_mean.tif / observed_sd.tif
    │       └── weight.tif                         # range×water×extent weight (built by 12A2)
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

1. ~~**Abiotic overlap / extrapolation risk**~~: **Addressed** by `10C_abiotic_extrapolation_diagnostics.R`. Per-subbasin KS and Mahalanobis diagnostics flag subbasins where BART extrapolates. Flags are annotated in 14B output. Remaining gap: diagnostics are post-hoc; they don't correct the extrapolation, only flag it.

2. ~~**Sector impacts are not additive**~~: **Addressed** by Shapley value attribution. 12B/12C now run all 2^8 = 256 sector coalitions. 14B computes exact Shapley values that sum to total HF impact. Remaining gap: Shapley values assume the coalition value function v(S) is well-estimated for all S; coalitions with large combined footprints may have more extrapolation risk.

3. **No spatial spillover**: The formula `cf = obs_on_non_coalition + backfilled_on_coalition` assumes removing a coalition's footprint only affects birds on those pixels. Edge effects, area sensitivity, and functional connectivity mean impacts extend beyond the footprint boundary (especially important for linear features like roads and seismic lines).

### Computational

4. ~~**Uncertainty combination is ad hoc**~~: **Addressed** by joint sampling in 12C. Each (bootstrap, scenario) pair draws a fresh BART posterior realization inside the BRT bootstrap loop, capturing the covariance between BART and BRT uncertainty naturally. No post-hoc variance combination.

5. **Backfill mosaic coverage is incomplete in some BCRs (UPSTREAM)**: 12B/12C drop footprint pixels whose backfilled design matrix is incomplete (any `_draw_*` covariate NA → `complete.cases` fails). In some BCRs this drops the overwhelming majority of coalition pixels — e.g. CAWA `can10` (verify run 2026-06-10) dropped **98.5–99.7%** of coalition pixels across cid 129/7/256, leaving `bf_on_coalition` based on a tiny remnant and emitting the `backfill mosaic likely degenerate or uncovered` warning. This does **not** affect obs/bf path equivalence (both paths drop the same pixels, so the restructure is still bit-identical), but it makes the counterfactual unreliable wherever coverage is this thin. Root cause is upstream in the backfill/mosaic stages (07 `train_and_backfill`, 08 `deploy_*bart`, 11 `premosaic`) — likely degenerate/flat BART output or uncovered subbasins, not a 12-series bug. Same family as the `12F negative-roads` finding. **Action before trusting per-BCR Shapley numbers**: audit `bart_models_mosaics/{year}/{bcr}_backfilled.tif` coverage (fraction of footprint pixels with a complete `_draw_*` set) per species×BCR; treat BCRs below some coverage floor as flagged. Currently only flagged via the per-BCR runtime warning, not corrected.

## Instructions from Masa
1.Always ignore the directory /Rscripts/misc when thinking. It's not immediately relevant to the project.
2.Note that population estimates made by boosted regression trees via `gbm` predict bird density at the hectare scale.
In order to match the scale of our covariate rasters (1km^2), we always multiply the pixel estimates by 100 to get the 1km^2 
density estimate. 
