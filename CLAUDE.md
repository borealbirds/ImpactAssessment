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

**Colleague's NationalModels directory** (`nm_root`): scripts 06 (on Fir), 11, 12A, 12C, 12D–12F and 13 reference
`/home/mannfred/projects/def-ecknight/NationalModels` for BRT bootstrap models, prediction tifs
and BCR covariate stacks. This path is cluster-only and must not be changed. **12A is an
exception**: it sets `cc` from `SLURM_JOB_ID`, and run locally it overrides `nm_root` to
`G:/Shared drives/BAM_NationalModels5`, so it reads Elly's bootstrap models and raw prediction
tifs from the local Google Drive.

Most prep scripts (01–05) run locally; 06 runs on Fir since 2026-10-14 (as 12A, it sets `cc` from `SLURM_JOB_ID`; run locally it still reads G:). Compute-heavy scripts (07, 11, 12A, 12C, 12D) run on the cluster via SLURM array jobs. 14B, 14C and 18 also run on the cluster, as two chained jobs after 12H (14B + 14C, then 18), so 12H's output never has to leave Fir; only their CSVs and figures are pulled.

## Pipeline

Scripts are numbered in execution order:

| Script | Where | Purpose |
|--------|-------|---------|
| `00_set_backfilling_order.R` | local | Define hierarchy order for biotic covariate backfilling; saves `biotic_variable_hierarchy.rds` |
| `01_reproject_soil_covariates.R` | local | Crop and reproject ISRIC soil data to BAM boundary |
| `02_reproject_disturbance_covariate.R` | local | Prepare CAfire time-since-disturbance layer |
| `03_reproject_humanfootprint.R` | local | Create low-HF (`CanHF_1km_lessthan1.tif`) and high-HF (`CanHF_1km_morethan1.tif`) masks from the AREA MEAN of Hirsh-Pearson's 300 m cumulative score (`raw_300m/cum_threat2020.02.18.tif`, copied from G: `CovariateRasters/Disturbance`) over each 1 km cell. Run its mask section only; the goodsectors / patch-summary sections are diagnostics |
| `04_reproject_and_crop_hydrobasins.R` | local | Crop Level 6 HydroBASINS to BAM study area |
| `05_merge_low_density_subbasins.R` | local | Merge data-sparse subbasins so every unit has ≥Q25 low-HF pixels; subset to subbasins with any high-HF pixels → `hydrobasins_masked_merged_subset.gpkg` (674 subbasins) |
| `06_build_covariate_stacks.R` + `.sh` | cluster (or local) | Mosaic BCR covariate stacks by year and add soil + CAfire layers → `covariates_mosaiced_{year}.tif`. `bash 06_build_covariate_stacks.sh 1995 2000` submits one array task per year (32G / 12 h; each year took ~1 h 53 min and 14 GB in 2026-10-15's run). On Fir it reads Elly's stacks from `{nm_root}/gis/stacks/` and `data/raw_data/{CAfire,isric_soil}/` (staged from G:). It stops if a BCR stack or the year's CAfire layer is missing, or if the mosaic's layer names (in order) or grid differ from `covariates_mosaiced_2020.tif`'s |
| `07_train_and_backfill.R` + `.sh` | cluster | Entry point: for a given SLURM array index (one year × subbasin), train BART models and backfill high-HF pixels. Years come from `YEARS` (comma-separated, default 2020); `bash 07_train_and_backfill.sh 2015 2020` submits itself, one array of 674 tasks per year, and chains that year's 11. Each year needs its own `covariates_mosaiced_{year}.tif` and footprint masks (`CanHF_1km_{lessthan1,morethan1}_{year}.tif`; the undated files serve 2020), or `HF_MASK_YEAR=<y>` to borrow year y's masks, else the task stops. Hirsh-Pearson exists only for 2020, so pre-2020 years run with `HF_MASK_YEAR=2020` (user, 2026-10-14): train on land with no footprint in 2020, backfill the whole 2020 footprint. 8h / 1 core; 72G for subbasins 57, 61, 62, 98, 107 and 24G for the rest (see **Submitting Cluster Jobs**) |
| `08A_train_and_backfill_subbasin_s.R` | sourced | Core `train_and_backfill_subbasin_s()` function; loops over biotic covariates in hierarchy order. Backfills with the footprint ("Disturbance") covariates at 0 (Open Limitation #9) |
| `08B_deploy_gbart.R` | sourced | `deploy_gbart()`: Gaussian BART for continuous biotic covariates (log1p-transformed, 90/10 train/holdout split) |
| `08B_deploy_mbart.R` | sourced | `deploy_mbart()`: Multinomial BART for categorical land-cover covariates. Writes the modal class (metrics only) and 100 drawn classes `<cov>_class_###`, each sampled from one posterior draw's normalised probit probabilities (see **Categorical covariates**) |
| `08B_deploy_pbart.R` | sourced | `deploy_pbart_blank()`: probit BART on all the training land for where a biotic covariate is blank (e.g. `SCANFIheightcv_*` where there is no canopy); 08A blanks the backfill where P > 0.5 and writes `<cov>_blank` / `<cov>_pblank`. A covariate blank on all the training land is blank everywhere |
| `09_collect_metrics_gbart/mbart.R` | sourced | Collect in-sample BART metrics |
| `09_collect_holdout_metrics_gbart/mbart.R` | sourced | Collect holdout BART metrics |
| `10A_process_backfill_metrics.R` | local | Aggregate per-subbasin metrics and confusion matrices into CSVs |
| `10B_inspect_backfill_metrics.R` | local | Inspect and visualize model accuracy |
| `10C_abiotic_extrapolation_diagnostics.R` | local | Area of applicability (Meyer & Pebesma 2021) of each subbasin's low-HF training pixels in the natural abiotic predictors 08A trains BART on (its footprint "Disturbance" class excluded; see Open Limitation #9); flags a subbasin when > 50% of its backfilled pixels fall outside → `extrapolation_flags.csv`. Needs `FNN` (installed locally, not on Fir) |
| `10D_BART_posterior_diagnostics.R` + `.sh` | cluster | Test whether the 100 stored BART posterior draws (subsampled from gbart()'s 700) cover the posterior's shape and tails — `12F` resamples 1 of the 100 per counterfactual scenario |
| `11_premosaic_backfilled_stacks.R` + `.sh` | cluster | Mosaic per-subbasin BART backfill rasters into BCR-wide stacks (run before 12). One array task per year × BCR (`YEARS` as for 07; 19 tasks per year); the .sh asks 32G / 4h / 1 core, but since C2d (~2400 layers per subbasin) can11, can12, can13 and can61 need 128G / 8h (at 32G terra falls to its block path and stalls), so 07's .sh submits those four separately (`BIG_BCRS`). Rebuilds any mosaic older than one of its subbasin backfills and writes under a `_partial` name renamed on completion. Each subbasin is resampled onto its own grid-aligned window of the BCR grid and the windows are `merge()`d per variable. `MOSAIC_DIR` redirects the output for tests. It caps terra at 30% of the Slurm allocation (`memmax`), because terra sizes its blocks from the node's free RAM, not the job's |
| `12A_observed.R` + `.sh` | cluster (or local) | Reads Elly's unclamped 32-bootstrap prediction tifs and bootstrap model `.Rdata` from `{nm_root}/output/{07_predictions,06_bootstraps}/{species}/`, applies V5's two-stage truncation (via `12B_v5_truncate.R`), and writes `observed_bootstraps.tif` (32-layer clamped stack, UNmasked), `observed_mean.tif`, `observed_sd.tif` and per-species `truncation_params.rds` to `data/derived_data/predictions/`. Canadian BCRs (`can*`) only. One array task per species (16G / 8 h; skips BCRs already done). Run locally (`SPECIES=OSFL Rscript Rscripts/12A_observed.R`) it reads G: instead, and the outputs must be Globus-transferred before 12D |
| `12B_v5_truncate.R` | sourced | `v5_truncate()`: line-for-line port of V5 `analysis/10.Truncate.R` (as of V5 `f082866`). Applies both upper caps (`densmax`, then the 99.9th-percentile `q99`) and the range/water/extent masks. Sourced by `12A`, so it is staged on the cluster; only `v5_load_masks()` (validation) needs `sf`, which it calls as `sf::`. `q99` can be passed in frozen and the legacy EPSG:3978 step skipped |
| `12C_build_prediction_weights.R` + `.sh` | cluster | Build per-species×BCR `weight.tif` (= range membership × not-water × inside-data-limit × **inside the BCR's own polygon**), replicating V5 `10.Truncate` range/water/extent/mosaic masking; all four terms use `touches = TRUE` (Open Limitation #7). Reads source masks from `data/raw_data/v5_gis/` (no G: access); each species needs `v5_gis/ranges/{species}.tif`, copied from G: `gis/ranges/` and staged. Grid template is the BCR stack. Run once per species before 12D. 12F multiplies BOTH observed and backfilled density by this weight (`w·bf − w·obs = w·(bf − obs)`). |
| `12D_repredict_all_coalitions.R` + `.sh` | cluster | Entry point: re-predict bird densities for ALL 255 coalitions, ONE SLURM array task per species × BCR (`coalition_task_table()` in 12F: every species' `06_bootstraps` files in `list.files()` order); its .sh sizes the array from `SPECIES` and chains 12H. Sources 12E, 12F, 12G (and compiles `12G_gbm_tree_walk.cpp`); writes one `density_tables/by_bcr/{species}_{year}_{bcr}.rds` per task (atomically; a BCR 12F skips is still written, with `result = NULL`). Requires `observed_bootstraps.tif` (hard error if missing). |
| `12G_gbm_tree_walk.R` + `.cpp` | sourced | Bit-identical replacement for gbm's compiled tree walk (`gbm_pred`), 3.6–4.4× faster: gbm makes three R API calls per node visited, this hoists them. `gbm_design()` builds predict.gbm's design matrix once per bootstrap; `gbm_check_fast()` `stop()`s in every 12F worker unless the fast path is `identical()` to `predict.gbm`. The `.cpp` also holds `coal_rowsum()`, 12F's copy-free, bit-identical coalition `rowsum` |
| `12H_merge_bcr_tables.R` + `.sh` | cluster | Run after 12D (`--dependency=afterany`; 12D's `.sh` chains it), with the same `SPECIES`, and writes only those species' tables. Checks every expected species × BCR has a per-BCR result (names the array indices to resubmit if not), that all came from one version of 12E/12F/12G (`code_md5`), then binds them in `06_bootstraps` order → `density_tables/{species}_{year}_coalition_{cid}.rds` (cid 2..256) and `arrays/`, plus `{species}_{year}_shapley_samples.rds` (per-sample subbasin Shapley values, stacked across BCRs in table row order) |
| `12E_shapley_utils.R` | sourced | Coalition enumeration, Shapley value computation utilities |
| `12F_predict_species_all_coalitions.R` | sourced | `predict_species_all_coalitions()`: builds the backfilled field ONCE per species×BCR over the all-8-sectors superset, then reduces all 255 coalitions as masked `rowsum`s (see **12F invariants**). Joint BRT×BART sampling; each bootstrap worker reduces its own field to per-coalition subbasin sums. Also predicts each bootstrap with observed vegetation and the footprint covariates at 0 (`d0`), for the direct/vegetation split (Open Limitation #9); backfilled vegetation goes only on a coalition's direct footprint, `d0` on the rest (Open Limitation #10; needs 14A's `{sector}_direct.tif`). Applies `shapley_weight_matrix()` (12E) to every (bootstrap, scenario) sample → per-sample subbasin Shapley values `[n_sub × 8 × 3200]` for 14B (Shapley is linear in v, so their mean equals the Shapley value of the tables' means). Writes `pixel_change/{species}_{year}_{bcr}.tif` (mean obs and cf of the all-sector coalition at each footprint pixel, for 18) after checking its pixels sum to the cid-256 table. Also holds `coalition_task_table()`, `coalition_array_ids()` and `combine_bcr_results()` (12H's cross-BCR bind). |
| `13_importance_of_covs_used_in_counterfactual.R` | cluster | Assess percentile importance of backfilled covariates in V5 bird models |
| `14A_reproject_hirshpearson.R` | local | Build each of the 8 sectors' 1 km layers on the CanHF mask grid from the raw 300 m Hirsh-Pearson rasters (`raw_300m/`, Borealis doi:10.5683/SP2/EVKAVL): `{sector}.tif` (max score) and `{sector}_direct.tif` (any direct 300 m cell; Open Limitation #10). Run after `03`, before `12D`; both are staged to the cluster. Optional args name a subset of sectors (~5 min each) |
| `14B_sector_attribution.R` + `.sh` | cluster | Its .sh runs 14B then 14C in one job (1 core, 16G / 4 h) and chains 18. Needs `data/raw_data/ModelReleaseDecisions.xlsx` and 10C's `rds_files/extrapolation_flags.csv` on Fir (10C runs locally). Reads 12H's per-sample subbasin Shapley values for every species in `density_tables/`, sums them bottom-up (subbasin → BCR → national) SAMPLE BY SAMPLE and summarises each level (mean, SD, 5th/95th percentiles) → `sector_effects/shapley_*.csv`; stops unless the sample means equal the Shapley values of the coalition tables' means. Splits every value into a direct and a vegetation part (Open Limitation #9). **This is where BAM's release filter lives** (`DROP_WITHHELD` / `withheld_models`): every Canadian row of the "remove" tab of BAM's `review/ModelReleaseDecisions.xlsx` (gitignored copy at `data/raw_data/ModelReleaseDecisions.xlsx`) is withheld, e.g. CAWA `can40` (AUC), but deliberately still produced upstream, so only the reported numbers are filtered. A species with every BCR withheld (8 of 146 as of 2026-10-14, e.g. BGGN) is skipped, so 14C and 18 leave it out too |
| `14C_results_tables.R` | cluster | Run after 14B, by 14B's `.sh`. One row per species, from the Shapley samples over the BCRs 14B kept (release filter): `sector_effects/results_summary.csv` (observed and no-industry populations, the change industry made, and that change as % of the no-industry population) and `results_sectors.csv` (each sector's change as % of the no-industry population), each with 5th/95th percentiles, and `results_bcr.csv` (both per species × BCR, % of the BCR's no-industry population). **Reported sign is the opposite of 14B's**: change = observed − counterfactual = −v(N), so negative = birds lost to industry, as the maps (18) show it. A % interval divides each sample by its own no-industry total. Stops unless its means equal 14B's national CSV to within 1 bird (14B rounds to whole birds) |
| `15A_plot_population_distributions.R` | local | Plot empirical population distributions, observed vs counterfactual, from the bootstrap × scenario arrays `12D` saves to `density_tables/arrays/` |
| `15B_preliminary_singletons.R` | local | Interim pre-Shapley single-sector standalone impacts at footprint / watershed / BCR scales → `logs/15B_singletons_summary.csv`. Superseded by `14B`'s exact Shapley values for attribution; retained for the scale-dependent reporting `15C` plots |
| `15C_singletons_plot.R` | local | Render `15B`'s summary CSV as sector impacts at the three geographic scales (map panel broken: its `predictions_coalitions/` input is deleted) |
| `16A_lowHF_isnot_a_proxy.R` | local | Test the premise that low-HF pixels are not a spatial proxy for the pre-industrial state |
| `16B_map_sampled_subbasins.R` | local | Map the 16 subbasins `16A` samples, reproducing the index from the same seed without re-running the RF pipeline |
| `17_detect_historical_hf_extent.R` | local | Detect whether the 2020 HF mask was smaller in earlier years, via PCA distance between observed and backfilled biotic vectors inside the HF mask |
| `18_map_population_change.R` + `.sh` | cluster | Chained after 14B's job, or alone with `SPECIES` (default: every species in 14B's tables) → `output_figures/`. Runs on Fir, where `pixel_change/` lives; terra only, no `sf`. Titles from `data/raw_data/species_names.csv` (BAMexploreR's `spp_tbl` plus 8 codes it lacks; a species not in it is titled by its code). A failing species is reported and skipped, and the job fails at the end. One map per species of what industry did to it (obs − cf; red = decline), clipped to Canada: the pixel-level change summed into 10 km cells (`change_grid_*`, from 12F's `pixel_change/` rasters; stops unless the cells add up to 14B's national change). The only figure type the user wants (2026-10-13): the per-subbasin and CAWA/OVEN % maps were dropped. Honours `DENS_MIN`/`DROP_Q99` like 14B and keeps production's colour scale for a sensitivity run |
| `misc/` | local | Downstream analysis (population summaries, visualization, vegetation vs. mines) plus one-off diagnostics/profiling not in the execution-order pipeline (`diag_*`, `analyze_seff`, `check_model_complexity`, `debug_*`, `figure_*`) |

## Transferring Files to the Cluster (Globus)

Use the Globus CLI at `C:\Users\mannf\AppData\Local\Python\pythoncore-3.14-64\Scripts\globus.exe`.

**IMPORTANT: never use `--batch` mode** — it always returns "permission denied". Submit one `globus transfer` call per file, or one filtered recursive transfer (`--include=<pattern>` ... `--exclude=*`) instead.

**Always pass `--notify off`.** Globus emails the user once per task by default, so a one-file-per-task
pull sends one email per file (500+ for the C2 pull). Check success locally instead: poll
`globus task list --filter-status ACTIVE` until empty, then confirm every file landed.

Endpoint IDs:
- Local: `a7878ccc-747b-11ef-b4b8-8fef73a45f39` (local paths must start with `/C/`)
- Cluster (Fir): `8dec4129-9ab4-451d-a45f-5b4b8471f7a3`

Example (one file at a time):
```powershell
$globus = "C:\Users\mannf\AppData\Local\Python\pythoncore-3.14-64\Scripts\globus.exe"
& $globus transfer `
  "a7878ccc-747b-11ef-b4b8-8fef73a45f39:/C/Users/mannf/Drive/.../Rscripts/12D_repredict_all_coalitions.R" `
  "8dec4129-9ab4-451d-a45f-5b4b8471f7a3:/home/mannfred/scratch/impact_assessment/Rscripts/12D_repredict_all_coalitions.R" `
  --label "my file" --sync-level checksum --notify off
```

More Globus gotchas (Windows glob expansion, Git Bash path rewriting, the scratchpad): `TODO.md`.

## Submitting Cluster Jobs

Before any cluster run, stage the entry point's **full `source()` closure**, not just the files
you edited — see `TODO.md` "Staging discipline".

Backfilling (SLURM array, one job per subbasin index):
```bash
# ONE script that submits itself: run with `bash` on a login node (no SLURM_JOB_ID) it submits;
# run by Slurm it trains/backfills its array index. (So never `bash` it inside an salloc
# session.) Years on the command line. Per year it exports YEARS=<year>, submits both memory
# tiers over 1..674 (indices must stay under Slurm's MaxArraySize of 1001, so years are never
# numbered in one array), and chains 11 after both: can11, can12, can13, can61 at 128G / 8 h,
# the other 15 BCRs at 32G / 4 h:
bash 07_train_and_backfill.sh 2020
HF_MASK_YEAR=2020 bash 07_train_and_backfill.sh 1995 2000 2005 2010 2015
# Each year needs covariates_mosaiced_{year}.tif first (06 on Fir):
bash 06_build_covariate_stacks.sh 1995 2000 2005 2010 2015
# By hand: YEARS=2015 sbatch --array=<indices> --mem=24G 07_train_and_backfill.sh
# (not --export=ALL,YEARS=...: sbatch splits --export on commas). The .sh's own
# --mem=128G is only a safe default for hand-submitted indices. 07 overwrites its outputs and
# 11 rebuilds any mosaic older than its inputs, so a rerun needs no wipe. 07 catches a
# subbasin's error and exits 0, so sacct shows it COMPLETED: grep the .out files for
# "Error in subbasin".

# Memory (`cluster_logs/07_2020_c2f/sacct_07_61937625.txt`; the C2d run fit the same tiers):
# peak MaxRSS median 3.8 GB, max 52.8 GB; only the five can11 prairie basins (57, 61, 62, 98,
# 107) exceed 20 GB, hence 72G for them and 24G for the rest. The 1995–2015 run (2026-10-20)
# peaked higher, at 67.3 GB (subbasin 98) and 18.7 GB, so 72G has little to spare. Memory is
# what Alliance bills.
# Single-threaded (BART::gbart), so --cpus-per-task=1.
#
# If one index OOMs, give that index (= subbasin) more room with the same YEARS, and add its
# subbasin to BIG_SUBS in 07_train_and_backfill.sh:
YEARS=<year> sbatch --array=<i> --mem=128G 07_train_and_backfill.sh

# 11 by hand (11 skips mosaics newer than their backfills); the index is the BCR's position in
# 11's sorted bcr_vec, so 2,3,4,12 = can11, can12, can13, can61:
YEARS=2020 sbatch --array=1-19 --mem=128G --time=08:00:00 11_premosaic_backfilled_stacks.sh
```

Re-predicting birds:
```bash
# Species are arguments: each of 12A, 12C and 12D's .sh, run with `bash` on a login node,
# exports SPECIES (comma-separated), sizes its own array and submits itself, as 07 does. The
# R scripts stop if SPECIES is unset. By hand, put it in the environment, never in --export
# (sbatch splits --export on commas): SPECIES=OSFL,GRSP sbatch --array=1-2 12C_build_prediction_weights.sh

# Phase 1: observed predictions, one array task per species. Reads Fir's 07_predictions and
# writes observed_bootstraps.tif + truncation_params.rds in place (12F stops if either is
# missing).
bash 12A_observed.sh OSFL GRSP

# Phase 1b: build prediction weights ONCE per species before 12D (independent of 12A, so the
# two can run together). Reads data/raw_data/v5_gis (stage ranges/{species}.tif first) +
# Regions/BAM_BCR_NationalModel_Unbuffered.shp; writes weight.tif per species x BCR.
bash 12C_build_prediction_weights.sh OSFL GRSP
# weight.tif is a HARD dependency and the only masking left (water / out-of-range /
# out-of-extent pixels no longer drop out of complete.cases(); an unmasked run over-counted
# CAWA can10 7.8x, quietly). 12F stops if it is missing, stale, or 0/NA over the whole BCR; a
# species whose range misses the footprint (BAWW, BBWA can41) runs on with impact 0.
# weight.tif carries its version in its band name ("weight_v3_touches"): 12C rebuilds a stale
# one and 12F refuses to run against it, so no manual `rm` is ever needed.

# Phase 2: ONE array task per species x BCR computes ALL 255 coalitions for that BCR in a
# single pass (see "12F invariants" below), then 12H merges.
# `bash 12D_repredict_all_coalitions.sh <species>` counts the species' can* models (honouring
# TEST_BCR), submits that array at 8 cores / 03:00:00 in two memory tiers (can11's tasks at
# 64G, the rest at 48G; BIG_BCRS in the .sh), and chains 12H (afterany on both, 32G / 3 h) with
# the same environment; the task table is printed at the top of every 12D log. An index past the
# table's end exits 0. Per-species outputs never overwrite other species', so a new batch needs
# no wipe; the wipe below is for re-running the SAME species.
# Sized from 807 tasks of 80 species (pre-C2d): MaxRSS median 16.9 GB, max 44.8 GB (can11; 20
# of 63 can11 tasks above 40 GB); median 5.7 min, max 66 min (AMRO can11). The six C2d species
# peaked at 44.5 GB in can11 and 40.9 GB elsewhere (can61). A task killed for memory is
# resubmitted alone with --mem=64G.
#
# For a FRESH full re-run of the same species, delete their stale outputs first. `rm -f` on
# files, never `rm -rf` on a directory: 15A reads arrays/. 12H refuses to merge per-BCR files
# from different code versions, but it cannot tell a stale file of the SAME version from a
# fresh one, so wipe by_bcr/ too.
rm -f ../data/derived_data/density_tables/{CAWA,OVEN}_*.rds
rm -f ../data/derived_data/density_tables/arrays/{CAWA,OVEN}_*.rds
rm -f ../data/derived_data/density_tables/by_bcr/{CAWA,OVEN}_*.rds
bash 12D_repredict_all_coalitions.sh CAWA OVEN          # 12D array + 12H; prints both job ids

# If 12H reports missing BCRs it prints the exact `SPECIES=... sbatch --array=<indices> ...`
# to resubmit.

# Many species: run_species_batches.sh reads species_todo.txt and runs 12A + 12C for all of
# them (`prep`), a can42 smoke (`smoke`), and 12D + 12H in waves of <= MAX tasks (default
# 400; MAX=900 fits under Slurm's MaxArraySize of 1001) (`plan`, then `wave 1`, `wave 2`, ...
# one at a time).

# Smoke test (2 bootstraps): TEST_BCR filters the task table to the listed BCRs of every
# species, and the chained 12H gets the same TEST_BCR and TEST_N_BOOT (it refuses to merge
# smoke files under production settings). The full run later overwrites the smoke's files.
TEST_BCR=can60 TEST_N_BOOT=2 bash 12D_repredict_all_coalitions.sh CAWA OVEN

# Phase 3: Shapley attribution and results tables (14B -> 14C), then maps (18), chained. 14B
# reads every species in density_tables/, so it waits for the LAST 12H (the last wave's, in a
# batch run). SPECIES in front limits only 18's maps:
bash 14B_sector_attribution.sh <12H job id>
# Maps alone, e.g. a test batch, or all species once that works:
SPECIES=CAWA,OVEN,OSFL,GRSP,LEYE,BOBO sbatch 18_map_population_change.sh
# Then pull data/derived_data/sector_effects/*.csv and output_figures/*.png (a few MB), never
# density_tables/.

# Sensitivity run without the frozen q99 cap (Open Limitation #6): put DROP_Q99=1 in front
# (DROP_Q99=1 bash 12D_repredict_all_coalitions.sh CAWA OVEN), and run 14B as
# `DROP_Q99=1 bash 14B_sector_attribution.sh`. They read and write density_tables_noq99/ and
# sector_effects_noq99/, never the production directories.
# Low-density floor (C4): DENS_MIN=0.005 does the same with every capped density below
# 0.005 birds/ha set to 0, per pixel and sample, before the weight, on the obs, d0 and bf sides
# (observed totals included). Directories: density_tables_dmin0.005/, sector_effects_dmin0.005/.
# 12E sensitivity_run() parses both switches and refuses them together.
```

## Fir Cluster Specifications

Fir (SFU, Digital Research Alliance of Canada). Login `fir.alliancecan.ca`; automation via
`robot.fir.alliancecan.ca` (no `crontab`); Globus collection `alliancecan#fir-globus`.
- CPU nodes: 192 cores (8 per CCD, so `--cpus-per-task=8` keeps threads cache-local), 750 GB RAM.
- Max walltime 7 days. No partition flag needed. Always `--account=def-bayne` (the allocation
  every job script uses — not `def-ecknight`, which is only the project *directory* path).
- Storage: `$HOME` is small; `$HOME/scratch` (where `ia_dir` lives) is large, **not backed up and
  purged automatically**, so copy outputs to keep to `$HOME/project/...` (backed up daily).

## Key Architecture Decisions

**Biotic variable hierarchy**: Continuous biotic covariates are backfilled in a fixed order (defined in `00_set_backfilling_order.R`). Each covariate's BART model uses all preceding biotic covariates as additional predictors, so backfilled values cascade down the hierarchy.

**Categorical covariates**: Six land-cover layers (`ABoVE_1km`, `NLCD_1km`, `MODISLCC_1km`, `MODISLCC_5x5`, `SCANFI_1km`, `VLCE_1km`) are treated as categorical throughout and use `mbart()` instead of `gbart()`. They are excluded from predicting continuous features.
- **No water.** `07`'s `water_classes` (SCANFI_1km 8, VLCE_1km 20, MODISLCC_1km/_5x5 17) are left out of each layer's training rows, so the backfill never predicts water on land.
- **Drawn classes, not the modal class.** 12F uses `<cov>_class_###` with the scenario's draw index, as it does the continuous `_draw_*` layers. The modal class gives the majority class to every pixel where it is likeliest and erases minority classes.
- **Not `mbart2`'s `prob.test`.** With `type = "pbart"` it is a softmax of probit scores, exp(z)/Σexp(z), which roughly doubles rare classes. 08B uses pnorm(z)/Σpnorm(z) from `yhat.test`.
- **Renumbered classes.** `mbart2` fits `y == h` for h in 1..K, so 08B renumbers the training split's classes 1..K.

**Blanks are modelled, not filled.** BART always returns a value, but some V5 covariates are missing by design (`SCANFIheightcv_*` where canopy height is 0), and gbm sends that NA down its missing branch.
- `08B_deploy_pbart.R` models where each biotic covariate is blank. A blanked backfill is NA to the later covariates in the hierarchy, so their Fix B `_isNA` flag is 1 there as on the training land.
- 12F passes blanks to gbm as `NA_real_` and counts them as complete.
- Without this, every backfilled pixel looked treed (`_isNA` = 0), which made GRSP's can11 decline an artefact (C2d).

**Coordinate scaling**: Within each subbasin, lat/lon are centered and scaled before being added as predictors to improve BART performance.

**Output per subbasin**: `data/derived_data/bart_models/{year}/subbasin_{i}/subbasin_{i}_backfill.tif`, `_metrics.rds`, `_confusion.rds`. Note **not every continuous biotic covariate has `_draw_*` layers**: where a covariate is constant across a subbasin's low-HF pixels, `08A` skips BART and writes `<cov>_mean`/`<cov>_sd` instead (all-NA if there were no valid training rows). The BCR mosaic unions layers across subbasins, so a covariate with draws in some subbasins is NA in the constant ones; `12F` fills those pixels from `<cov>_mean` (raw scale — no `expm1`) for every draw, or its complete-case gate would drop them (CAWA can60: 43% of weight > 0 footprint pixels). A covariate with no draws anywhere in the BCR is carried as a one-draw covariate from `<cov>_mean`, so it too gets the low-HF constant rather than its observed value. `BART_MODELS_DIR` redirects 07's output for tests; 11 always reads `bart_models/`, so a production 07 runs with it unset.

**Re-prediction**: `11_premosaic` mosaics backfilled subbasin rasters into BCR-wide stacks. `12A_observed.R` (on Fir) writes the truncated `observed_bootstraps.tif` per species×BCR, a hard dependency for `12D`/`12F`. `12D_repredict_all_coalitions.R` runs ONE array task per species × BCR and computes all 255 coalitions in a single pass: it sources `12F_predict_species_all_coalitions.R`, which builds the backfilled field ONCE over the BCR's all-8-sectors superset and reduces every coalition as a cheap masked `rowsum`; `12H` merges the per-BCR files. For a given coalition S of sectors, pixels where any sector in S has footprint (AND CanHF ≥ 1) get the footprint covariates at 0, and those in the DIRECT footprint of a sector in S also get backfilled vegetation (Open Limitation #10); all other pixels use observed. Joint BART×BRT sampling nests BART posterior draws inside BRT bootstrap iterations.

**12F invariants**: The joint-sampling seed depends only on species, BCR, bootstrap `i`, and scenario `k` — never the coalition. So for fixed (species, BCR, i, k) the backfilled density field over pixels is identical across all 255 coalitions; the coalition only selects which pixels are masked in, never the backfilled value at a pixel. This is what licenses computing the backfilled field ONCE over the all-8-sectors superset and reducing each coalition as a masked `rowsum`. Constraints that must hold:
- The complete-case mask is identical between obs and bf (a single `keep` drives both).
- BART draws are log1p-scaled → `expm1` before use; non-finite → NA, never 0.
- Categorical `var.levels` are indexed by `match()` against `var.names`, never by name.
- BOTH V5 upper caps are applied per pixel in the gbm step and DO feed the density table: `pmin(pmin(pred, qsp), q99)`, where `qsp` is the densmax cap (`q.out$densmax`; there is no `$q` column). `q99` is the per-BCR 99.9th percentile **frozen from the observed landscape** by 12A and read from `predictions/{species}/truncation_params.rds`. It must never be re-derived from a counterfactual, or the cap enters the obs/bf contrast. The observed side is pre-clamped by 12A, so both sides carry identical caps.
- Grouped `_draw_*` reads use INTERLEAVE=BAND (once per BCR).
- **Missing covariates reach gbm as `NA_real_`, never `NaN`.** `terra::values()` returns missing
  cells as `NaN`, and gbm's tree walk sends ONLY `NA_real_` down a split's missing branch (a
  `NaN` fails the `<` test and goes right). V5's observed predictions used the missing branch,
  so a `NaN` on the bf side makes a fake impact (CAWA can11: ≈ +10% of the observed birds on
  that footprint, before the fix). 12F recodes `NaN → NA_real_` where `X_obs_super`
  is built. `SurfaceWater_1km` is gbm-categorical (`var.type` 2, levels "0","1") but not in
  `categorical_responses`, so it is passed raw as V5 did; that is correct only because its
  levels are exactly "0","1", and a raw `NaN` there is undefined behaviour in `gbm_pred`.
- **Identity gate.** Before sampling, 12F predicts EVERY bootstrap from the OBSERVED design
  (nothing backfilled) through the workers' exact path and `stop()`s unless it reproduces
  `observed_bootstraps.tif` (float32 tolerance), over-sampling pixels with a missing covariate or
  with a class some bootstrap never saw. It catches any obs/bf rule divergence, and a
  `b.list[[i]]` ↔ observed layer `i` mis-pairing.
- **Per-bootstrap reduction.** Workers return per-coalition `[n_sub × n_scen]` sums, not their
  pixel fields. `rowsum()` adds each column independently in row order, so this is bit-identical
  to reducing the assembled `[n_super × (n_boot·n_scen)]` matrix (up to ~17 GB for can11). The
  sum itself is `coal_rowsum()` (12G `.cpp`): R's `rowsum` loop over the kept rows in place, so
  it never copies `sc[kr, ]` 255 times.
- **Categorical levels are PER BOOTSTRAP.** Each bootstrap is fitted to a different resample, so
  `var.levels` differ between bootstraps of one model (CAWA can14: 28 of 32 know VLCE_1km classes
  32/40 that bootstrap 1 never saw). `as_model_factors()` converts with each model's own levels,
  as `predict.gbm` and V5 do. The complete-case gate tests the RAW backfilled class for NA,
  never a factor.
- **Worker memory.** A forked worker inherits the parent's GC trigger, which the BART-draw reads
  push to ~10–15 GB on the big BCRs, so 8 workers each piling up garbage OOM-ed C1 at 64G.
  Workers allocate only per-draw vectors, collect them with `gc(full = FALSE)` (young objects
  only: cheap, and it leaves the parent's pages alone), and each fork runs ONE bootstrap
  (`mc.preschedule = FALSE`). A worker the kernel kills returns `NULL`; 12F names it.

**Shapley attribution**: 8 sectors → 256 coalitions (2^8; cid 1 = empty is skipped → 255 computed). All 255 are produced in one pass per species × BCR (12D), merged by 12H. `14B_sector_attribution.R` computes exact Shapley values from the coalition density tables. Shapley values sum exactly to the total HF impact. `12E_shapley_utils.R` provides coalition enumeration and the Shapley formula.

**Bottom-up aggregation**: The subbasin is the atomic spatial unit. BCR totals = sum of subbasin values within the BCR. National totals = sum of BCR values. Uncertainty is aggregated per sample, never by adding variances: 14B sums each (bootstrap, scenario) sample's subbasin Shapley values within a BCR, then across BCRs (bootstrap i paired across BCRs, as V5's national mosaic pairs them), and takes the SD and percentiles of those sums. The SD columns of the coalition tables are per-subbasin only and must not be combined as if independent: every subbasin of a BCR shares its 32 bird models. The same holds for ratios: a "% of observed" interval divides each sample by its own bootstrap's observed total (`obs_total` in the Shapley samples file).

**Abiotic extrapolation diagnostics**: `10C` computes each subbasin's area of applicability: a backfilled pixel is outside it when its standardised distance to the nearest training pixel exceeds what training pixels show to each other across spatial blocks. `frac_outside_aoa` and the flag (> 0.5) are annotated in `14B`'s subbasin output, and `14B` reports the impact carried by flagged subbasins per BCR and nationally.

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
│   ├── CAfire/CAfire_{year}_masked.tif      # 06's input (02's output, staged from G:)
│   ├── isric_soil/isric_soil_covariates_masked.tif   # 06's input (01's output, staged from G:)
│   ├── covariates_mosaiced/
│   │   └── covariates_mosaiced_{year}.tif
│   ├── hirshpearson/
│   │   ├── CanHF_1km_lessthan1.tif
│   │   ├── CanHF_1km_morethan1.tif
│   │   ├── {sector}.tif / {sector}_direct.tif   # 1 km, built by 14A (8 sectors)
│   │   └── raw_300m/                 # Hirsh-Pearson 300 m originals + cum_threat (deleted locally 2026-10-08; re-download to rebuild 03/14A)
│   ├── hydrobasins_masked_merged_subset.gpkg
│   ├── species_names.csv                # species code → common name (18's titles; versioned)
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
    │   ├── {species}_{year}_coalition_{id}.rds   # subbasin-level population arrays per coalition
    │   ├── {species}_{year}_shapley_samples.rds  # per-sample subbasin Shapley values (12H → 14B)
    │   ├── arrays/                               # per-pixel bootstrap × scenario arrays (12D → 15A)
    │   ├── pixel_change/                         # {species}_{year}_{bcr}.tif: mean obs and cf (all-sector coalition) at footprint pixels, birds/km² weighted (12F → 18)
    │   └── by_bcr/                               # per-task 12D output, merged by 12H
    ├── predictions/
    │   └── {species}/
    │       ├── truncation_params.rds              # frozen {densmax, q99} per BCR (12A; 12F stops if absent)
    │       └── {bcr_code}/{year}/
    │           ├── observed_bootstraps.tif        # canonical 32-layer bootstrap stack (UNweighted)
    │           ├── observed_mean.tif / observed_sd.tif
    │           └── weight.tif                     # range×water×extent×BCR weight (built by 12C)
    ├── sector_effects/
    │   ├── shapley_subbasin.csv
    │   ├── shapley_bcr.csv
    │   ├── shapley_national.csv
    │   └── results_summary.csv / results_sectors.csv   # species tables (14C)
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

Large spatial files (`.tif`, `.gpkg`, `.shp`) and most `.rds` files are gitignored. The CSV outputs in `data/derived_data/rds_files/` (accuracy/confusion matrices) and `sector_effects*/` (14B) are kept local, not versioned (user decision, 2026-10-04).

## Open Limitations (updated 2026-10-13)

Numbers are cited from code and `TODO.md`, so closed items keep their number.

### Conceptual

1. ~~**Abiotic extrapolation risk**~~: **Addressed** by `10C`'s per-subbasin area of applicability, annotated in 14B. It flags extrapolation; it does not correct it.

2. ~~**Sector impacts are not additive**~~: **Addressed** by exact Shapley values over all 2^8 coalitions. They assume every v(S) is well estimated; coalitions with large combined footprints may carry more extrapolation risk.

3. **No spatial spillover**: The formula `cf = obs_on_non_coalition + backfilled_on_coalition` assumes removing a coalition's footprint only affects birds on those pixels. Edge effects, area sensitivity, and functional connectivity mean impacts extend beyond the footprint boundary (especially important for linear features like roads and seismic lines).

9. **Footprint covariates in the counterfactual (DECIDED 2026-09-25, C2e; record in `git show 3fcef6d:TODO.md`)**. V5's
   "Disturbance" predictor class (CanHF_1km/_5x5, canroad_1km/_5x5, CCNL_1km night lights) is
   human footprint. The counterfactual has no industry, so a backfilled pixel carries no
   footprint: a backfilled forest cannot sit on top of a road. Both uses therefore set the
   class to 0.
   (a) **Bird model.** 12F predicts with it at 0 at every backfilled pixel. It also predicts
   each bootstrap once with the OBSERVED vegetation and the class at 0 (`d0`), which splits
   every Shapley value exactly into a direct part (`d0 − obs`: the bird model's response to the
   footprint itself) and a vegetation part (`bf − d0`). Positive responses to modest
   disturbance can be real (Mahon et al. 2019, Ecol. Appl. 29:e01895: deciduous-associated
   species, Ovenbird included, increase with wide linear features), so the direct part is
   reported, not removed. It can be large and of opposite sign (OVEN roads, C2d run: −1.27 M
   direct, +1.74 M vegetation).
   The part rests on thin survey support: in the can11/12/13 models at least 90%
   of survey locations have CanHF_1km ≥ 6–7 and canroad_5x5 ≥ 0.3–0.5 (gbm `var.levels`
   deciles), so the response at 0 comes from the low tail. The 5×5 covariates are exact only
   for the all-sector coalition: for a partial one, zeroing them also removes other sectors'
   footprint within the neighbourhood.
   (b) **BART.** 07's `abiotic_vars` includes the class. 08A trains on observed values and
   backfills with the class at 0, which is inside the low-HF training range (each pixel's
   observed values, beyond that range, held vegetation at the training edge: as if the road
   were still there). `10C` excludes the class from its extrapolation diagnostic: its observed values are
   outside the training range by construction.

10. **What counts as a sector's footprint (DECIDED 2026-09-28, C2f; same record): two masks.**
    Hirsh-Pearson scores roads, rail, mines, reservoirs and oil and gas by distance band, so a
    sector's pressure layer includes influence zones around the mapped feature. The split
    follows their own: direct pressure is "the actual footprint such as physical removal of
    viable habitat"; indirect pressure is access, edge effects and other pressures the feature
    makes possible. Vegetation is backfilled only on the direct footprint. Whether vegetation
    in the influence zones is unaltered was never measured (the C2f test hints it differs:
    CAWA can60's vegetation part halved when they stopped being backfilled). The cumulative
    score is the plain SUM of the sector layers, so `CanHF ≥ 1` does not separate the two.
    Per sector, `14A` builds from the raw 300 m layers:
    - the **footprint** `{j}.tif > 0` (and CanHF ≥ 1): the counterfactual sets the footprint
      covariates to 0 there, on observed vegetation (`d0`);
    - the **direct footprint** `{j}_direct.tif == 1`, a subset: only there is the vegetation
      backfilled. Direct = any 300 m cell with built/crop/pasture > 0, roads/rail/dams ≥ 6
      (0–300 m band of secondary roads, rail and reservoirs, plus inner highway bands), mines ≥ 6
      (0–600 m, or 0–1500 m for large mines, which most are), oil and gas = 10 (a 900 m disc
      around the site, not one cell).
    - **Mines and oil and gas understate converted land.** HP mapped them as points, so the
      direct area is a disc around each point, whatever the real pit or well field covers. It is
      4.8% (mines) and 7.7% (oil and gas) of their 1 km footprint, and only the disc is backfilled.
      The 1 km "any" rule widens each disc to whole cells (×2.2 and ×2.9), a margin of a few
      hundred metres. A known limitation: there are no polygon data to fix it.
    A coalition's counterfactual on its footprint is `bf` over its direct pixels plus `d0` over
    the rest; the direct part (`d0 − obs`, whole footprint) and vegetation part (`bf − d0`,
    direct pixels) still sum exactly to the Shapley value. `bf_on_coalition_*` keeps its name
    but holds that combined counterfactual. `14A` and `03` aggregate 300 m → 1 km with max /
    area mean, never bilinear (which smears each 300 m cell over 2×2 1 km cells).

### Computational

4. ~~**Uncertainty combination is ad hoc**~~: **Addressed** by joint sampling in 12F: each (bootstrap, scenario) pair draws a fresh BART realization inside the BRT bootstrap loop, so no post-hoc variance combination is needed.

5. ~~**Backfill coverage collapse**~~: **Fixed.** `12F` drops footprint pixels whose
   backfilled design is incomplete; this once dropped 98.5–99.7% in some BCRs. `02` recodes
   CAfire's "never burned" `NA → 0`; **Fix B**: `08A` median-imputes partial-NA V5 covariates
   (phenology) with an `_isNA` flag; **Fix A**: `12F` gates only on backfilled covariates. Every
   12D run since keeps ≥ 99.9% of weight > 0 footprint pixels in every BCR except can3 (82.7%
   before C2d, which should fix it; `TODO.md` Loose ends).

6. **The frozen `q99.9` cap biases impact downward, slightly**: `T(S) = clamp(clamp(S, densmax),
   q99) × weight` is applied identically to observed and every counterfactual, with **`q99`
   frozen from the observed landscape** (`predictions/{species}/truncation_params.rds`).
   Freezing is mandatory: a cap re-derived per counterfactual would adapt to the landscape being
   differenced. It is not neutral, because counterfactual densities hit the frozen ceiling more
   often. **Measured 2026-10-04 (C3, `DROP_Q99=1`, pre-C2d): small.** Uncapped, v(N) rises
   1.16% (CAWA) and 0.05% (OVEN), and no sector share moves; q99 binds at only 0–1% of observed
   footprint pixels. With `DROP_Q99`, 12F re-predicts the observed footprint (the tif is
   q99-clamped) and stops unless re-capping it reproduces the tif.
   `q.out`'s schema is `spp, thresh, countmax, densmax` — **there is no `$q` column**, and
   `pmin(x, NULL)` returns `numeric(0)` *silently*, so `12F` guards the schema and `stop()`s.
   `q99` is the binding cap almost everywhere (1.2×–35× below `densmax` for CAWA/OVEN).

7. ~~**Two masking defects in `weight.tif`**~~: **Fixed 2026-09-11.** Kept because both are easy
   to reintroduce:
   (a) **V5 prediction grids are BUFFERED 100 km past their subunit** (for V5's border blend;
   V5 crops after mosaicking, `10.Truncate.R:146`), and we read the pre-mosaic
   `07_predictions`. Without `12C`'s `inbcr` crop (`Regions/BAM_BCR_NationalModel_Unbuffered.shp`,
   identical to V5's `Subregions_Mosaics_EPSG3978.shp`), each subbasin straddling BCRs was
   summed under every BCR it touches (populations 2.2×–14.6× too high).
   (b) **`terra::rasterize()` defaults to `touches = FALSE`**, while V5's `mask()`/`crop()`
   behave as `touches = TRUE` (the centre rule let through 4.9–7.7% of abundance on water). All
   terms pass `touches = TRUE`. Harness: `Rscripts/misc/verify_weight_vs_v5_masking.R`.
   **Load-bearing invariant when adding a mask term**: `weight.tif` must be non-NA everywhere.
   That is how we reproduce V5's `mask.i[is.na(mask.i)] <- 0` (`10.Truncate.R:141`) under a sum —
   range NA→0 via `classify`, and every `rasterize()` called with `background = 0`. A new term
   without `background = 0` would silently *drop* pixels instead of zeroing them. `12C` `stop()`s
   if any weight cell is NA.

8. **V5's published abundances read ~3% low (UPSTREAM)**: `13.Summarize.R:227` sums
   `density × 100` over the `10_truncated` product in EPSG:3978. `× 100` (birds/ha → birds/km²)
   yields birds only if every pixel is 1 km² of ground, which holds in equal-area projections.
   EPSG:5072 (ours) is Albers **Equal Area**; EPSG:3978 is Lambert **Conformal** Conic, whose
   1000 m cells hold ~1.027 km² between its standard parallels. So V5's totals read low and ours
   are the correct side; weighting by `terra::cellSize()` makes the two agree to 0.05–0.23%. The
   fix is projection-agnostic: `global(t * 100 * cellSize(t, unit = "km"), "sum")`. It is a bias
   inside the ±25% bootstrap intervals and cancels in the obs/bf contrast, so no Shapley number
   moves. **Do not deliver in 3978 like V5**: `10.Truncate.R:19` calls that step legacy, and
   `clamp()` does not commute with projection, which would break 12F's superset → masked-rowsum
   decomposition. Harness: `Rscripts/misc/verify_projection_area_units.R`.


## Instructions from Masa
1.Always ignore the directory /Rscripts/misc when thinking. It's not immediately relevant to the project.
2.Note that population estimates made by boosted regression trees via `gbm` predict bird density at the hectare scale.
In order to match the scale of our covariate rasters (1km^2), we always multiply the pixel estimates by 100 to get the 1km^2 
density estimate.
