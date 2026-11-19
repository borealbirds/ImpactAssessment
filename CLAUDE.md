# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.
Open work lives in `TODO.md`. Fuller history of any decision below: `git log` and the
`git show <commit>:CLAUDE.md` / `TODO.md` snapshots that TODO.md lists.

## Project Overview

A boreal bird impact assessment pipeline: the counterfactual effect of industrial human footprint
(HF) on bird populations across Canada.

1. In each hydrological subbasin, train BART models on low-HF pixels to predict biotic
   (vegetation) covariates from abiotic ones.
2. Backfill (predict) the biotic covariates on high-HF pixels, removing the imprint of industry.
3. Re-predict bird densities with BAM's V5 models on the backfilled landscape, and attribute the
   difference to 8 industrial sectors with exact Shapley values.

**Only 2020 is reported.** BAM's 1995–2015 bird models are not validated, so every pre-2020
output was deleted (2026-11-02). Kept for a later rerun: Fir's `covariates_mosaiced_{year}.tif`
and `CAfire_{year}_masked.tif` for 1995–2015, and the `YEARS` code. Until then 17, 19A–19G and 20
cannot run.

## Instructions from Masa
1.Always ignore the directory /Rscripts/misc when thinking. It's not immediately relevant to the project.
2.Note that population estimates made by boosted regression trees via `gbm` predict bird density at the hectare scale.
In order to match the scale of our covariate rasters (1km^2), we always multiply the pixel estimates by 100 to get the 1km^2
density estimate.
3.Do not put subtitles or footnotes (captions) on figures unless asked for. A title, axis labels and a legend are enough;
the method notes belong in the script's header and in chat.

## Execution Contexts

```r
ia_dir <- if (nzchar(Sys.getenv("SLURM_JOB_ID"))) "/home/mannfred/scratch/impact_assessment" else getwd()
```

The cluster and the local RProject share the layout under **Data Directory Structure**. Scripts
that run in both places (07, 11, 12A, 14B, 14C, 18, 19A) use the line above; Fir-only scripts
(06, 12C, 12D, 12H) set the cluster path; local-only ones use `getwd()` (run from the repo root).
00–05 write to `data/raw_data/` and read only raw inputs from `G:/Shared drives/BAM_NationalModels5`.
Scripts 06, 11, 12A, 12C, 12D and 12F read Elly's BRT bootstraps, prediction tifs and BCR stacks
from `nm_root = /home/mannfred/projects/def-ecknight/NationalModels` (cluster-only; never change
it). 12A run locally reads `G:/Shared drives/BAM_NationalModels5` instead.

## Pipeline

Scripts are numbered in execution order. OL = Open Limitation (below).

| Script | Where | Purpose |
|--------|-------|---------|
| `00_set_backfilling_order.R` | local | Backfilling order of the biotic covariates → `biotic_variable_hierarchy.rds` |
| `01_reproject_soil_covariates.R` | local | Crop and reproject ISRIC soil data to the BAM boundary |
| `02_reproject_disturbance_covariate.R` | local | CAfire time-since-disturbance layer (never burned = 0, not NA) |
| `03_reproject_humanfootprint.R` | local | Low-/high-HF masks `CanHF_1km_{lessthan1,morethan1}.tif` from the AREA MEAN of Hirsh-Pearson's 300 m cumulative score over each 1 km cell. Its mask section runs alone (raw_300m/ is deleted locally) |
| `04_reproject_and_crop_hydrobasins.R` | local | Crop Level 6 HydroBASINS to the study area |
| `05_merge_low_density_subbasins.R` | local | Merge data-sparse subbasins within BCRs (≥ Q25 low-HF pixels each); keep those with high-HF pixels → `hydrobasins_masked_merged_subset.gpkg` (674). The production file (Feb 2026) predates 03's 2026-09-28 mask fix, so a rerun would not reproduce it; 05 refuses to overwrite it |
| `06_build_covariate_stacks.R` + `.sh` | cluster | Mosaic BCR stacks + soil + CAfire → `covariates_mosaiced_{year}.tif`. `bash 06_build_covariate_stacks.sh 1995 2000` (32G / 12 h per year). Stops if an input is missing or the layer names or grid differ from 2020's |
| `07_train_and_backfill.R` + `.sh` | cluster | One array task per subbasin (one array per year): train BART and backfill high-HF pixels; chains 11 (see **Submitting Cluster Jobs**). Needs that year's mosaic and footprint masks; Hirsh-Pearson is 2020-only, so earlier years run with `HF_MASK_YEAR=2020`. Locally: `Rscript Rscripts/07_train_and_backfill.R <subbasin>`. Sources 08A, 08B_*, 09 |
| `08A_train_and_backfill_subbasin_s.R` | sourced | `train_and_backfill_subbasin_s()`: loops over the biotic covariates in hierarchy order; backfills with the footprint covariates at 0 (OL #9) |
| `08B_deploy_gbart.R` | sourced | Gaussian BART for continuous covariates (log1p, 90/10 train/holdout) |
| `08B_deploy_mbart.R` | sourced | Multinomial BART for land-cover classes: modal class (metrics only) and 100 drawn classes `<cov>_class_###` |
| `08B_deploy_pbart.R` | sourced | Probit BART for where a covariate is blank; 08A blanks the backfill where P > 0.5 (`<cov>_blank`, `<cov>_pblank`) |
| `09_collect_metrics.R` | sourced | In-sample and holdout BART metrics (gbart, mbart) |
| `10A_process_backfill_metrics.R` / `10B_inspect_backfill_metrics.R` | local | Aggregate and inspect per-subbasin metrics and confusion matrices |
| `10C_abiotic_extrapolation_diagnostics.R` | local | Area of applicability (Meyer & Pebesma 2021) of each subbasin's training pixels in the natural abiotic predictors (footprint class excluded); flags a subbasin when > 50% of its backfilled pixels fall outside → `rds_files/extrapolation_flags.csv`. Needs `FNN` (not on Fir) |
| `11_premosaic_backfilled_stacks.R` + `.sh` | cluster | Mosaic the subbasin backfills into `bart_models_mosaics/{year}/{bcr}_backfilled.tif`; one task per BCR (one array per year). Rebuilds any mosaic older than its inputs (writes `_partial`, renames on completion). `MOSAIC_DIR` redirects for tests. Caps terra at 30% of the Slurm allocation (`memmax`): terra sizes blocks from the node's RAM |
| `12A_observed.R` + `.sh` | cluster (or local) | Elly's unclamped 32-bootstrap prediction tifs → V5's two-stage truncation (12B) → `observed_bootstraps.tif` (32 layers, unmasked), `observed_mean/sd.tif`, `truncation_params.rds` (frozen densmax and q99, keyed `{bcr}_{year}`). Canadian BCRs only. One task per species (16G / 8 h; 32G for several years). Skips BCR-years already done. Run locally (`SPECIES=OSFL Rscript Rscripts/12A_observed.R`) it reads G:, and its outputs must be Globus-transferred before 12D |
| `12B_v5_truncate.R` | sourced | `v5_truncate()`: port of V5 `analysis/10.Truncate.R` (V5 `f082866`): `densmax` cap, then `q99` (99.9th percentile), then range/water/extent masks |
| `12C_build_prediction_weights.R` + `.sh` | cluster | Per species × BCR `weight.tif` = range × not water × inside data limit × inside the BCR's own polygon, all `touches = TRUE` (OL #7). Reads `data/raw_data/v5_gis/` (each species needs `ranges/{species}.tif`, copied from G: `gis/ranges/`). Once per species before 12D. 12F multiplies obs and bf by it |
| `12D_repredict_all_coalitions.R` + `.sh` | cluster | One task per species × BCR, all 255 coalitions in one pass → `density_tables/by_bcr/{species}_{year}_{bcr}.rds` (written even when 12F skips the BCR, with `result = NULL`). The .sh chains 12H. Sources 12E–12G |
| `12E_shapley_utils.R` | sourced | Coalitions, Shapley weights (matrix form, and the textbook formula 14B checks it against) |
| `12F_predict_species_all_coalitions.R` | sourced | `predict_species_all_coalitions()`, one species × BCR file per call (see **12F invariants**): obs, `d0` (OL #9) and bf predictions, coalition sums, per-sample subbasin Shapley values `[n_sub × 8 × 3200]`, `pixel_change/{species}_{year}_{bcr}.tif` for 18 (checked against the cid-256 table). Also `coalition_task_table()`, `combine_bcr_results()` |
| `12G_gbm_tree_walk.R` + `.cpp` | sourced | Bit-identical, 3.6–4.4× faster gbm tree walk; `gbm_check_fast()` stops every worker unless it is `identical()` to `predict.gbm`. Also `coal_rowsum()` |
| `12H_merge_bcr_tables.R` + `.sh` | cluster | After 12D (afterany), same `SPECIES`. Checks every species × BCR is present (prints the indices to resubmit) and one `code_md5`, then binds in `06_bootstraps` order → `density_tables/{species}_{year}_coalition_{cid}.rds` (cid 2..256), `arrays/`, `{species}_{year}_shapley_samples.rds` |
| `14A_reproject_hirshpearson.R` | local | 8 sectors' 1 km layers from the raw 300 m Hirsh-Pearson rasters: `{sector}.tif` (max score) and `{sector}_direct.tif` (OL #10). After 03, before 12D; stage both |
| `14B_sector_attribution.R` + `.sh` | cluster | The .sh runs 14B then 14C (1 core, 16G / 4 h) and chains 18, one year per job. Sums the per-sample subbasin Shapley values bottom-up, sample by sample → `sector_effects/shapley_{subbasin,bcr,national}.csv` (mean, SD, 5/95%), each split into direct and vegetation parts. **BAM's release filter lives here**: every Canadian row of the "remove" tab of `data/raw_data/ModelReleaseDecisions.xlsx` (gitignored) is in `withheld_models` and dropped from the reported numbers only; a species with every BCR withheld (8) is skipped downstream too. Needs that xlsx and `extrapolation_flags.csv` on Fir |
| `14C_results_tables.R` | cluster | `results_summary.csv`, `results_sectors.csv`, `results_bcr.csv`, with 5/95%. **Sign opposite to 14B's**: change = obs − cf (negative = birds lost to industry); % is of the no-industry population, each sample over its own total. Stops unless its means equal 14B's to within 1 bird |
| `14D_master_species_list.R` | local | `data/raw_data/anna_data/master_species_list.csv` (gitignored): habitat group and `include` (113 reported). Excluded: BAM-withheld, Birds Canada's wetland and grassland lists, Anna's open list. Labels from Anna's classes, each BCR voting with the species' observed 2020 birds |
| `14E_sector_effects_birds.R` | local | `results_sectors_birds.csv`: 14C's sector table in birds |
| `15A_plot_population_distributions.R` | local | Population distributions from `density_tables/arrays/` (pull them from Fir first) |
| `17_detect_historical_hf_extent.R` | local | How much of the 2020 footprint existed in earlier years (needs the pre-2020 mosaics locally) |
| `18_map_population_change.R` + `.sh` | cluster | Chained after 14B, or alone with `SPECIES` → `output_figures/change_maps/change_grid_{species}_{year}.png`: obs − cf summed into 10 km cells, clipped to Canada (red = decline); stops unless the cells sum to 14B's national change. Titles from `data/raw_data/species_names.csv`. The only map type the user wants |
| `19A`–`19G` | 19A cluster, rest local | Rosenberg et al. (2019) Fig. 1-style panels over 1995–2020 → `output_figures/figure1/`, laid out in `figure1.pptx` (slide 1 birds, slide 2 % via `_pct.png`). Need pre-2020 results |
| `20_export_review_app.R` | local | Data bundle for `impactexploreR`, the co-authors' Shiny app (private sibling repo `mannfred/impactexploreR`; `APP_DIR` redirects). Reads 1995–2020 |
| `misc/` | local | Outside the pipeline, run from the project root: `harnesses/` (`verify_*`: 12F regression, weight vs V5 masking, truncation port, area units), `figures/`, and one folder per open investigation: `v5_upstream/` (PR1/PR2), `bart_mixing/`, `bc_valley_bottoms/` |

## Transferring Files to the Cluster (Globus)

CLI: `C:\Users\mannf\AppData\Local\Python\pythoncore-3.14-64\Scripts\globus.exe`.
Endpoints: local `a7878ccc-747b-11ef-b4b8-8fef73a45f39` (paths start `/C/`), Fir
`8dec4129-9ab4-451d-a45f-5b4b8471f7a3`.

- **Never use `--batch`** (always "permission denied"): one `globus transfer` per file, or one
  recursive transfer filtered with `--include=<pattern>` ... `--exclude=*`.
- **Always pass `--notify off`** (Globus emails once per task otherwise). Poll
  `globus task list --filter-status ACTIVE` until empty, then check the files landed.
- More gotchas (Windows globbing, Git Bash path rewriting, the scratchpad): `TODO.md`.

```powershell
$globus = "C:\Users\mannf\AppData\Local\Python\pythoncore-3.14-64\Scripts\globus.exe"
& $globus transfer `
  "a7878ccc-747b-11ef-b4b8-8fef73a45f39:/C/Users/mannf/Drive/.../Rscripts/12D_repredict_all_coalitions.R" `
  "8dec4129-9ab4-451d-a45f-5b4b8471f7a3:/home/mannfred/scratch/impact_assessment/Rscripts/12D_repredict_all_coalitions.R" `
  --label "my file" --sync-level checksum --notify off
```

## Submitting Cluster Jobs

Stage the entry point's **full `source()` closure** first (`TODO.md` "Staging discipline"). Each
`.sh` run with `bash` on a login node submits itself (never inside `salloc`). Put `YEARS` and
`SPECIES` in the environment, never in `--export` (sbatch splits it on commas).

```bash
# Backfill: per year, one array per memory tier over 1..674 (72G for subbasins 57, 61, 62, 98,
# 107 = BIG_SUBS, peak 67.3 GB; 24G otherwise; 8 h, 1 core), then 11 (can11/12/13/61 = BIG_BCRS
# at 128G / 8 h, the rest 32G / 4 h). Array indices must stay under MaxArraySize 1001.
bash 07_train_and_backfill.sh 2020
HF_MASK_YEAR=2020 bash 07_train_and_backfill.sh 1995 2000   # after 06 for those years
YEARS=2020 sbatch --array=<i> --mem=128G 07_train_and_backfill.sh   # an OOM-ed subbasin; add it to BIG_SUBS
YEARS=2020 sbatch --array=1-19 --mem=128G --time=08:00:00 11_premosaic_backfilled_stacks.sh
# 07 catches a subbasin's error and exits 0: grep the .out files for "Error in subbasin".
# Reruns need no wipe: 07 overwrites, 11 rebuilds mosaics older than their inputs.
# 11's index = position in its sorted bcr_vec (2,3,4,12 = can11, can12, can13, can61).

# Birds: species as arguments.
bash 12A_observed.sh OSFL GRSP
bash 12C_build_prediction_weights.sh OSFL GRSP   # independent of 12A
bash 12D_repredict_all_coalitions.sh OSFL GRSP   # 12D array + chained 12H; prints both job ids
# 12D: 8 cores / 3 h; can11 tasks 64G, the rest 48G (BIG_BCRS). MaxRSS median ~17 GB, peak ~45 GB;
# runtime median ~6 min, max ~66 min. A task killed for memory: resubmit it alone at --mem=64G.
# 12H prints the exact sbatch line for any missing BCR. An index past the task table exits 0.
# Fresh rerun of the SAME species: rm -f its files first (never rm -rf a directory). 12H refuses
# mixed code versions but cannot spot a stale file of the same version:
rm -f ../data/derived_data/density_tables/{CAWA,OVEN}_*.rds ../data/derived_data/density_tables/{arrays,by_bcr}/{CAWA,OVEN}_*.rds
# Many species: run_species_batches.sh (species_113.txt; LIST=<file> for another):
#   prep (12A + 12C), smoke (can42), plan, then wave 1, wave 2, ... (MAX tasks, default 400, <= 900)
TEST_BCR=can60 TEST_N_BOOT=2 bash 12D_repredict_all_coalitions.sh CAWA OVEN   # smoke; 12H gets the same
YEARS=1995,2000 MAX=900 bash run_species_batches.sh wave 1                    # other years

# Attribution, tables, maps. 14B reads every species in density_tables/: wait for the LAST 12H.
bash 14B_sector_attribution.sh <12H job id>
SPECIES=CAWA,OVEN sbatch 18_map_population_change.sh   # maps alone
# Pull sector_effects/*.csv and output_figures/change_maps/*.png, never density_tables/.
```

The DROP_Q99 / DENS_MIN sensitivity switches were removed 2026-11-13 (their runs predate C2d;
code at `git show e5024dc:Rscripts/12F_predict_species_all_coalitions.R`).

**`weight.tif` is a hard dependency** and the only masking: 12F stops if it is missing, stale, or
0/NA over the whole BCR (an unmasked run once over-counted CAWA can10 7.8×). Its band name
carries its version (`weight_v3_touches`); 12C rebuilds a stale one. A species whose range misses
the footprint (BAWW, BBWA can41) runs with impact 0.

## Fir Cluster

SFU, Digital Research Alliance of Canada. Login `fir.alliancecan.ca` (SSH is 2FA: the user runs
every cluster command). CPU nodes: 192 cores (8 per CCD, so `--cpus-per-task=8`), 750 GB RAM.
Max walltime 7 days, no partition flag. Always `--account=def-bayne` (`def-ecknight` is only the
project directory). `$HOME/scratch` (where `ia_dir` lives) is not backed up and is purged; keep
outputs in `$HOME/project/...`.

## Key Architecture Decisions

**Biotic hierarchy**: continuous biotic covariates are backfilled in a fixed order (00); each
BART model uses all preceding biotic covariates as predictors, so backfilled values cascade.

**Categorical covariates**: `ABoVE_1km`, `NLCD_1km`, `MODISLCC_1km`, `MODISLCC_5x5`, `SCANFI_1km`,
`VLCE_1km` use `mbart()` and are not predictors of continuous covariates.
- No water: 07's `water_classes` (SCANFI_1km 8, VLCE_1km 20, MODISLCC 17) are left out of
  training, so the backfill never predicts water on land.
- 12F uses the drawn class `<cov>_class_###` with the scenario's draw index, never the modal
  class (which erases minority classes).
- Class probabilities are `pnorm(z)/Σpnorm(z)` from `yhat.test`, not `mbart2`'s `prob.test` (a
  softmax of probit scores that roughly doubles rare classes). Classes are renumbered 1..K.

**Blanks are modelled, not filled.** Some V5 covariates are missing by design
(`SCANFIheightcv_*` where canopy height is 0) and gbm sends NA down its missing branch. 08B's
pbart models where a covariate is blank; a blanked backfill is NA to later covariates (their
`_isNA` flag is 1), and 12F passes it to gbm as `NA_real_` and counts it as complete.

**Coordinates**: lat/lon are centred and scaled within each subbasin before entering BART.

**Output per subbasin**: `bart_models/{year}/subbasin_{i}/subbasin_{i}_{backfill.tif,metrics.rds,confusion.rds}`.
A covariate constant across a subbasin's training pixels skips BART and gets `<cov>_mean` /
`<cov>_sd`; 12F fills those pixels from `<cov>_mean` (raw scale, no `expm1`) for every draw, and
a covariate with no draws anywhere in a BCR is carried as one draw from `<cov>_mean`.
`BART_MODELS_DIR` redirects 07 for tests (11 always reads `bart_models/`).

**Re-prediction**: for coalition S, pixels where any sector in S has footprint (and CanHF ≥ 1) get
the footprint covariates at 0; those in the direct footprint of a sector in S also get backfilled
vegetation (OL #10); all other pixels keep observed values. BART draws are nested inside BRT
bootstraps (joint sampling).

**12F invariants**: the joint-sampling seed depends only on species, BCR, bootstrap `i` and
scenario `k`, never the coalition, so the backfilled field is computed once over the
all-sector superset and each coalition is a masked `rowsum`. Must hold:
- One complete-case mask (`keep`) drives obs and bf.
- BART draws are log1p-scaled → `expm1`; non-finite → NA, never 0.
- Categorical `var.levels` are indexed by `match()` against `var.names`, and converted with
  EACH bootstrap's own levels (`as_model_factors()`): bootstraps of one model know different
  classes. The complete-case gate tests the raw backfilled class, never a factor.
- Both V5 caps per pixel, feeding the tables: `pmin(pmin(pred, densmax), q99)`, with `q99`
  frozen from the observed landscape (12A's `truncation_params.rds`), never re-derived from a
  counterfactual. `q.out` has no `$q` column (`pmin(x, NULL)` is silently `numeric(0)`); 12F
  guards the schema.
- Missing covariates reach gbm as `NA_real_`, never `NaN` (terra returns `NaN`, which gbm sends
  right instead of down the missing branch); 12F recodes where `X_obs_super` is built.
  `SurfaceWater_1km` (gbm levels "0","1") is passed raw, as V5 did.
- Identity gate: 12F first predicts every bootstrap from the observed design through the
  workers' path and stops unless it reproduces `observed_bootstraps.tif` (float32 tolerance).
- Workers return per-coalition `[n_sub × n_scen]` sums (`coal_rowsum()`), not pixel fields;
  bit-identical to reducing the full matrix.
- Forked workers inherit the parent's GC trigger: they allocate only per-draw vectors, call
  `gc(full = FALSE)`, and run one bootstrap each (`mc.preschedule = FALSE`). A killed worker
  returns `NULL` and 12F names it.
- Grouped `_draw_*` reads use INTERLEAVE=BAND.

**Shapley attribution**: 8 sectors → 255 non-empty coalitions; exact Shapley values (12E, 14B)
sum to the total impact.

**Bottom-up aggregation, per sample**: subbasin → BCR → national. 14B sums each
(bootstrap, scenario) sample's subbasin Shapley values (bootstrap `i` paired across BCRs) and
takes SDs and percentiles of the sums; never add variances (a BCR's subbasins share 32 bird
models). A % interval divides each sample by its own total.

## Data Directory Structure

```
data/
├── raw_data/
│   ├── biotic_variable_hierarchy.rds
│   ├── CAfire/CAfire_{year}_masked.tif        # 06's input. G:'s 1995–2015 copies are pre-fix: stage from here
│   ├── isric_soil/isric_soil_covariates_masked.tif
│   ├── covariates_mosaiced/covariates_mosaiced_{year}.tif   # 2020 locally; 1995–2020 on Fir
│   ├── hirshpearson/
│   │   ├── CanHF_1km_lessthan1.tif / CanHF_1km_morethan1.tif
│   │   └── {sector}.tif / {sector}_direct.tif   # 14A (raw_300m/ deleted locally; re-download to rebuild 03/14A)
│   ├── hydrobasins_masked_merged_subset.gpkg
│   ├── species_names.csv
│   ├── ModelReleaseDecisions.xlsx               # BAM's; gitignored
│   ├── v5_gis/{ranges/{species}.tif, WaterMask_Canada.*, DataLimitationsMask.*}
│   └── Regions/BAM_BCR_NationalModel_Unbuffered.shp
└── derived_data/
    ├── bart_models/{year}/subbasin_{i}/
    ├── bart_models_mosaics/{year}/{bcr}_backfilled.tif
    ├── density_tables/
    │   ├── {species}_{year}_coalition_{cid}.rds    # subbasin rows per coalition
    │   ├── {species}_{year}_shapley_samples.rds    # 12H → 14B
    │   ├── arrays/                                 # 12H → 15A
    │   ├── pixel_change/{species}_{year}_{bcr}.tif # 12F → 18
    │   └── by_bcr/                                 # 12D → 12H
    ├── predictions/{species}/
    │   ├── truncation_params.rds
    │   └── {bcr}/{year}/{observed_bootstraps,observed_mean,observed_sd,weight}.tif
    ├── sector_effects/                             # 14B, 14C, 14E CSVs (local, never committed)
    └── rds_files/                                  # metrics, extrapolation_flags.csv (local, never committed)
```

Coalition table columns: `species`, `subbasin`, `bcr`, `coalition_id`, `obs_total_mean/sd`,
`obs_on_coalition_mean/sd`, `bf_on_coalition_mean/sd`. A subbasin with no kept pixels has
`obs_on_coalition_mean = NaN`, `_sd = NA`, `bf_on_coalition_mean = 0`: key on the NaN mean.
Large spatial files and most `.rds` are gitignored; `writing/` is private and never committed.

## Open Limitations

Numbers are kept stable because code and notes cite them.

1. ~~Abiotic extrapolation~~: flagged (not corrected) by 10C's area of applicability, annotated
   in 14B.
2. ~~Non-additive sector impacts~~: exact Shapley values over all 2^8 coalitions. They assume
   every v(S) is well estimated; large combined footprints carry more extrapolation risk.
3. **No spatial spillover**: `cf = obs off the coalition + backfilled on it` assumes removing a
   footprint affects only its own pixels. Edge effects, area sensitivity and connectivity extend
   impacts beyond it (especially linear features).
4. ~~Ad hoc uncertainty combination~~: joint BRT × BART sampling in 12F.
5. ~~Backfill coverage collapse~~: fixed by 02's CAfire `NA → 0`, 08A's median imputation of
   partial-NA covariates with an `_isNA` flag (Fix B), and 12F gating only on backfilled
   covariates (Fix A). Every BCR keeps ≥ 99.9% of weight > 0 footprint pixels.
6. **The frozen q99.9 cap biases impact slightly downward**: counterfactual densities hit the
   observed-landscape ceiling more often. Measured (`DROP_Q99=1`, pre-C2d): v(N) +1.16% (CAWA),
   +0.05% (OVEN), no sector share moves; q99 binds at 0–1% of observed footprint pixels. Freezing
   is mandatory: a per-counterfactual cap would enter the contrast.
7. ~~Two `weight.tif` masking defects~~, fixed; easy to reintroduce:
   (a) V5's pre-mosaic prediction grids are buffered 100 km past their subunit, so 12C crops to
   the unbuffered BCR polygon (else straddling subbasins count under every BCR, 2.2–14.6×).
   (b) `terra::rasterize()` defaults to `touches = FALSE`; V5's masks behave as `TRUE`.
   `weight.tif` must be non-NA everywhere (V5's `mask.i[is.na(mask.i)] <- 0`): every new mask
   term needs `background = 0`, and 12C stops on any NA. Harness:
   `misc/harnesses/verify_weight_vs_v5_masking.R`.
8. **V5's published abundances read ~3% low (upstream)**: V5 sums `density × 100` in EPSG:3978,
   whose 1000 m cells hold ~1.027 km²; our EPSG:5072 totals are correct. It cancels in obs − cf.
   Never deliver in 3978: `clamp()` does not commute with projection.
9. **Footprint covariates in the counterfactual are set to 0.** V5's "Disturbance" class
   (CanHF_1km/_5x5, canroad_1km/_5x5, CCNL_1km) is human footprint, which the counterfactual
   lacks. (a) Bird model: 12F also predicts each bootstrap with observed vegetation and the class
   at 0 (`d0`), splitting every Shapley value exactly into a direct part (`d0 − obs`, the bird
   model's response to footprint itself) and a vegetation part (`bf − d0`). Positive responses to
   modest disturbance can be real (Mahon et al. 2019, Ecol. Appl. 29:e01895), so the direct part
   is reported, not removed; it can be large and opposite in sign (OVEN roads: −1.27 M direct,
   +1.74 M vegetation). Its support is thin: in the can11/12/13 models ≥ 90% of surveys have
   CanHF_1km ≥ 6–7 and canroad_5x5 ≥ 0.3–0.5. The 5×5 covariates are exact only for the
   all-sector coalition. (b) BART trains on observed values and backfills with the class at 0
   (inside the training range); 10C excludes the class.
10. **Two masks per sector.** Hirsh-Pearson scores roads, rail, mines, reservoirs and oil and gas
    by distance band, so a sector's layer includes influence zones. Footprint = `{j}.tif > 0`
    (and CanHF ≥ 1): footprint covariates at 0 on observed vegetation (`d0`). Direct footprint =
    `{j}_direct.tif == 1`, a subset where vegetation is also backfilled: any 300 m cell with
    built/crop/pasture > 0, roads/rail/dams ≥ 6, mines ≥ 6, oil and gas = 10. Mines and oil and
    gas are mapped as points, so their direct area is a disc (4.8% and 7.7% of their footprint)
    and converted land is understated; no polygon data exist to fix it. Whether influence-zone
    vegetation is unaltered was never measured. `bf_on_coalition_*` holds the combined
    counterfactual. 03 and 14A aggregate 300 m → 1 km by max / area mean, never bilinear.
