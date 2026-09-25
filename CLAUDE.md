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

**Colleague's NationalModels directory** (`nm_root`): scripts 12D–12F and 13 reference
`/home/mannfred/projects/def-ecknight/NationalModels` for BRT bootstrap models and BCR
covariate stacks. This path is cluster-only and must not be changed. **12A is an exception**:
when `cc=FALSE`, `nm_root` is overridden to `G:/Shared drives/BAM_NationalModels5` so the
script reads Elly's bootstrap models and raw prediction tifs from the local Google Drive.

Most prep scripts (01–06) run locally. Compute-heavy scripts (07, 12D) run on the cluster via SLURM array jobs.

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
| `07_train_and_backfill.R` + `.sh` | cluster | Entry point: for a given SLURM array index (subbasin), train BART models and backfill high-HF pixels. ONE array job over all 674 subbasins at 128G / 8h / 1 core. The former 64G+750G tiering (and its second script, `07_train_and_backfill_larger.sh`) was an artifact of the `08A` assembly OOM fixed in `73e8b68`, not of subbasin size |
| `08A_train_and_backfill_subbasin_s.R` | sourced | Core `train_and_backfill_subbasin_s()` function; loops over biotic covariates in hierarchy order |
| `08B_deploy_gbart.R` | sourced | `deploy_gbart()`: Gaussian BART for continuous biotic covariates (log1p-transformed, 90/10 train/holdout split) |
| `08B_deploy_mbart.R` | sourced | `deploy_mbart()`: Multinomial BART for categorical land-cover covariates |
| `09_collect_metrics_gbart/mbart.R` | sourced | Collect in-sample BART metrics |
| `09_collect_holdout_metrics_gbart/mbart.R` | sourced | Collect holdout BART metrics |
| `10A_process_backfill_metrics.R` | local | Aggregate per-subbasin metrics and confusion matrices into CSVs |
| `10B_inspect_backfill_metrics.R` | local | Inspect and visualize model accuracy |
| `10C_abiotic_extrapolation_diagnostics.R` | local | Area of applicability (Meyer & Pebesma 2021) of each subbasin's low-HF training pixels in the natural abiotic predictors 08A trains BART on (its footprint "Disturbance" class excluded; see Open Limitation #9); flags a subbasin when > 50% of its backfilled pixels fall outside → `extrapolation_flags.csv`. Needs `FNN` (installed locally, not on Fir) |
| `10D_BART_posterior_diagnostics.R` + `.sh` | cluster | Test whether the 100 stored BART posterior draws (subsampled from gbart()'s 700) cover the posterior's shape and tails — `12F` resamples 1 of the 100 per counterfactual scenario |
| `11_premosaic_backfilled_stacks.R` + `.sh` | cluster | Mosaic per-subbasin BART backfill rasters into BCR-wide stacks (run before 12) |
| `12A_observed.R` | local | Reads Elly's unclamped 32-bootstrap prediction tifs and bootstrap model `.Rdata` from `G:/Shared drives/BAM_NationalModels5/output/{07_predictions,06_bootstraps}/{species}/`, applies V5's two-stage truncation (via `12B_v5_truncate.R`), and writes `observed_bootstraps.tif` (32-layer clamped stack, UNmasked), `observed_mean.tif`, `observed_sd.tif` and per-species `truncation_params.rds` to `data/derived_data/predictions/`. Canadian BCRs (`can*`) only. Globus-transfer the results to the cluster before running 12D — they are a hard dependency. |
| `12B_v5_truncate.R` | sourced | `v5_truncate()`: line-for-line port of V5 `analysis/10.Truncate.R` (as of V5 `f082866`). Applies both upper caps (`densmax`, then the 99.9th-percentile `q99`) and the range/water/extent masks. Sourced by `12A` only, so it is local-only by design and is deliberately NOT staged on the cluster; `q99` can be passed in frozen and the legacy EPSG:3978 step skipped |
| `12C_build_prediction_weights.R` + `.sh` | cluster | Build per-species×BCR `weight.tif` (= range membership × not-water × inside-data-limit × **inside the BCR's own polygon**), replicating V5 `10.Truncate` range/water/extent/mosaic masking; all four terms use `touches = TRUE` (Open Limitation #7). Reads source masks from `data/raw_data/v5_gis/` (no G: access); grid template is the BCR stack. Run ONCE before 12D. 12F multiplies BOTH observed and backfilled density by this weight, preserving obs/bf symmetry (`w·bf − w·obs = w·(bf − obs)`). |
| `12D_repredict_all_coalitions.R` + `.sh` | cluster | Entry point: re-predict bird densities for ALL 255 coalitions, ONE SLURM array task per species × BCR (`coalition_task_table()` in 12F: every species' `06_bootstraps` files in `list.files()` order; 25 tasks for CAWA + OVEN). Sources 12E, 12F, 12G (and compiles `12G_gbm_tree_walk.cpp`); writes one `density_tables/by_bcr/{species}_{year}_{bcr}.rds` per task (atomically; a BCR 12F skips is still written, with `result = NULL`). Requires `observed_bootstraps.tif` present (hard error if missing — observed rasters are now always staged on the cluster). |
| `12G_gbm_tree_walk.R` + `.cpp` | sourced | Bit-identical replacement for gbm's compiled tree walk (`gbm_pred`), 3.6–4.4× faster: gbm makes three R API calls per node visited, this hoists them. `gbm_design()` builds predict.gbm's design matrix once per bootstrap; `gbm_check_fast()` `stop()`s in every 12F worker unless the fast path is `identical()` to `predict.gbm`. The `.cpp` also holds `coal_rowsum()`, 12F's copy-free, bit-identical coalition `rowsum` |
| `12H_merge_bcr_tables.R` + `.sh` | cluster | Run after 12D (`--dependency=afterany`). Checks every expected species × BCR has a per-BCR result (names the array indices to resubmit if not), that all came from one version of 12E/12F/12G (`code_md5`), then binds them in `06_bootstraps` order → `density_tables/{species}_{year}_coalition_{cid}.rds` (cid 2..256) and `arrays/`, bit-identical to the old one-job-per-species output, plus `{species}_{year}_shapley_samples.rds` (per-sample subbasin Shapley values, stacked across BCRs in table row order) |
| `12E_shapley_utils.R` | sourced | Coalition enumeration, Shapley value computation utilities |
| `12F_predict_species_all_coalitions.R` | sourced | `predict_species_all_coalitions()`: builds the backfilled field ONCE per species×BCR over the all-8-sectors superset, then reduces all 255 coalitions as cheap masked `rowsum`s (verified bit-identical to the retired per-coalition path). Runs joint BRT×BART sampling; each bootstrap worker reduces its own field to per-coalition subbasin sums, so the full pixels × 3200 matrix is never built. `rdata_files` / `return_per_bcr` let 12D run one BCR; `combine_bcr_results()` (used by 12H) does the cross-BCR bind. Also applies `shapley_weight_matrix()` (12E) to every (bootstrap, scenario) sample, giving per-sample subbasin Shapley values `[n_sub × 8 × 3200]` for 14B's uncertainty (Shapley is linear in v, so the samples' mean equals the Shapley value of the tables' means). Also holds `coalition_task_table()` and `coalition_array_ids()`. |
| `13_importance_of_covs_used_in_counterfactual.R` | cluster | Assess percentile importance of backfilled covariates in V5 bird models |
| `14A_reproject_hirshpearson.R` | local | Reproject the per-sector Hirsh-Pearson footprint rasters (built, crop, mines, …) to EPSG:5072 on the hydrobasins grid at 1000 m. Run once before `14B`; `CanHF*` left untouched |
| `14B_sector_attribution.R` | local | Reads 12H's per-sample subbasin Shapley values, sums them bottom-up (subbasin → BCR → national) SAMPLE BY SAMPLE and summarises each level (mean, SD, 5th/95th percentiles) → `sector_effects/shapley_*.csv`; stops unless the sample means equal the Shapley values of the coalition tables' means. **This is where BAM's release filter lives** (`DROP_WITHHELD` / `withheld_models`): CAWA `can40` is withheld by BAM (`review/ModelReleaseDecisions.xlsx`, "remove" tab, AUC) but is deliberately still produced upstream, so the products stay a complete record of what we ran and only the reported numbers are filtered. Until 2026-09-25 it propagated the tables' SDs as if subbasins (which share one BCR's 32 bird models) and nested coalitions were independent, which understated v(S) SDs 1.1–2.7× and inflated small sectors' SDs |
| `15A_plot_population_distributions.R` | local | Plot empirical population distributions, observed vs counterfactual, from the bootstrap × scenario arrays `12D` saves to `density_tables/arrays/` |
| `15B_preliminary_singletons.R` | local | Interim pre-Shapley single-sector standalone impacts at footprint / watershed / BCR scales → `logs/15B_singletons_summary.csv`. Superseded by `14B`'s exact Shapley values for attribution; retained for the scale-dependent reporting `15C` plots |
| `15C_singletons_plot.R` | local | Render `15B`'s summary CSV as sector impacts at the three geographic scales |
| `16A_lowHF_isnot_a_proxy.R` | local | Test the premise that low-HF pixels are not a spatial proxy for the pre-industrial state |
| `16B_map_sampled_subbasins.R` | local | Map the 16 subbasins `16A` samples, reproducing the index from the same seed without re-running the RF pipeline |
| `17_detect_historical_hf_extent.R` | local | Detect whether the 2020 HF mask was smaller in earlier years, via PCA distance between observed and backfilled biotic vectors inside the HF mask |
| `misc/` | local | Downstream analysis (population summaries, visualization, vegetation vs. mines) plus one-off diagnostics/profiling not in the execution-order pipeline (`diag_*`, `analyze_seff`, `check_model_complexity`, `debug_*`, `figure_*`) |

## Transferring Files to the Cluster (Globus)

Use the Globus CLI at `C:\Users\mannf\AppData\Local\Python\pythoncore-3.14-64\Scripts\globus.exe`.

**IMPORTANT: never use `--batch` mode** — it always returns "permission denied". Submit one `globus transfer` call per file instead.

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

## Submitting Cluster Jobs

Before any cluster run, stage the entry point's **full `source()` closure**, not just the files
you edited — see `TODO.md` "Staging discipline".

Backfilling (SLURM array, one job per subbasin index):
```bash
# ONE script, ONE tier: --array=1-674%30 at 128G / 8h / 1 core, all baked in.
sbatch 07_train_and_backfill.sh

# The old two-script split (64G + 750G, with complementary --array lists that had to partition
# 1-674 exactly or race each other writing the same subbasin_{i}_backfill.tif) is GONE as of
# 2026-09-19. It existed only because 08A's layer-by-layer assembly churned ~nlyr^2 * ncell * 8
# bytes through the allocator; after 73e8b68 the worst subbasin in the set peaks at 50.8 GB, so
# a single 128G tier covers everything with ~2.5x margin — and costs less in total than the two
# tiers did. Single-threaded by construction (BART::gbart, not mc.gbart), so --cpus-per-task=1.
#
# If one index OOMs, give that index more room instead of re-tiering the file:
sbatch --array=<i> --mem=256G 07_train_and_backfill.sh
```

Re-predicting birds:
```bash
# Phase 1 (LOCAL): Rscript Rscripts/12A_observed.R, array task 1 then 2. Reads G:, writes
# observed_bootstraps.tif + truncation_params.rds to data/derived_data/predictions/.
# Globus-transfer both to the cluster BEFORE 12D — 12F stops with an error if either is missing.

# Phase 1b: build prediction weights ONCE before 12D.
# Reads data/raw_data/v5_gis + Regions/BAM_BCR_NationalModel_Unbuffered.shp; writes
# weight.tif per species x BCR.
sbatch 12C_build_prediction_weights.sh   # --array=1-<n_species>
# weight.tif is a HARD dependency: 12F stops if it is missing, and rejects an all-zero/NA one.
# Since Fix B median-imputes partial-NA covariates, water / out-of-range / out-of-extent pixels
# no longer drop out via complete.cases() on their own, so weight.tif is the only masking left.
# An unmasked run is quietly wrong, not obviously broken (7.8x over-count on CAWA can10).
#
# VERSIONED: weight.tif carries its version in its band name ("weight_v3_touches"). 12C rebuilds
# any weight whose stamp does not match instead of skipping it, and 12F refuses to run against a
# stale one, so no manual `rm` is ever needed — just re-run 12C.

# Phase 2: ONE array task per species x BCR computes ALL 255 coalitions for that BCR in a
# single pass (superset restructure — see "12F restructure invariants" below), then 12H merges.
# 12D_repredict_all_coalitions.sh fixes --array=1-25 (CAWA + OVEN) at 8 cores / 48G / 03:00:00;
# the task table is printed at the top of every 12D log. An index past the table's end exits 0.
# Sized from C1 (2026-09-24): peak MaxRSS 28.8 GiB (OVEN can61), longest task 42 min (CAWA can11).
# Memory is what Alliance bills here (48G on 8 cores = 12.3 core-equivalents), so it is not
# padded further; a task killed for memory is resubmitted alone with --mem=64G.
#
# For a FRESH full re-run, delete stale outputs first. `rm -f` on files, never `rm -rf` on a
# directory: 15A reads arrays/. 12H refuses to merge per-BCR files from different code versions,
# but it cannot tell a stale file of the SAME version from a fresh one, so wipe by_bcr/ too.
rm -f ../data/derived_data/density_tables/*.rds
rm -f ../data/derived_data/density_tables/arrays/*.rds
rm -f ../data/derived_data/density_tables/by_bcr/*.rds
sbatch 12D_repredict_all_coalitions.sh                                  # note the job id
sbatch --dependency=afterany:<12D job id> 12H_merge_bcr_tables.sh      # tables + arrays

# If 12H reports missing BCRs it prints the exact `sbatch --array=<indices> ...` to resubmit.

# Smoke test (can60, 2 bootstraps). TEST_BCR filters the task table to every species' can60,
# so it has TWO tasks (1 = CAWA, 2 = OVEN) and 12H expects both. Give 12H the same TEST_BCR
# AND TEST_N_BOOT: it refuses to merge smoke files under production settings.
sbatch --array=1-2 --time=01:00:00 --export=ALL,TEST_BCR=can60,TEST_N_BOOT=2 12D_repredict_all_coalitions.sh
sbatch --dependency=afterany:<smoke job id> --export=ALL,TEST_BCR=can60,TEST_N_BOOT=2 12H_merge_bcr_tables.sh
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
- Always specify `--account=def-bayne` (the allocation every job script in `Rscripts/` uses — not `def-ecknight`, which is only the project *directory* path)
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

**Output per subbasin**: `data/derived_data/bart_models/{year}/subbasin_{i}/subbasin_{i}_backfill.tif`, `_metrics.rds`, `_confusion.rds`. Note **not every continuous biotic covariate has `_draw_*` layers**: where a covariate is constant across a subbasin's low-HF pixels, `08A` skips BART and writes `<cov>_mean`/`<cov>_sd` instead (all-NA if there were no valid training rows). This is NOT harmless by default: the BCR mosaic unions layers across subbasins, so a covariate with draws in some subbasins is NA in the constant ones, and `12F`'s complete-case gate dropped them (CAWA can60: 43% of weight > 0 footprint pixels, all from `SCANFIBalsamFir_5x5`). `12F` now fills draw-less pixels from `<cov>_mean` (raw scale — no `expm1`) for every draw. A covariate with no draws anywhere in the BCR (constant in every subbasin) is carried as a one-draw covariate from `<cov>_mean`, so it too gets the low-HF constant rather than its observed value.

**Re-prediction**: `11_premosaic` mosaics backfilled subbasin rasters into BCR-wide stacks. `12A_observed.R` runs locally and writes the truncated `observed_bootstraps.tif` per species×BCR (see the pipeline table); these are Globus-transferred to the cluster, where they are a hard dependency for `12D`/`12F`. `12D_repredict_all_coalitions.R` runs ONE array task per species × BCR and computes all 255 coalitions in a single pass: it sources `12F_predict_species_all_coalitions.R`, which builds the backfilled field ONCE over the BCR's all-8-sectors superset and reduces every coalition as a cheap masked `rowsum`; `12H` merges the per-BCR files. For a given coalition S of sectors, pixels where any sector in S has footprint (AND CanHF ≥ 1) use backfilled covariates; all other pixels use observed. Joint BART×BRT sampling nests BART posterior draws inside BRT bootstrap iterations.

**12F restructure invariants**: The joint-sampling seed depends only on species, BCR, bootstrap `i`, and scenario `k` — never the coalition. So for fixed (species, BCR, i, k) the backfilled density field over pixels is identical across all 255 coalitions; the coalition only selects which pixels are masked in, never the backfilled value at a pixel. This is what licenses computing the backfilled field ONCE over the all-8-sectors superset and reducing each coalition as a cheap masked `rowsum`. Correctness constraints that must hold for the restructure to stay bit-identical to the (retired) per-coalition path:
- The complete-case mask is identical between obs and bf (a single `keep` drives both).
- BART draws are log1p-scaled → `expm1` before use; non-finite → NA, never 0.
- Categorical `var.levels` are indexed by `match()` against `var.names`, never by name.
- BOTH V5 upper caps are applied per pixel in the gbm step and DO feed the density table: `pmin(pmin(pred_vec, qsp), q99)` (12F:350, where `qsp` is the densmax cap). `densmax` is `q.out$densmax` (the old `$q` column no longer exists); `q99` is the per-BCR 99.9th percentile **frozen from the observed landscape** by 12A and read from `predictions/{species}/truncation_params.rds`. It must never be re-derived from a counterfactual, or the cap enters the obs/bf contrast. The observed side is pre-clamped by 12A, so both sides carry identical caps.
- Grouped `_draw_*` reads use INTERLEAVE=BAND (now once per BCR).
- **Missing covariates reach gbm as `NA_real_`, never `NaN`.** `terra::values()` returns missing
  cells as `NaN`, and gbm's tree walk sends ONLY `NA_real_` down a split's missing branch (a
  `NaN` fails the `<` test and goes right). V5's observed predictions used the missing branch,
  so until 2026-09-24 the bf side was predicted under a different rule wherever an
  observed-only covariate was missing: CAWA can11, 3.2% of weight > 0 footprint pixels
  (`StandardGreenup`/`StandardDormancy`), ~35× V5's density there, ≈ +10% of the observed
  birds on that footprint as a fake impact. 12F recodes `NaN → NA_real_` where `X_obs_super`
  is built. `SurfaceWater_1km` is gbm-categorical (`var.type` 2, levels "0","1") but not in
  `categorical_responses`, so it is passed raw as V5 did; that is correct only because its
  levels are exactly "0","1", and a raw `NaN` there is undefined behaviour in `gbm_pred`.
- **Identity gate.** Before sampling, 12F predicts EVERY bootstrap from the OBSERVED design
  (nothing backfilled) through the workers' exact path and `stop()`s unless it reproduces
  `observed_bootstraps.tif` (float32 tolerance), over-sampling pixels with a missing covariate or
  with a class some bootstrap never saw. It catches any obs/bf rule divergence, and a
  `b.list[[i]]` ↔ observed layer `i` mis-pairing. (It tested only bootstraps 1 and 32 until
  2026-09-24; the categorical-levels bug below showed that was too thin.)
- **Per-bootstrap reduction.** Workers return per-coalition `[n_sub × n_scen]` sums, not their
  pixel fields. `rowsum()` adds each column independently in row order, so this is bit-identical
  to reducing the assembled `[n_super × (n_boot·n_scen)]` matrix the old code built (up to
  ~17 GB for can11, held 2–3× over at peak). The sum itself is `coal_rowsum()` (12G `.cpp`):
  R's `rowsum` loop over the kept rows in place, so it never copies `sc[kr, ]` 255 times.
- **Categorical levels are PER BOOTSTRAP.** Each bootstrap is fitted to a different resample, so
  `var.levels` differ between bootstraps of one model (CAWA can14: 28 of 32 know VLCE_1km classes
  32/40 that bootstrap 1 never saw). `as_model_factors()` converts with each model's own levels,
  as `predict.gbm` and V5 do. Until 2026-09-24 12F used bootstrap 1's for all 32, so those
  classes went down the missing branch in other bootstraps (identity gate, can14 bootstrap 32:
  155 cells off V5) and a backfilled class unknown to bootstrap 1 failed the complete-case gate.
  The gate now tests the RAW backfilled class for NA, never a factor.
- **Worker memory.** A forked worker inherits the parent's GC trigger, which the BART-draw reads
  push to ~10–15 GB on the big BCRs, so 8 workers each piling up garbage OOM-ed C1 at 64G
  (can12/13/80). Workers allocate only per-draw vectors, collect them with `gc(full = FALSE)`
  (young objects only: cheap, and it leaves the parent's pages alone), and each fork runs ONE
  bootstrap (`mc.preschedule = FALSE`). A worker the kernel kills returns `NULL`; 12F names it.

**Shapley attribution**: 8 sectors → 256 coalitions (2^8; cid 1 = empty is skipped → 255 computed). All 255 are produced in one pass per species × BCR (12D), merged by 12H. `14B_sector_attribution.R` computes exact Shapley values from the coalition density tables. Shapley values sum exactly to the total HF impact. `12E_shapley_utils.R` provides coalition enumeration and the Shapley formula.

**Bottom-up aggregation**: The subbasin is the atomic spatial unit. BCR totals = sum of subbasin values within the BCR. National totals = sum of BCR values. Uncertainty is aggregated per sample, never by adding variances: 14B sums each (bootstrap, scenario) sample's subbasin Shapley values within a BCR, then across BCRs (bootstrap i paired across BCRs, as V5's national mosaic pairs them), and takes the SD and percentiles of those sums. The SD columns of the coalition tables are per-subbasin only and must not be combined as if independent: every subbasin of a BCR shares its 32 bird models.

**Abiotic extrapolation diagnostics**: `10C` computes each subbasin's area of applicability: a backfilled pixel is outside it when its standardised distance to the nearest training pixel exceeds what training pixels show to each other across spatial blocks. `frac_outside_aoa` and the flag (> 0.5) are annotated in `14B`'s subbasin output, and `14B` reports the impact carried by flagged subbasins per BCR and nationally. (Until 2026-09-25 it used KS + Mahalanobis rules, which flagged 667 of 667 subbasins: KS measures shift, not extrapolation, and Mahalanobis over ~40 collinear climate covariates inverts a near-singular covariance.)

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
    │   ├── {species}_{year}_coalition_{id}.rds   # subbasin-level population arrays per coalition
    │   ├── {species}_{year}_shapley_samples.rds  # per-sample subbasin Shapley values (12H → 14B)
    │   └── arrays/                               # per-pixel bootstrap × scenario arrays (12D → 15A)
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

## Open Limitations (updated 2026-09-15)

### Conceptual

1. ~~**Abiotic overlap / extrapolation risk**~~: **Addressed** by `10C_abiotic_extrapolation_diagnostics.R`. A per-subbasin area-of-applicability diagnostic flags subbasins where BART extrapolates in the natural abiotic predictors. Flags are annotated in 14B output. Remaining gap: diagnostics are post-hoc; they don't correct the extrapolation, only flag it.

2. ~~**Sector impacts are not additive**~~: **Addressed** by Shapley value attribution. 12D/12F run all 2^8 = 256 sector coalitions; 14B computes exact Shapley values that sum to total HF impact. Remaining gap: Shapley values assume the coalition value function v(S) is well-estimated for all S; coalitions with large combined footprints may carry more extrapolation risk.

3. **No spatial spillover**: The formula `cf = obs_on_non_coalition + backfilled_on_coalition` assumes removing a coalition's footprint only affects birds on those pixels. Edge effects, area sensitivity, and functional connectivity mean impacts extend beyond the footprint boundary (especially important for linear features like roads and seismic lines).

9. **Footprint covariates in the counterfactual (OPEN, design decision; `TODO.md` C2e)**. V5's
   "Disturbance" predictor class (CanHF_1km/_5x5, canroad_1km/_5x5, CCNL_1km night lights) is
   human footprint, and it enters the counterfactual twice, inconsistently.
   (a) **Bird model.** 12F sets it to 0 at every backfilled pixel. That term alone can dominate an
   impact: OVEN can12's −828k is −1,041k from zeroing it (observed vegetation kept) plus +213k
   from the backfilled vegetation, although the class holds only 2.8% of the model's relative
   influence. Zero footprint is outside the data on footprint pixels.
   (b) **BART.** 07's `abiotic_vars` includes the class, so 08A trains on it, then backfills with
   each pixel's OBSERVED values, which lie beyond the low-HF training range (subbasin 57:
   CanHF_1km 0–10 in training, median 12 at backfilled pixels). Trees hold predictions at the
   training edge, so the backfill is conditioned on the most-disturbed training pixels, not on
   no industry. A local test (2026-09-25) fitted BART as 08A does (without the preceding biotics) in subbasins 57 (can11) and 98 (can12), for deciduous %, canopy height and biomass. Backfilling with the observed footprint values, instead of dropping the class, moved the backfilled means by −15% to +32%; setting the values to 0 moved them by −6% to +73%. The class took 4–10% of BART's splits (`Rscripts/misc/diag_extreme_bcr_bart_footprint.R`).
   `10C` excludes the class from its extrapolation diagnostic for this reason.

10. **What counts as a sector's footprint (OPEN, design decision; `TODO.md` C2f)**. A pixel is in
    sector j's footprint when j's Hirsh-Pearson pressure is > 0 after `14A` and CanHF ≥ 1. The
    pressure layers are continuous and include indirect-influence zones: roads.tif is > 0 on 19%
    of Canada's cells, median 4.55, and roads-only footprint pixels have a median of 2.7–3.4
    against a direct score of 8. In boreal BCRs 46–61% of the footprint has no sector above
    pressure 4. The 1 km vegetation there is largely intact, so backfilling it measures how
    roaded land differs from unroaded land, not habitat the road converted.

### Computational

4. ~~**Uncertainty combination is ad hoc**~~: **Addressed** by joint sampling in 12F. Each (bootstrap, scenario) pair draws a fresh BART posterior realization inside the BRT bootstrap loop, capturing the covariance between BART and BRT uncertainty naturally. No post-hoc variance combination.

5. **Backfill coverage collapse — fixed in code, cluster re-run still pending**: `12F` drops
   footprint pixels whose backfilled design matrix is incomplete (any `_draw_*` covariate NA →
   `complete.cases` fails). Before the fix this dropped **98.5–99.7%** of coalition pixels in some
   BCRs, leaving the counterfactual based on a tiny remnant. Two causes, both in our own pipeline
   and neither carrying V5 parity risk (neither covariate class can enter the bird BRT — see
   `12F:169` `model_vars_shared`):

   (a) **`CAfire`** (time-since-disturbance, an IA-only backfill predictor) encodes "unburned in
   the 1985–2020 record" as `NA` — true of ~99% of high-HF pixels. Partial-NA columns survive the
   all-NA drop at `08A:155`, BART returns `NaN` for any row with an NA predictor, and
   `08A:124–132` cascades that NaN down the whole biotic hierarchy. Fixed by recoding `NA → 0`
   in `02` (with `1/(ysf+1)`, never-burned is the `ysf → ∞` limit = 0, so 0 is semantically
   correct).

   (b) **Partial-NA V5 covariates** (`StandardGreenup`/`StandardDormancy` phenology on water,
   soil) were dropped twice — as BART `NaN` draws, and again at `12F`'s `complete.cases` gate —
   even though V5 predicted at 100% of those pixels, because `gbm` sends `NA_real_` down a
   missing-value branch and BART does not. Fixed by **Fix B** (`08A` median-imputes with an
   `_isNA` flag) and **Fix A** (`12F` gates only on backfilled covariates). Recoding phenology was rejected: unlike
   CAfire, these ARE V5 covariates.

   **Validation state.** Fix B is confirmed at scale — the cluster `07` smoke gives 100% complete
   `_draw_*` coverage on 1784/1784 high-HF pixels (was 0.8%). **Fix A passed** in the A7 smokes
   (100% of weight > 0 footprint pixels complete). Fix A exposed a second defect: it let pixels with a missing observed-only covariate into the
   prediction, where terra's `NaN` (not `NA_real_`) sent gbm down the wrong branch - see "12F
   restructure invariants". Fixed 2026-09-24; every 12F run now checks obs/bf parity itself
   (the identity gate).

6. **The frozen `q99.9` cap biases impact downward**: density truncation now conforms to current
   V5 packaging (`12B_v5_truncate.R` ports `10.Truncate.R`; both caps applied; gates G1–G4 passed
   — see `TODO.md`). The transform is `T(S) = clamp(clamp(S, densmax), q99) × weight`, applied
   identically to observed and every counterfactual, with **`q99` frozen from the observed
   landscape** (`predictions/{species}/truncation_params.rds`). Freezing is mandatory — a cap
   re-derived per counterfactual would adapt to the landscape being differenced and contaminate
   the contrast — but it is not neutral: counterfactual densities are higher, so they hit the
   frozen ceiling more often than observed, **under-estimating impact in the highest-density
   pixels**. The cap removes 4–12% of abundance, so this is not negligible. `TODO.md` C3 is the
   uncapped sensitivity pass.

   Two durable facts from that conformance work. `q.out`'s schema is
   `spp, thresh, countmax, densmax` — **there is no `$q` column**, and `pmin(x, NULL)` returns
   `numeric(0)` *silently*, so `12F` guards the schema and `stop()`s. And `q99` binds **1.2×–35.4×
   lower** than `densmax` across all 25 species×BCR pairs (median 6.9× CAWA / 2.1× OVEN), so it is
   the dominant cap everywhere, not a refinement.

7. **Two masking defects in `weight.tif` — both fixed 2026-09-11, weights rebuilt 25/25**. Kept
   here because both are easy to reintroduce:

   (a) **V5 prediction grids are BUFFERED well past their subunit.** The buffer is deliberate and
   load-bearing: `04.Stratify.R` buffers every BCR by 100 km "so that we can feather predictions
   from adjacent regions together", and `08.MosaicPredictions.R` consumes that overlap in a
   distance-weighted border blend (`p * w`, `mosaic(fun = "sum")`, then divide by the summed
   weights). The per-subunit crop at `10.Truncate.R:146` runs **after** mosaicking, on the
   delivery products — it is not a pre-stitch discard of the overlap. We stage the pre-mosaic
   `07_predictions` stacks, so the buffer is fully present in what `12C` sees; `12C` originally
   replicated V5's range / water / data-limit masks but not that crop. On our staged stacks
   **59–71% of non-NA pixels** lie outside the subunit polygon. Because `12D:38` assigns a subbasin to *every* BCR it intersects, and **329 of 674 subbasins
   (59% of area) straddle more than one Canadian BCR**, each straddling subbasin was summed in
   full under each of its BCRs — inflating populations **2.2×–14.6×**. `12C` now multiplies in an
   `inbcr` term from `Regions/BAM_BCR_NationalModel_Unbuffered.shp` (verified identical geometry
   to V5's `Subregions_Mosaics_EPSG3978.shp`: per-BCR areas agree to <1 km², IoU = 1.0000).
   (Wording corrected 2026-09-21: this used to say V5 crops back "before mosaicking", which
   wrongly implied the buffered overlap is thrown away rather than feathered. The defect and the
   `inbcr` fix are unaffected — only V5's order of operations was misdescribed.)

   (b) **`terra::rasterize()` defaults to `touches = FALSE`** (cell-centre rule), while V5's
   `mask(vect, inverse = TRUE)` and `crop(vect, mask = TRUE)` both behave as `touches = TRUE`.
   The gap scales with fragmentation: `WaterMask_Canada` has **29,545 polygons**, mostly thin
   rivers and small lakes that touch a cell without covering its centre, and the centre rule let
   through **4.9–7.7% of total abundance**; the single-blob BCR polygon differs by 0.7–2.2% and
   the data-limit mask by 0.02%. All three terms now pass `touches = TRUE`.

   These moved **absolute populations only** — the weight multiplies both sides, so
   `w·bf − w·obs = w·(bf − obs)` held throughout and Shapley *shares* were far less distorted
   than totals. Any per-BCR or national **total** produced before the fix is wrong. `weight.tif`
   now carries a version stamp in its band name; `12C` rebuilds a stale weight rather than
   skipping it, and `12F` refuses to run against one. Harness:
   `Rscripts/misc/verify_weight_vs_v5_masking.R`.

   **Load-bearing invariant when adding a mask term**: `weight.tif` must be non-NA everywhere.
   That is how we reproduce V5's `mask.i[is.na(mask.i)] <- 0` (`10.Truncate.R:141`) under a sum —
   range NA→0 via `classify`, and every `rasterize()` called with `background = 0`. A new term
   without `background = 0` would silently *drop* pixels instead of zeroing them. `12C` `stop()`s
   if any weight cell is NA.

8. **V5 sums density over a CONFORMAL grid, so published abundances read ~3% low (UPSTREAM)**:
   `13.Summarize.R:227` computes population as `global(t * 100, sum)` over the `10_truncated`
   product, which lives in EPSG:3978. The `× 100` step converts birds/ha to birds/km² (see
   "Instructions from Masa"), and summing that over pixels yields birds **only if every pixel is
   one square kilometre of ground** — a property of equal-area projections specifically.

   EPSG:5072 (our production CRS) is Albers **Equal Area**. EPSG:3978 (V5's delivery CRS) is
   Lambert **Conformal** Conic: between its standard parallels (49°N, 77°N) the scale factor is
   below 1, so a 1000 m × 1000 m cell holds **~1.027 km² of ground**. Each 3978 cell therefore
   carries ~2.7% more ground than it is credited with, the region needs ~2.7% fewer cells, and
   the naive sum falls short by about that much (measured naive 3978/5072 ratio: 1.0284–1.0346).
   Weight every cell by `terra::cellSize()` and the two projections agree to **0.05–0.23%** — so
   essentially the whole gap is the area-units mismatch, and bilinear interpolation contributes
   only that ~0.1% remainder.

   **Consequences.** (a) Our 5072 totals satisfy the `× 100` assumption and V5's 3978 totals do
   not, so the ~3% is V5 reading low, not us reading high. (b) The fix is one multiplication on
   either side: `global(t * 100 * cellSize(t, unit = "km"), "sum")`, which is projection-agnostic.
   (c) It is a bias, not an invalidating error — 3% sits well inside the published 5–95% bootstrap
   intervals (roughly ±25%) — and it cancels in the obs/bf contrast, so no Shapley number moves.
   (d) Gate G2 still stands precisely *because* it reproduced V5's own convention in V5's own CRS:
   it tested transform parity, not unit correctness.

   **Do not "fix" this by delivering in 3978 like V5 does.** `10.Truncate.R:19` documents the
   3978 reprojection as legacy ("Future versions will not require this step"), and `clamp()` does
   not commute with projection — reprojecting would break the exact superset→masked-rowsum
   decomposition that the whole 12F restructure rests on. Production stays in 5072; 3978 is a
   validation harness only (`Rscripts/misc/verify_projection_area_units.R`).

   **Generalizable lesson**: `× 100` silently encodes "one pixel = 1 km²", which only equal-area
   projections honour, so any per-area quantity summed over a grid is re-scaled by a mid-analysis
   projection change unless cell areas are re-derived. (Until 2026-09-14 this gap was documented
   throughout the repo as "bilinear resampling does not conserve sums" — wrong mechanism and
   wrong direction; the affected script headers and memories have been corrected.)


## Instructions from Masa
1.Always ignore the directory /Rscripts/misc when thinking. It's not immediately relevant to the project.
2.Note that population estimates made by boosted regression trees via `gbm` predict bird density at the hectare scale.
In order to match the scale of our covariate rasters (1km^2), we always multiply the pixel estimates by 100 to get the 1km^2 
density estimate. 
