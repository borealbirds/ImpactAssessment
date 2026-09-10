# HANDOFF — Fix CAfire NA-by-design that collapses backfill coverage (Open Limitation #5)

**Author:** Claude (handoff written 2026-06-13; status block updated 2026-09-10)
**Status:** root cause CONFIRMED; fix BUILT, locally validated, and STAGED on the cluster.
**Cluster compute has never been launched — that is the only thing left.** See §0 below.
**Scope:** backfilling instruments only (`02` → `06` → `07`/`08` → `11` → `12`). Does NOT touch the
V5 bird model or methodology parity (see "Why this is safe" below).

---

## 0. STATE AS OF 2026-09-10 — PICK UP HERE

Sections 1–4 below are the original diagnosis and remain accurate as history. Sections 5–8 are
now largely DONE — read this block first, it supersedes their "TODO" framing.

### What is finished

| Item | State |
|---|---|
| CAfire `NA → 0` recode in `02` + `CAfire/` path reconcile + `06` year=2020 | **committed & pushed** `0bbf2ea` (2026-09-10) |
| Fix B — `08A` median-impute + `_isNA` flag for partial-NA BART predictors | **committed** (in `4d23436`) |
| Fix A — `12C:247` gate `complete.cases` on `backfilled_vars` only | **written, NOT committed** |
| `12B` weight.tif preflight (hard-stop if any species×BCR weight missing) | **written, NOT committed** |
| Local de-risk of CAfire recode | **PASSED** — mosaic CAfire NA 91.7% → 0.1%; `frac_backfillable` ~0.01 → 0.68–0.87 |
| Local validation of Fix B | **PASSED** — subbasin 1: all continuous-biotic `_draw_001` 100% finite, 0 NaN; `frac_backfillable` 0.8% → ~80% (post-CAfire) → **100%** (post-Fix-B) |
| `06` 2020 mosaic rebuild (local) + Globus to cluster | **DONE 2026-07-02**, 1636027925 bytes, byte-exact verified |
| Edited `02`/`08A`/`12C` Globus'd to cluster `Rscripts/` | **DONE 2026-07-02** |
| Cluster compute (`07` → `11` → `12B`) | **NOT STARTED** |

### The phenology question (§8) is CLOSED — do not re-litigate it

`StandardGreenup_1km`/`StandardDormancy_1km` were the residual NA driver after CAfire. Resolved
2026-06-13:
1. They **are** V5 bird covariates (in `model_vars_shared`), so recoding them would corrupt the V5
   bird model. **Recoding is OFF the table** (unlike CAfire, which is IA-only and had a principled 0).
2. Their NA is **genuine V5 source data** (raw `gis/stacks/{bcr}_2020.tif`: 8.6% can11, 18.8% can60,
   30.8% can80 within footprint), concentrated on water (VLCE 20) — not a `06` artifact.
3. **V5 predicted at 100% of phenology-NA pixels.** `gbm::predict.gbm` tolerates NA via surrogate
   splits; BART does not. So the residual gap was a **double over-restriction in our own pipeline**,
   not a real limitation — hence Fix A + Fix B, which are complementary and BOTH required (such a
   pixel is dropped twice: once as a BART NaN draw, once at the `12C` gate).

### Why the `12B` preflight matters now (do not skip it)

`12C` silently falls back to an UNMASKED run (`w == 1`) with only a warning when `weight.tif` is
missing. That got materially more dangerous once Fix B started median-imputing: water pixels used to
fail `complete.cases` in `12C` and drop out of BOTH obs and bf on their own, so **`weight.tif` is now
the only thing keeping water / out-of-range / out-of-extent pixels out of the density tables.** The
preflight fails fast instead of 20 h into a job whose totals are quietly wrong.

### Next session, in order

1. Commit Fix A (`12C`) + the `12B` preflight. (`02`/`06`/`08A` are already in.)
2. Re-Globus `12C` and `12B` to the cluster if either changed after the 2026-07-02 staging.
3. Confirm `weight.tif` exists for every species×BCR — the preflight now enforces this, but `12A2`
   may never have been run. If absent: `sbatch 12A2_build_prediction_weights.sh`.
4. Smoke `07` on a few subbasins (Fix B is untested at scale; subbasin 3 OOM'd on the *local*
   workstation only — expected fine at cluster 384G).
5. Full run: `sbatch --array=1-674 07_train_and_backfill.sh` → `sbatch 11_premosaic_backfilled_stacks.sh`
   → `rm -f ../data/derived_data/density_tables/*.rds` → `sbatch 12B_repredict_all_coalitions.sh`
   (`--array=1-2`). Do **NOT** rerun `12A`/`12A2` for CAfire/phenology reasons — neither touches
   observed density or the weight.
6. Validate Fix A on the cluster: watch the `12C` runtime line `complete superset pixels: N / M (X%)`
   — it must be ≫ the old 1–3%. This is the ONLY untested fix.
7. `14B_sector_attribution.R` locally for corrected Shapley CSVs. Then update `CLAUDE.md` Open
   Limitation #5 and the memory files.

### Known blocker

Cluster automation SSH fails with **"Host key verification failed"** — compute needs the user's own
authenticated session. Globus transfers work fine (one file per call, never `--batch`).

### Caveat on old outputs

Bit-identity with pre-fix runs is **intentionally moot**: Fix A changes the gate, so post-fix numbers
legitimately differ. Existing `density_tables/*.rds` are stale on two counts (pre-weighting AND
pre-gate-change) — delete before the `12` run, per step 5.

---

## 1. The problem, in one paragraph

In many BCRs, `12C_predict_species_all_coalitions.R:238`'s `complete.cases` gate keeps only ~1–3%
of high-HF coalition pixels (CAWA `can10` ≈ 3.5%, `can61` ≈ 0.4%, OVEN `can80` ≈ 0%), so the
counterfactual — and the per-BCR Shapley sector numbers — are unreliable there. The dropout is NOT
a V5 problem and NOT the V5 climate sparsity (FFP etc., which is a harmless `06` union artifact,
all-NA columns dropped per subbasin at zero pixel cost — see `memory/project_backfill_coverage_audit.md`
V5-exoneration block). **The real culprit is `CAfire` (time-since-disturbance), an ImpactAssessment-
added abiotic backfilling predictor, whose `NA` means "unburned in the 1985–2020 fire record" — true
of ~99% of high-HF/industrial pixels.** It is being treated as *missing data* when it is actually a
*meaningful "no recent fire" state*.

## 2. Root cause, verified in code (read these to confirm before changing anything)

- **`02_reproject_disturbance_covariate.R:37–52`** — CAfire = normalized years-since-fire `1/(ysf+1)`
  (recent fire → 1, old fire → 0). `NA` is assigned for: no fire in the 1985–2020 record (`:39`),
  and `ysf < 0` future-fire pixels (`:47`). All of these mean **"unburned as of this year."**
- **`08A_train_and_backfill_subbasin_s.R`** cascade:
  - `:155` drops only *all-NA* backfill columns. CAfire is NA at ~99% of `backfill_idx` pixels but
    NOT all-NA (the ~1% that burned), so it **survives** as a live predictor.
  - `08B_deploy_gbart.R` / `08B_deploy_mbart.R`: `BART::gbart`/`mbart` return `NaN` for any test row
    with an NA predictor → that `NaN` is written into `df_backfill[[b]]` for ~99% of pixels.
  - `08A:124–132`: each later continuous biotic covariate uses already-backfilled earlier ones as
    `b_before` predictors → **one NaN cascades to NaN for every downstream covariate**. By the last
    hierarchy level ~99% of pixels are NaN.
- **`12C:182–191, 237–238`** — those NaNs become NA `_draw_*` layers; `complete.cases(X_rep[,
  model_vars_shared])` then drops the pixel. That is the 98–99% dropout.
- **Quantitative confirmation** (already on disk): `logs/rootcause_abiotic_na_by_cov.csv` shows CAfire
  ("Time Since Disturbance") `na_frac ≈ 0.987`, `all_na = FALSE` = THE seed; soil ≈ 0.018 and the
  phenology bands (StandardGreenup ≈ 0.26, StandardDormancy ≈ 0.25) are secondary. `logs/
  rootcause_cascade_by_subbasin.csv` shows `frac_any_abiotic_na ≈ dominant_na_frac ≈ 0.99` with
  `dominant_na_class = "Time Since Disturbance"` for nearly every subbasin, and `frac_backfillable_last`
  collapsing to ≈ 0.035 (can10).

## 3. Why this is safe (methodology parity is NOT at risk)

Verified in `12C`: the bird BRT is predicted with `gbm::predict.gbm(model, X_k[, model_vars_shared])`
(`12C:269–271`), `model_vars_shared <- b.list[[1]]$var.names` (`12C:123`) = the exact V5 covariate
list. **CAfire and soil are never in `model_vars_shared`** — they are not V5 covariates and cannot
enter the bird model. They exist only as backfilling predictors in `08`. The counterfactual differs
from observed ONLY by (a) biotic covariates swapped to backfilled values and (b) disturbance vars
zeroed (`12C:254–266`); the same weight + `qsp` clamp apply to both. So fixing CAfire changes only the
**coverage/quality of the backfilled biotic layers** — zero risk to V5 parity. (Full parity write-up
is in the conversation that produced this handoff.)

---

## 4. THE FIX (recommended path)

**Recode CAfire `NA → 0` at the source in `02`.** Rationale: with `1/(ysf+1)`, "never burned" is the
limit `ysf → ∞ → 0`, so **0 is the semantically correct value, not a hack**. This fills the cascade
seed, recovers ~99% of pixels, and *keeps* the genuine burned/unburned signal. Do it at the raster
source so both training (low-HF) and backfill (high-HF) pixels see a fully-populated band, and do it
BEFORE `06`'s bilinear reproject (so 0s blend correctly instead of NA-contaminating edges).

### 4a. Edit `02_reproject_disturbance_covariate.R`, function `find_year_since_fire` (`:43–52`)

After `ysf_norm` is computed and named, add the recode and re-apply the BAM boundary mask so
off-study-area stays NA (within-subbasin pixels are all in-study, so they get a real 0):

```r
    ysf_norm <- 1 / (ysf + 1)                 # existing
    ysf_norm[is.na(ysf_norm)] <- 0            # NEW: unburned-in-record = 0 (real value, not missing)
    ysf_norm <- terra::mask(ysf_norm, bam_boundary)  # NEW: keep off-study-area as NA (harmless; masked downstream anyway)
    names(ysf_norm) <- paste0("CAfire_", current_year)  # existing
```

(`bam_boundary` is already in scope from `02:29`.) Leave `02:38–40` as-is — those still produce the
pre-recode NAs that the new line converts.

### 4b. CRITICAL path-reconciliation gotcha — DO NOT SKIP

`02` currently writes `CAfire_<year>_masked.tif` to the **sandbox root**
(`02:61–64`, `ia_dir = G:/Shared drives/BAM_NationalModels5/data/Extras/sandbox_data/impactassessment_sandbox`),
but `06` reads CAfire from a **`CAfire/` subdirectory**:
`06:160` `caf_path <- file.path(ia_dir, "CAfire", paste0("CAfire_", y, "_masked.tif"))`.
Before rebuilding, **confirm where `06`'s `ia_dir` points (read `06`'s header) and place the recoded
`CAfire_<year>_masked.tif` files into the exact `CAfire/` dir `06` reads.** Otherwise `06` will hit the
`:177–179` "CAfire not found … (skipping)" branch and silently rebuild the mosaic WITHOUT CAfire.

### 4c. Confirm the production `year`

There is a known year ambiguity: `06:22` currently `years <- seq(1990, 2015, by=5)` (2020 line
commented at `:23`); memory notes 07 has used `year <- 2015` while the 12-series and observed
bootstraps are **2020** (`density_tables/{species}_{year}_...`, `predictions/.../2020/`). **Determine
the year the current `12B` run consumes** (read `07_train_and_backfill.R` `year <-` and
`12B_repredict_all_coalitions.R`), and rebuild the mosaic for THAT year (almost certainly 2020). Set
`06:22` accordingly before rerunning.

---

## 5. De-risk LOCALLY before any cluster run (do this first)

The local machine has `data/raw_data/covariates_mosaiced/covariates_mosaiced_2020.tif` and G: source
data. Prove the fix recovers pixels before committing to a multi-day cluster re-run:

1. Apply the 4a edit.
2. Re-run `02` locally to regenerate the recoded `CAfire_<year>_masked.tif` (note `02:25` aggregate is
   fast; the expensive `project` at `:32–35` is ~hours — if a prior aggregated/reprojected CAfire exists
   you may be able to recode it directly with terra rather than re-projecting from 30 m source).
   Place outputs per 4b.
3. Re-run `06` locally for the production year only to rebuild `covariates_mosaiced_<year>.tif`
   (set `06:22` to a single year). Confirm the run printed "successfully added soil stack" AND did NOT
   print "CAfire not found".
4. **Coverage check:** run a lightweight version of the existing cascade audit
   (`Rscripts/misc/diag_backfill_rootcause.R`, or `diag_backfill_coverage_audit.R`) over 2–3 subbasins
   (one known-bad e.g. a `can10` subbasin, one healthy `can14`) against the rebuilt mosaic with
   `cc=FALSE, local=TRUE`. **Expected result: `frac_any_abiotic_na` collapses from ≈ 0.99 to ≈ 0**
   (CAfire no longer NA), and `frac_backfillable_last` jumps toward ~1.0. If it does not, STOP and
   re-diagnose before spending cluster time.

## 6. Full re-run sequence (cluster) — only after §5 passes

The recode propagates: `02` (regen CAfire) → `06` (rebuild mosaic) → `07`/`08` (re-backfill, sources
08) → `11_premosaic` (re-mosaic backfill) → `12B`/`12C` (re-predict). **`12A_observed` and `12A2`
weights do NOT change** (CAfire never touches observed density or the range/water/extent weight) — do
not rerun them.

Because `06` reads G: source BCR stacks + CAfire + soil, rebuild the mosaic LOCALLY (§5 step 3), then
Globus-transfer the rebuilt mosaic to the cluster scratch path
(`/home/mannfred/scratch/impact_assessment/data/raw_data/covariates_mosaiced/`). Then on the cluster:

```bash
sbatch 07_train_and_backfill.sh   # --array=1-674 full re-backfill (sources 08); multi-hour per task
```
```bash
sbatch 11_premosaic_backfilled_stacks.sh   # re-mosaic per-subbasin backfill into BCR stacks
```
```bash
rm -f ../data/derived_data/density_tables/*.rds   # clear stale tables before a fresh 12 run
```
```bash
sbatch 12B_repredict_all_coalitions.sh   # --array=1-2 (1=CAWA, 2=OVEN); writes {species}_{year}_coalition_{cid}.rds
```

Smoke test one species/BCR/2 bootstraps first:
```bash
sbatch --array=1 --time=01:00:00 --mem=192G --export=ALL,TEST_BCR=can60,TEST_N_BOOT=2 12B_repredict_all_coalitions.sh
```

After `12B`, re-run `14B_sector_attribution.R` (local) for the corrected Shapley CSVs.

## 7. Verify the fix worked

- During `12C`, the runtime message `complete superset pixels: N / M (X%)` (`12C:240–242`) should now
  report HIGH percentages (≫ the previous 1–3%) in the previously-degenerate BCRs (can10, can61,
  can70, can80). The `backfill mosaic likely degenerate or uncovered` warning should largely disappear.
- Optionally re-run the cluster cascade audit (`diag_backfill_rootcause.sh`) and confirm BCR-rollup
  `frac_backfillable_last` is near 1.0 where it was ~0.03–0.06.

## 8. Secondary follow-ups (after the CAfire fix lands)

- **Phenology bands** `StandardGreenup` (~0.26 NA) and `StandardDormancy` (~0.25 NA) are the next
  contributors. Check whether their NA is also a meaningful state (e.g. water / non-vegetated) and
  whether to similarly recode or exclude them. Soil (~0.018) is minor; likely ignorable.
- Update `CLAUDE.md` Open Limitation #5 and `memory/project_backfill_coverage_audit.md` /
  `project_backfill_coverage_degeneracy.md` once coverage is confirmed restored.

## 9. House rules (carry forward)

- **Globus:** one `transfer` call per file, **never `--batch`** (always "permission denied"). CLI at
  `C:\Users\mannf\AppData\Local\Python\pythoncore-3.14-64\Scripts\globus.exe`. Endpoints — Local
  `a7878ccc-747b-11ef-b4b8-8fef73a45f39` (local paths start `/C/`), Cluster Fir
  `8dec4129-9ab4-451d-a45f-5b4b8471f7a3`. Use `--sync-level checksum`.
- **`nm_root` = `/home/mannfred/projects/def-ecknight/NationalModels`** — cluster-only, never change.
- Bash for the user must be **single-line** (no `\` continuations, no multiline blocks).
- Bird density: gbm predicts birds/ha; ×100 for birds/km² (1 km² = 100 ha).
- Ignore `Rscripts/misc` for conceptual reasoning (it's diagnostics/one-offs).
- SLURM account `--account=def-ecknight`; keep threads within one CCD (`--cpus-per-task=8`); use
  `$SCRATCH` for job I/O, copy keepers to `$HOME/project`.
