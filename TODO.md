# TODO — consolidated plan (written 2026-09-11)

Single entry point for the two open workstreams. The diagnosis behind **A** is in
`CLAUDE.md` Open Limitation #5; refuted hypotheses are in the `project_pipeline_history`
memory. (This file replaced `HANDOFF_cafire_backfill_fix.md`, deleted 2026-09-11.)

| | workstream | where | state |
|---|---|---|---|
| **A** | CAfire / phenology backfill fix (Open Limitation #5) | cluster compute | fixes built + staged; **cleanup done 2026-09-11**; compute never launched |
| **B** | Conform observed + counterfactual density to current V5 packaging (Open Limitation #6) | local | **B1 done (G1 passed)**; B2 next |

**Both edit `12C`. Both must land before the single `12B` run** (2 × 24 h @ 384 G per
species — doing them in separate passes pays for it twice).

---

## Critical path

A's cluster compute (`07`: 674 array tasks) is the long pole and depends on nothing in B.
B is entirely local. **Run them in parallel and converge at `12B`.**

```
cluster:   [cleanup DONE] ──► sbatch 07 ──► sbatch 11 ────────────────────┐
                                                                          ├──► sbatch 12B ──► validate ──► 14B
local:     [B1 DONE] ──► B2 rewrite 12A ──► B3 re-run 12A ──► Globus ──────┘
```

Blocker for every cluster step: SSH host key verification fails, so cluster commands
must be run from your own authenticated session. Globus works (one file per call,
**never** `--batch`).

---

## 0. Commit what's already written

- [ ] Commit Fix A (`12C` — gate `complete.cases` on `backfilled_vars` only) and the
      `12B` weight.tif preflight. Both are uncommitted in the working tree.
      `02`/`06`/`08A` are already in (`0bbf2ea`, `4d23436`).

---

## A. CAfire backfill re-run (cluster)

- [x] **A1. DONE 2026-09-11.** Cleanup Tier A — reclaim ~14.4 G, safe any time:
      six pre-fix `covariates_mosaiced_{1990..2015}.tif` (2026-02-06, pre-CAfire-fix and
      therefore defective; only consumer `17` skips missing years) + `predictions_coalitions/`
      (5.3 G, output of the retired per-coalition path).
- [ ] **A2.** Confirm `weight.tif` exists for all **25** species×BCR pairs
      (11 CAWA + 14 OVEN `can*` models). If short, `sbatch 12A2_build_prediction_weights.sh`.
      **Do not skip.** `12C` falls back to an UNMASKED run (`w == 1`) with only a warning when
      `weight.tif` is missing, and that got materially more dangerous once Fix B started
      median-imputing: water pixels used to fail `complete.cases` and drop out of BOTH obs and
      bf on their own, so `weight.tif` is now the **only** thing keeping water, out-of-range
      and out-of-extent pixels out of the density tables. The preflight fails fast instead of
      20 h into a job whose totals are quietly wrong.
- [ ] **A3.** Re-Globus `12B`/`12C` if either changed since the 2026-07-02 staging (they have —
      see step 0 and B4).
- [x] **A4. DONE 2026-09-11** (run early, not just-before-`07` — strictly safer, nothing left to mix). Cleanup Tier B — ~27 G, run **immediately before** `sbatch 07`:
      `bart_models/2020` + `bart_models_mosaics/2020`, and
      `rm -f density_tables/*.rds density_tables/arrays/*.rds`.
      Timing matters: a partial `07` would otherwise interleave pre- and post-fix subbasins
      with nothing but mtime to tell them apart. The `_metrics.rds`/`_confusion.rds` baseline
      is already safe locally at 674/674.
- [ ] **A5.** Smoke `07` on a few subbasins — Fix B is untested at scale.
- [ ] **A6.** `sbatch --array=1-674 07_train_and_backfill.sh` → `sbatch 11_premosaic_backfilled_stacks.sh`.
- [ ] **A7.** Validate Fix A: the `12C` runtime line `complete superset pixels: N / M (X%)`
      must be ≫ the old 1–3%. This is the only fix never tested.

## B. V5 packaging conformance (local)

V5 split `10.Package.R` into `10.Truncate.R` + `11.Package.R` (commit `f082866`, 2026-06-04)
and rewrote the truncation values (`4e7fc83`). Measured gap on `can10` 2020:

| stage | CAWA | OVEN |
|---|---|---|
| `densmax` cap removes | 1.29 % | 0.32 % |
| frozen `q99.9` cap removes a further | **12.19 %** | **4.13 %** |
| projection 5072→3978 costs a further | ~2.3 % | — |

`q99.9` (≈0.133 for CAWA) binds 6.3× lower than `densmax` (0.833) and is doing nearly all
the work — and our density tables have **never** applied it.

**Design rule:** one transform `T`, applied identically to observed and every counterfactual,
with all data-dependent parameters **frozen from the observed landscape**.
`T(S) = clamp(clamp(S, densmax), q99) × range × notwater × inlimit`.
Recomputing `q99.9` per counterfactual would let the cap adapt to the landscape being
differenced and contaminate the contrast.

**Stay in EPSG:5072.** `10.Truncate.R:19` states the 3978 reprojection is legacy
("Future versions will not require this step"). It costs 2.3 % of abundance (bilinear
isn't conservative) and would break the exact superset→masked-rowsum decomposition,
since `clamp` doesn't commute with projection. The *parameter* is CRS-invariant
(q99.9 = 0.132897 in 5072 vs 0.132940 in 3978, 0.03 % apart), so the 5072-derived cap
**is** V5's cap. Use 3978 only as a validation harness.

- [x] **B1. DONE — Gate G1 PASSED 2026-09-11.** `Rscripts/12A0_v5_truncate.R` ports
      `10.Truncate.R`; `Rscripts/misc/verify_v5_truncate_port.R` is the harness.
      Reproduces `10_truncated` on CAWA can10/can71 and OVEN can10: dims, extent and CRS
      identical, sum ratios within 4.3e-5, and **flat terrain bit-identical** (mean |diff| in the
      smoothest decile 1.0e-09 / 5.6e-07 / 7.6e-08). Residual differences are confined to steep
      terrain and are reprojection noise from a PROJ/GDAL version difference — see the
      `project_terra_gotchas` memory. The port generalises the transform in two ways the
      counterfactual needs: `q99` can be passed in (frozen) rather than derived, and
      `project_to = NULL` skips the legacy 3978 step.
- [ ] **B2.** Rewrite `12A_observed.R`:
      - `q.out$q` → `q.out$densmax` (`:55`). Schema changed; `$q` no longer exists.
      - drop the `q0`/`denshthresh` step (`:56`, `:117-119`) — V5 deleted it.
      - apply densmax **and** q99.9 to the saved stack, so `observed_bootstraps.tif`
        becomes our 5072 analogue of `10_truncated`.
      - write `truncation_params.rds` per species×BCR: `{densmax, q99, spp, bcr, year}`.
- [ ] **B3.** Re-run `12A` locally for CAWA + OVEN (25 stacks). **Gate G2** — check observed
      BCR/national totals against `output/13_summary/BAMV5-abundance.RData` /
      `BAMV5-results.xlsx`. `13.Summarize.R:227` is `global(t * 100, sum)` per bootstrap
      layer, median + 5–95 % — directly comparable to `obs_total_mean/sd`. First exact
      external check we've ever had on the observed side.
      **Gate G3** — record the 5072-vs-3978 total ratio per species×BCR as the published
      crosswalk to BAM's numbers.
- [ ] **B4.** Update `12C_predict_species_all_coalitions.R`:
      - `$q` → `$densmax` (`:34`). **Currently `pmin(pred_vec, NULL)` returns `numeric(0)`
        silently** if the new Rdata is dropped in without this.
      - delete `q0` (`:35`).
      - read `q99` from `truncation_params.rds`; insert the second clamp at `:282`
        *before* the weight multiply (V5 order: truncate → mask; `:293` already has this right):
        `sc[, k] <- pmin(pmin(pred_vec, densmax), q99)`
- [ ] **B5.** Restage `data/raw_data/SpeciesPredictionTruncationValues.Rdata` from
      `G:/Shared drives/BAM_NationalModels5/data/` (2026-06-04 version).
- [ ] **B6.** `12A2` — no change needed. `weight = range × notwater × inlimit` already matches
      `10.Truncate.R:137-148`, and the source GIS layers are unchanged since staging
      (WaterMask_Canada Apr 23, DataLimitationsMask Jan 16, ranges Apr 22 — all predate
      our 2026-06-05 copy, byte sizes identical). Check only that our NA→0 handling
      matches `10.Truncate.R:141`.
- [ ] **B7.** Globus the 25 regenerated `observed_bootstraps.tif` + `truncation_params.rds`
      to the cluster. One `globus transfer` call per file.
      **Note:** this overrides the old "`predictions/` — never regenerate" rule. Keep the
      `weight.tif` files that live in the same directories.
- [ ] **B8.** `14B_sector_attribution.R` — exclude or flag **CAWA `can40`**.
      `review/ModelReleaseDecisions.xlsx` "remove" tab row 47 withholds it (AUC).
      Our BCR discovery reads `06_bootstraps/{spp}/can*.Rdata`, which yields exactly the
      release set **plus can40** for CAWA; OVEN matches all 14. One anti-join closes it.

## C. Converge

- [ ] **C1.** `sbatch 12B_repredict_all_coalitions.sh` (`--array=1-2`). The usual pre-run wipe of
      `density_tables/*.rds` **and** `arrays/*.rds` already happened in A4; redo it only if
      anything writes there first. The old tables were stale on three counts: pre-weighting,
      pre-gate-change (A), and pre-truncation-conformance (B).
- [ ] **C2.** `14B_sector_attribution.R` locally → corrected Shapley CSVs.
- [ ] **C3.** Sensitivity pass with the `q99.9` stage disabled, and report the spread.
      The frozen cap has a known-direction bias: counterfactual densities are higher, so
      they hit the frozen ceiling more often than observed, which systematically
      **under-estimates impact in the highest-density pixels**. With the cap removing
      4–12 % of abundance, that is not negligible. Headline = conformed; sensitivity = uncapped.
- [ ] **C4.** Update `CLAUDE.md` (Open Limitation #5; the pre-`12B` wipe instruction should
      also clear `density_tables/arrays/`; V5 script renumbering) and the memory files.

---

## Loose ends

- [ ] Decide the fate of local `covariates_mosaiced_2020_PREORIG.tif` (1.5 G) — the only
      surviving pre-CAfire-fix 2020 mosaic. The G: sandbox copy was overwritten by the rebuild.
- [ ] `12F_singletons_plot.R:150` reads `predictions_coalitions/` (local, coalition 9 / CAWA /
      mines). That directory was deleted in the local cleanup, so the map panel is already
      broken; A1 removes the cluster copy too. Regenerate via `save_arrays_ids` if wanted.

## V5 reference (for `CLAUDE.md`)

Renumbering in `f082866`: `10.Package.R` → `10.Truncate.R` + `11.Package.R`;
`11.Validate` → `12.Validate`; `12.Summarize` → `13.Summarize`;
`output/10_packaged/` → `output/11_packaged/`, new `output/10_truncated/`.
**`06_bootstraps` and `07_predictions` were NOT renumbered** — they are the only V5 output
folders our scripts read, so no path in our pipeline is broken. The CAWA/OVEN 2020
prediction tifs date May 10–11, before our 2026-05-15 `12A` run; inputs are unchanged.
