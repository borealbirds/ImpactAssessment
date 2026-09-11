# TODO — consolidated plan (written 2026-09-11)

Single entry point for the two open workstreams. The diagnosis behind **A** is in
`CLAUDE.md` Open Limitation #5; refuted hypotheses are in the `project_pipeline_history`
memory. (This file replaced `HANDOFF_cafire_backfill_fix.md`, deleted 2026-09-11.)

| | workstream | where | state |
|---|---|---|---|
| **A** | CAfire / phenology backfill fix (Open Limitation #5) | cluster compute | fixes built; **all scripts staged + cleanup done 2026-09-11 (A3, A4)**; compute never launched |
| **B** | Conform observed + counterfactual density to current V5 packaging (Open Limitation #6) | local | **DONE - B1-B8, C4 and gates G1-G4 all passed.** G3/G4 found + fixed two `weight.tif` defects (Open Limitation #7); only residual is the deliberate +2.9-3.4% 5072-vs-3978 choice |

**Both edit `12C`. Both must land before the single `12B` run** (2 × 24 h @ 384 G per
species — doing them in separate passes pays for it twice).

---

## Critical path (revised 2026-09-11, after B6/B8/C4)

All staging is done except `weight.tif` (A2). `07` (674 array tasks) is the long pole and
depends on nothing in B, so it should start first and everything else runs in its shadow.

```
cluster:  [cleanup DONE] ─► A5 smoke 07 ─► A6 sbatch 07 ─► sbatch 11 ─┐
          [A2 weight.tif, AFTER re-Globus of 12A2] ─────────┐         ├─► 12B smoke ─► C1 12B ─► C2 14B
local:    [B1-B8, C4, G1-G4 ALL DONE] ─────────────┴─────────┘
```

**A2 changed and now has a prerequisite.** `12A2_build_prediction_weights.R` and
`12C_predict_species_all_coalitions.R` were both edited 2026-09-11 (BCR cut + version
stamp + `touches = TRUE` + NA audit) and must be re-Globus'd before A2 runs, or A2
rebuilds the same defective weights. No `rm` of the old `weight.tif` is needed — the version stamp
(`weight_v3_touches`) forces a rebuild, and 12C refuses to run against anything older.

Verified 2026-09-11 and NOT a blocker any more:
- `covariates_mosaiced_2020.tif` is the CAfire-fixed build on **both** ends (local partial-NA
  0.06 %, and a checksum-sync Globus transfer moved 0 bytes). No `06` re-run needed.
- All seven changed cluster scripts are staged (A3) — `08A` carrying Fix B among them, which
  was stale and would have made the `07` re-run reproduce the bug it exists to fix.
- All 25 observed stacks + both `truncation_params.rds` are staged (B7), byte-exact.

Blocker for every cluster step: SSH is keyboard-interactive (2FA), so cluster commands must be
run from your own authenticated session. Globus works (one file per call, **never** `--batch`).

---

## 0. Commit what's already written

- [x] **DONE 2026-09-11.** Fix A + the `12B` weight preflight committed (`6455756`); the V5
      port, truncation conformance and docs followed in `d4b803e`, `7057c91`, `cb573cb`.

---


## A. CAfire backfill re-run (cluster)

- [x] **A1. DONE 2026-09-11.** Cleanup Tier A — reclaim ~14.4 G, safe any time:
      six pre-fix `covariates_mosaiced_{1990..2015}.tif` (2026-02-06, pre-CAfire-fix and
      therefore defective; only consumer `17` skips missing years) + `predictions_coalitions/`
      (5.3 G, output of the retired per-coalition path).
- [ ] **A2.** Confirm `weight.tif` exists for all **25** species×BCR pairs
      (11 CAWA + 14 OVEN `can*` models). If short, `sbatch 12A2_build_prediction_weights.sh`.
      **Do not skip.** `12C` now `stop()`s outright if `weight.tif` is missing (changed
      2026-09-11; it used to fall back to `w == 1` with only a message). Fix B made that
      fallback dangerous: 08A median-imputes partial-NA covariates, so water, out-of-range
      and out-of-extent pixels no longer fail `complete.cases()` and drop out of BOTH sides
      on their own. `weight.tif` is now the only masking left, and an unmasked run would be
      quietly wrong rather than obviously broken — a 7.8x over-count on CAWA can10. 12C also
      rejects an all-zero/NA weight, which would zero every density in the BCR.
- [x] **A3. DONE 2026-09-11.** Re-Globus'd every cluster-side script changed since the
      2026-07-02 staging — **seven files, not the two this item originally named.** All seven
      transferred non-zero bytes, i.e. all seven were stale on the cluster:

      | file | why it had to go | in `07` path? |
      |---|---|---|
      | `08A_train_and_backfill_subbasin_s.R` | **Fix B** (median-impute + `_isNA`) | **yes** — sourced by `07` |
      | `09_collect_metrics_mbart.R` | sourced by `08B_deploy_mbart` | **yes** |
      | `11_premosaic_backfilled_stacks.R` | `cc = TRUE`, changed | after `07` |
      | `12B_repredict_all_coalitions.R` | weight preflight | no |
      | `12C_predict_species_all_coalitions.R` | Fix A + B4 caps + weight hard-stop | no |
      | `10C_abiotic_extrapolation_diagnostics.R` | `cc = TRUE`, changed | no |
      | `10D_BART_posterior_diagnostics.R` | hardcoded cluster path, changed | no |

      **Had A3 been executed as scoped (`12B`/`12C` only), the 674-task `07` run would have
      used the pre-Fix-B `08A` and reproduced the coverage collapse we are re-running to fix.**
      Lesson for future stagings: derive the file set from
      `git log --since=<staging date> --name-only -- Rscripts/` intersected with the sourced-from
      entry points, rather than from memory of which scripts "are cluster scripts".
      `12A0_v5_truncate.R` was deliberately NOT staged — `12A` runs locally and nothing on the
      cluster sources it.

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
(CAWA can10 q99.9 = 0.132457 in 5072 vs 0.132881 in 3978, 0.32 % apart), so the 5072-derived cap
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
- [x] **B2. DONE 2026-09-11.** Rewrote `12A_observed.R`:
      - `q.out$q` → `q.out$densmax` (`:55`). Schema changed; `$q` no longer exists.
      - drop the `q0`/`denshthresh` step (`:56`, `:117-119`) — V5 deleted it.
      - apply densmax **and** q99.9 to the saved stack, so `observed_bootstraps.tif`
        becomes our 5072 analogue of `10_truncated`.
      - write `truncation_params.rds` per species×BCR: `{densmax, q99, spp, bcr, year}`.
- [x] **B3. DONE 2026-09-11.** Re-ran `12A` locally for CAWA + OVEN. All 25 stacks rebuilt
      (11 CAWA incl. the withheld `can40`, 14 OVEN); both exit 0; `truncation_params.rds`
      holds 11 and 14 entries; 32 layers and `INTERLEAVE=BAND` verified on spot-checks.
      Stack maxima sit +4e-8 *relative* above the recorded `q99` — that is FLT4S round-to-
      nearest on write (float32 eps = 1.2e-7), not a clamp failure.

      **The `densmax/q99` ratio at scale confirms the finding was not local to can10.**
      `densmax` is essentially never the binding cap:

      | | range | median |
      |---|---|---|
      | CAWA (11 BCRs) | 3.4x - 24.8x | 6.9x |
      | OVEN (14 BCRs) | 1.2x - 35.4x | 2.1x |

      Extremes: CAWA `can71` q99=0.0336 vs densmax 0.833 (24.8x); OVEN `can82` q99=0.0393 vs
      1.394 (35.4x). Every one of the 25 pairs was previously truncated at `densmax` only.

- [x] **Gate G2 PASSED 2026-09-11.** Our population estimates reproduce BAM's published
      numbers to full published precision, point estimate *and* both bootstrap bounds:

      | | ours (M) | BAM (M) | ours 5-95% | BAM 5-95% |
      |---|---|---|---|---|
      | CAWA can10 | 0.0194 | 0.019 | 0.0146-0.0239 | 0.015-0.024 |
      | OVEN can10 | 0.0541 | 0.054 | 0.0486-0.0588 | 0.049-0.059 |
      | CAWA can71 | 0.1619 | 0.162 | 0.1472-0.2085 | 0.147-0.208 |

      Confirms the transform, the frozen parameters, the x100 hectare->km2 convention, and the
      bootstrap-interval construction, against `output/13_summary/BAMV5-abundance.RData`.
      Caveat: this is the *harness* config (3978 + V5 vector masks), not production — see G3/G4.
- [~] **G3/G4 — masking: ONE BIG DEFECT FOUND AND FIXED 2026-09-11; residual still open.**
      Ran the check locally (rebuilding weights on this machine rather than waiting for A2)
      and it immediately failed — for a reason that was not on the list.

      **`weight.tif` was missing V5's mosaic cut.** V5 predicts each subunit on a grid
      buffered well past the subunit and crops it back to the subunit polygon at
      `10.Truncate.R:146` before mosaicking. `12A2` had V5's range / water / data-limit
      masks but not that crop. On our own staged stacks **59-71% of non-NA pixels** lie
      outside the subunit's own polygon, carrying **41-84% of the raw density sum**; and
      because `12B:38` assigns a subbasin to every BCR it intersects (**329 of 674, 59% of
      area**), every straddling subbasin was summed in full under each of its BCRs.
      Measured inflation vs the corrected weight:

      | pair | old weight | BCR-cut weight | inflation |
      |---|---|---|---|
      | CAWA can10 | 125,000 | 21,195 | **5.90x** |
      | CAWA can71 | 369,570 | 168,887 | **2.19x** |
      | OVEN can10 | 918,402 | 62,965 | **14.59x** |

      Fix: `12A2` multiplies in an `inbcr` term from `Regions/BAM_BCR_NationalModel_
      Unbuffered.shp`, verified to be the same geometry as V5's
      `Subregions_Mosaics_EPSG3978.shp` (areas agree to <1 km2, IoU = 1.0000 on
      can10/can11/can60) — so no new Globus staging. `weight.tif` now carries a version
      stamp in its band name (`weight_v2_bcrcut`); `12A2` REBUILDS a stale weight instead
      of skipping it, and `12C` refuses to run against one. See Open Limitation #7.

      **Residual, decomposed 2026-09-11.** Ran three configurations on the same
      `observed_bootstraps.tif` to separate the two candidate causes:

      | pair | A 5072+weight.tif | B 5072+V5 vectors | C 3978+V5 vectors | published |
      |---|---|---|---|---|
      | CAWA can10 | 21,195 | 19,956 | **19,389** | 0.019 M |
      | CAWA can71 | 168,887 | 166,549 | **161,922** | 0.162 M |
      | OVEN can10 | 62,965 | 55,891 | **54,062** | 0.054 M |

      **C reproduces the published numbers on all three**, and matches the G2 harness
      figures (0.0194 / 0.1619 / 0.0541 M) — so everything from `observed_bootstraps.tif`
      onward is sound, and the whole residual lives in the two deliberate production
      departures.

      - **B->C, the 5072-vs-3978 crosswalk: +2.9%, +2.9%, +3.4%.** Consistent, close to the
        ~2.3% previously estimated, and an accepted cost of staying in 5072 (documented in
        `12A0_v5_truncate.R`: `clamp()` does not commute with projection, so projecting
        would break the superset -> masked-rowsum decomposition).
      - **A->B, raster-vs-vector masking: +6.2%, +1.4%, +12.7%.** Not uniform, so not a
        simple edge effect. Isolated per term, and it was a **single root cause**:
        `terra::rasterize()` defaults to `touches = FALSE` (cell-centre rule), while
        V5's `mask(vect, inverse = TRUE)` and `crop(vect, mask = TRUE)` both behave as
        `touches = TRUE`. Per-term ratios (ours / V5): range 1.0000 (control - it is a
        raster in both paths), data-limit 0.9998, BCR polygon 0.978-0.993, **water
        1.049-1.077**. Water dominates because it is by far the most fragmented layer -
        29,545 polygons, mostly thin rivers and small lakes that touch a cell without
        covering its centre. Setting `touches = TRUE` reproduced V5 **to the digit** on
        every term (water 874105 vs 874105, 370971 vs 370971).

        A hypothesis that proved WRONG and is worth not re-testing: that `crop_to_grid`
        was dropping polygons by cropping in the source CRS against a straight-edged
        projected extent rectangle. It does drop some (36 on OVEN can10, 165 the other
        way on CAWA can71), but they are outside the grid - the sums agree to 0.01%.

- [x] **G3/G4 PASSED 2026-09-11.** With the BCR cut and `touches = TRUE` both in,
      `sum(observed x weight x 100)` reproduces V5's own vector masking at
      **A/B = 1.000011 / 1.000002 / 1.000001** (CAWA can10, CAWA can71, OVEN can10) -
      float32 noise. The masking half of the conformance work is now validated against
      the PRODUCTION config (5072 + `weight.tif`), not just the G2 harness.

      Only departure left is the deliberate one: staying in EPSG:5072 costs **+2.9% to
      +3.4%** vs V5's legacy 3978 delivery projection (previously estimated at 2.3%).
      That is the documented trade for keeping `clamp()` commutable with the superset ->
      masked-rowsum decomposition; it applies equally to obs and bf, so it cancels in the
      contrast. Harness: `Rscripts/misc/verify_weight_vs_v5_masking.R`, which prints all
      three configurations and is the thing to re-run after any weight change.

      Note this defect moved **absolute populations only**. The weight multiplies both
      sides, so `w*bf - w*obs = w*(bf - obs)` held throughout and Shapley *shares* were far
      less distorted than totals — but every per-BCR and national total produced before the
      fix is wrong.
- [x] **B4. DONE 2026-09-11.** Updated `12C_predict_species_all_coalitions.R`:
      - `$q` → `$densmax`, behind a schema guard that `stop()`s rather than letting
        `pmin(x, NULL)` return `numeric(0)` silently.
      - dropped `q0`. Note `l.out` **still ships** in the `.Rdata` (151 rows, `denshthresh`
        intact) — V5 deleted the *step*, not the object, so this was dead weight rather
        than a latent error.
      - loads `truncation_params.rds` once per species; per BCR it looks up the frozen
        `q99`, `stop()`s if absent, and cross-checks `densmax` against `q.out` so 12A and
        12C can never read different `.Rdata` versions without failing loudly.
      - `:329` now `pmin(pmin(pred_vec, qsp), q99)`, still ahead of the weight multiply
        (V5 order: truncate → mask).
      - logs `caps: densmax=... q99=... (q99 binds Nx lower)` per BCR.
      Observed side needs no change: it is read pre-clamped from `observed_bootstraps.tif`,
      so both sides of the contrast now carry identical caps.
      `CLAUDE.md`'s "12C restructure invariants" bullet was rewritten — it previously
      asserted the opposite ("q99/q0 caps only ever touched inspection rasters").

- [x] **B5. DONE 2026-09-11.** Restaged `data/raw_data/SpeciesPredictionTruncationValues.Rdata`
      from `G:/Shared drives/BAM_NationalModels5/data/` (2026-06-04 version). The old file is
      recoverable from git history if the pre-rewrite `$q` values are ever needed.
- [x] **B6. DONE 2026-09-11 — and the premise was wrong.** The task assumed "no change
      needed, check only NA→0". The NA→0 half checked out; the "no change needed" half did not.

      *NA→0 (`10.Truncate.R:141`) — equivalent, no change.* V5 does `mask.i <- truncate2.i *
      range.i; mask.i[is.na(mask.i)] <- 0`. Ours reaches the same totals by two routes:
      `weight.tif` is non-NA everywhere by construction (range NA→0 via `classify`, every
      `rasterize` has `background = 0`), so mask-induced NA becomes a hard 0 exactly as in V5;
      and prediction-origin NAs are zeroed at aggregation instead (`zonal(..., na.rm = TRUE)`
      at 12C:363, `Ok[is.na(Ok)] <- 0` at 12C:429). V5's step-8 crops leave NA where we leave
      literal 0, which is identical under a sum. Added an explicit `stop()` in 12A2 if any
      weight cell is NA, since that assumption is now load-bearing.

      *What the check actually turned up:* the missing BCR-polygon crop — see G3/G4 above and
      Open Limitation #7.

      *Third finding, guarded not fixed:* the obs side zeroes its NAs but the bf side's
      `rowsum(M[keep, ])` has no `na.rm`, so one NA would NA out a whole subbasin on the
      backfilled side only. And `complete_mask` is built from **draw column 1** as a proxy for
      all 100 draws, while each scenario samples a different draw — so a pixel complete in
      draw 1 can be NA in the draw actually used. 12C now audits M column-by-column once per
      BCR and `stop()`s with the count. If it ever fires, widen the gate to all draws.
- [x] **B7. DONE 2026-09-11.** Globus'd all 25 regenerated `observed_bootstraps.tif` +
      both `truncation_params.rds` to the cluster, one `globus transfer` call per file
      (never `--batch`). 27/27 SUCCEEDED, 0 faults, and every transferred byte count
      matches the local file size exactly (1,981,925,510 B ≈ 1.85 GiB). No file reported
      0 bytes, i.e. none was a checksum no-op — consistent with all 25 having been
      rewritten in the B3 re-run. `weight.tif` alongside them was untouched, as intended.

- [x] **B8. DONE 2026-09-11.** `14B_sector_attribution.R` now drops withheld models at read
      time via `DROP_WITHHELD` + a `withheld_models` table, with a message reporting the row
      count (and a distinct warning if nothing matched, which would mean the BCR code changed).
      Verified against the workbook rather than trusting the note: the "remove" tab has 662
      rows, of which exactly one touches our species — CAWA can40, AUC = 1; OVEN has none.
      Kept upstream production of can40 deliberately, so the products stay a complete record
      of what we ran and only the reported numbers are filtered. Commented the read-from-xlsx
      one-liner to swap in at the planned ~60-species scale.

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
- [x] **C4. DONE 2026-09-11.** `CLAUDE.md`: pre-`12B` wipe now also clears
      `density_tables/arrays/*.rds` (with the `rm -f` not `rm -rf` caveat); the two
      `10.Package.R` references retargeted to `10.Truncate.R` via `12A0_v5_truncate.R`; the
      `12A2` row and Phase 1b block rewritten for the four-term, version-stamped weight; new
      **Open Limitation #7** for the BCR-cut double count; #6 updated to record that the
      `14B` release filter is now in place. Memory files updated.

---

## Loose ends

- [ ] `12A` loads each full `b.list` `.Rdata` (up to 654 MB, off Google Drive) purely to read
      `attr(b.list[[1]], "bcr")` — ~4.3 GB of I/O per species to extract 11 strings. The BCR code
      is in the filename, and the `12B` preflight already derives it that way. Verify the filename
      always equals the attribute, then drop the load. Matters at the planned ~60-species scale.

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
