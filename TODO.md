# TODO — consolidated plan (last updated 2026-09-23)

Open work is the cluster backfill re-run (**A**) and the converge steps (**C**). Workstream
**B** (conform density to current V5 packaging) is finished and validated — gates G1–G4 all
passed. The diagnosis behind A is `CLAUDE.md` Open Limitation #5; refuted hypotheses live in
the `project_pipeline_history` memory.

## Critical path

```
cluster:  [A6 07 + 11 DONE] ─► 12D smoke = A7 ─► wipe ─► C1 12D ─┐
local:    [B, G1–G4 DONE]                                        ├─► C2 14B ─► C3
```

A6 closed 2026-09-22, so **the `12D` smoke is the next action** — and it doubles as A7.

**Blocker for every cluster step**: SSH is keyboard-interactive (2FA), so cluster commands must
be run from your own authenticated session. Globus works (one file per call, **never** `--batch`).

---

## Staging discipline — read before any cluster run

Derive the set of files to stage from the **`source()` closure of the entry point**, never from
"what did I change". Then checksum-sync the whole closure: a 0-byte transfer is a free proof of
equality, so there is no reason to sync only the suspects. Three separate times a file that had
never been edited — and so appeared in no diff-derived list — turned out to be stale or absent:

- `08A` (carrying Fix B) was nearly left stale before the 674-task `07` run, which would have
  reproduced the exact bug that run exists to fix.
- 7 of the 9 files in `07`'s source chain were pre-2026-07-02 versions, including
  `08B_deploy_gbart.R`, which would have written a different backfill layer set.
- `12E_shapley_utils.R` had never been on the cluster at all; `12D:73` sources it, so C1 would
  have died on its first `source()` after queueing for a 384 G node.
- The closure includes **data files**, not just scripts: `SpeciesPredictionTruncationValues.Rdata`
  was recorded as "restaged" in B but Fir still had the old-schema copy, and the A7 smoke died
  on `12F`'s `densmax` guard.

Two mechanical gotchas found the same way:

- `.sh` files MUST be transferred as LF. A CRLF script dies on Linux with
  `/bin/bash^M: bad interpreter`. Verify with `tr -dc '\r' < f | wc -c` (0 = LF), **not**
  `grep -c $'\r'`, which silently matches the letter `r` in some shells.
- Globus refuses sources outside the local endpoint's configured root, so a temp-dir staging
  copy fails with `Path not allowed`. Stage from inside the repo tree.

**Current state (2026-09-16)**: `/Rscripts` on Fir matches HEAD with no known divergence.
The four files held back while A6 was queued — `08A_train_and_backfill_subbasin_s.R`,
`08B_deploy_gbart.R`, `12D_repredict_all_coalitions.R` and `.sh` — were all re-staged once the
queue drained. `/Rscripts/12*` is exactly six files (old-named duplicates removed after the
renumbering; `12A` is local-only by design).

---

## A. CAfire / phenology backfill re-run (cluster)

All fixes (CAfire recode, Fix A, Fix B, the `12D` weight preflight) are committed and staged.
Backfill and premosaic (A6) are done for 2020 — see **Completed**. Only A7 remains.

- [x] **A7.** Validate Fix A — the only fix never tested. Run it on the `12D` smoke, not the
      full run: the `12F` runtime line `complete superset pixels: N / M (X%)` must be ≫ the old
      1–3%. If it is not, stop — do not burn the 2 × 24 h C1 run.
      `sbatch --array=1 --time=01:00:00 --mem=192G --export=ALL,TEST_BCR=can60,TEST_N_BOOT=2 12D_repredict_all_coalitions.sh`
      The smoke writes real files: `CAWA_2020_coalition_*.rds` into `density_tables/` (can60
      only, 2 bootstraps) and per-pixel arrays into `arrays/`. **Wipe both again before C1.**
      **Smoke 1 (job 61147703, 2026-09-23): 13097 / 113017 (11.6%)** — up from 1–3% but not ≫,
      so A7 is NOT closed. The denominator includes the stack grid's 100 km buffer, which `11`
      masks to NA and weight zeroes anyway; `12F` now also logs coverage among weight > 0
      pixels plus per-var / per-subbasin NA counts.
      **Smoke 2 (job 61160614): 11869 / 20931 weight > 0 (56.7%).** Every one of the 9062 lost
      pixels is NA in one var, `SCANFIBalsamFir_5x5`, across 24 subbasins (most lost in full).
      Cause: `08A:101–116` writes `<cov>_mean` instead of `_draw_*` where a covariate is constant
      in a subbasin; the BCR mosaic then has draws elsewhere but NA there, and the gate drops it.
      Fixed in `12F`: draw-less pixels take the `_mean` constant for every draw (logged as
      `filled N draw-less superset pixels`). Staged.
      Smoke 3 (61162614) hung on node `fc30537` (`ALLOCATED+NOT_RESPONDING`) — node fault, not code.
      **Smoke 4 (job 61167742): 20931 / 20931 weight > 0 (100%). A7 PASSES.** Fills: BalsamFir
      11929, DouglasFir 24599, LodgepolePine 131/225 (1km/5x5). M audit silent; 255 tables written.
      Sampling took 90 s vs 52 s (1.9× the complete pixels) — expect C1's gbm stage to scale alike.
      **Follow-up before C1**: a covariate constant in *every* subbasin of a BCR had no draws, so
      it kept its observed value on the bf side. `12F` now carries it as a one-draw covariate from
      `<cov>_mean` (logged `no draws in BCR, using <cov>_mean`); it also enters the gate.
      **Smoke 5 (job 61170330, all CAWA BCRs, 2 boots, 3 h): hit the time limit after 3 BCRs.**
      can10 100294/100294, can11 427690/427690, can12 132604/132604 — all 100%. The partial fill
      is large there (can11: PonderosaPine/WhiteRedPine 329k pixels each). No `no draws in BCR`
      line yet, so the mean-only path is still unexercised. Sampling time per BCR (≈ one C1 wave):
      can60 1.5 min, can10 17 min, can11 **2 h 10 min**, can12 > 13 min. C1 runs 32 boots on 16
      cores = 2 waves, so can11 alone is ~4.5 h; the old 24 h envelope predates the fill.
      C1 overwrites every coalition it writes, but `12D:131` skips an empty coalition without
      writing, so a smoke table could survive into C2 unnoticed.
      **Speedups before C1 (2026-09-23, bit-identical):** `12F` now predicts once per *distinct*
      BART draw per bootstrap (~64 of 100 scenario picks are distinct under our seeds → ~36% fewer
      gbm calls), skips weight-0 pixels (2–15% of complete pixels), and reads draws only for model
      covariates (~17 of ~50; the draw-load phase was 5–7 min/BCR). Tested old vs new on CAWA can11's
      real models: `M * weight` identical, 1.87× faster. Rejected: a static/dynamic gbm tree split
      (67–97% of trees split on a backfilled covariate; 1.11×). Not done: summing pixels by
      (subbasin, sector-signature) instead of `M[keep,]` per coalition — 11× faster reduction and
      no `M`, but only equal to ~1e-15, not bit-identical.
      **`arrays/` was never written** by the superset path — `save_arrays_ids` was accepted but
      unused since `a428ee7`, so `15A` had no input. Restored: `12F` returns national
      `[n_boot × n_scen]` `obs_total/obs_on_coal/bf_on_coal` matrices for the 9 target coalitions
      and `12D` writes `arrays/{species}_{year}_coalition_{cid}_arrays.rds`. Unlike the tables
      (and the retired code, which dropped the BCR), a BCR with no footprint for that coalition
      contributes its real `obs_total` to the arrays, so 15A's national observed total stays whole.
      **Smoke 6 PASSED** (job `61306338`, CAWA can60, 2 boots; 2026-09-24): 20931 / 20931 weight > 0
      complete (as smoke 4), `predicting 20931 … skipping 4095 with weight 0`, M NA audit silent,
      `wrote 255 coalition tables`, `wrote 9 array files`, `nice.`. Sampling 55 s vs smoke 4's 91 s
      (1.65×); observed+draw load 66 s vs 98 s. Arrays are 2 × 100 and their means equal the tables'
      subbasin sums exactly (cid 256, 2, 9); cid 9 has no can60 footprint and correctly carries the
      real `obs_total` with zero on/bf. The mean-only-covariate path did not fire in can60 — first
      real exercise is C1. (Smoke 61304676 before it was an `--mem=192M` typo, not a code fault.)

## C. Converge

- [ ] **C1.** `sbatch 12D_repredict_all_coalitions.sh` (`--array=1-2`). Redo the wipe of
      `density_tables/*.rds` **and** `arrays/*.rds` first — the A7 smoke writes to both.
      Old tables are stale on three counts: pre-weighting, pre-gate-change
      (A), pre-truncation-conformance (B).

      **Watch for one guarded-not-fixed condition.** `complete_mask` is built from draw column 1
      as a proxy for all 100 draws, but each scenario samples a different draw — so a pixel
      complete in draw 1 can be NA in the draw actually used, and the bf side's
      `rowsum(M[keep, ])` has no `na.rm`, so a single NA would NA out a whole subbasin on the
      backfilled side only. `12F` audits `M` column-by-column once per BCR and `stop()`s with a
      count. If it fires, widen the gate to all draws.
- [ ] **C2.** `14B_sector_attribution.R` locally → corrected Shapley CSVs.
- [ ] **C3.** Sensitivity pass with the `q99.9` stage disabled; report the spread. The frozen cap
      has a known-direction bias — counterfactual densities are higher, so they hit the frozen
      ceiling more often than observed, systematically **under-estimating impact in the
      highest-density pixels**. With the cap removing 4–12% of abundance that is not negligible.
      Headline = conformed; sensitivity = uncapped.

## Loose ends

- [ ] **Move `12A` to the cluster** (prerequisite for the planned 100+ species). Evidence
      gathered 2026-09-23 that Elly's `def-ecknight/NationalModels/output/` is a sound source:
      `06_bootstraps` and the 2020 `07_predictions` pair 1:1 (2,868 each, 1,686 Canadian); a
      name+size comparison against G: over 32 of 151 species matched 4,504/4,505 common files
      (the odd one is `AMPI_can3_1990`, not 2020); `q.out` covers exactly the same 151
      species. `12A` reads only raw `07_predictions` + `06_bootstraps` filenames + our own
      `SpeciesPredictionTruncationValues.Rdata`, so it does not depend on whether Elly ran
      V5's truncation revision there. Remaining proof of byte-equivalence: run `12A` on the
      cluster for CAWA into a scratch dir and diff against the G:-built
      `observed_bootstraps.tif` / `truncation_params.rds`. Then set `cc <- TRUE`, add a `.sh`,
      stage `12B_v5_truncate.R` (production path uses no G: resource: `apply_masks = FALSE`,
      `project_to = NULL`), and drop the Globus step from CLAUDE.md.
      **Done 2026-09-23**: `12A` no longer loads each `b.list` `.Rdata` (up to 654 MB, ~4.3 GB
      per species) just to read `attr(b.list[[1]], "bcr")`; it parses the BCR from the
      filename. `12F`, which loads `b.list` anyway, now `stop()`s if the attribute and the
      filename disagree. Verified equal for all 25 CAWA/OVEN pairs.
- [ ] Decide the fate of local `covariates_mosaiced_2020_PREORIG.tif` (1.5 G) — the only
      surviving pre-CAfire-fix 2020 mosaic.
- [ ] `15C_singletons_plot.R:150` reads `predictions_coalitions/`, deleted on both ends, so the
      map panel is broken. Regenerate via `save_arrays_ids` if wanted.

---

## Completed

**A6 (backfill + premosaic), 2026-09-16 → 09-22**
- `07`: 674/674. Run 1 left six `OUT_OF_MEMORY` gaps (57, 62, 98 at 750G; 454, 474, 583 at
  64G) caused by allocation churn in `08A`'s assembly, not subbasin size. Fixed in `73e8b68`;
  the six reruns peaked at 7–51 GB. The 64G/750G array tiering is retired (`12f1ad6`): one
  `--array=1-674%30` at 128G. Oracle for `07` is `_confusion.rds` + mtime — `tryCatch` makes a
  failed subbasin exit `COMPLETED` and the logs append across runs.
- `11`: 19/19 BCR mosaics, job `60806856`, written 2026-09-21 16:20 → 09-22 11:00, every
  `.out` ending `done.`. A first submission (`60564008`, 09-19) was a no-op: the 19 stale
  2026-05-24 mosaics had never actually been deleted, so `11:187`'s skip-if-exists guard
  `quit(status = 0)`'d every task in seconds and `sacct` reported all 19 `COMPLETED 0:0`.
  **Oracle for `11` is the `.out`, not `sacct`**: three code paths exit 0 without writing
  (`:189` exists, `:202` no subbasins, `:216` no backfills), and `terra::mask(filename=)`
  streams, so a killed task leaves a readable truncated `.tif` that the guard would later
  accept. Only `:241–242`'s `masked and written to` / `done.` prove the write returned.

**A (cluster prep), 2026-09-11 → 09-14**
- Cleanup: pre-fix `covariates_mosaiced_{1990..2015}.tif`, `predictions_coalitions/`,
  `bart_models/2020`, `bart_models_mosaics/2020`, `density_tables/*.rds` + `arrays/*.rds`.
- `covariates_mosaiced_2020.tif` is the CAfire-fixed build on both ends (checksum sync moved
  0 bytes). No `06` re-run needed.
- **weight.tif rebuilt 25/25** — every one logged as stale against the `weight_v3_touches`
  stamp, i.e. all were the old range-only rasters, and the BCR cut is the dominant masking term.
- **`07` smoke (array 1–3) validated Fix B at scale**: subbasin 1 carries values
  at 1784/1784 high-HF pixels, **100.00%** complete across all 17 `_draw_*` covariates — the
  direct successor to the 0.8% figure that opened the investigation. `np` constant down the
  whole hierarchy, no NaN cascade.

**B (V5 packaging conformance), 2026-09-11.** `12B_v5_truncate.R` ports `10.Truncate.R`;
`12A` rewritten and re-run (all 25 stacks carry both caps + a frozen per-BCR `q99` in
`truncation_params.rds`); `12F` reads `$densmax` behind a schema guard and applies
`pmin(pmin(pred_vec, densmax), q99)`; `SpeciesPredictionTruncationValues.Rdata` restaged
locally (the Fir copy was NOT — it was still the 2026-05-04 old-schema file until the first A7
smoke hit `12F`'s schema guard on 2026-09-23; Globus'd then, 6387 bytes);
`14B` drops withheld models (`DROP_WITHHELD`; the workbook's "remove" tab flags exactly one of
our 25 pairs — CAWA can40); all 25 observed stacks Globus'd to the cluster (27/27, byte-exact).
Gates: **G1** reproduces V5's `10_truncated` (flat terrain bit-identical). **G2** reproduces
BAM's published abundances to full published precision, point estimate and both 5–95% bounds.
**G3/G4** reproduce V5's own vector masking in the *production* config to A/B = 1.000011 /
1.000002 / 1.000001. Harnesses: `Rscripts/misc/verify_v5_truncate_port.R` and
`verify_weight_vs_v5_masking.R` — re-run the latter after any weight change.

---

## V5 reference

Renumbering in `f082866`: `10.Package.R` → `10.Truncate.R` + `11.Package.R`;
`11.Validate` → `12.Validate`; `12.Summarize` → `13.Summarize`;
`output/10_packaged/` → `output/11_packaged/`, new `output/10_truncated/`.
**`06_bootstraps` and `07_predictions` were NOT renumbered** — they are the only V5 output
folders our scripts read, so no path in our pipeline is broken.
