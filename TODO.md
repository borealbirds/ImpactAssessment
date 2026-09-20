# TODO — consolidated plan (last updated 2026-09-15)

Open work is the cluster backfill re-run (**A**) and the converge steps (**C**). Workstream
**B** (conform density to current V5 packaging) is finished and validated — gates G1–G4 all
passed. The diagnosis behind A is `CLAUDE.md` Open Limitation #5; refuted hypotheses live in
the `project_pipeline_history` memory.

## Critical path

```
cluster:  A6 sbatch 07 (674 tasks) ─► sbatch 11 ─┐
                                                 ├─► 12D smoke ─► C1 12D ─► C2 14B ─► C3
local:    [B, G1–G4 all DONE] ───────────────────┘
```

`07` is the long pole and depends on nothing in B, so **A6 is the next action**.

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

All fixes (CAfire recode, Fix A, Fix B, the `12D` weight preflight) are committed as of `6455756` and staged; the compute has never been launched. Prerequisites
(cleanup, weight rebuild 25/25, `07` smoke) are all done — see **Completed**.

- [ ] **A6.** ~~`sbatch 07_train_and_backfill.sh`~~ **DONE — 674/674 on 2026-09-19.** Remaining
      half is `sbatch 11_premosaic_backfilled_stacks.sh`, which must NOT run until the
      `_confusion.rds` count reads 674 (`11:71` silently drops non-existent paths and `11:189`
      skips BCRs whose mosaic already exists, so a hole would be baked in permanently).
      **Run 1 (jobs `60009763` large / `60010888` small) finished 668/674 on 2026-09-16.** The
      six gaps were all `OUT_OF_MEMORY` with MaxRSS pegged exactly at the request: 57, 62, 98 at
      750G and 454, 474, 583 at 64G. Cause was **not** subbasin size — S454 (`ncell=45288`) died
      at 64G while S5 (`ncell=44896`) succeeded — but allocation churn in `08A`'s raster
      assembly, where `result_raster[[j]] <- v` copied the whole ~2000-layer stack on each of
      ~2000 assignments. Patched 2026-09-16 (`73e8b68`) to fill one `ncell x nlyr` matrix and
      call `terra::values()` once. Mixing patched and unpatched outputs is safe — the patch
      changes assembly only, no modelled value.
      **Run 2 (job `60134338`, the six reruns) all COMPLETED `0:0` on 2026-09-17**, confirming
      the diagnosis quantitatively: S57 50.8 GB / S62 44.4 GB / S98 50.2 GB (all three had
      pegged a full 750G node) and S454 7.4 GB / S474 9.5 GB / S583 7.3 GB (all three had OOM'd
      at 64G). Roughly a 15× drop; none of these subbasins was ever big.
      **Consequently the array tiering is retired (2026-09-19).**
      `07_train_and_backfill_larger.sh` is deleted and `07_train_and_backfill.sh` now carries a
      single `--array=1-674%30` at 128G / 8h / `--cpus-per-task=1` (single-threaded:
      `BART::gbart`, not `mc.gbart`). No more complementary index lists to keep in sync, and no
      more race risk from overlapping arrays. One index OOMing is now handled by
      `sbatch --array=<i> --mem=256G 07_train_and_backfill.sh`, not by re-tiering.
      **Diagnostic note**: neither `.out`/`sacct` state nor `done` in `logs/Y2020_S*.log` proves
      success. `07_train_and_backfill.R`'s `tryCatch` makes a failed subbasin exit `COMPLETED`,
      and `08A:24`'s `on.exit(logp("done"))` fires on error unwind. Logs also *append* across
      runs, so they cannot isolate one run. The oracle is the filesystem: `_confusion.rds` is
      written last, plus an mtime check so pre-CAfire-fix leftovers cannot pass as fresh.
- [ ] **A7.** Validate Fix A — the only fix never tested. The `12F` runtime line
      `complete superset pixels: N / M (X%)` must be ≫ the old 1–3%.

## C. Converge

- [ ] **C1.** `sbatch 12D_repredict_all_coalitions.sh` (`--array=1-2`). The pre-run wipe of
      `density_tables/*.rds` **and** `arrays/*.rds` already happened; redo it only if anything
      writes there first. Old tables are stale on three counts: pre-weighting, pre-gate-change
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

- [ ] `12A` loads each full `b.list` `.Rdata` (up to 654 MB, off Google Drive) purely to read
      `attr(b.list[[1]], "bcr")` — ~4.3 GB of I/O per species to extract 11 strings. The BCR code
      is in the filename, and the `12D` preflight already derives it that way. Verify the filename
      always equals the attribute, then drop the load. Matters at the planned ~60-species scale.
- [ ] Decide the fate of local `covariates_mosaiced_2020_PREORIG.tif` (1.5 G) — the only
      surviving pre-CAfire-fix 2020 mosaic.
- [ ] `15C_singletons_plot.R:150` reads `predictions_coalitions/`, deleted on both ends, so the
      map panel is broken. Regenerate via `save_arrays_ids` if wanted.

---

## Completed

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
`pmin(pmin(pred_vec, densmax), q99)`; `SpeciesPredictionTruncationValues.Rdata` restaged;
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
