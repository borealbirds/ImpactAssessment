# TODO (last updated 2026-10-04)

> **Where we left off (2026-10-04).** The C2f rerun is done end to end for 2020: 07 + 11 on the
> new masks, C1 (12D + 12H) and 14B. The headline numbers are under Next action 5. C3 is done
> too (Next action 6): without the frozen `q99.9` cap v(N) moves +1.2% (CAWA) and +0.05% (OVEN).
> The first new-species batch (OSFL, GRSP, LEYE, BOBO) is under way (Next action 7).
> The user wants the project finished: log small new problems under Loose ends and do not
> reopen finished stages for them. Keep `writing/` off the internet: never commit or push it.
>
> **Uncommitted as of 2026-10-02:**
> - the 11 window/merge rewrite: `11.R`, both `.sh` files, and `misc/compare_mosaics.{R,sh}`;
> - the C3 `DROP_Q99` switch in 12D, 12F, 12H and 14B (staged on Fir);
> - species as arguments (`SPECIES`) in 12A, 12C, 12D, 12H and their `.sh`; 12A on the cluster;
>   12B without `library(sf)`; 14B's withheld list from the workbook; `.gitignore` (staged on Fir,
>   14B local only);
> - `CLAUDE.md` and this file.
>
> Not this session's to commit:
> - `14A_reproject_hirshpearson.R` has a comment-only change that nobody in this session made.
>   Ask the user before committing it.
> - `data/derived_data/rds_files/`, `sector_effects/` (now the C2f results) and
>   `sector_effects_noq99/` (C3) stay out of git until the user decides.
>
> The full record of finished work is in `git show 3fcef6d:TODO.md` (trimmed 2026-09-28). The
> 07/11 checking recipe that used to open this file is in `git show 994d90b:TODO.md`.

## Globus from this laptop

- Write filters as `--opt=value` (`"--include=slurm-<id>_*.out" "--exclude=*"`): a bare `*`
  is glob-expanded on Windows. `--exclude=by_bcr` did not keep that directory out of a
  recursive pull (2026-10-02).
- From Git Bash, commands that take a Fir path as its own argument (`globus rename`) need
  `MSYS_NO_PATHCONV=1`, or the path is rewritten to `C:/Program Files/Git/...`.
- Transfers stuck ACTIVE with `GC_NOT_CONNECTED` mean Globus Connect Personal is not running.
  Start `C:\Program Files (x86)\Globus Connect Personal\bin\globus_connect_personal.exe`.
- A log still being written fails checksum verification over and over; pass
  `--no-verify-checksum` to pull it.
- Globus Connect Personal cannot read Claude's scratchpad under `AppData` (`PERMISSION_DENIED`,
  task stuck ACTIVE); stage outgoing files from the project folder.
- Fir's `sacct -S/-E` times are Pacific, not UTC (Globus listings are UTC). Logs land in
  `Rscripts/` on Fir; sort a listing by its time column before concluding nothing ran.

## Critical path

```
cluster:  [07 + 11: done] ─► [C1: done] ─► [14B: done] ─► [C3: done]
```

**Next actions, in order:**

1. ~~**Check 07 + 11**~~ (2020, new masks), done 2026-10-02.
   - 07: all 674 `.out` end `$ok TRUE` with no error/OOM/timeout line. Pull them in ONE
     filtered Globus transfer: `--recursive "--include=slurm-<07 job>_*.out" "--exclude=*"`.
     `sacct` is not enough: `tryCatch` makes a failed subbasin exit `COMPLETED`, and the logs
     append across runs. So also check that every `subbasin_{i}_confusion.rds` postdates
     2026-09-28 19:00 UTC and that the last run in every `logs/Y2020_S{i}.log` has
     `footprint covariates set to 0`.
   - **07 = job `61937625`, checked 2026-10-01** (`cluster_logs/07_2020_c2f/`): 672 `.out` end
     `$ok TRUE` with no error line (tasks 1 and 2 have no `.out` on Fir, but their logs show
     complete runs). In all 674 logs the last run sets the footprint covariates to 0 and ends
     `done`, and 673 of 674 changed their pixel counts from the old-mask run (423 has one
     backfill pixel both times). Log times (Fir local) run 09-28 12:11 → 09-29 18:44; the
     "19:36 UTC" above is wrong, so run the `find` check with `-newermt '2026-09-28 12:00'`.
     Runtime: median 15 min, max 4.9 h, 302 task-hours, driven by training pixels (r = 0.76).
     The `%30` cap bound only 2% of the 30.6 h span; for 16% nothing ran, so the scheduler,
     not the cap, set the pace. `STALE_CONFUSION=0` on Fir.
     `sacct` (`cluster_logs/07_2020_c2f/sacct_07_61937625.txt`; `Rscripts/misc/sacct_summary.R`):
     all 674 COMPLETED, CPU efficiency 0.98. Peak memory median 3.8 GB, 99th percentile
     15.6 GB, max 52.8 GB. Only subbasins 57, 61, 62, 98 and 107 exceed 20 GB (25–53 GB; the
     can11 prairie basins). Peak memory is predictable from subbasin size (ncell, backfill
     rows, training pixels; residual SD 0.6 GB). At 128G every task bills 33 core-equivalents:
     10,037 core-eq-h for this run. 24G for all but those five, and 72G for them, would have
     billed 1,926 core-eq-h (see Loose ends).
   - **11: done 2026-10-02, all 19 BCRs** (`cluster_logs/07_2020_c2f/11/`; every `.out` ends
     `done.`). Job `62431589` (64G / 12 h) built 16 of them with the `278bc85` `11.R`, which
     caps terra at 30% of the allocation (`memmax`). They peaked at 5.6–33.5 GB and ran 0.57–1.0×
     A6's wall time (`sacct_11_62431589.txt`). can41 and can5 stalled at 64G: one modest
     subbasin each, hours at 100% CPU, flat ~5 GB, no disk I/O. They were rebuilt at 512G (job
     `62468313`, 52 and 67 min). Both stalled inside the whole-grid resample, which the fix
     removes.
     can3 (71 subbasins, ~9M-cell grid) projected at 11–15 h, and its 24 h job pended a day.
     **The fix:** `11.R` resamples each subbasin onto its own grid-aligned window and
     `merge()`s per variable. The test build (job `62557048`, 32G / 3 h, `MOSAIC_DIR`) matched
     the production can11, can42 and can81 in every layer (`misc/compare_mosaics.R`, job
     `62557069`). It built can3 in 1.4 h and can81 in 2.0 h (was 5.0 h), peaking at 3.1–13.0 GB
     of 32G. can3's mosaic was moved
     into `bart_models_mosaics/2020/` and the test folder deleted. `.sh` defaults are now
     32G / 4 h, and `07_train_and_backfill.sh` chains one 11 array (no 24 h can3 split).
     Fir has this `11.R`, the new `.sh` and `07_train_and_backfill.sh`; none committed yet.
     If a task OOMs: `YEARS=2020 sbatch --array=<i> --mem=64G 11_premosaic_backfilled_stacks.sh`.
2. ~~10C~~, done 2026-09-28: 4 of 674 flagged (50, 72, 430, 481; was 72, 98, 481). Subbasin 98
   fell from 0.51 to 0.14 outside the AOA. `frac_outside_aoa` correlates 0.89 with the old run;
   median 1.6%, 90th percentile 14.3%. Log: `cluster_logs/10C_2026-09-28.log`.
3. ~~**Wipe, then rerun C1**~~, done 2026-10-02: tables wiped; 12D = job `62596963`, 12H =
   `62597111` (logs: `cluster_logs/C1_2026-10-02_c2f/`). All 25 tasks `nice.`. In each, the
   identity gate passes for all 32 bootstraps at 1,123–3,117 observed-design pixels.
   Complete weight > 0 superset pixels: 100% everywhere (CAWA can81: 75,908 of 75,919).
   The direct footprint is 61–95% of the superset (OVEN can41 61%, can11 95%). Tasks ran 1–44 min.
   12H: 255 tables per species, 9 arrays each, and the Shapley samples (CAWA 734 and OVEN 831
   subbasin rows × 3,200).
4. ~~**Pull**~~, done 2026-10-02: 512 tables + Shapley samples and 18 arrays, all dated today
   locally (one recursive transfer; `--exclude=by_bcr` did not exclude, so `by_bcr/` came too).
5. ~~**14B**~~ locally, done 2026-10-02 (`logs/14B_2026-10-02_c2f.log`; CAWA can40 withheld).
   Sample means match the tables' Shapley values (max |diff| 4.7e-10), and sum(phi) = v(N) in
   every sample. **CAWA v(N) = +757k** on 4.47 M observed (+16.9%; 5–95%: 649k–870k). Shares
   of v(N): roads 42%, pasture 22%, crop 21%, built 11%, rail 2.2%; oil and gas, mines and
   dams < 0.4% each. Direct +138k, vegetation +619k. **OVEN v(N) = +4.59 M** on 37.65 M
   (+12.2%; 3.70–5.18 M). Shares: crop 35%, pasture 32%, roads 19%, built 13%, rail 1.4%;
   mines −19k, dams −15k, oil and gas ≈ 0. Direct −1.68 M, vegetation +6.27 M; roads is
   −1.27 M direct and +2.14 M vegetation, so roads carries the widest interval
   (325k–1.21 M). The 4 AOA-flagged subbasins carry 23k (CAWA) and −3k (OVEN).
   The first C2 run (2026-09-24: old masks, pre-C2e, can40 withheld) gave CAWA +574k (roads
   34%, pasture 22%, crop 21%, built 18%) and OVEN +3.53 M (crop 33%, pasture 32%, built 21%,
   roads 11%).
6. **C3:** the uncapped `q99.9` sensitivity pass. Built 2026-10-02 as `DROP_Q99=1`, read by
   12F, 12D, 12H and 14B. Only densmax caps the observed, d0 and bf sides. 12F re-predicts the
   observed footprint, because `observed_bootstraps.tif` is already q99-clamped, and stops
   unless re-capping that prediction at q99 reproduces the tif at every predicted pixel and
   bootstrap. Output goes to `density_tables_noq99/` and `sector_effects_noq99/`. The %
   denominators stay the capped observed totals, because only the footprint is re-predicted.
   Submitted as a can60 2-bootstrap smoke, then the full 12D (`afterok`), then 12H:
   `S=$(sbatch --parsable --array=1-2 --time=01:00:00 --export=ALL,TEST_BCR=can60,TEST_N_BOOT=2,DROP_Q99=1 12D_repredict_all_coalitions.sh)`,
   then `--dependency=afterok:$S --export=ALL,DROP_Q99=1` for 12D and `afterany` for 12H.
   **Done 2026-10-02/04.** Smoke = job `62636061`, 12D = `62636072`, 12H = `62636104` (logs:
   `cluster_logs/C3_2026-10-02_noq99/`). All 25 tasks `nice.`, with the identity gate passing for
   all 32 bootstraps, and the observed side reproducing the tif once re-capped at every predicted
   pixel × bootstrap. q99 binds at only 0–1.02% of those (OVEN can12 the most). 14B
   (`logs/14B_noq99_2026-10-04.log`; `sector_effects_noq99/`): means match the tables
   (max |diff| 5.2e-10), additivity exact. **Uncapped, v(N) moves +1.16% for CAWA** (756,594 →
   765,389; 5–95% 650k–886k) **and +0.05% for OVEN** (4,591,911 → 4,594,430). By sector: CAWA
   roads +6.1k (+1.9%), every other sector < 1k; OVEN roads −3.1k and pasture +3.7k, every other
   < 1.4k. Shares do not move. By BCR, the biggest shift is CAWA can61 +5.5k (+2.1%) and OVEN
   can14 +10.7k (+1.6%). Large % changes occur only where an impact is near 0 (CAWA can80
   −1,945 → −1,262). Shapley SDs rise 0–7%. So the frozen cap's downward bias is real for CAWA
   but small: 1% of v(N), far inside its 5–95% interval. The cap binds on the high-density
   landscape, not on the footprint, where densities are low on both sides.
7. **First new-species batch: OSFL, GRSP, LEYE, BOBO** (started 2026-10-04). They are species at
   risk, chosen to differ from CAWA/OVEN: widespread boreal, grassland, wetland/northern,
   farmland. (The rest of the at-risk list once in 12C: BANS, BARS, EAWP, EVGR, GCTH, GWWA.)
   42 species × BCR tasks: OSFL 17, LEYE 11, BOBO 9, GRSP 5. They reach five BCRs CAWA/OVEN
   never ran through 12F: can3 (LEYE), can5 (OSFL, LEYE), can9 (OSFL, BOBO), can72 (OSFL).
   Withheld by BAM (14B drops them): GRSP can10 and can61, BOBO can10, LEYE can10.
   Species are now arguments (`bash 12A_observed.sh OSFL GRSP LEYE BOBO`, same for 12C and 12D;
   12D chains 12H). 12A runs on Fir (Fir's `07_predictions` has every species' 2020 tifs).
   Ranges for the four are staged on Fir. Steps:
   - [ ] 12A + 12C on Fir (independent; run together).
   - [ ] Smoke: `TEST_BCR=can3,can5,can9,can72,can13 TEST_N_BOOT=2 bash 12D_repredict_all_coalitions.sh OSFL GRSP LEYE BOBO`
     (9 tasks: every species and every new BCR).
   - [ ] Full 12D + 12H, pull, 14B. Watch for 12F's "weight.tif is all-zero/NA over the
     superset" stop (narrow ranges), and for tasks over 48G.

**Multi-year:** 07 and 11 take years (`bash 07_train_and_backfill.sh 2010 2015 2020`), but
each year needs `covariates_mosaiced_{year}.tif` (06; only 2020 exists) and that year's
footprint masks (`CanHF_1km_{lessthan1,morethan1}_{year}.tif`, or `HF_MASK_YEAR=2020` to borrow
2020's). 12A, 12C, 12D and 12H still fix `year <- 2020`.

**SSH is keyboard-interactive (2FA)**, so the user runs every cluster command from their own
session. Globus works: one file per call (or one filtered recursive transfer), **never**
`--batch`, always `--notify off`.

---

## Staging discipline: read before any cluster run

Stage the entry point's whole **`source()` closure**, including data files, never just what you
edited, and checksum-sync all of it (a 0-byte transfer proves equality for free). Unedited files
have been stale or missing on Fir four times: `08A` with Fix B, 7 of 9 files in 07's chain,
`12E` (never staged), and the old-schema `SpeciesPredictionTruncationValues.Rdata`.

- **Never restage an `.R` file while a job that runs it is active.** Rscript buffers its file
  and reads again at the end; Globus rewrites the file in place, so a running job can read the
  new file's tail as code (a 2026-10-01 comment edit to `11.R` would have left a stray
  ` | done.")` after every running task's `done.`). Restored in time; see the 11 entry.
- `.sh` files must be LF (a CRLF script dies with `/bin/bash^M: bad interpreter`). Check with
  `tr -dc '\r' < f | wc -c` (0 = LF), not `grep -c $'\r'`.
- Globus refuses sources outside the local endpoint's root, so stage from inside the repo tree.
- Before staging a 12F change, run `Rscript Rscripts/misc/verify_12f_can60_local.R` (~4 min) on
  Fir's real CAWA can60 inputs (`cluster_logs/localtest/`, `cluster_logs/smoke7/`, gitignored).
  Tables and arrays must be `identical()` to smoke 7 unless predictions or masks changed on
  purpose. Arrays may differ from Fir by ~3e-16 relative (Fir vs Windows floating point).
  The C2f variant is `cluster_logs/c2f_test/verify_c2f.R <all|none|real>`.

**Current state (2026-09-28):** `/Rscripts` on Fir matches `3fcef6d` for both closures:
- 07/11: `07` `.R`+`.sh`, `08A`, `08B_*`, `09_*`, `11` `.R`+`.sh`;
- 12D/12H: `12D` `.R`+`.sh`, `12E`, `12F`, `12G` `.R`+`.cpp`, `12H` `.R`+`.sh`.

`data/raw_data/hirshpearson/` on Fir holds the 2026-09-28 CanHF masks and the 16 sector layers.
`density_tables/` on Fir still holds C1's 2026-09-24 output; wipe it (Next action 3).

---

## C. Converge

- [ ] **C2a. Shapley SDs: fixed in code (`e747504`); needs the C1 rerun.** 14B used to
      propagate the tables' SDs as if subbasins and nested coalitions were independent. It now
      sums 12F's per-sample subbasin Shapley values (`[n_sub × 8 × 3200]`) within and across
      BCRs, sample by sample. Checked locally on CAWA can60: sample means equal the table
      Shapley values to 7e-13, and 14B's v(N) SD equals the arrays' joint SD.
- [x] **C2b. BART draws per bootstrap: kept as is** (user, 2026-09-25): 100 scenarios per
      bootstrap, each one of the 100 stored draws with replacement. Straddling subbasins get
      independent draws per BCR, which slightly understates their BART spread.
- [x] **C2c. Extrapolation flags: 10C rewritten as an area of applicability** (`686c16a`;
      Meyer & Pebesma 2021). The old KS + Mahalanobis rule flagged 667 of 667. The 0.5 cut is a
      judgement call; `frac_outside_aoa` is carried into 14B's subbasin table for any other cut.
- [ ] **C2d. Ecological check of the extreme BCR impacts (for the write-up).** Harnesses:
      `Rscripts/misc/diag_extreme_bcr_*.R`, inputs in `cluster_logs/sanity/`. All numbers are
      from the old masks and the pre-C2e backfill.
      - can11 (CAWA +227%, OVEN +275%): the footprint is 95% of the BCR, 93% of it crop/pasture,
        and low-HF land is 4.3%. The backfill turns cropland to grassland/shrub, which is right
        for prairie, but is more treed than its own training land. Percentages are large because
        the baseline is tiny. Driving subbasins (57, 61, 62) learn from ~3k pixels to backfill
        68k–112k; they are inside the AOA.
      - can13 (CAWA +301%, OVEN +89%): the footprint is 97%, low-HF land 0.4%. Cropland becomes
        mixed forest, consistent with the pre-settlement forest, but training is very thin.
      - OVEN negatives (can12 −8.5%, can80 −1.7%, can81 −6.3%) come from roads through the bird
        model's footprint response: can12's −828k was −1,041k direct and +213k vegetation. That
        response can be real (Mahon et al. 2019, Ecol. Appl. 29:e01895) but has thin survey
        support (CLAUDE.md Open Limitation #9).
      - Recheck these BCRs on the rerun's direct/vegetation split before reporting.
- [x] **C2e. Footprint covariates are 0 wherever the counterfactual removes footprint**
      (decided 2026-09-25, `8487982`, `bcfe665`; CLAUDE.md Open Limitation #9). 12F splits every
      Shapley value into direct (`d0 − obs`) and vegetation (`bf − d0`) parts; 08A backfills with
      the class at 0.
- [x] **C2f. Two masks per sector** (decided and built 2026-09-28, `3fcef6d`; CLAUDE.md Open
      Limitation #10). Footprint = `{sector}.tif > 0` and CanHF ≥ 1 (covariates zeroed on
      observed vegetation); direct = `{sector}_direct.tif == 1` (vegetation backfilled).
      - Masks: high-HF 1.80 M → 1.87 M px (97% of old kept); every subbasin keeps ≥ 911 low-HF
        px; the 674-subbasin set stays.
      - Footprint with HF ≥ 1, old → new (direct share): built 653k → 354k (100%), crop 677k →
        605k (100%), pasture 584k → 425k (100%), rail 153k → 119k (82%), roads 1,497k → 1,496k
        (85%), dams 65k → 85k (82%), mines 44k → 44k (5%), oil_gas 28k → 24k (8%); any sector
        1,569k → 1,598k (87%).
      - 12F harness (CAWA can60, 2 boots, pre-C2e backfill, so indicative): all-direct is
        `identical()` to smoke 7; none-direct gives a vegetation part of exactly 0; the real
        layers move v(N) −3,592 → −1,595 (roads −3,037 → −1,582).
- [x] **C3.** Sensitivity pass with the `q99.9` cap disabled (done 2026-10-04, Next action 6).
      Uncapped v(N): CAWA +1.16%, OVEN +0.05%; no sector share moves. Headline = conformed;
      sensitivity = uncapped.

## Loose ends

- [x] **Right-size 07's memory: done 2026-10-01 (`278bc85`), staged on Fir.**
      `07_train_and_backfill.sh` (run with `bash`; it absorbed the old submit helper on
      2026-10-01) submits two complementary arrays built from one list:
      subbasins 57, 61, 62, 98, 107 (in every year block) at 72G, everything else at 24G; 11
      waits on both. Tested with a stub `sbatch` for 1, 2 and 3 years: the two lists cover
      every index exactly once. Headroom: 1.25× (24G) and 1.36× (72G) over the 2026-09-28 peaks.
      If the `%30` cap starts binding on the next run, raise it.

- [ ] **Move `12A` to the cluster** (needed for 100+ species). Elly's
      `def-ecknight/NationalModels/output/` matches G: (06_bootstraps and 2020 07_predictions
      pair 1:1; 4,504/4,505 files match by name and size over 32 species; `q.out` covers the same
      151 species). Remaining: run `12A` on Fir for CAWA into a scratch dir and diff against the
      G:-built `observed_bootstraps.tif` / `truncation_params.rds`. Then set `cc <- TRUE`, add a
      `.sh`, stage `12B_v5_truncate.R` (`apply_masks = FALSE`, `project_to = NULL`) and drop the
      Globus step from CLAUDE.md.
- [ ] **Cheaper C1 for many species (optional):** 16 distinct BART draws per bootstrap instead of
      100 picks with replacement would cut a C1 from ~64 to ~25 billed core-equivalent-hours and
      moves the all-sector mean by 0.08% (CAWA) / 0.18% (OVEN) of the impact.
- [ ] Local files to delete once C1 is through (all gitignored):
      - C2f: `data/raw_data/hirshpearson/raw_300m/` (~4.6 GB; 03 and 14A need it only to
        rebuild), `data/raw_data/hirshpearson/old_bilinear/`,
        `cluster_logs/localtest/ia/data/raw_data/hirshpearson_c2f_{all,none,real}/`,
        `cluster_logs/c2f_test/`, `cluster_logs/extrapolation_flags_2026-09-25_bilinear_masks.csv`.
      - C2: `cluster_logs/sanity/` (~7 GB; read by `diag_extreme_bcr_*.R`),
        `cluster_logs/sub1_zeroed/`, `cluster_logs/extrapolation_flags_KS_mahal_2026-09-24.csv`,
        `cluster_logs/c1_old/`, `cluster_logs/c1_new/`.
- [ ] Decide the fate of local `covariates_mosaiced_2020_PREORIG.tif` (1.5 G), the only
      surviving pre-CAfire-fix 2020 mosaic.
- [ ] `15C_singletons_plot.R:150` reads `predictions_coalitions/`, deleted on both ends, so the
      map panel is broken. Regenerate via `save_arrays_ids` if wanted.
- [ ] Pre-2020 years need historical footprint layers, not the 2020 masks, and a per-year 07.

---

## Completed

Full detail: `git show 3fcef6d:TODO.md`.

- **C1 passed** 2026-09-24 (12D `61367890`, 12H `61367891`): 25/25 tasks, identity gate
  passed on all 32 bootstraps everywhere. Attempt 1 (`61344103`) exposed the fork GC OOM and
  the per-bootstrap categorical levels, both fixed. **First C2 run** the same day gave the
  Shapley means quoted in Next action 5; its SDs were wrong (fixed in C2a).
- **D. 12D rework**, 2026-09-24: NaN → `NA_real_` fix with the identity gate, the 12G C++ tree
  walk (3.6–4.4×), per-bootstrap reduction, one task per species × BCR at 8 cores / 48G / 3 h.
  Smokes 6–7 and 12H passed.
- **A7. Fix A validated**, 2026-09-23 (smoke 4: 100% of weight > 0 footprint pixels complete,
  after the `<cov>_mean` fill for draw-less subbasins).
- **A6. Backfill + premosaic**, 2026-09-16 → 09-22: 07 674/674 after the 08A assembly OOM fix
  (`73e8b68`); 11 19/19.
- **A. Cluster prep**, 2026-09-11 → 09-14: weight.tif rebuilt 25/25; the 07 smoke validated
  Fix B (100% complete draws on subbasin 1).
- **B. V5 packaging conformance**, 2026-09-11: `12B_v5_truncate.R` ports `10.Truncate.R`; gates
  G1–G4 pass. Re-run `Rscripts/misc/verify_weight_vs_v5_masking.R` after any weight change;
  the truncation harness is `verify_v5_truncate_port.R`.

## V5 reference

Renumbering in `f082866`: `10.Package.R` → `10.Truncate.R` + `11.Package.R`;
`11.Validate` → `12.Validate`; `12.Summarize` → `13.Summarize`;
`output/10_packaged/` → `output/11_packaged/`, new `output/10_truncated/`.
`06_bootstraps` and `07_predictions` were NOT renumbered, and they are the only V5 output
folders our scripts read.
