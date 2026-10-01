# TODO (last updated 2026-09-28)

> **Where we left off (2026-09-28).** C2f is decided, built and committed (`3fcef6d`; CLAUDE.md
> Open Limitation #10). The HF masks (03) and the 16 sector layers (14A) were rebuilt from the
> raw 300 m Hirsh-Pearson data, so 07 + 11 are rerunning for 2020 on Fir (submitted 2026-09-28).
> 10C is done on the new masks. Masks, sector layers and both closures are staged on Fir.
> The user wants the project finished: log small new problems under Loose ends and do not
> reopen finished stages for them. Keep `writing/` off the internet: never commit or push it.
>
> This file was trimmed on 2026-09-28. The full record of finished work (every smoke, job id
> and measurement) is `git show 3fcef6d:TODO.md`.

## Handoff to a new chat (written 2026-09-28)

**Open the new chat with:** "Continue from the Handoff section of TODO.md. Here is my
`sacct` output: <paste>". CLAUDE.md and the memory index load on their own.

**State at hand-off.** All code and CLAUDE.md are committed and pushed (`e290d00`). Nothing is running
locally. On Fir, 07 then 11 (2020, new masks) were submitted at 2026-09-28 19:36 UTC with
`bash 07_submit_backfill_years.sh 2020` (since 2026-10-01 merged into `bash 07_train_and_backfill.sh 2020`).

**What the new chat does not know:**
- **The new 07 and 11 job ids.** The submit helper printed them as
  `BACKFILL years=2020 07=<id> ... 11=<id>`. Otherwise run this on Fir and paste the output:
  `sacct -X -S 2026-09-28 --name=backfill,premosaic -o JobID%24,JobName%12 | awk '{split($1,a,"_"); print a[1], $2}' | sort -u`
- **The older run is still on Fir.** Don't mix it up with the new one. 07 `61546382` ran on
  the old masks; its 674 `.out` files are still in `Rscripts/`, with copies in
  `cluster_logs/07_2020_c2e/`. Its 11, `61546385`, was cancelled before it ran.
- **Two directories are untracked on purpose.** `data/derived_data/rds_files/` (98 MB) and
  `data/derived_data/sector_effects/` (the first C2 run's CSVs, which are stale and will be
  overwritten by 14B) are not committed. The user was asked whether to commit them and has not
  answered, so leave them out of git.

**Checking 07 and 11 without SSH (Next action 1).** In PowerShell, with
`$g = "C:\Users\mannf\AppData\Local\Python\pythoncore-3.14-64\Scripts\globus.exe"`,
`$fir = "8dec4129-9ab4-451d-a45f-5b4b8471f7a3:/home/mannfred/scratch/impact_assessment"` and
`$loc = "a7878ccc-747b-11ef-b4b8-8fef73a45f39:/C/Users/mannf/Drive/boreal_avian_modelling_project/ImpactAssessment/cluster_logs/07_2020_c2f"`,
each pull is ONE filtered recursive task. Write the filters as `--opt=value`, or click
glob-expands the `*` on Windows.
- `& $g transfer --recursive "--include=slurm-<07 id>_*.out" "--exclude=*" "$fir/Rscripts/" "$loc/" --label "07 c2f outs" --notify off`
- the same with `"--include=Y2020_S*.log"` from `"$fir/logs/"` into `"$loc/subbasin_logs/"`;
- the same with `"--include=slurm-<11 id>_*.out"` into `"$loc/11/"`.

Then:
- Wait with `& $g task wait <task id>`.
- Grep for the pass lines listed under Next action 1.
- Ask the user to run this on Fir; it must print 0:
  `find /home/mannfred/scratch/impact_assessment/data/derived_data/bart_models/2020 -name '*_confusion.rds' ! -newermt '2026-09-28 19:00 UTC' | wc -l`

On the old masks, 7 subbasins had a single backfill pixel (119, 160, 221, 246, 247, 304, 423).
08A skips BART for those, so their `.out` is tiny and their metrics file is 45 bytes. That is
expected, and the set may shift on the new masks.

**After that,** follow Next actions 3–6 in order. The user runs every cluster command. Give
each command as a single line.

## Critical path

```
cluster:  [07 + 11 rerun, 2020, new masks] ─► [wipe tables] ─► [C1 rerun: 12D + 12H] ─► [14B] ─► C3
local:    [10C rerun: done] ──────────────────────────────────────────────────────────┘
```

**Next actions, in order:**

1. **Check 07 + 11** (2020, new masks, submitted 2026-09-28).
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
   - **11 = job `61937626` (128G, 24 h): pending 2+ days on `(Priority)`; to be cancelled and
     resubmitted (2026-10-01).** A6's 11 (`60806856`, `--mem=512G`) peaked at 6.6–316 GB
     (`cluster_logs/07_2020_c2f/11/sacct_11_60806856.txt`): can3 (task 6) 316 GB / 17.4 h;
     can61, can80, can81 (12, 16, 17) 138–153 GB; can60 115 GB; all others ≤ 94 GB and ≤ 5.3 h
     wall. terra sizes its blocks from the node's `MemAvailable` (`/proc/meminfo`; checked in
     terra's `src/ram.cpp`) and cannot see the Slurm cap, so those peaks are what terra chose
     to hold, not what 11 needs, and 128G would likely have been killed on those four BCRs.
     `11.R` now caps terra at 30% of `SLURM_MEM_PER_NODE` (`memmax`; capped and uncapped runs
     were `identical()` locally for resample/cover/mask); `.sh` defaults 64G / 12 h, and
     `07_train_and_backfill.sh` submits can3 (task 6 of each year block) separately at 24 h.
     Resubmit: `--array=1-5,7-19 --time=12:00:00` and `--array=6 --time=24:00:00`. Check the
     first short tasks (can42 = 9, can9 = 19, can82 = 18) with `sacct` for MaxRSS well under 64G.
     **Resubmitted: the 12 h array is job `62431589`** (tasks 1–3 started 10:50 PDT 10-01; every
     `.out` logs `terra memmax = 19.2 GB` and processes in blocks). **Fir's `11.R` is the
     `278bc85` version on purpose** (12,154 bytes): HEAD differs only by a comment on line 190.
     Stage HEAD's `11.R` only after both 11 arrays have finished (see Staging discipline).
   - 11: every `.out` ends with `done.`. Here too `sacct` is not enough: three code paths exit 0
     without writing. If a task OOMs:
     `YEARS=2020 sbatch --array=<i> --mem=256G 11_premosaic_backfilled_stacks.sh`.
2. ~~10C~~, done 2026-09-28: 4 of 674 flagged (50, 72, 430, 481; was 72, 98, 481). Subbasin 98
   fell from 0.51 to 0.14 outside the AOA. `frac_outside_aoa` correlates 0.89 with the old run;
   median 1.6%, 90th percentile 14.3%. Log: `cluster_logs/10C_2026-09-28.log`.
3. **Wipe, then rerun C1** (~1 h wall, ~64 billed core-equivalent-hours):
   `cd /home/mannfred/scratch/impact_assessment/Rscripts && rm -f ../data/derived_data/density_tables/*.rds ../data/derived_data/density_tables/arrays/*.rds ../data/derived_data/density_tables/by_bcr/*.rds && sbatch 12D_repredict_all_coalitions.sh`
   then `sbatch --dependency=afterany:<12D job id> 12H_merge_bcr_tables.sh`.
   Pass: all 25 tasks `nice.`; the identity gate passes; each `superset pixels` line reports the
   direct footprint; 12H logs `wrote Shapley samples` for both species.
4. **Pull** the 510 tables, 18 arrays and the two `*_shapley_samples.rds` (`--notify off`). The
   masks changed, so `obs_on_coalition` changes too; there is no `identical()` gate.
5. **14B** locally (C2 rerun). For comparison, the first C2 run (2026-09-24, old masks, pre-C2e,
   CAWA can40 withheld) gave CAWA v(N) = +574k on 4.47 M observed (roads 34%, pasture 22%,
   crop 21%, built 18%) and OVEN v(N) = +3.53 M on 37.65 M (crop 33%, pasture 32%, built 21%,
   roads 11%).
6. **C3:** the uncapped `q99.9` sensitivity pass.

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
- [ ] **C3.** Sensitivity pass with the `q99.9` cap disabled; report the spread. The frozen cap
      under-estimates impact in the highest-density pixels (counterfactual densities hit it more
      often than observed), and it removes 4–12% of abundance. Headline = conformed;
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
