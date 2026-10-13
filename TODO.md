# TODO (last updated 2026-10-06)

> **Where we left off (2026-10-12): the C2d rerun is under way.** The C2d fix (blank model, no
> water, drawn classes) passed its 14-subbasin tests, and the user approved the full rerun.
> 07 for all 674 subbasins was submitted 2026-10-12 ~22:20 (11 chained), with `BART_MODELS_DIR`
> unset; it overwrites production's backfills. 12F (blanks, drawn classes) was staged the
> same evening. The first finished backfill (subbasin 6) checks out: 400 drawn-class layers, no
> water class drawn, and its four all-blank covariates are blank on all its training land too.
> 07 finished 2026-10-13 ~05:40 (Fir time): 664 of 674 tasks ok. Tasks 23–30 failed to load
> `dplyr` (transient), and subbasins 247 and 423 (one backfilled pixel each) failed in the new
> blank model, since BART drops a one-row `x.test` to a vector (fixed in `08B_deploy_pbart.R`:
> P(blank) = the training share, as gbart/mbart take the mean; staged 2026-10-13). 11 rebuilt 15
> BCRs; can11, can12, can13 and can61 timed out at 4 h, each stuck pre-resampling one subbasin
> (62, 98, 107). Those stacks now exceed what terra holds in memory at 32G (memmax 9.6 GB:
> 6.6–14.6 GB as doubles, up from 5.5 GB before the class and blank layers), and terra's block
> path stalls on ~2400-layer pixel-interleaved stacks. Resubmitted: 07 for 23–30, 247, 423, then
> 11 for all 19 at 128G / 8 h (it skips mosaics newer than their backfills).
> Next:
> 1. Run a 12D smoke on can11 (`TEST_BCR=can11 TEST_N_BOOT=2`), since 12F's new paths have
>    only run through 16C's mirror.
> 2. Wipe the six species' tables, then run 12D + 12H for CAWA, OVEN, OSFL, LEYE, GRSP, BOBO.
> 3. Run 14B, 14C and 18 locally.
> 4. Redo the C3 (`DROP_Q99`) and C4 (`DENS_MIN`) sensitivity runs if they are to be reported.
> 5. The PDF report is the very last step.
> The Results table below is pre-C2d.
> The user wants the project finished: log small new problems under Loose ends and do not
> reopen finished stages for them. Keep `writing/` off the internet: never commit or push it.
>
> Committed and pushed through `e619087` (2026-10-04), and Fir's `Rscripts/` matches that
> commit (12D closure checksum-synced 2026-10-05; retired files removed 2026-10-06). The C4
> `DENS_MIN` edits to 12D/12E/12F/12H/14B (2026-10-06), C5's 12F pixel output (staged on Fir
> 2026-10-07) and the new script 18 are uncommitted. `rds_files/`,
> `sector_effects/` and `sector_effects_noq99/` stay local, never committed (user, 2026-10-04).
> SSH is 2FA, so the user runs every cluster command; Globus works without it.
>
> Full records of finished work: `git show 3fcef6d:TODO.md` (to 2026-09-28) and
> `git show e619087:TODO.md` (to 2026-10-04: the 07/11 checks, C1, C3 run details).

## Results (2020, BAM release filter applied)

Total impact v(N) = counterfactual − observed: positive means fewer birds with industry.
Per-sector values: `data/derived_data/sector_effects/shapley_{national,bcr,subbasin}.csv`.

| Species | Observed | v(N) (5–95%) | % of observed | Largest sectors |
|---|---|---|---|---|
| CAWA | 4.47 M | +757k (649k, 870k) | +16.9% | roads, pasture, crop, built |
| OVEN | 37.65 M | +4.59 M (3.70 M, 5.18 M) | +12.2% | crop, pasture, roads, built |
| OSFL | 10.74 M | +417k (288k, 582k) | +3.9% | roads, crop, pasture |
| LEYE | 5.79 M | +38k (21k, 57k) | +0.7% | roads |
| GRSP | 0.96 M | −463k (−695k, −173k) | −48% | roads, crop, pasture |
| BOBO | 6.76 M | −3.42 M (−3.63 M, −3.26 M) | −51% | roads, crop, pasture |

- **Withheld by BAM** (14B drops them): CAWA can40; GRSP can10, can61; BOBO can10; LEYE can10.
- **Direct vs vegetation** (Open Limitation #9): CAWA direct +138k, vegetation +619k. OVEN
  direct −1.68 M, vegetation +6.27 M; roads is −1.27 M direct, +2.14 M vegetation, hence its
  wide interval.
- **C3, without the frozen q99 cap:** v(N) +1.16% (CAWA) and +0.05% (OVEN); no sector share moves.
  Headline = capped; sensitivity = uncapped.
- **C4, densities < 0.005 birds/ha set to 0** (`sector_effects_dmin0.005/`): CAWA v(N) −2.7%
  (757k → 736k) on observed −6.0% (4.47 M → 4.20 M), so +16.9% → +17.5% of observed; OVEN
  v(N) −0.5%, observed −0.6%. Sector shares move < 1 point. By BCR it is not a uniform cut:
  CAWA can11 58.5k → 10.9k (−81%, the prairie floor), but can13 228k → 247k and can61
  258k → 267k. There the floor removes thin observed density (can13 obs 59k → 41k) while the
  backfilled forest stays above 0.005, so the impact grows by what the observed side lost.
  OVEN can11 697k → 673k (−3.5%). Undecided whether it becomes the headline.
- **10C:** 4 of 674 subbasins flagged (50, 72, 430, 481); they carry 23k (CAWA) and −3k (OVEN).

Runs and logs:
- 07 `61937625` + 11 `62431589` (can41/can5 `62468313`): `cluster_logs/07_2020_c2f/`.
- 10C: `cluster_logs/10C_2026-09-28.log`.
- C1 (CAWA, OVEN): 12D `62596963`, 12H `62597111`, `cluster_logs/C1_2026-10-02_c2f/`;
  14B `logs/14B_2026-10-02_c2f.log`.
- C3: 12D `62636072`, 12H `62636104`, `cluster_logs/C3_2026-10-02_noq99/`;
  14B `logs/14B_noq99_2026-10-04.log`.
- C4: 12D `63402677`, 12H `63402678`, `cluster_logs/C4_2026-10-06_dmin/`; 14B
  `logs/14B_dmin0.005_2026-10-07.log`. 25/25 `nice.`, identity gate 25/25, complete weight > 0
  footprint ≥ 99.9% everywhere; 14B means match the tables (max |diff| 4.1e-10), additivity exact.
- C5 (six species, production, pixel_change): 12D `63502004`, 12H `63502005`,
  `cluster_logs/C5_2026-10-07_pixel/`. 67/67 `nice.` and pixel rasters; all 1,536 tables
  byte-identical to the pre-run copies, so 14B was not rerun.
- New species:
  - 12A `62902543`, 12C `62902545`: `cluster_logs/newspp_2026-10-04/`.
  - Smoke 12D `63130780`.
  - Full run 12D `63133307`, 12H `63133308`: `full/` in the same folder.
  - 14B: `logs/14B_newspp_2026-10-06.log`.
  - All 42 tasks passed the identity gate. The longest took 44 min (can11), and 100% of
    weight > 0 footprint pixels were complete everywhere except LEYE can3 (see Loose ends).

## Adding a species batch

Remaining species-at-risk candidates: BANS, BARS, EAWP, EVGR, GCTH, GWWA. Check BAM's "remove"
tab for their withheld BCRs.

1. Stage `data/raw_data/v5_gis/ranges/{SP}.tif` (from G: `gis/ranges/`) on Fir.
2. Run `bash 12A_observed.sh SP1 SP2` and `bash 12C_build_prediction_weights.sh SP1 SP2`. They
   can run together.
3. Smoke-test any BCR no earlier species ran through 12F:
   `TEST_BCR=<bcrs> TEST_N_BOOT=2 bash 12D_repredict_all_coalitions.sh SP1 SP2`.
4. Full run: `bash 12D_repredict_all_coalitions.sh SP1 SP2` (12H is chained). Check every task
   for `nice.` and the identity gate.
5. Pull with one filtered recursive transfer of `density_tables/`
   (`--include=SP1_*`, one per species, then `--exclude=*`). It also descends into `arrays/`
   and `by_bcr/`.
6. Run 14B locally. It reads every species in `density_tables/`.

**All remaining species (started 2026-10-08):** `Rscripts/run_species_batches.sh` drives steps
2–4 for `Rscripts/species_todo.txt`: 140 species, 1,557 species × BCR tasks (Fir listing:
`cluster_logs/species_listing/spp_counts.txt`). Left out: DEJU, DOWO, EUST, FOSP (bootstrap
models but no 2020 `07_predictions` on Fir) and ROSA (no Canadian models).
`bash run_species_batches.sh prep` (12A + 12C, 140 each), then `smoke` (can42, the one BCR 12F
has never run: 3 species, 2 bootstraps), `plan`, and `wave 1` … `wave 4` (≤ 400 tasks each, one
at a time; `MAX=` resizes). Tested locally against a mock of Fir's `06_bootstraps`.
- [x] Range rasters for the 140 copied from G: `gis/ranges/{SP}.tif` by exact name (2026-10-08,
      after a Google Drive restart: it had stopped reading files when the disk filled). All open,
      same CRS and 10 km grid as OSFL's, values 0.13–1. Staged to Fir with
      `run_species_batches.sh` and `species_todo.txt`.
- [x] `prep` (2026-10-08): 12A `63740250`, 12C `63740251`, logs in
      `cluster_logs/allspp_2026-10-09/prep/`. 140/140 complete each; 1,557 observed stacks and
      1,557 weights (the planned task count), nothing skipped. Lowest mean weight RBGR can41
      (0.001); densmax/q99 median 4.4×, max 4,243×. Globus CLI re-login: run `globus session
      update computecanada.ca --no-local-server` in a PowerShell window (not via `!`).
- [x] `smoke` (2026-10-09): 12D `63884817`, 12H `63884818`, logs in
      `cluster_logs/allspp_2026-10-09/smoke/`. ALFL, AMRO, ATSP can42: 3/3 `nice.`, identity
      gate passes, 12H merged all three. can42 carries only 9 weight > 0 footprint pixels (all
      complete) of its 2,438 superset pixels, so it adds next to nothing.
- [x] `wave 1` (2026-10-09): 12D `63888611`, 12H `63888612` (stopped), logs in
      `cluster_logs/allspp_2026-10-09/wave1/`. 387/389 `nice.`, nothing skipped; can3 is 82.7%
      complete for every species (the LEYE can3 loose end), every other BCR ≥ 99.9%. BAWW and
      BBWA can41 stopped on "weight.tif is all-zero/NA over the superset": their ranges miss
      can41's footprint. Fixed in 12F (stops only if the weight is 0 over the whole BCR; a
      zero-weight footprint runs on with impact 0 and keeps obs_total). That changes 12F's md5,
      so 12H for the other 33 species merges the old-version files, and BAWW + BBWA rerun in
      full (26 tasks) on the new 12F. 12H raised from 8G / 1 h to 32G / 3 h: it holds a whole
      wave's per-BCR files at once (~1.6 GB of rds for wave 1). MaxRSS 7.20 GB for 363 files,
      so 8G would have been at the edge; keep 32G.
      Done: 12H `63924302` merged the 33 (363 results, 3 min); 12D `63924303` reran BAWW + BBWA
      (26/26 `nice.`; can41 logged "weight is 0 at every footprint pixel", predicted 0 of 4,228
      pixels, pixel_change sums to the table) and 12H `63924304` merged them. 35/35 species
      have 255 tables, 9 array files and Shapley samples. Logs in `.../wave1fix/`.
- [x] `wave 2` (2026-10-09): 12D `63948940`, 12H `63948941`, logs in `.../wave2/`. 394/394
      `nice.`, 12H `nice.` (39 species), nothing skipped, no zero-weight BCR; can3 82.7% again,
      every other BCR ≥ 99.9%. All 80 species on Fir (74 + the 6 earlier) have 255 tables and
      Shapley samples.
- [ ] **On hold: 07 must be re-run (backfill defect found 2026-10-12)**, so every density
      product is stale and waves 3–4 wait for 07 → 11. Kept for the rerun: per-task runtimes
      (`cluster_logs/allspp_2026-10-09/task_timing.tsv`, 807 tasks: median 5.7 min, max 66 min
      AMRO can11, 128 task-h; each wave took ~1.5 h from first start to last end) and the
      old-backfill baseline per species × BCR (`.../baseline_old_backfill_species_bcr.csv`,
      from the 80 cid-256 tables in `.../c256/`) to measure how far the fix moves each species.
      12A/12C outputs do not depend on 07 and stay; `prep` and `smoke` need no rerun.
      Deleted 2026-10-12: every file in Fir's `density_tables/` (22,900, all 80 species), in
      `density_tables_noq99/` and `density_tables_dmin0.005/` (1,110; folders kept) and
      the local `sector_effects*/` CSVs. `species_todo.txt` now holds all 146 species (the 6
      earlier ones added), 1,624 tasks, for the rerun.
      Resources (`sacct_12D_allspp.txt` → `task_resources.csv`, `sacct_summary.R`): 12D MaxRSS
      median 16.9 GB, max 44.8 GB (BAWW can11; can11 median 37.4, can61 max 41.6, every other
      BCR ≤ 33 GB), so 48G is needed and can11 runs at 93% of it: a task the fixed backfill
      pushes over is resubmitted at 64G (12H names it). CPU efficiency 0.75 overall (0.1–0.4
      on the small BCRs, where single-threaded loading dominates). Queue wait median 18–22 min,
      max 40 min; each ~390-task wave was done 1.5 h after submission. Queue time is not the
      bottleneck, the gaps between waves are: QOS `normal` sets no MaxSubmitPU/MaxJobsPU, so rerun
      with `MAX=900` (2 waves; under Slurm's default MaxArraySize of 1001). 12H MaxRSS 10.5 GB for 394 files
      (~0.027 GB/file), so 32G holds a wave of up to ~1,000 files.
- [ ] Pulling all 140 needs ~10–12 GB locally (Shapley samples ~50–70 MB per species).

## Globus from this laptop

- Write filters as `--opt=value` (`"--include=slurm-<id>_*.out" "--exclude=*"`): a bare `*`
  is glob-expanded on Windows. Filters do not keep subdirectories out of a recursive pull.
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

## Staging discipline: read before any cluster run

Stage the entry point's whole **`source()` closure**, including data files, never just what you
edited, and checksum-sync all of it (a 0-byte transfer proves equality for free). Unedited files
have been stale or missing on Fir four times: `08A` with Fix B, 7 of 9 files in 07's chain,
`12E` (never staged), and the old-schema `SpeciesPredictionTruncationValues.Rdata`.

- **Never restage an `.R` file while a job that runs it is active.** Rscript buffers its file
  and reads again at the end; Globus rewrites the file in place, so a running job can read the
  new file's tail as code.
- `.sh` files must be LF (a CRLF script dies with `/bin/bash^M: bad interpreter`). Check with
  `tr -dc '\r' < f | wc -c` (0 = LF), not `grep -c $'\r'`.
- Globus refuses sources outside the local endpoint's root, so stage from inside the repo tree.
- Before staging a 12F change, run `Rscript Rscripts/misc/verify_12f_can60_local.R` (~4 min) on
  Fir's real CAWA can60 inputs (`cluster_logs/localtest/`, `cluster_logs/smoke7/`, gitignored).
  Tables and arrays must be `identical()` to smoke 7 unless predictions or masks changed on
  purpose. Arrays may differ from Fir by ~3e-16 relative (Fir vs Windows floating point).
- Closures: 07/11 = `07` `.R`+`.sh`, `08A`, `08B_*`, `09_*`, `11` `.R`+`.sh`.
  12D/12H = `12D` `.R`+`.sh`, `12E`, `12F`, `12G` `.R`+`.cpp`, `12H` `.R`+`.sh`.
  12A = `12A` `.R`+`.sh`, `12B`. 12C = `12C` `.R`+`.sh`.

---

## Open before the write-up

- [ ] **C2d. Ecological check of the farmland BCRs (can11–13).** Harnesses:
      `Rscripts/misc/diag_extreme_bcr_*.R`, inputs in `cluster_logs/sanity/`. The CAWA/OVEN
      numbers below are from the old masks and the backfill made before C2e.
      - can11 (CAWA +227%, OVEN +275%): the footprint is 95% of the BCR, 93% of it crop/pasture,
        and low-HF land is 4.3%. The backfill turns cropland to grassland/shrub, which is right
        for prairie, but it is more treed than its own training land. Subbasins 57, 61 and 62
        learn from ~3k pixels to backfill 68k–112k; they are inside the AOA.
      - can13 (CAWA +301%, OVEN +89%): the footprint is 97%, low-HF land 0.4%. Cropland becomes
        mixed forest, consistent with the pre-settlement forest, but training is very thin.
      - OVEN negatives (can12, can80, can81) come from roads through the bird model's footprint
        response (Open Limitation #9).
      - New species (2026-10-06):
        - OSFL, a forest species, gains more birds than it has today in can13 (+34k on 2.9k
          observed) and can11 (+16k on 12k). can11–13 carry 28% of its national +417k. This
          looks like CAWA's prairie density floor: a low BRT density on backfilled forest × a
          large footprint.
        - **GRSP and BOBO (16C, 2026-10-09; report `writing/grassland_decline_report_2026-10-09.pdf`).**
          `16C_grassland_decline_diagnostics.R` on the subbasins carrying 82–99% of each BCR's change
          (GRSP can11; BOBO can11, can13), with 12F's models and design rules (identity and cid-256
          checks pass everywhere). Output `rds_files/grassland_diagnostics/`, log
          `logs/16C_2026-10-09.log`. Share of today's birds kept on the direct footprint:

          | | BART (12F) | kNN, all training pixels | in / out of AOA |
          |---|---|---|---|
          | GRSP can11 | 0.37 | 1.41 | 0.37 / 0.43 |
          | BOBO can11 | 0.63 | 0.65 | 0.64 / 0.63 |
          | BOBO can13 | 0.13 | 0.23 | 0.13 / 0.13 |

          - **GRSP can11 is a BART artefact.** `SCANFIheightcv_*` is NA wherever SCANFI height is 0
            (no trees): 91–98% of the training land in subbasins 56, 63 and 60. V5 sends NA down
            gbm's missing branch (open land). BART never returns NA or an exact 0, so every
            backfilled pixel gets a canopy CV and a little broadleaf cover; GRSP's partial dependence
            on `SCANFIprcD_5x5` falls to a third between 0.03% and 0.3%. One at a time: prcD_5x5
            −382k, prcC_5x5 −212k, heightcv_1km −179k (relative influence 0.1%); all vegetation −365k.
            Quantile-mapping BART onto the training values (restores zeros, not NAs) gives 0.92.
            Production GRSP can11 (−440k, −49%) should not be reported as it stands.
          - BOBO can11 is what the training land implies (V5: 6.8 vs 14.5 birds/km² on training land
            vs footprint in subbasin 57). So is BOBO can13 (0.5–3 vs 19–28; the training land is can12
            Shield forest), but BART overstates it by ~11% (loss 1.97M vs kNN 1.76M on the studied
            direct footprint), mostly through 08B's modal class (Loose ends).
          - Not extrapolation in any case: the same inside and outside the AOA; ≤ 0.1% of backfilled
            values outside `var.levels`. T4: GRSP can11 subbasins with more footprint per training
            pixel lose *less* (Spearman 0.59); BOBO can13's seven large subbasins lose 85–90% across
            ratios 0.95–34 (its ρ −0.62 comes from small subbasins with ample training land).
          - Probably wider: any species whose model splits near zero canopy or on a missing heightcv,
            wherever the training land is treeless (much of can11). For forest species it would
            inflate the counterfactual (CAWA, OVEN, OSFL can11 gains). Untested.
          - ~~kNN hot-deck backfill~~ (the report's proposal): **rejected by the user 2026-10-12**,
            "a new method on top of BART that is a bit arbitrary"; make BART itself robust.
          - **Cause, found 2026-10-12: 08A's `_isNA` flag, not BART's zeros.** BART never returns NA,
            so a backfilled heightcv is never blank and its `_isNA` flag is 0 on every backfilled
            pixel, while it is 1 on ~95% of the training land. prcC/prcD (after heightcv in the
            hierarchy) then read every backfilled pixel as one of the few treed training pixels.
            Refit test (`Rscripts/misc/diag_isna_flag_refit.R`, subbasin 56; replays production's
            draws to 2.4e-7): modelling where heightcv is blank (probit BART, P > 0.5) gives it 97%
            blank (training 95%) and moves prcD_5x5 from median 1.08% / 0.2% zeros to 0.005% / 30%
            (training 0.001% / 45%). Height, backfilled before the flag, already has 46% zeros.
          - **Fix (2026-10-12), 07 rerun + 12D for all species:**
            (a) `08B_deploy_pbart.R`: a probit BART models where each biotic covariate is blank on
            its training land; 08A blanks the backfill where P > 0.5 (cascade) and writes
            `<cov>_blank` / `<cov>_pblank`; 12F sends blanks to gbm as NA and counts them complete.
            A covariate blank on ALL the training land is blank everywhere (fixes the can3 Loose end).
            (b) 08A leaves each land-cover layer's water class (SCANFI 8, VLCE 20, MODIS 17) out of
            its training rows: the backfill never predicts water (user: reviewers). Production
            labelled 94% of can13 subbasin 72's direct footprint water (SCANFI_1km; training 7%).
            (c) 08B writes 100 drawn classes `<cov>_class_###` (each from one posterior draw's class
            probabilities); 12F picks them with the scenario's draw index (fixes the modal-class
            Loose end). (d) 08B renumbers the training split's classes 1..K before `mbart2`, which
            fits `y == h` for h in 1..K (a holdout-only class shifted every class above it).
            (e) `prep_bart_design()` keeps a predictor that is blank on the whole backfill if it is
            partly blank on the training land (its flag carries the information).
            12F on an old mosaic is `identical()` to Fir's production 12F (can60, 255/255 tables,
            9/9 arrays, Shapley samples; `Rscripts/misc/verify_12f_can60_local.R`'s smoke-7 baseline
            predates C2f, so the comparison is old vs new 12F on the same inputs).
          - **14-subbasin test 1, 2026-10-12** (jobs 64595258/9, `BART_MODELS_DIR=
            bart_models_blank_test`; backfills in `cluster_logs/grassland_diag/bart_models_blank_test/`;
            `BF_NEW_DIR=bart_models_blank_test Rscript Rscripts/16C_...`, log
            `logs/16C_2026-10-12_blanktest.log`, output `grassland_diagnostics/bart_models_blank_test/`).
            Share of observed birds kept on the direct footprint (16C's 8-bootstrap × 10-draw check,
            which reads production ~10% low for GRSP):

            | | production BART | new BART | kNN |
            |---|---|---|---|
            | GRSP can11 | 0.34 | 0.97 | 1.41 |
            | BOBO can11 | 0.63 | 0.60 | 0.65 |
            | BOBO can13 | 0.13 | 0.29 | 0.23 |

            GRSP's artefact is mostly gone (heightcv blank shares now match kNN's in 5 of 7
            subbasins; prcD_5x5 back to the training land's level in 56, 63; still above kNN in 55,
            60). BOBO can11 unchanged, as the earlier tests implied. No backfilled NA other than
            modelled blanks. BOBO can13 went ABOVE kNN (subbasin 98: 0.37 vs 0.23) because the
            drawn classes over-drew rare classes (herbs 6.2%, shrub 10.3% vs training 1.1, 2.2%):
            **mbart2's `prob.test` is a softmax of PROBIT scores**, exp(z)/Σexp(z), which flattens
            the probabilities. On held-out training pixels it roughly doubles rare classes (98: herbs
            3.1% observed, 6.1% softmax, 3.2% with pnorm(z)/Σpnorm(z)); normalised probit has the
            lower log loss in subbasins 98, 72, 56. 08B now draws classes from normalised pnorm(z)
            (the modal class is unaffected). Test 1 ran without it and without the all-blank rule.
          - **14-subbasin test 2, 2026-10-12** (both changes; jobs 64633705/6, `BART_MODELS_DIR=
            bart_models_blank_test2`; backfills in `cluster_logs/grassland_diag/bart_models_blank_test2/`,
            log `logs/16C_2026-10-12_blanktest2.log`, output `grassland_diagnostics/bart_models_blank_test2/`).
            Same seeds as test 1, so only the class conversion differs. Production → test 1 → test 2
            (kNN): GRSP can11 0.34 → 0.97 → 0.99 (1.41); BOBO can11 0.63 → 0.60 → 0.61 (0.65);
            **BOBO can13 0.13 → 0.29 → 0.28 (0.23)**. The softmax explained only a fifth of the can13
            gap. All of the rest is subbasin 98 (0.36 vs kNN 0.23; +177k birds; 92 goes the other way,
            0.23 vs 0.40). There SCANFI herbs are 10.6% (test 1 12.5%, kNN 6.1%, training 2.8%), shrub
            18.8% (kNN 15.0, training 3.6); VLCE and MODIS classes are close to kNN. On 98's held-out
            training pixels the classes are calibrated (herbs 3.1 vs 3.2%), so the footprint's excess
            is BART extending abiotic trends beyond its training land: 98 trains on can12 land and 59%
            of can13's observed birds are outside the AOA, where kNN copies the nearest training pixel
            instead. Not a coding defect; Loose ends. No backfilled NA other than modelled blanks
            (`new_frac_na_not_blank` 0 everywhere); no water. 11 resamples every layer with "near", so
            the class layers mosaic as they are.
          - Smoke on subbasin 72 locally: no water in 875,500 drawn classes; SCANFI_1km mixed 53 /
            broadleaf 27 / conifer 15%, like its training land. Fir and local agree exactly on every
            continuous draw; mbart's probit chains differ by platform (VLCE, SCANFI), so land-cover
            replays only match on Fir.
          - The 2026-10-08 run took complete-case donors only, which dropped every treeless training
            pixel and read can11 as wooded (GRSP kNN 0.20). Its outputs are kept in
            `grassland_diagnostics/run_2026-10-08_completecase/`; don't reuse them.
- [x] **C4. Low-density floor sensitivity run (`DENS_MIN=0.005`).** Done 2026-10-07; results
      under Results, runs under Runs and logs. Colleagues (2026-10-05):
      there should be no CAWA in the prairies. V5's CAWA range raster is 1.0 over all of can11,
      so the weight cannot remove them, and can11's CAWA is a ~0.002 birds/ha BRT floor that the
      backfill spreads over 440k km² of footprint (+58k). They proposed zeroing densities below
      0.005 birds/ha (confirmed /ha, not /km²). Estimated on the 2020 observed means: −6.8% of
      CAWA's observed total (can10 96%, can11 82%, can71 73%, can13 32%), −0.9% of OVEN's.
      Built 2026-10-06: 12F floors capped density per pixel and sample, before the weight, on
      obs, d0 and bf; 12D/12H/14B use `density_tables_dmin0.005/`, `sector_effects_dmin0.005/`.
      Local can60 harness: production path `identical()` to the C2f run (255/255 tables, 9/9
      arrays, Shapley samples); floored obs totals match an independent recompute (1e-15).
      Run: `DENS_MIN=0.005 bash 12D_repredict_all_coalitions.sh CAWA OVEN`, pull, then
      `DENS_MIN=0.005 Rscript Rscripts/14B_sector_attribution.R`. Compare with production.
      Not committed yet.
- [x] **C5. Pixel-level change maps** (user, 2026-10-07: subbasin choropleths jump at watershed
      lines). 12F now writes `pixel_change/{species}_2020_{bcr}.tif` (mean obs and cf of the
      all-sector coalition at each footprint pixel, birds/km² weighted) after checking they sum
      to the cid-256 table; script 18 sums obs − cf into 10 km cells over the species' modelled
      Canadian area → `output_figures/change_grid_*`. Local can60 harness: tables, arrays and
      Shapley samples `identical()` to the C2f run; pixels sum to the table (1,594.5).
      Done 2026-10-07 for the six species, production only (user dropped the `DENS_MIN=0.005`
      rerun); run under Runs and logs. Every species' 10 km cells add up to 14B's national
      change within 0.5 birds. 18 counts a cell's area as its 1 km pixels touching the modelled
      BCRs (12C's rule; `cover = TRUE` dropped 4 BOBO birds at BCR edges) and checks against
      the national CSV (14B rounds each subbasin to a bird, so their sum drifts).
      Script 18 now also clips subbasins to Canada (41 clipped; areas are Canadian only).
- [x] **Results tables** (2026-10-07): `14C_results_tables.R` →
      `sector_effects/results_summary.csv` and `results_sectors.csv` (sectors as % of observed).
      Its means match 14B's national CSV for all six species after the C5 re-run.
      2026-10-08, at a colleague's request (positive = loss confused them): 14C and the 18 map
      titles now report change = observed − cf (negative = birds lost to industry) as % of the
      no-industry population (CAWA −14.5%, was +16.9% of observed). Flipping the sign alone would
      have left GRSP's near doubling (+94%) reading as +48%. Report: `writing/results_report.typ`.

## Loose ends

- [ ] **08B's modal land-cover class distorts class shares (16C, BOBO can13).** mbart's per-pixel
      `max.col` class, used for every draw, erases minority classes. On the studied can13 farmland,
      SCANFI_1km after BART: mixed forest 65% (kNN 29%, training land 38%), broadleaf 21% (kNN 43%),
      herbs 0% (kNN 5%), shrub 3% (kNN 11%). It carries most of BOBO can13's BART − kNN gap (~210k
      birds); forest species (OVEN, CAWA) could feel it more. Fixed in code 2026-10-12 (drawn
      classes, C2d fix (c)); takes effect with the 07 rerun.
- [ ] 16C's BART reconstruction (8 bootstraps × 10 draws) is 10% below 12F's cf for GRSP (27%
      in subbasins 56/63/60; BOBO 1.00): GRSP's prcD draws straddle the 0.03% split, so few
      draws are noisy. T1 uses 12F's own cf, so only T2/T5's magnitudes are approximate.
- [ ] **can3 (Arctic): VLCE_1km has no backfill in four subbasins.** In LEYE can3, 126 of the
      730 weight > 0 footprint pixels were dropped, all from subbasins 8060193680 (90),
      8060011800 (21), 8060018840 (9) and 8060189340 (3). Likely cause: VLCE maps forested
      ecozones only, so these tundra subbasins have no valid training rows and 08A writes VLCE
      all-NA. If observed VLCE is NA there too, V5 predicted those pixels down gbm's missing
      branch, and 12F could pass the backfilled VLCE as NA instead of dropping the pixel. The
      effect is small for LEYE; check it before running Arctic species. Fixed in code
      2026-10-12: a covariate blank on all the training land is a modelled blank everywhere
      (C2d fix (a)); takes effect with the 07 rerun.
- [ ] **After the C2d fix, BART and kNN still differ where training land is thin or far (test 2).**
      BOBO can13 subbasin 98 keeps 0.36 of its observed birds vs kNN 0.23 (+177k; more SCANFI
      herbs and shrub than kNN; trains on can12 land, 59% of can13's birds outside the AOA).
      GRSP can11 subbasins 55 and 60 keep 0.73 and 0.55 vs kNN 1.12 and 1.18 (prcD_5x5 above kNN;
      likely the Gaussian amount model's missing point mass at 0). Neither is a coding defect, and
      both methods agree on the direction. Report it as method sensitivity, with 10C's AOA flags.
- [ ] **11 needs 128G for can11, can12, can13 and can61 since the C2d backfills** (~2400 layers
      per subbasin). At the .sh's 32G, subbasins 62, 98 and 107 fall to terra's block path, which
      stalled > 4 h on one subbasin (2026-10-13). The .sh still says 32G / 4 h; a later rerun
      passes `--mem=128G --time=08:00:00` (or 07's .sh could pass it when it chains 11).
- [ ] 12A on Fir (optional check): run it for CAWA into a scratch dir and diff against the
      G:-built `observed_bootstraps.tif` / `truncation_params.rds`.
- [ ] **Cheaper 12D for many species (optional):** 16 distinct BART draws per bootstrap instead of
      100 picks with replacement would cut a CAWA + OVEN run from ~64 to ~25 billed
      core-equivalent-hours. It moves the all-sector mean by 0.08% (CAWA) / 0.18% (OVEN) of the
      impact.
- [ ] **The laptop disk is small (236 GB) and filled to 0 bytes on 2026-10-08**, which failed the
      16C backfill pulls (`NO_SPACE_LEFT`). Deleted that day (user): `hirshpearson/raw_300m/`
      (re-download from Borealis doi:10.5683/SP2/EVKAVL + `cum_threat` from G: to rebuild 03/14A),
      `covariates_mosaiced_2020_PREORIG.tif`, `old_bilinear/`, the localtest
      `hirshpearson_c2f_*` layers, `c2f_test/`, `sub1_zeroed/`, `c1_old/`, `c1_new/` and the two
      old extrapolation-flag CSVs (~6 GB). Earlier (2026-10-07): local `by_bcr/` folders, two
      pre-v3 `weight.tif`, the 16A/17 outputs. Still deletable once C2d is done:
      `cluster_logs/sanity/` (6.7 GB). Check `Get-PSDrive C` before any large pull.
- [ ] `15C_singletons_plot.R:150` reads `predictions_coalitions/`, deleted on both ends, so the
      map panel is broken. Regenerate via `save_arrays_ids` if wanted.
- [ ] **Pre-2020 years** need historical footprint layers (not the 2020 masks), each year's
      `covariates_mosaiced_{year}.tif` (06; only 2020 exists) and a per-year 07. 07 and 11 take
      years (`bash 07_train_and_backfill.sh 2015 2020`); 12A, 12C, 12D and 12H fix 2020.

---

## Completed

Full detail: the two `git show` commands at the top.

- **New species batch 1** (OSFL, GRSP, LEYE, BOBO), 2026-10-04 → 10-06: 12A moved to Fir; 12A,
  12C, 12D and 12H take species as arguments.
- **C3** (no-q99 sensitivity), 2026-10-02 → 10-04: `DROP_Q99=1`. 12F stops unless the
  re-predicted observed footprint, re-capped, reproduces `observed_bootstraps.tif`.
- **C2f rerun**, 2026-09-28 → 10-02: two masks per sector (Open Limitation #10); 07 674/674
  and 11 19/19 on the new masks. 11 now mosaics each subbasin onto its own window (can3 17.4 h
  → 1.4 h), and 07's memory is right-sized (24G / 72G).
- **C2a–C2e**, 2026-09-24 → 09-25:
  - C2a: 14B sums per-sample Shapley values.
  - C2b: kept 100 BART draws with replacement per bootstrap (user decision).
  - C2c: 10C rewritten as an area of applicability.
  - C2e: footprint covariates are 0 in the counterfactual (Open Limitation #9).
- **C1 passed** 2026-09-24. The **D. 12D rework** the same day added the NaN → `NA_real_`
  identity gate, the 12G C++ tree walk, the per-bootstrap reduction, and one task per
  species × BCR.
- **A7** (Fix A validated, 2026-09-23); **A6** (backfill + premosaic, 2026-09-22); **A**
  (cluster prep, 2026-09-14); **B** (V5 packaging conformance, gates G1–G4, 2026-09-11).
  - Re-run `Rscripts/misc/verify_weight_vs_v5_masking.R` after any weight change.
  - The truncation harness is `verify_v5_truncate_port.R`.

## V5 reference

Renumbering in `f082866`: `10.Package.R` → `10.Truncate.R` + `11.Package.R`;
`11.Validate` → `12.Validate`; `12.Summarize` → `13.Summarize`;
`output/10_packaged/` → `output/11_packaged/`, new `output/10_truncated/`.
`06_bootstraps` and `07_predictions` were NOT renumbered, and they are the only V5 output
folders our scripts read.
