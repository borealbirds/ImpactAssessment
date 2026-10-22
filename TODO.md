# TODO (last updated 2026-10-20)

> **Where we left off (2026-10-14).** The C2d backfill (modelled blanks, no water, drawn
> land-cover classes; CLAUDE.md "Categorical covariates" and "Blanks are modelled") is in
> production: 07 674/674 and 11 19/19. 12D + 12H ran on it for all 146 species, and 14B + 14C +
> 18 for all of them on Fir (2026-10-14): 138 species reported (8 have every BCR withheld by
> BAM), 138 maps; CSVs and maps pulled to `sector_effects/` and `output_figures/`.
> Next:
> 1. ~~Rerun the other 140 species on the C2d backfill~~: done 2026-10-14, all clean ("All species").
> 2. ~~14B → 14C → 18 for all species~~: done 2026-10-14 ("All species", "Results").
> 3. ~~Redo C3 and C4~~: not to be redone (user, 2026-10-14). Their pre-C2d results below stay,
>    for reporting; production (no floor, frozen q99) is the headline.
> 4. Supplementary results for the Google Doc "Supp. Mat. for: Declining boreal birds", from page
>    25 on, instead of a PDF (138 maps would make one huge): `writing/supp_results_docx.R` builds a
>    .docx in `results_report.typ`'s format plus a BCR sector table after each species' map
>    (14C's new `results_bcr.csv`, 14C rerun alone as `65037894`, 2026-10-14; 1,462 species × BCR
>    rows). Built: `writing/supp_results_2026-10-14.docx` (22 MB, 239 Letter pages, Tables
>    S1–S140, Figures S1–S138; the build takes ~28 min). Left: paste it into the Doc after the
>    Supplementary Methods. The Drive connector reads the Doc (shared "anyone with the link" since
>    2026-10-14) but has no tool that edits a Doc's body.
> 5. ~~Backfill 1995, 2000, 2005, 2010, 2015~~: done 2026-10-20, all clean ("Pre-2020
>    backfill" below). Birds are not re-predicted for those years (12A–14B fix 2020).
>
> The user wants the project finished: log small new problems under Loose ends and do not
> reopen finished stages for them. Keep `writing/` off the internet: never commit or push it.
> `rds_files/` and `sector_effects*/` CSVs stay local, never committed (user, 2026-10-04).
> SSH is 2FA, so the user runs every cluster command; Globus works without it.
>
> Full records of finished work: `git show c91a863:TODO.md` (to 2026-10-13: the C2d
> investigation and its tests, waves 1–2, C4, C5), `git show e619087:TODO.md` (to 2026-10-04)
> and `git show 3fcef6d:TODO.md` (to 2026-09-28).

## Results (2020, C2d backfill, BAM release filter applied)

From 14C (`sector_effects/results_summary.csv`, `results_sectors.csv`). Change = observed −
no-industry (negative = birds lost to industry); % is of the no-industry population.

| Species | Observed | No industry | Change (5–95%) | % | Largest sectors (%) |
|---|---|---|---|---|---|
| CAWA | 4.47 M | 5.15 M | −683k (−782k, −564k) | −13.3 | roads −5.7, pasture −3.0, crop −2.7 |
| OVEN | 37.65 M | 41.27 M | −3.62 M (−4.14 M, −2.75 M) | −8.8 | crop −3.3, pasture −3.0, built −1.3 |
| OSFL | 10.74 M | 11.11 M | −372k (−537k, −217k) | −3.3 | roads −2.1 |
| LEYE | 5.79 M | 5.83 M | −41k (−61k, −25k) | −0.7 | roads −0.6 |
| GRSP | 0.96 M | 0.90 M | +59k (−293k, +408k) | +6.6 | roads +7.1, pasture +2.6, crop −2.5 |
| BOBO | 6.76 M | 3.63 M | +3.14 M (+2.96 M, +3.33 M) | +86.4 | roads +30.3, crop +26.6, pasture +20.2 |

- **All 138 species** (`results_summary.csv`, sorted by %): 88 decline (72 with the whole 5–95%
  interval below 0), 50 increase (38 above 0), 28 intervals span 0; median −1.7%. Largest
  declines: EATO −62.2%, FISP −43.3%, GCFL −30.7%, RBGR −23.7%, PIWA −21.7%. Largest increases,
  all open-country or human-associated species: BARS +156.9%, ROPI +151.7%, YHBL +112.7%,
  COGR +104.4%, HOSP +104.1%. The sector with the largest |%| is roads for 98 species, crop 33,
  pasture 5, built 1, reservoirs 1. EATO's and FISP's no-industry populations are 2.7× and 1.8×
  today's, mostly in can13 (southern Ontario; EATO also 91k in can11); not examined.
- **Withheld by BAM** (14B drops them): CAWA can40; GRSP can10, can61; BOBO can10; LEYE can10.
  Every BCR of BGGN, BWWA, DUNL, GWWA, NESP, SPGR, YBCU and YTVI, so 14B skips those 8.
- **Direct vs vegetation** (Open Limitation #9; 14B's sign, cf − obs): CAWA direct +138k,
  vegetation +545k. OVEN direct −1.68 M, vegetation +5.30 M; its roads are −1.27 M direct and
  +1.74 M vegetation, hence their interval spanning 0.
- **What C2d changed** (share of today's birds the no-industry landscape keeps on the all-sector
  footprint, old → C2d): GRSP can11 0.39 → 0.94, so industry adds +41k there (SD 210k) instead
  of +440k; GRSP nationally +463k → +59k. BOBO can13 0.14 → 0.27 (a kNN check gave 0.23 on its main
  subbasins), can11 0.65 → 0.64; BOBO nationally +3.42 M → +3.14 M. (The can11 smoke's 1.32
  came from 2 bootstraps; production's 32 do not bear it out.) Old-backfill baseline per
  species × BCR: `cluster_logs/allspp_2026-10-09/baseline_old_backfill_species_bcr.csv`.
- **10C:** 4 of 674 subbasins flagged (50, 72, 430, 481); they carry 22k for CAWA, −71k for
  BOBO and < 4k for the others (14B's sign).
- **Pre-C2d sensitivity runs (CAWA, OVEN):** C3, without the frozen q99 cap: v(N) +1.16% (CAWA),
  +0.05% (OVEN), no sector share moves. C4, densities < 0.005 birds/ha set to 0: national v(N)
  −2.7% (CAWA) and −0.5% (OVEN), but CAWA can11 −81% (the prairie floor) while can13 and can61
  rise (the floor removes thin observed density while the backfilled forest stays above it).
  Neither is rerun on C2d (user, 2026-10-14): production is the headline, these stay as reported.

Runs and logs (C2d):
- 07 `64664419` (reruns of 23–30, 247, 423) and 11 `64788561` (128G), plus the GRSP can11 12D
  smoke (`64818278`, 12H `64818279`): `cluster_logs/07_2020_c2d/`.
- Six species: 12D `64874085` (67 tasks), 12H `64874086`: `cluster_logs/c2d_6spp/`; 14B and 14C
  ran locally (`logs/14{B,C}_c2d_6spp_2026-10-13.log`), then on Fir: 14B + 14C `64902653`, 18
  `64902654` (logs and Fir's CSVs in `cluster_logs/c2d_6spp/`). 18's 10 km grids sum to 14B's
  national change for every species.
- All species: 12D/12H waves 1–2, 14B + 14C `65014117`, 18 `65014118`:
  `cluster_logs/c2d_allspp_2026-10-14/` ("All species").
- Pre-C2d runs (C1, C3, C4, C5, new species, waves 1–2): `git show c91a863:TODO.md`.

## All species

`Rscripts/run_species_batches.sh` runs 12A + 12C (`prep`), a can42 smoke (`smoke`), then 12D +
12H in waves (`plan`, then `wave 1`, `wave 2`, ... one at a time; `MAX=` tasks per wave) for
`Rscripts/species_todo.txt`: 146 species, 1,624 species × BCR tasks. Left out: DEJU, DOWO,
EUST, FOSP (bootstrap models but no 2020 `07_predictions` on Fir) and ROSA (no Canadian models).
Commands for a hand-run batch: CLAUDE.md "Re-predicting birds".

- [x] `prep` (12A `63740250`, 12C `63740251`, 2026-10-08) and `smoke` (2026-10-09), logs in
      `cluster_logs/allspp_2026-10-09/`. Neither depends on 07, so neither reruns.
- [x] Waves 1–2 (80 species) ran on the old backfill; every density product was deleted
      2026-10-12. Kept: per-task runtimes (`cluster_logs/allspp_2026-10-09/task_timing.tsv`:
      median 5.7 min, max 66 min AMRO can11) and the old-backfill baseline above.
- [x] **Rerun on the C2d backfill: done 2026-10-14** (logs `cluster_logs/c2d_allspp_2026-10-14/`).
      Wave 1: 12D `64916943` (48G) + `64916944` (64G), 12H `64916945`; wave 2: 12D `64916949` +
      `64916950`, 12H `64916952`. 1,557/1,557 tasks `nice.`, no errors; identity gate 1,554/1,554
      (the other 3, BAWW/BBWA/RBGR can41, have weight 0 on the whole footprint: impact 0); complete
      weight > 0 footprint ≥ 99.9% in every BCR, can3 included; both 12H `nice.`.
      14B `64916954` stopped at BGGN: BAM withheld every BCR of 8 species (BGGN, BWWA, DUNL, GWWA,
      NESP, SPGR, YBCU, YTVI), leaving 0 rows. 14B now skips such a species (staged); 18 never ran.
      Rerun `65008610`: 14B completed all 146 (138 reported, the 8 skipped), then 14C stopped at
      EATO: its observed-population check was relative (1e-6) against 14B's whole-bird rounding
      (198,987.4 vs 198,987); everything else matched. Now within 1 bird, like 18's (staged).
      Rerun `65014117` (14B + 14C) and `65014118` (18): 14B skipped the 8 and, for the other
      138, the sample means matched the tables' Shapley values (max diff 2.1e-9) with exact
      additivity; 14C wrote both tables (138 rows); 18 mapped all 138 (every pixel raster
      present, each 10 km grid within 0.5 bird of 14B's national change, ≤ 21 s a species).
      The six done earlier reproduce their 2026-10-13 numbers exactly.
- Setup of that rerun (record): `species_todo.txt` lists all 146 with the six already done
      commented out (140 left; staged, and the 12D/12H and 14B closures checksum-synced,
      2026-10-13: Fir's copies equal the code that ran the six). Use `MAX=900` (2 waves): QOS `normal` sets
      no MaxSubmitPU/MaxJobsPU and Slurm's MaxArraySize is 1001, and the gaps between waves, not
      queue time, were the bottleneck. Sizing from waves 1–2 (`task_resources.csv`): 12D MaxRSS
      median 16.9 GB, max 44.8 GB (BAWW can11; 20 of 63 can11 tasks above 40 GB). The six C2d
      species reached 44.5 GB in can11 (`sacct` 64874085), so 12D's .sh now gives can11 tasks
      64G and the rest 48G (2026-10-13, staged); a task still killed for memory is resubmitted
      at 64G (12H names it). 12H needs ~0.027 GB per per-BCR file, so its 32G holds a wave of
      ~1,000.
- [x] Species-at-risk candidates BANS, BARS, EAWP, EVGR, GCTH and GWWA are in the list and ran
      (GWWA has every BCR withheld, so it is not reported).

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
- Re-login: run `globus session update computecanada.ca --no-local-server` in a PowerShell
  window (not via `!`).
- Fir's `sacct -S/-E` times are Pacific, not UTC (Globus listings are UTC). Logs land in
  `Rscripts/` on Fir; sort a listing by its time column before concluding nothing ran.
- The laptop's C: drive is 236 GB and has filled to 0 bytes (`NO_SPACE_LEFT`): check
  `Get-PSDrive C` before any large pull.

## Staging discipline: read before any cluster run

Stage the entry point's whole **`source()` closure**, including data files, never just what you
edited, and checksum-sync all of it (a 0-byte transfer proves equality for free). Unedited files
have been stale or missing on Fir four times.

- **Never restage an `.R` file while a job that runs it is active.** Rscript buffers its file
  and reads again at the end; Globus rewrites the file in place, so a running job can read the
  new file's tail as code.
- `.sh` files must be LF (a CRLF script dies with `/bin/bash^M: bad interpreter`). Check with
  `tr -dc '\r' < f | wc -c` (0 = LF), not `grep -c $'\r'`.
- Globus refuses sources outside the local endpoint's root, so stage from inside the repo tree.
- Before staging a 12F change, run `Rscript Rscripts/misc/verify_12f_can60_local.R` (~4 min) on
  Fir's real CAWA can60 inputs (`cluster_logs/localtest/`, `cluster_logs/smoke7/`, gitignored).
  Its smoke-7 baseline predates C2f, so compare the new 12F against the previous 12F on the same
  inputs: tables and arrays `identical()` unless predictions or masks changed on purpose. Arrays
  may differ from Fir by ~3e-16 relative (Fir vs Windows floating point).
- Closures: 06 = `06` `.R`+`.sh`, `data/raw_data/CAfire/CAfire_{year}_masked.tif`,
  `data/raw_data/isric_soil/isric_soil_covariates_masked.tif`.
  07/11 = `07` `.R`+`.sh`, `08A`, `08B_*`, `09_*`, `11` `.R`+`.sh`.
  12D/12H = `12D` `.R`+`.sh`, `12E`, `12F`, `12G` `.R`+`.cpp`, `12H` `.R`+`.sh`.
  12A = `12A` `.R`+`.sh`, `12B`. 12C = `12C` `.R`+`.sh`.
  14B = `14B` `.R`+`.sh`, `14C`, `18` `.R`+`.sh`, `12E`, `data/raw_data/species_names.csv`,
  `data/raw_data/ModelReleaseDecisions.xlsx`, `data/derived_data/rds_files/extrapolation_flags.csv`.

---

## Open before the write-up

- [ ] **Ecological check of the farmland BCRs (can11–13) on the C2d results.** Before C2d, can11
      and can13 carried implausible changes: forest species (CAWA, OVEN, OSFL) gained more birds
      than they have today, and GRSP's can11 decline was a BART artefact: BART never returned a
      blank canopy CV, so every backfilled pixel looked treed. C2d models blanks and draws
      classes; on 14 test subbasins (share of today's birds kept on the direct footprint) it moved
      GRSP can11 0.34 → 0.99 (kNN 1.41), BOBO can11 0.63 → 0.61 (0.65), BOBO can13 0.13 → 0.28
      (0.23). The user rejected a kNN hot-deck backfill (2026-10-12: "a new method on top of BART
      that is a bit arbitrary"). OVEN's negatives in can12, can80 and can81 come from roads
      through the bird model's footprint response (Open Limitation #9).
      C2d barely moved the forest species' farmland gains (no-industry birds on the all-sector
      footprint per bird today, old → C2d): OVEN can11 4.7 → 4.9, can13 2.5 → 2.1; CAWA can11
      3.7 → 2.6, can13 5.1 → 4.3; OSFL can13 12.9 → 5.9. The user (2026-10-13) accepted them for
      the all-species run, to be handled in this check before the write-up.

## Loose ends

- [ ] **BART and kNN still differ where training land is thin or far (C2d test 2).** BOBO can13
      subbasin 98 keeps 0.36 of its observed birds vs kNN 0.23 (+177k; more SCANFI herbs and
      shrub; it trains on can12 land and 59% of can13's birds are outside the AOA). GRSP can11
      subbasins 55 and 60 keep 0.73 and 0.55 vs kNN 1.12 and 1.18 (prcD_5x5 above kNN; likely
      the Gaussian amount model's missing point mass at 0). Neither is a coding defect, and both
      methods agree on the direction. Report it as method sensitivity, with 10C's AOA flags.
- [x] **can3 (Arctic) completeness:** before C2d, LEYE can3 (and every species' can3) dropped
      17.3% of weight > 0 footprint pixels (604 of 730 kept), because VLCE_1km had no training
      rows in four tundra subbasins. C2d makes a covariate blank on all the training land a
      modelled blank everywhere. Checked 2026-10-13 on LEYE's C2d `pixel_change/` raster (the 12D
      logs were cleared from Fir): 723 footprint pixels with observed density > 0, so ≥ 99% kept.
- [x] **11 needs 128G / 8 h for can11, can12, can13 and can61 since C2d** (~2400 layers per
      subbasin; at 32G terra's block path stalled > 4 h on one subbasin, at 128G all of can11
      took 30 min). Fixed 2026-10-14: 07's .sh submits those four (`BIG_BCRS`) at 128G / 8 h.
- [ ] 12A on Fir (optional check): run it for CAWA into a scratch dir and diff against the
      G:-built `observed_bootstraps.tif` / `truncation_params.rds`.
- [ ] **Cheaper 12D (optional):** 16 distinct BART draws per bootstrap instead of 100 picks with
      replacement would cut a CAWA + OVEN run from ~64 to ~25 billed core-equivalent-hours, and
      move the all-sector mean by 0.08% (CAWA) / 0.18% (OVEN) of the impact.
- [ ] `15C_singletons_plot.R:150` reads `predictions_coalitions/`, deleted on both ends, so the
      map panel is broken. Regenerate via `save_arrays_ids` if wanted.
- [x] **Pre-2020 backfill (1995–2015), 2026-10-14 → 10-20.**
      - **Masks: the 2020 ones** (`HF_MASK_YEAR=2020`, user 2026-10-14). Hirsh-Pearson exists
        only for 2020. Every year trains on land with no footprint in 2020 and backfills the
        whole 2020 footprint; where a 2020 feature didn't exist yet, observed and backfilled
        should agree (17 measures this).
      - **06 moved to Fir.** Elly's stacks for 1985–2020 are in `{nm_root}/gis/stacks/` (Fir's
        `can60_2020.tif` is byte-identical to G:'s). Building the mosaics on the laptop would
        have streamed ~3 GB of stacks per year from G: onto a disk with 23 GB free.
      - **The CAfire layers for 1995–2015 on G: were pre-fix.** Unburned land was NA, with no
        zeros; only 2020 had been rebuilt after `0bbf2ea`. Instead of a 6 h rerun of 02, the fix's
        two steps (NA → 0, then `mask(bam_boundary)`) were applied to the existing files, into
        `data/raw_data/CAfire/` (G: left as is). Checked: the NA cells are identical to the fixed
        2020 file's (20,700,502), the burned values are unchanged, and every burned cell is
        burned in 2020 too. A rerun of 02 would rebuild the same files.
      - **07.sh submits one array per year.** The old year-block numbering needed indices up to
        674 × years, above Slurm's MaxArraySize of 1001 for any 2+ years.
      - Staged 2026-10-14 (06 `.R`+`.sh`, 07 `.sh`, CAfire 1995–2020, soil; the rest of the
        07/11 closure was already identical on Fir).
      - **06 done**: `65094895`, tasks 1–5, finished 2026-10-15 ~05:50 UTC (logs in
        `cluster_logs/06_covmosaic/`). Each year read all 19 stacks plus CAfire and soil, and has
        88 layers, with the same names, order and grid as 2020. Each file is 1.65 GB.
        `sacct`: 1 h 52–53 min and 13.8–14.2 GB MaxRSS per year, so 32G / 12 h has room.
      - **Smoke passed** 2026-10-19: `66003517` (1995, subbasin 1, 24G; 4 min 36 s, 2.7 GB;
        logs in `cluster_logs/07_1995_smoke/`). It read the undated 2020 masks, set the 5
        footprint covariates to 0, wrote 2146 layers (11.9 MB vs 2020's 12.0 MB) and returned
        `ok = TRUE`. Its only warning is the dplyr `across(...)` deprecation at `07:53`, which
        the 2020 run printed too.
      - **Full run submitted 2026-10-19** (`HF_MASK_YEAR=2020`). Job ids, in the order
        07 small (24G) / 07 big (72G) / 11 small (32G) / 11 big (128G):
        - 1995: `66010503` / `66010505` / `66010566` / `66010567`
        - 2000: `66010568` / `66010569` / `66010570` / `66010572`
        - 2005: `66010573` / `66010574` / `66010575` / `66010576`
        - 2010: `66010587` / `66010588` / `66010589` / `66010590`
        - 2015: `66010591` / `66010592` / `66010593` / `66010594`

        Smoke command:
        `YEARS=1995 HF_MASK_YEAR=2020 sbatch --array=1 --mem=24G 07_train_and_backfill.sh`; then
        `HF_MASK_YEAR=2020 bash 07_train_and_backfill.sh 1995 2000 2005 2010 2015`.
      - **Done, all clean** (checked 2026-10-20 on Fir):
        - Every year has 674/674 subbasins with a `_backfill.tif` and 19/19 BCR mosaics, with no
          `_partial` left.
        - `sacct` shows all 3,370 07 tasks and 95 11 tasks COMPLETED.
        - No log has `Error in subbasin`, any capitalised `Error`, `Execution halted`, a time
          limit or an OOM kill. Every 07 log matches lowercase "error" only in the tidyverse
          banner ("force all conflicts to become errors"), as the 2020 logs do.
        - Peak MaxRSS per tier: 07 small 17.6–18.7 GB of 24G; 07 big 65.8–67.3 GB of 72G, all
          from subbasin 98 (2020's max was 52.8 GB, so 72G has only ~5 GB to spare); 11 small
          11.5 GB of 32G; 11 big 23.4 GB of 128G (can61). The 128G is for time, not memory:
          terra's block path stalls at 32G.
        - Longest task: 07 2 h 17 min; 11 2 h 00 min.
        - Every mosaic is within 4% of its 2020 size (can13 is the largest gap, 879 vs 880–916 MB).
      - **Fir scratch cleaned 2026-10-20:**
        - the Slurm logs older than this run, already pulled to `cluster_logs/`;
        - the empty `Rscripts/{misc,deviation_from_expected}/` and
          `density_tables_{noq99,dmin0.005}/` folders, which 12D and 12H recreate if needed.

        Left for the user to decide: `/home/mannfred/scratch/v4_bootstraps/` (196 GB of V4
        models from 2025-08, read by no script here).
      - Downstream, 12A, 12C, 12D, 12H, 14A and 14B still fix 2020. Elly's
        `07_predictions/{sp}/{sp}_{bcr}_{1990..2020}.tif` are on Fir.
- Deleted locally to free disk, needed to rebuild: `hirshpearson/raw_300m/` (03, 14A;
  re-download from Borealis doi:10.5683/SP2/EVKAVL + `cum_threat` from G:).

---

## Completed

Full detail: the `git show` commands at the top.

- **C2d** (backfill defect), 2026-10-09 → 10-13: modelled blanks (`08B_deploy_pbart.R`), no water
  class in training, 100 drawn land-cover classes from normalised probit probabilities (not
  `mbart2`'s softmax, not the modal class), classes renumbered 1..K; 07 + 11 rerun; six species
  through 12D/14C. 14B, 14C and 18 run on Fir (first run 2026-10-13); 18 dropped `sf`.
- **Results tables and maps** (14C, 18), 2026-10-07 → 10-08: change = observed − cf as % of the
  no-industry population (a colleague found positive = loss confusing).
- **C5** (pixel-level change maps, `pixel_change/` + 18's 10 km grid) and **C4** (`DENS_MIN`
  floor), 2026-10-06 → 10-07.
- **New species batch 1** (OSFL, GRSP, LEYE, BOBO), 2026-10-04 → 10-06: 12A moved to Fir; 12A,
  12C, 12D and 12H take species as arguments.
- **C3** (no-q99 sensitivity), 2026-10-02 → 10-04: `DROP_Q99=1`.
- **C2f rerun**, 2026-09-28 → 10-02: two masks per sector (Open Limitation #10); 11 mosaics each
  subbasin onto its own window; 07's memory right-sized (24G / 72G).
- **C2a–C2e**, 2026-09-24 → 09-25: per-sample Shapley sums (C2a), 100 BART draws with
  replacement (C2b), 10C as an area of applicability (C2c), footprint covariates at 0 (C2e).
- **C1 passed** and the **12D rework** (identity gate, C++ tree walk, per-bootstrap reduction,
  one task per species × BCR), 2026-09-24.
- **A7, A6, A** (Fix A, backfill, cluster prep) and **B** (V5 packaging conformance, gates
  G1–G4), 2026-09-11 → 09-23. Harnesses: `Rscripts/misc/verify_weight_vs_v5_masking.R` (rerun
  after any weight change), `verify_v5_truncate_port.R`.

## V5 reference

Renumbering in `f082866`: `10.Package.R` → `10.Truncate.R` + `11.Package.R`;
`11.Validate` → `12.Validate`; `12.Summarize` → `13.Summarize`;
`output/10_packaged/` → `output/11_packaged/`, new `output/10_truncated/`.
`06_bootstraps` and `07_predictions` were NOT renumbered, and they are the only V5 output
folders our scripts read.
