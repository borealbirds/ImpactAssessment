# TODO (last updated 2026-10-13)

> **Where we left off (2026-10-13).** The C2d backfill (modelled blanks, no water, drawn
> land-cover classes; CLAUDE.md "Categorical covariates" and "Blanks are modelled") is in
> production: 07 674/674 and 11 19/19. 12D + 12H ran on it for CAWA, OVEN, OSFL, LEYE, GRSP and
> BOBO (all `nice.`), and 14B + 14C ran locally on their pulled outputs, so the local
> `sector_effects/` CSVs are current for these six.
> Next:
> 1. Run 14B, 14C and 18 on Fir: `bash 14B_sector_attribution.sh` (14B + 14C, then 18). It
>    should reproduce the local CSVs, and it is the first Fir run of the terra-only 18. Pull only
>    `sector_effects/*.csv` and `output_figures/*.png`.
> 2. Rerun the other 140 species on the C2d backfill (see "All species"). Checksum-sync the 12D
>    closure first ("Staging discipline"): Fir's `Rscripts/` was staged file by file during C2d.
> 3. 14B → 14C → 18 for all species, chained after the last 12H.
> 4. Redo the C3 (`DROP_Q99`) and C4 (`DENS_MIN`) sensitivity runs if they are to be reported;
>    their numbers below are pre-C2d.
> 5. The PDF report (`writing/results_report.typ`) is the very last step.
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

- **Withheld by BAM** (14B drops them): CAWA can40; GRSP can10, can61; BOBO can10; LEYE can10.
- **Direct vs vegetation** (Open Limitation #9; 14B's sign, cf − obs): CAWA direct +138k,
  vegetation +545k. OVEN direct −1.68 M, vegetation +5.30 M; its roads are −1.27 M direct and
  +1.74 M vegetation, hence their interval spanning 0.
- **What C2d changed:** GRSP can11 keeps 1.32 of today's birds on the all-sector footprint
  (old backfill 0.39; 16C's kNN 1.41), so industry now costs GRSP there ~220k birds instead of
  adding 440k. Old-backfill baseline per species × BCR:
  `cluster_logs/allspp_2026-10-09/baseline_old_backfill_species_bcr.csv`.
- **10C:** 4 of 674 subbasins flagged (50, 72, 430, 481); they carry 22k for CAWA, −71k for
  BOBO and < 4k for the others (14B's sign).
- **Pre-C2d sensitivity runs (CAWA, OVEN):** C3, without the frozen q99 cap: v(N) +1.16% (CAWA),
  +0.05% (OVEN), no sector share moves. C4, densities < 0.005 birds/ha set to 0: national v(N)
  −2.7% (CAWA) and −0.5% (OVEN), but CAWA can11 −81% (the prairie floor) while can13 and can61
  rise (the floor removes thin observed density while the backfilled forest stays above it).
  Undecided whether C4 becomes the headline.

Runs and logs (C2d):
- 07 `64664419` (reruns of 23–30, 247, 423) and 11 `64788561` (128G), plus the GRSP can11 12D
  smoke (`64818278`, 12H `64818279`): `cluster_logs/07_2020_c2d/`.
- Six species: 12D `64874085` (67 tasks), 12H `64874086`: `cluster_logs/c2d_6spp/`; 14B and 14C
  ran locally (`logs/14{B,C}_c2d_6spp_2026-10-13.log`).
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
- [ ] **Rerun on the C2d backfill.** `species_todo.txt` lists all 146, including the six already
      done, so drop those first or accept their rerun. Use `MAX=900` (2 waves): QOS `normal` sets
      no MaxSubmitPU/MaxJobsPU and Slurm's MaxArraySize is 1001, and the gaps between waves, not
      queue time, were the bottleneck. Sizing from waves 1–2 (`task_resources.csv`): 12D MaxRSS
      median 16.9 GB, max 44.8 GB (BAWW can11; can11 runs at 93% of 48G), so a task the C2d
      backfill pushes over is resubmitted at 64G (12H names it). 12H needs ~0.027 GB per per-BCR
      file, so its 32G holds a wave of ~1,000.
- [ ] Remaining species-at-risk candidates, not in the list: BANS, BARS, EAWP, EVGR, GCTH, GWWA.
      Check BAM's "remove" tab for their withheld BCRs, and stage `v5_gis/ranges/{SP}.tif` (from
      G: `gis/ranges/`) on Fir first.

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
- Closures: 07/11 = `07` `.R`+`.sh`, `08A`, `08B_*`, `09_*`, `11` `.R`+`.sh`.
  12D/12H = `12D` `.R`+`.sh`, `12E`, `12F`, `12G` `.R`+`.cpp`, `12H` `.R`+`.sh`.
  12A = `12A` `.R`+`.sh`, `12B`. 12C = `12C` `.R`+`.sh`.
  14B = `14B` `.R`+`.sh`, `14C`, `18` `.R`+`.sh`, `12E`, `data/raw_data/species_names.csv`,
  `data/raw_data/ModelReleaseDecisions.xlsx`, `data/derived_data/rds_files/extrapolation_flags.csv`.

---

## Open before the write-up

- [ ] **Ecological check of the farmland BCRs (can11–13) on the C2d results.** Before C2d, can11
      and can13 carried implausible changes: forest species (CAWA, OVEN, OSFL) gained more birds
      than they have today, and GRSP's can11 decline was a BART artefact (16C; report
      `writing/grassland_decline_report_2026-10-09.pdf`): BART never returned a blank canopy CV,
      so every backfilled pixel looked treed. C2d models blanks and draws classes; on 14 test
      subbasins (16C's share of today's birds kept on the direct footprint) it moved GRSP can11
      0.34 → 0.99 (kNN 1.41), BOBO can11 0.63 → 0.61 (0.65), BOBO can13 0.13 → 0.28 (0.23).
      The user rejected a kNN hot-deck backfill (2026-10-12: "a new method on top of BART that is
      a bit arbitrary"). OVEN's negatives in can12, can80 and can81 come from roads through the
      bird model's footprint response (Open Limitation #9).

## Loose ends

- [ ] **BART and kNN still differ where training land is thin or far (C2d test 2).** BOBO can13
      subbasin 98 keeps 0.36 of its observed birds vs kNN 0.23 (+177k; more SCANFI herbs and
      shrub; it trains on can12 land and 59% of can13's birds are outside the AOA). GRSP can11
      subbasins 55 and 60 keep 0.73 and 0.55 vs kNN 1.12 and 1.18 (prcD_5x5 above kNN; likely
      the Gaussian amount model's missing point mass at 0). Neither is a coding defect, and both
      methods agree on the direction. Report it as method sensitivity, with 10C's AOA flags.
- [ ] **can3 (Arctic) completeness:** before C2d, LEYE can3 (and every species' can3) dropped
      17.3% of weight > 0 footprint pixels, because VLCE_1km had no training rows in four tundra
      subbasins. C2d makes a covariate blank on all the training land a modelled blank
      everywhere; confirm 100% complete in the rerun's 12D logs (not pulled for the six).
- [ ] **11 needs 128G / 8 h for can11, can12, can13 and can61 since C2d** (~2400 layers per
      subbasin). At the .sh's 32G, terra's block path stalled > 4 h on one subbasin; at 128G all
      of can11 took 30 min. The .sh still says 32G / 4 h, and 07's .sh chains 11 without
      overriding it.
- [ ] 16C's BART reconstruction (8 bootstraps × 10 draws) reads GRSP ~10% below 12F's cf, so
      only its T2/T5 magnitudes are approximate. Its T1 reads `cluster_logs/grassland_diag/backfill/`,
      deleted 2026-10-13: pull the C2d production backfills there before rerunning 16C.
- [ ] 12A on Fir (optional check): run it for CAWA into a scratch dir and diff against the
      G:-built `observed_bootstraps.tif` / `truncation_params.rds`.
- [ ] **Cheaper 12D (optional):** 16 distinct BART draws per bootstrap instead of 100 picks with
      replacement would cut a CAWA + OVEN run from ~64 to ~25 billed core-equivalent-hours, and
      move the all-sector mean by 0.08% (CAWA) / 0.18% (OVEN) of the impact.
- [ ] `15C_singletons_plot.R:150` reads `predictions_coalitions/`, deleted on both ends, so the
      map panel is broken. Regenerate via `save_arrays_ids` if wanted.
- [ ] **Pre-2020 years** need historical footprint layers (not the 2020 masks), each year's
      `covariates_mosaiced_{year}.tif` (06; only 2020 exists) and a per-year 07. 07 and 11 take
      years; 12A, 12C, 12D and 12H fix 2020.
- Deleted locally to free disk, needed to rebuild: `hirshpearson/raw_300m/` (03, 14A;
  re-download from Borealis doi:10.5683/SP2/EVKAVL + `cum_threat` from G:).

---

## Completed

Full detail: the `git show` commands at the top.

- **C2d** (backfill defect), 2026-10-09 → 10-13: modelled blanks (`08B_deploy_pbart.R`), no water
  class in training, 100 drawn land-cover classes from normalised probit probabilities (not
  `mbart2`'s softmax, not the modal class), classes renumbered 1..K; 07 + 11 rerun; six species
  through 12D/14C. 14B, 14C and 18 can now run on Fir (not yet run there); 18 dropped `sf`.
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
