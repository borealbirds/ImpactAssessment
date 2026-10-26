# TODO (last updated 2026-10-22)

> **Where we left off (2026-10-22).** 2020 is done for all 146 species on the C2d backfill
> (07 + 11, then 12D/12H and 14B/14C/18 on Fir, 2026-10-14): 138 species reported (BAM withholds
> every BCR of the other 8), 138 maps; CSVs and maps pulled to `sector_effects/` and
> `output_figures/`. The 1995–2015 backfill (07 + 11) is done and clean (2026-10-20).
> Open:
> 1. **Supplementary results** for the Google Doc "Supp. Mat. for: Declining boreal birds", from
>    page 25 on, instead of a PDF: `writing/supp_results_docx.R` builds
>    `writing/supp_results_2026-10-14.docx` (each species' map, then its BCR sector table from
>    14C's `results_bcr.csv`; the build takes ~28 min). Left: paste it into the Doc after the
>    Supplementary Methods. The Drive connector reads the Doc but cannot edit its body.
> 2. **Pre-2020 densities (1995–2015) for the 113 reported species** (`Rscripts/species_113.txt`:
>    species_todo.txt's 146 minus BAM's 8 fully withheld and the 25 wetland/grassland/open
>    species). 12A, 12C, 12D and 12H read `YEARS` (staged 2026-10-22, whole closure
>    checksum-synced). Each year uses its own observed stack (`{bcr}_{year}.tif`), Elly's year-y
>    predictions, a q99 frozen from year y's observed landscape (as V5 truncates each year),
>    year y's backfill, and the 2020 footprint (as 07 did). 12E, 12F and 12G are unchanged, so
>    `code_md5` matches the 2020 run.
>    - Done: `prep`: 12A `66813898` (39 tasks OOM-killed at 16G, resumed at 32G as `66820731`),
>      12C `66813899`; logs and Fir's `truncation_params.rds` in `cluster_logs/pre2020_prep/`.
>      1,316 species × BCR per year (6,580 in all), each with Elly's prediction tif, a weight and
>      a finite q99. q99 / 2020's q99: median 1.00, 5–95% 0.61–1.47; the 284 BCR-years beyond ½×
>      or 2× are irruptive or budworm species (EVGR can80 1995 42.6×, BBWA can80 2000–05 ~5×),
>      likely real change in Elly's year models, not a defect.
>    - Done: `smoke` (can42; AMRO, ATSP, CAJA; 2 bootstraps; 12D / 12H: 1995 `66834279` /
>      `66834280`, 2000 `66834281` / `66834282`, 2005 `66834283` / `66834284`, 2010 `66834285` /
>      `66834286`, 2015 `66834287` / `66834288`; logs in `cluster_logs/pre2020_smoke/`). All 15
>      tasks: identity gate 9/9, complete weight > 0 footprint 100%, `nice.`; all five 12H `nice.`.
>      can42 has only 9 footprint pixels, but their observed birds differ between years (CAJA
>      58.7 in 1995 vs 102.3 in 2005, same 2 bootstraps), so passing the gate every year shows
>      each year read its own stack. `code_md5` equals the 2020 production files'.
>    - Done: 12D + 12H for all 113 species × 5 years (2026-10-22 to 10-24). Wave 1 (78 species AMCR–PUFI,
>      891 tasks per year; 12H 1995 `66836901`, 2000 `66836904`, 2005 `66836907`, 2010
>      `66836910`, 2015 `66836913`); wave 2 (35 species RBGR–YRWA, 425 per year; Fir's job IDs
>      had wrapped: 12H `74263`, `74266`, `74302`, `74305`, `74308`). Logs in
>      `cluster_logs/pre2020_full/`. All 6,580 tasks `nice.`, each species × BCR × year once, no
>      OOM or timeout. Identity gate all 32 bootstraps at n/n pixels in 6,565; the other 15 are
>      BAWW, BBWA and RBGR can41 × 5 years (no weight > 0 footprint pixels, impact 0, as in
>      2020). Every other task kept 100% of its weight > 0 footprint pixels. All ten 12H `nice.`,
>      no warnings: 255 tables + Shapley samples for 113 species in every year.
>    - Done: 14B, 14C and 18 run one year each (`YEARS`, default 2020; staged 2026-10-24, 12E
>      unchanged so `code_md5` holds). 14B's .sh submits one 14B + 14C and one 18 per year.
>      2020 keeps its file names; other years' CSVs and maps carry `_{year}`, so 14D and
>      `writing/` still read 2020. Each year's map uses its species' 2020 colour scale. Tested
>      locally on the six species' 2020 tables relabelled 1995: all six CSVs identical to
>      production 2020, maps identical but for the title.
>    - Done (2026-10-24): 14B + 14C / 18 for 1995 `289752` / `289753`, 2000 `289754` /
>      `289755`, 2005 `289756` / `289757`, 2010 `289758` / `289759`, 2015 `289760` / `289761`;
>      logs in `cluster_logs/pre2020_report/`, the 30 CSVs in `sector_effects/` and 565 maps in
>      `output_figures/` pulled. Every year: 113 species, sample means = tables' Shapley values
>      (max diff 2.2e-9), per-sample additivity residual 0, 14C's checks passed, 113 maps with
>      every pixel raster present and each grid within 0.5 bird of 14B's national change.
>      Results: each year's species % changes correlate 0.96 with 2020's; median % −3.4 (1995),
>      −2.7, −2.3, −3.0, −2.8, −2.6 (2020); 61–66 species with the whole interval below 0 each
>      year. CAWA −15.7% (1995) to −13.3% (2020), OVEN −10.5% to −8.8%, OSFL −4.9% to −3.4%.
>      Outliers, not examined: EVGR's 1995 observed population is 5.2× 2020's (irruptive; see
>      the q99 note above); WITU (41–66k birds) is +1 to +21% before 2020, +93% in 2020.
> 3. **Rosenberg et al. (2019) Fig. 1 style figures, one file per panel** (2026-10-26). Their
>    panels A–C sum ALL species of a group (increasing and declining); D is the share declining.
>    `19A` (Fir) writes `sector_effects/population_samples.rds` (per-sample observed and
>    no-industry population, species × year), since 14C's percentiles cannot be summed or
>    differenced; tested locally on the six 2020 species (means and every 5/95 percentile equal
>    14C's). `19B` (local) plots panel A, net change since 1995, observed and without industry,
>    each from its own 1995 total. Done: 19A job `603212` (7 min, every check passed; log in
>    `cluster_logs/rosenberg_19A/`), rds pulled, 19B → `output_figures/net_change_1995_2020.png`.
>    2020 vs 1995: observed +58 M (5–95% +12 to +99 M), no industry +66 M (+12 to +113 M). Per
>    sample, no industry − observed is −20 M (1995) to −12 M (2020), every interval below 0;
>    the trends differ by +8 M by 2020 (+0.1 to +15 M), within 0 in every earlier year.
>    Done: panel B, `19C` → `net_change_by_habitat_1995_2020.png` (Rosenberg's B is one
>    trajectory per biome; S1 in the manuscript PDF shows the layout). Birds industry costs
>    (no industry − observed), 1995 → 2020: conifer +54 → +45 M, deciduous +28 → +27 M,
>    mixedwood +6 → +7 M, forest no type −2 M, non-forest −106 → −90 M; so the near-zero total is
>    forest losses cancelled by non-forest gains. Trend differences by 2020: conifer −9 M (−14 to
>    −6), non-forest +17 M (+12 to +21), the rest within 0.
>    Done: panel C, `19D` → `percent_change_by_habitat_1995_2020.png`. 1995→2020 observed / no
>    industry: conifer −4.7 / −5.1%, mixedwood +14.7 / +14.7%, deciduous +8.1 / +7.0%, forest no
>    type −7.4 / −7.1%, non-forest +9.0 / +12.0%. Per sample (no industry − observed): conifer
>    −0.4 points (−0.8 to −0.2), deciduous −1.1 (−1.6 to −0.7), non-forest +3.0 (+2.2 to +3.7),
>    mixedwood and forest no type within 0.
>    Done: panel D, `19E` → `species_declining_by_habitat_1995_2020.png` (declining = mean change
>    < 0; solid = whole interval < 0). Observed / no industry: conifer 12 / 13 of 27 (8 / 7 sure),
>    mixedwood 2 / 2 of 10 (none sure), deciduous 4 / 4 of 22, forest no type 2 / 2 of 6,
>    non-forest 18 / 15 of 48 (12 / 10 sure). Decline only as observed: AMRO, EAKI, EAPH, ROPI,
>    TRES; only without industry: ATTW, HOWR, WEWP. All four panels done.
>    From 14C's means: the 113-species total is 3.02 bn (1995) to 3.08 bn (2020) observed, and
>    the no-industry total is LOWER every year (by 20 M in 1995, 12 M in 2020: AMRO, BHCO, AMGO,
>    SOSP, COGR, BARS, HOSP gain 120 M on the footprint, the rest lose 108 M), so the two lines
>    nearly overlap. The +104 M observed jump 2015→2020 is 45 M PISI (irruptive; Anna and Diana
>    drop the irruptives, Loose ends).
> 4. **Maps moved to `output_figures/change_maps/`** (2026-10-26): all 703 local maps, the 138
>    2020 maps renamed `change_grid_{species}_2020.png`. 18 now writes there with the year on
>    every map, and `writing/` reads the new 2020 names. Not staged to Fir yet: stage 18 before
>    its next run. Fir's `output_figures/` still holds the 2020 maps under their old names.
>
> The user wants the project finished: log small new problems under Loose ends and do not
> reopen finished stages for them. Keep `writing/` off the internet: never commit or push it.
> `rds_files/` and `sector_effects*/` CSVs stay local, never committed (user, 2026-10-04).
> SSH is 2FA, so the user runs every cluster command; Globus works without it.
>
> Full records of finished work: `git show f4b0c36:TODO.md` (to 2026-10-20: the all-species
> C2d run, the pre-2020 backfill), `git show c91a863:TODO.md` (to 2026-10-13: the C2d
> investigation and its tests, waves 1–2, C4, C5), `git show e619087:TODO.md` (to 2026-10-04)
> and `git show 3fcef6d:TODO.md` (to 2026-09-28).

## Results (2020, C2d backfill, BAM release filter applied)

From 14C (`sector_effects/results_summary.csv`, `results_sectors.csv`, `results_bcr.csv`).
Change = observed − no-industry (negative = birds lost to industry); % is of the no-industry
population.

| Species | Observed | No industry | Change (5–95%) | % | Largest sectors (%) |
|---|---|---|---|---|---|
| CAWA | 4.47 M | 5.15 M | −683k (−782k, −564k) | −13.3 | roads −5.7, pasture −3.0, crop −2.7 |
| OVEN | 37.65 M | 41.27 M | −3.62 M (−4.14 M, −2.75 M) | −8.8 | crop −3.3, pasture −3.0, built −1.3 |
| OSFL | 10.74 M | 11.11 M | −372k (−537k, −217k) | −3.3 | roads −2.1 |
| LEYE | 5.79 M | 5.83 M | −41k (−61k, −25k) | −0.7 | roads −0.6 |
| GRSP | 0.96 M | 0.90 M | +59k (−293k, +408k) | +6.6 | roads +7.1, pasture +2.6, crop −2.5 |
| BOBO | 6.76 M | 3.63 M | +3.14 M (+2.96 M, +3.33 M) | +86.4 | roads +30.3, crop +26.6, pasture +20.2 |

- **All 138 species** (sorted by %): 88 decline (72 with the whole 5–95% interval below 0), 50
  increase (38 above 0), 28 intervals span 0; median −1.7%. Largest declines: EATO −62.2%,
  FISP −43.3%, GCFL −30.7%, RBGR −23.7%, PIWA −21.7%. Largest increases, all open-country or
  human-associated species: BARS +156.9%, ROPI +151.7%, YHBL +112.7%, COGR +104.4%, HOSP
  +104.1%. The sector with the largest |%| is roads for 98 species, crop 33, pasture 5, built 1,
  reservoirs 1. EATO's and FISP's no-industry populations are 2.7× and 1.8× today's, mostly in
  can13 (southern Ontario; EATO also 91k in can11); not examined.
- **Withheld by BAM** (14B drops them): CAWA can40; GRSP can10, can61; BOBO can10; LEYE can10.
  Every BCR of BGGN, BWWA, DUNL, GWWA, NESP, SPGR, YBCU and YTVI, so 14B skips those 8.
- **Direct vs vegetation** (Open Limitation #9; 14B's sign, cf − obs): CAWA direct +138k,
  vegetation +545k. OVEN direct −1.68 M, vegetation +5.30 M; its roads are −1.27 M direct and
  +1.74 M vegetation, hence their interval spanning 0.
- **10C:** 4 of 674 subbasins flagged (50, 72, 430, 481); they carry 22k for CAWA, −71k for
  BOBO and < 4k for the others (14B's sign).
- **Sensitivity runs (pre-C2d, CAWA and OVEN; not rerun, user 2026-10-14):** C3, without the
  frozen q99 cap: v(N) +1.16% (CAWA), +0.05% (OVEN), no sector share moves. C4, densities
  < 0.005 birds/ha set to 0: national v(N) −2.7% (CAWA) and −0.5% (OVEN), but CAWA can11 −81%
  (the prairie floor) while can13 and can61 rise. Production is the headline.
- Logs: `cluster_logs/07_2020_c2d/`, `c2d_6spp/`, `c2d_allspp_2026-10-14/`. Old-backfill
  baseline per species × BCR: `cluster_logs/allspp_2026-10-09/baseline_old_backfill_species_bcr.csv`.

## Species lists

- `Rscripts/species_todo.txt`: 146 species, 1,624 species × BCR tasks in 2020. **BOBO, CAWA,
  GRSP, LEYE, OSFL and OVEN are commented out** (they ran first), so uncomment them before
  using it for a full run. Left out: DEJU, DOWO, EUST, FOSP (bootstrap models but no 2020
  `07_predictions` on Fir) and ROSA (no Canadian models).
- `Rscripts/species_113.txt`: the 113 reported species (14D's `include`).
- 12H needs ~0.027 GB per per-BCR file, so its 32G holds a wave of ~1,000.

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
- **Fir's job IDs wrapped** between 2026-10-22 and 10-23 (Slurm's ceiling is 67,043,328): wave
  1 ran as 66836899–913, wave 2 as 74261–74308. Pull logs by the job IDs sbatch printed, never
  by "IDs above the last one". A filter's `[3-9]` is not a range: `slurm-668[3-9]*` missed
  6684…–6688….
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
  19A = `19A` `.R`+`.sh` (reads 14B's and 14C's CSVs and 12H's Shapley samples, already on Fir).

---

## Open before the write-up

- [ ] **Ecological check of the farmland BCRs (can11–13).** Forest species gain more birds there
      without industry than they have today (no-industry birds on the all-sector footprint per
      bird today, pre-C2d → C2d): OVEN can11 4.7 → 4.9, can13 2.5 → 2.1; CAWA can11 3.7 → 2.6,
      can13 5.1 → 4.3; OSFL can13 12.9 → 5.9. The user accepted them for the all-species run
      (2026-10-13), to be handled here. C2d removed GRSP's can11 decline, a BART artefact (no
      blank canopy CV, so every backfilled pixel looked treed). The user rejected a kNN hot-deck
      backfill (2026-10-12: "a new method on top of BART that is a bit arbitrary"). OVEN's
      negatives in can12, can80 and can81 come from roads through the bird model's footprint
      response (Open Limitation #9).

## Loose ends

- [ ] **12A's memory grows with every BCR-year a task does** (2026-10-22, `66813898`). With
      `YEARS` set, one task runs up to 19 BCRs × 5 years; 39 of 113 were OOM-killed at 16G after
      58–69 BCR-years (the 2020 run never did more than 19), whatever the BCR. Resumable (each
      BCR-year is saved and skipped on a rerun). Likely GDAL's block cache, which defaults to 5%
      of the node's RAM (~37 GB on Fir) and ignores terra's `memmax`; untested. If 12A is ever
      run for many years again, cap it (`terra::setGDALconfig("GDAL_CACHEMAX", "1024")`) or give
      it 32G.
- [ ] **Anna's extreme densities (email 2026-10-22) do not reach our tables.** Her land-change-only
      stacks (CLSW can14 max 1.9e7, LISP can14 1.9e20, WAVI can5 flat at 26.7 birds/ha) are new
      predictions that never passed V5's `10.Truncate.R`, and a 99% cap taken from the stack
      itself fails once more than 1% of it is extreme. We apply that script's two caps (densmax,
      then the frozen observed q99.9) to obs, d0 and bf; V5 `origin/main` `f50a6b8` (2026-10-08)
      has not changed them since `f082866`. C2d caps (birds/ha): CLSW can14 0.134, LISP can14
      0.212, WAVI can5 1.16; q99 < densmax in all three, so < 0.1% of Elly's observed pixels
      exceed them. Open (reporting only): Anna and Diana drop GWWA can13 (already withheld),
      CLSW can14, LISP can14, WAVI can5 and the irruptives PISI, EVGR, PIGR, WWCR, RECR; all of
      those are `include = TRUE` in 14D. Dropping pairs = add them to 14B's `withheld_models`.
- [ ] **BART and kNN still differ where training land is thin or far (C2d test 2).** BOBO can13
      subbasin 98 keeps 0.36 of its observed birds vs kNN 0.23 (+177k; more SCANFI herbs and
      shrub; it trains on can12 land and 59% of can13's birds are outside the AOA). GRSP can11
      subbasins 55 and 60 keep 0.73 and 0.55 vs kNN 1.12 and 1.18 (prcD_5x5 above kNN; likely
      the Gaussian amount model's missing point mass at 0). Neither is a coding defect, and both
      methods agree on the direction. Report it as method sensitivity, with 10C's AOA flags.
- [ ] Fir scratch, for the user to decide: `/home/mannfred/scratch/v4_bootstraps/` (196 GB of V4
      models from 2025-08, read by no script here).
- [ ] 12A on Fir (optional check): run it for CAWA into a scratch dir and diff against the
      G:-built `observed_bootstraps.tif` / `truncation_params.rds`.
- [ ] **Cheaper 12D (optional):** 16 distinct BART draws per bootstrap instead of 100 picks with
      replacement would cut a CAWA + OVEN run from ~64 to ~25 billed core-equivalent-hours, and
      move the all-sector mean by 0.08% (CAWA) / 0.18% (OVEN) of the impact.
- [ ] `15C_singletons_plot.R:150` reads `predictions_coalitions/`, deleted on both ends, so the
      map panel is broken. Regenerate via `save_arrays_ids` if wanted.

---

## Workstream codes

Code comments, logs and older notes name finished work by these codes. Full records: the
`git show` commands at the top.

- **A, A6, A7** (2026-09-11 → 09-23): the CAfire fix (Fix A, Fix B; Open Limitation #5),
  backfill and cluster prep. **B**: V5 packaging conformance, gates G1–G4. Harnesses:
  `Rscripts/misc/verify_weight_vs_v5_masking.R` (rerun after any weight change),
  `verify_v5_truncate_port.R`.
- **C1** passed, with the **12D rework** (identity gate, C++ tree walk, per-bootstrap
  reduction, one task per species × BCR), 2026-09-24.
- **C2a–C2e** (2026-09-24 → 09-25): per-sample Shapley sums (C2a), 100 BART draws with
  replacement (C2b), 10C as an area of applicability (C2c), footprint covariates at 0 (C2e).
- **C2f** (2026-09-28 → 10-02): two masks per sector (Open Limitation #10).
- **C2d** (2026-10-09 → 10-13): modelled blanks, no water class in training, drawn land-cover
  classes; the production backfill.
- **C3**: no-q99 sensitivity (`DROP_Q99=1`). **C4**: low-density floor (`DENS_MIN`). **C5**:
  pixel-level change maps (`pixel_change/`, 18's 10 km grid).

## V5 reference

Renumbering in `f082866`: `10.Package.R` → `10.Truncate.R` + `11.Package.R`;
`11.Validate` → `12.Validate`; `12.Summarize` → `13.Summarize`;
`output/10_packaged/` → `output/11_packaged/`, new `output/10_truncated/`.
`06_bootstraps` and `07_predictions` were NOT renumbered, and they are the only V5 output
folders our scripts read.
