# TODO (last updated 2026-11-02)

> **Where we left off (2026-11-02).** 2020 is done for all 146 species on the C2d backfill
> (07 + 11, then 12D/12H and 14B/14C/18 on Fir, 2026-10-14): 113 reported since 2026-10-31
> (14D), CSVs and maps pulled to `sector_effects/` and `output_figures/change_maps/`. 2020 is
> the only year reported: the 1995–2015 outputs are withdrawn and deleted (item 2).
> Open:
> 1. **Supplementary results** for the Google Doc "Supp. Mat. for: Declining boreal birds", from
>    page 25 on, instead of a PDF: `writing/supp_results_docx.R` builds
>    `writing/supp_results_2026-10-14.docx` (each species' map, then its BCR sector table from
>    14C's `results_bcr.csv`; the build takes ~28 min). Left: paste it into the Doc after the
>    Supplementary Methods. The Drive connector reads the Doc but cannot edit its body.
>    2026-10-31: rebuilt for the 113 reported species (`supp_results_2026-10-31.docx`), with the
>    Doc's own Table S1/S2 captions, to replace the 138-species tables in "Extended Data for:
>    Declining boreal birds".
> 2. **Pre-2020 results withdrawn (2026-11-02).** BAM's 1995–2015 bird density models are not
>    validated, so only 2020 is presented. Every 1995–2015 output is deleted: BART backfills and
>    mosaics, 12A/12C predictions, 12D/12H tables, 14B/14C/14E CSVs, 19A's
>    `population_samples.rds`, 18's maps, the 19B–19G figure 1 panels and 17's outputs
>    (`figure1.pptx` kept), and every non-2020 key of each `truncation_params.rds` (12F reads only
>    `{bcr}_2020`; 12A rebuilds any BCR-year whose stack is missing). Done on both ends. Fir: a
>    one-off script deleted 5 + 5 BART year dirs (145 GB), 6,580 prediction dirs, 162,885
>    density-table files, 31 CSVs/rds and 565 maps, and 6,580 q99 keys from 113 files; the 2020
>    counts were identical before and after (674 subbasins, 19 mosaics, 1,316 prediction dirs,
>    32,577 density-table files, 6 CSVs, 113 maps).
>    Kept: the inputs (`covariates_mosaiced_{year}.tif`, `CAfire_{year}_masked.tif`), the logs,
>    and the `YEARS` code, so validated year models can be run later. Until then 17, 19A–19G and
>    20 (which reads 1995–2020) cannot run. The record of the pre-2020 work is under "Withdrawn"
>    at the end.
> 3. **Maps moved to `output_figures/change_maps/`** (2026-10-26): all 703 local maps, the 138
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

- `Rscripts/species_todo.txt`: 146 species, 1,624 species × BCR tasks in 2020; now only 14D's
  candidate list. **BOBO, CAWA, GRSP, LEYE, OSFL and OVEN are commented out** (they ran first).
  Left out: DEJU, DOWO, EUST, FOSP (bootstrap models but no 2020 `07_predictions` on Fir) and
  ROSA (no Canadian models).
- `Rscripts/species_113.txt`: the 113 reported species (14D's `include`), and since 2026-10-31
  `run_species_batches.sh`'s default list. That day the other 33 species' results were deleted
  locally and on Fir (`predictions/`, density tables, `arrays/`, `by_bcr/`, `pixel_change/`,
  maps) and their rows dropped from the 2020 `sector_effects/` CSVs, so every year's tables
  hold the same 113. Their range rasters (`v5_gis/ranges/`) and logs were kept.
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
- [ ] **VATH's industry-positive result is the bird model's direct part, not the backfill**
      (2026-11-03). 14C: +332k (+2.9%), roads, can10 +287k and can5 +73k. 14B split (obs − cf):
      direct +409k, vegetation −77k (199 of 535 subbasins lose birds to the backfilled
      vegetation, none 10C-flagged among the 128 that gain), so BART is not the cause. Local
      reproduction (8 bootstraps, `gbm_tree_walk`) gives obs − d0 = +298k can10 / +93k can5 vs
      14B's +309k / +83k. Zeroing each footprint covariate alone: can10 CanHF_1km +198k,
      CanHF_5x5 +96k, canroad +15–26k, CCNL 0; can5 canroad_5x5 +54k, CanHF_1km +27k. obs/d0 is
      a near-flat multiplier (can10 1.28–1.51, can5 1.12–1.25) across MAT, PAS and mTPI
      quintiles, so it is not tied to elevation or valley floors; it is largest at moderate
      footprint (CanHF 4–7, canroad_5x5 0.2–0.6) and gone at the most intense (can5
      canroad_5x5 > 1.1: 0.96). d0 sits below the 10th percentile of can10's surveys (CanHF_1km
      3.07, canroad_5x5 0.21). Open Limitation #9, not a bug; 58 of 113 species have a positive
      direct road part. For the write-up: report VATH with its split. Scripts:
      `vath_pixel.R`, `vath_env.R` (session scratchpad, not kept). Figure for colleagues:
      `Rscripts/misc/figure_vath_footprint_pdp.R` → `output_figures/direct_part/VATH_footprint_pdp.png`
      (16 bootstraps; direct part +325k can10 / +95k can5 vs 14B's +309k / +83k). In can10 the
      gain is the step from no footprint to a little: scaling every footprint covariate from 0
      to 10% of its observed value gives ~93% of it, a range below the lowest survey decile. In
      can5 it builds gradually up to half the observed values. The format going forward (user,
      2026-11-03): that 0 → observed scaling curve only, `Rscripts/misc/figure_footprint_scaling.R`
      → `output_figures/direct_part/footprint_scaling_{can10,can5,can9}.png`, six western species
      (RBNU, SWTH, RCKI, OCWA, TOWA, VATH), 8 bootstraps, 5,000 footprint pixels per BCR; 14B's
      direct part falls inside the 8-bootstrap 5–95% range for 15 of 18 species × BCR
      (`footprint_scaling_direct.csv`). Positive for most; OCWA can9 negative (−150k).
      **Who the zero-footprint surveys are (VATH can10, 2026-11-05).** From V5's
      `G:/.../data/04_NM5.0_data_stratify.Rdata` and its `bootlist` (the survey each of the 32
      bootstraps drew per 2.5 km cell-year; each trains on 28,517 of can10's 63,540 surveys). The
      model's zero bin (every split puts the survey on 0's side) is essentially exact 0. Rows per
      bootstrap in it: CanHF_1km 482 (1.7%), CanHF_5x5 129 (0.45%), canroad_5x5 1,821 (6.4%),
      canroad_1km 4,164 (14.6%), CCNL_1km 11,457 (40%; 0 does mean no lights, negatives are
      small noise); all five at once 23 (0.08%, 63 locations, VATH detected at 8). All five:
      74% BC Breeding Bird Atlas 2008–12, the rest ABMI, Banff, BU; ARUs 26% of them vs 1.6% of
      all rows. CanHF_1km zero: Atlas 58%, BBS 19%, national parks ~11%, eBird 5%, CANFOR 4%.
      **19% of the CanHF_1km zero bin is 22 roadside BBS stops** (canroad_5x5 > 0 at 98% of
      them, so CanHF_1km misses a road V5's road layer sees), with 0 VATH detections against
      10–12% at other BBS stops: a lead for the zero-footprint drop, untested. Raw (unadjusted)
      Atlas detection: 22.5% at CanHF_1km zero vs 18.0% elsewhere. Scripts (session scratchpad,
      not kept): `zero_fp_surveys.R`, `zero_fp_profile.R`, `zero_fp_bbs.R`. **Correction:** those
      counts cover only `bcrlist` code 1 (inside can10, 28,517 rows per bootstrap); the model also
      trains on code 2 (the 100 km buffer, 28,518 rows), 57,035 in all. With the buffer, bootstrap
      1's CanHF_1km zero bin is 717 rows (479 inside, 238 buffer); the BBS part is still the same
      91 rows at the 22 locations, all inside, VATH detected in none (7.0% at other BBS stops).
      **BBS lead refuted** (Fir job 2683226): bootstrap 1 refitted exactly as V5, 5 seeds per set
      (`cluster_logs/refit_vath/`). V5's own pairing gives direct +317k (seed sd ~9k); without the
      22 stops, +4k more (paired by seed, sd 33k); without 22 random other BBS stops, −2k.
      **The colleague's check: confounding among the bird model's surveys (2026-11-05).**
      Bootstrap 1, all 57,035 rows. Z = CanHF_1km zero bin (717 rows, 632 locations), D = the
      rest, P = the 5,000 can10 footprint pixels. Elevation from `data/raw_data/dem_bc_1km.tif`
      (elevatr, `Rscripts/misc/bc_dem.R`; covers every row and pixel); temperature = DD18
      (ClimateNA degree-days above 18 °C). V5's ERA covariates are left out (item below).
      - *Premise holds.* Z sits higher and colder than D: elevation median 1,146 vs 923 m (SMD
        +0.71), DD18 19 vs 32 (−0.63), precipitation as snow 252 vs 146 mm (+1.18, the largest
        gap; r = 0.43 with elevation); also mTPI +0.72 and conifer cover 5×5 +0.83 (habitat, not
        elevation: r = 0.19). P lies between (1,012 m, PAS 209 mm).
      - *The model's drop is flat across elevation:* 26–31% of observed density in every pixel
        elevation band, from < 600 m (10 Z surveys, 0.14% of rows there) to > 1,500 m (155, 3.5%).
      - *The deficit is in the data, model-free, weaker once elevation is held.* Each Z survey vs
        its 5 nearest D surveys (importance-weighted distance over 29 continuous non-footprint,
        non-ERA covariates, 78% of influence; |SMD| ≤ 0.12 after matching, elevation ≤ 0.04):
        offset-corrected VATH density (count / exp(offset matched by id)) Z / D = 0.54 (90%
        location bootstrap 0.44–0.64); also within ±100 m of elevation 0.63 (0.51–0.75). By Z's
        elevation thirds: low (94–998 m) 0.81 (0.60–1.09), middle 0.44 (0.26–0.67), high 0.62
        (0.47–0.80). Within ±100 m: without the 22 BBS stops 0.63, inside can10 only 0.74
        (0.59–0.90), all five at 0 0.37 (72 surveys; 0.19–0.62), within the BC Atlas only (one
        protocol) 0.71 (0.57–0.86). The model's d0 / obs over P is 0.71. Unmatched, Z detects VATH
        more often than D (13.5% vs 8.6%): Z sits in better thrush habitat.
      - *Reading:* the model's drop matches a deficit the surveys show at middle and high
        elevation. Below ~1,000 m the surveys do not show one clearly, yet the model applies the
        same ~30% there; pixels below 900 m carry ~26% of the figure's direct part. Left: a real
        response at 1 km, or a difference none of the covariates (nor elevation) carries.
        Figure: `output_figures/direct_part/VATH_zero_footprint_surveys_can10.png` (+
        `VATH_zero_footprint_matched_ratio.csv`, `VATH_zero_footprint_drop_by_elevation.csv`).
        Scripts (session scratchpad, not kept): `elev_extract.R`, `confounding_vath.R`,
        `confounding_vath_elev.R`.
      - *Anna's email points (2026-11-05).* Season: every can10 training survey inside the BCR is
        1 June–15 July (jday 152–196), Z and D alike (buffer rows not in the cached table); raw Z
        detection is lower only in 1–15 June (6.5% vs 9.5%), higher later (late arrival up high?
        unmatched). Her predicted order of the drop (can9 > can10 > can5) holds for none of the six
        (`footprint_scaling_curves.csv`, x = 1 vs 0): VATH can10 27% > can5 17% > can9 0.7%; SWTH
        can9 27.5%, can5 19.6%, can10 19.4%. Water layers (SurfaceWater, WetOccur/Recur/Season,
        LakeEdge) are abiotic, held at observed values; only Peatland is backfilled. VATH can10 has
        no river covariate (mTPI 1.3%, WetOccur_5x5 1.2%, Peatland 0.9%, LakeEdge 0.9% of
        influence). Untested: navigable waterways, which HP scores into CanHF and which is not one
        of our 8 sectors (nor is forestry). Scripts: `email_checks.R`, `email_dates.R` (scratchpad).
- [ ] **V5's bootstrap models pair offsets with the wrong surveys (UPSTREAM, found 2026-11-05;
      not yet reported to BAM; email drafted 2026-11-05).** Line numbers: V5 `main` a53c6ab.
      **Cause (proven 2026-11-05):** `02:220` leaves `offsets` in `visit` order; `03:826-832`
      (`st_intersection(can)` then `full_join(visit)`) moves non-Canadian surveys to the end of
      `visit` and `offsets` is not re-sorted (saved `03:926`); `04:128-139` filters both; `06:105`,
      `05.Tune:145` and `12.Validate:127` take offsets in their own order and pair by position.
      In `04_NM5.0_data_stratify.Rdata`, `visit$id` is `identical()` to `offsets$id` with the
      non-Canadian rows moved to the end; `cov` and `bcrlist` match `visit`. (An earlier note here
      blamed 02's `rbind`: wrong, 02:220 undoes it.) A few border surveys suffice, each shifting
      every later row. **Scope (bootstrap 1, share of rows paired with another survey's id):** 13
      of 19 Canadian BCRs, 56.4% of Canadian training rows: can5 99.9%, can81 98.6, can14 98.4,
      can12 98.1, can82 95.9, can71 88.8, can70 84.6, can13 81.8, can40 81.7, can3 80.6, can72 79.1,
      can9 26.4, can10 25.3 (25.3–25.4 over all 32); 0 in can11, 41, 42, 60, 61, 80. 13 of 14 US
      BCRs too. cor(assigned, own) for VATH 0.15 (can5) to 0.80 (can10), mean |diff| up to 0.33 on
      the log scale, so can10 is one of the mildest. Also reaches tuning and the validation AUCs
      (`12.Validate:152-155`). Fix: one line at `04:364`,
      `offsets <- offsets.use[match(visit$id, offsets.use$id), ]`.
      VATH can10 bootstrap 1: 13,260 of 57,035 rows (23%) carry another survey's offset,
      from row 42,607 on; eBird 96%, WildTrax 19%; cor(assigned, own) 0.80, mean |diff| 0.077 on
      the log scale. **Confirmed that V5 fitted this way:** refitting with V5's exact construction
      and `set.seed(1)` (gbm 2.2.2) reproduces the stored model's structure for the first 450 trees
      (leaf values differ by ~1e-18, then a near-tie splits differently); with the offsets matched by
      id, tree 1's leaves differ by 7e-4. Effect, one fit each, on the figure's 5,000 can10 footprint
      pixels: observed 9.53 → 10.66 birds/km² (+12%), direct part +303k → +387k, CanHF_1km alone
      +250k → +260k. On Fir, 5 seeds each (job 2683226; the V5-pairing seed-1 fit reproduces
      V5's stored bootstrap 1 exactly on this evaluation): observed 9.56 → 10.64 (+11%), direct
      +317k → +396k, paired by seed +79k (+68k to +88k in every seed; seed sd ~9k). The VATH drop
      is not an artefact of the misalignment: matched offsets enlarge it. Effect measured only in
      VATH can10; not yet in a heavily affected BCR. Scripts: `refit_bbs.R`, `compare_trees.R`,
      `s5_offsets.R` (session scratchpads); fits cached in `scratchpad/refit/`.
      **Observed-only test, can81 (2026-11-11):** `Rscripts/misc/PR1_offset_refit_observed.{R,sh}`
      (Fir; `04_NM5.0_data_stratify.Rdata` staged to `data/raw_data/`, gitignored). OVEN, SWTH,
      REVI; 8 bootstraps × {V5 pairing, matched by id}, V5's tuning kept; predicted on V5's cells,
      12B truncation (each arm's own q99), × weight.tif, summed. No counterfactual. User chose
      boreal; can61 would show nothing (0% mismatch: no US surveys in its table). Out:
      `offset_test/summary/PR1_can81_*.csv` (copies in `cluster_logs/pr1/`).
      **Result (jobs 3628277/8, 2026-11-12):** the V5-pairing arm reproduced V5's stored trees
      in all 24 bootstraps and its raw predictions to float32 (max rel diff 6e-8), so the arms
      differ only in offsets. cor(V5 offset, own): OVEN 0.40, SWTH 0.35, REVI 0.09. Matched −
      V5, can81 population (8 paired bootstraps): OVEN +31k (+0.97%, 8/8 up), SWTH +187k
      (+0.84%, 7/8), REVI +121k (+2.35%, 8/8); each about a third of the bootstrap 5–95% width
      or less. Maps: cor 0.987–0.999, but OVEN's median cell falls 8% (41.5% of cells down
      > 10%) while the total rises: the low-density cells drop. Each arm's own q99 is 4–8%
      higher when matched. Not tested: the effect on the obs − cf contrast.
      **can5 (jobs 3651813/4, VATH v5 rerun 3655832/3, 2026-11-12), VATH, RCKI, TOWA, same
      design:** 99.9% of rows carry another survey's offset in every bootstrap; cor(V5 offset,
      own) 0.14–0.24, mean |error| 0.28–0.32 ≈ the own offsets' SD (0.29–0.40), i.e. near
      noise. V5-pairing trees identical to V5's stored in 24/24, raw predictions to 6e-8.
      Matched − V5: VATH +141k (+4.9%, 5/8 up, paired 5–95% −226k to +473k), RCKI −27.6k
      (−5.7%, 3/8 up, −99k to +23k), TOWA +730k (+11.0%, 8/8 up, min +163k). Maps cor
      0.975/0.961/0.983; 69–74% of weighted cells move > 10% either way; q99 +4% (VATH), +11%
      (RCKI), +2% (TOWA). VATH's v5 arm hit CVMFS I/O errors twice (nodes fc30564, fc30613)
      before a clean rerun on another node.
- [ ] **Ten V5 training covariates are attached to the wrong surveys (UPSTREAM, found 2026-11-05;
      not yet reported to BAM; email drafted 2026-11-05).** The eight ERA layers, CCNL_1km and
      MODISLCC_1km: every temporally matched GEE layer with `Extraction == "point"`.
      **Cause (proven 2026-11-05):** `03.ExtractCovariates.R` extracts by year (`:676`), stacks the
      year blocks (`:769`), then `cbind()`s them by position (`:770`) onto `loc.gee`, which is not in
      year order. Radius layers join by id (`:796-801`) and are fine. `03_AK:564-566` has the same
      code (untested). Committed 03 does not parse: the j loop's `}` is missing after line 766.
      **Evidence:** all 50 covariates of VATH can10 vs the stack at the survey (2005/10/15): these
      ten r −0.20 to 0.08, the other 40 r 0.30–1.00. Undoing the sort (loc.n order from
      `03_NM5.0_data_covariates_GD.csv`), 2000–2015: ERA r 0.96–1.00, ~98% exact bar ERATavesmt;
      CCNL r 0.98–0.99; MODISLCC agreement 14–15% → 70–74%. ERAMAT recovers only on the GEE-match
      CSV's own row order (r 1.000, 98% within 0.05 K), so it was bound in a separate pass.
      ERATavesmt keeps a median −1.3 K vs the stack (lookup has June–September, its siblings
      June–August; cause unconfirmed). Training table = that CSV on 100% of 57,035 rows. 3,312,270 of
      3,312,336 rows moved; an ERA value belongs to a survey a median 1,251 km away (0.2% within one
      ERA cell). 07 predicts on true values. Relative influence of the ten: median 10.7% (5–95%
      1.6–24.3) over 867 Canadian models in `bam_predictor_importance_v5.rds` (70 species, Feb 2026,
      partial); VATH can10 7.0%; extremes PUFI can61 95%, RUBL can12 51%, PISI can82 41% (CCNL).
      Effect on predictions not measured (needs a refit).
      **Invalidates, in the VATH item above:** every training-side CCNL_1km figure (its 40% zero
      bin, "all five at once" 23 rows, "all five at 0" 0.37). CanHF and canroad figures stand.
      Scripts: `s3_unscramble.R`, `s4_audit.R`, `s6_final.R`, `s8_influence_can.R` (scratchpad).
      Earlier symptoms: VATH can10 bootstrap 1, survey years with a local
      stack (2000, 2005, 2010, 2015, 2020; 430–3,902 rows each): training ERAMAT_1km,
      ERATavesmt_1km and ERAMAP_1km correlate with `covariates_mosaiced_{year}.tif` at the
      survey's V5 `visit` coordinates at r = −0.15 to 0.25 (mostly ≈ 0), while DD18, PAS and
      SCANFIprcC_5x5 match at 0.98–0.995. Training ERAMAT is uncorrelated with elevation, latitude
      and survey date (the stack: −0.52 with elevation, −0.40 with latitude), varies by 3.4 K (sd)
      between years at one location, and is wider and warmer than the stack (2010 median 277.0 vs
      275.6 K); 1.4% of 2010 rows are within 0.05 K of their own stack value. ERA covariates carry
      ~4.7% of VATH can10's relative influence: noise in training, real values at prediction.
      Scripts: `era_check.R`, `era_check2.R` (session scratchpad).
- [ ] **BC valley bottoms: backfill vs d0 and vs intact land, six western species** (2026-11-05).
      Colleagues: valley bottoms are backfilled by models trained on high intact land, reducing
      birds. BC = 149 (subbasin, BCR) units in can5/9/10/40/60/61, 149,717 km² of direct footprint
      (99.5% inside BC). **Premise holds:** direct footprint median 848 m vs training 1,175 m (can5
      252 vs 1,071); 92% below its own subbasin's training median, 41% below the 5th percentile,
      5.2% below every training pixel. **bf vs d0** (= 14B's vegetation part): RBNU +6.7%, SWTH
      +10.3%, RCKI +6.0%, TOWA +17.4%, VATH +5.8%, OCWA −5.9%; all intervals exclude 0. This
      cannot clear the backfill: d0 is cleared land. **bf vs intact land in the same subbasin and
      100 m band, both at footprint 0** (birds/km²; 91% of pixels matched): OCWA +8.6%, RBNU +6.9%,
      RCKI +10.6%, SWTH +3.6%, TOWA −9.5%, VATH −23.3% (bf 7.4, d0 7.1, intact 9.6). VATH's gap is
      widest in valley bottoms (−28% below the 5th percentile, −16 to −18% above the median; can5
      −36% to −3%); TOWA's is not (−5% to −23%, widest highest up). In birds, intact minus bf is
      +255k VATH and +275k TOWA, against national changes of +332k and +236k, the two
      industry-positive species; −98k to −234k for the other four, which is small beside their
      changes (OCWA's +127k excepted). A what-if, not a correction: intact valley land is a
      selected sample (RBNU can5's d0 is already 51% above it). 9% of pixels (19% of the lowest
      band) have no intact land in their band; they carry 32% of VATH's bf − d0.
      **Trap:** observed density on intact land is not comparable with bf. Only 16–32% of low-HF
      land has CanHF_1km == 0, and the footprint step inflates observed density 8–21% for five of
      the six (OCWA 2%), which made VATH look −32%. Predict intact land at footprint 0.
      `any_elev` is context only: density up the slope reflects climate as well as vegetation.
      **Next:** compare backfilled canopy height, closure and conifer share with intact land at
      the same elevation, where VATH falls short; it needs the BART mosaics, so on Fir.
      Scripts `Rscripts/misc/bc_{dem,d0_pixels,d0_tables,intact_d0,intact_benchmark}.R`; tables
      `output_figures/bc_backfill/`. Re-predicted d0 reproduces 12F (0 pixels off in 35 species ×
      BCR; bf − d0 sums match 14B's samples to 1.5e-9).
- [x] **BAM releases V5 for 2020 only; 1995–2015 "by request"** (2026-11-02). **Why (meeting
      minutes, user, 2026-11-02):** "comparisons with independent population trend estimates
      suggest that the models underestimate the magnitude of both positive and negative trends.
      Although the reasons remain uncertain, this limits confidence in using model predictions to
      estimate absolute population change over time." The group judged counterfactual comparisons
      that hold other conditions constant more defensible. **Decision:** 2020 observed vs
      no-industry only; historical trajectories are not a primary analysis. Pre-2020 outputs
      deleted ("Withdrawn" below). The minutes are the only record: V5's code, its release
      workbook (`release` tab: 2020 only), BAMexploreR's 2026-07 commits and the V5 drive state
      no reason.
      **For the write-up:** if the trends are attenuated because the bird models' habitat
      responses are (BRT shrinkage, or covariate error at 1 km), our counterfactual differences
      go through the same responses and are attenuated too: the industry effects would be
      conservative. Attenuation from drivers outside the covariates, or from year-to-year
      covariate inconsistency, cancels in a one-year obs − cf contrast. One sentence in the
      Limitations.

## Loose ends

- [ ] **BART does not mix at 08B's settings (2026-11-04; decision open).** 8 production pairs
      (subbasin × continuous covariate, n 1,328–11,107, 62–81 predictors), 4 chains of 300 + 1,700
      sweeps at 50 and 200 trees: `Rscripts/misc/bart_production_settings.R` and `_summary.R`, CSVs
      in `data/derived_data/rds_files/bart_production_settings/`. Replaces the 5-predictor benchmark
      (`bart_vs_xgboost_benchmark.R`, subbasin 388), whose "200 trees: R-hat 1.01" and "about 10
      effective draws" do not carry over.
      - **Mixing.** Median pixel R-hat 1.92 (50 trees), 1.69 (200); none of 20 arms is under 1.05.
        Sigma is still falling at sweep 2,000. By between-chain ESS, one chain's 100 stored draws
        are worth 1.3 (50 trees) to 1.6 (200) draws, and thinning, longer chains and a longer
        burn-in do not change that. Their SD is 61–72% of the posterior's.
      - **Trees.** 200 vs 50: held-out RMSE −9% (7 of 8 pairs) for 2.4× the time, but −5% on a
        footprint-shaped holdout and −3% with cascaded inputs (2 of 7 pairs worse, one by 66%).
        The backfill mean moves by less than chain noise in 7 of 8 pairs.
      - **Chains.** The only setting that helped all 8 pairs: 4 pooled chains, RMSE −7%. Effective
        draws add per chain. 100 kept sweeps after burn-in 500 score the same as 700 after 300, and
        are worth 1.12 effective draws per chain.
      - **`mbart2` is fine.** 3 chains at 08B's settings, SCANFI_1km, subbasins 9 and 645: pixel
        R-hat 1.02–1.05, class-share R-hat 1.04–1.06, 13–16 effective draws. Leave it. `pbart`
        (the blank model) untested.
      - **Validation.** Footprint-shaped holdout vs 08B's random one: RMSE ×1.5 with observed
        inputs, ×1.7 cascaded; R² 0.89 → 0.77 → 0.47; 90% interval covers 90% → 81% → 63%. (The
        cascaded column uses a rough 100 + 100-sweep cascade, without blanks, so it overstates.)
        PSIS-LOO from one fit's `yhat.train` tracks the real holdout (median ratio 0.97; 8 s), but
        shares the random holdout's blind spot. 100% vs 90% of pixels: RMSE −2% (median).
      - **Effect on results.** BART is 1.4% of the variance of national change (median species;
        q90 22%, NOCA 56%). A subbasin's chains differ by 4% of its backfill mean: 7–13% of the
        contrast where contrast is large, more where it is near zero. That error is missing from
        the reported intervals. Rough bound (vegetation part, independent across subbasins): SDs
        widen 2–13% for the median species, ≥ 10% for 13–62 of 113. Not measured at bird level.
      - **What fixes mixing (2026-11-05).** `Rscripts/misc/bart_mixing_fixes.R` and `_summary.R`,
        4 of the 8 fits, 4 chains; CSVs `mixing_fixes`, `footprint_learning_curves`, `chain_length`.
        No tweak of `gbart` does (median R-hat; 1.84 today): sparse off 1.83, one predictor per
        |r| > 0.9 cluster 1.70, no coordinates 1.70, `stochtree` warm start 1.94, 1,000 pixels 1.51
        (RMSE +16%). `gbart` for 23,000 sweeps stays at 1.67 (2 fits): random-holdout RMSE falls
        22% but footprint-shaped RMSE rises 4% (10% with cascaded inputs), and sigma keeps falling
        (0.150 to 0.036 on the Peatland fit, where x and y take 47% of the splits against
        `dbarts`' 9%). It is memorising the training land. So the sparse prior's and long runs'
        random-holdout gains are local detail that a backfill cannot use.
        `dbarts::bart2` at 200 trees, no sparse prior, burn-in 3,000 (3 fits, n ≤ 3,762): R-hat
        1.10–1.39 at 5,000 kept sweeps, 1.06–1.31 at 18,000; per chain 2.4–8 effective draws at
        5,000, plateauing near 3.5 at n = 3,762. Footprint-shaped RMSE about −5% with observed
        inputs, −1 to +1% cascaded, and flat from sweep 1,000 on. 0.43× `gbart`'s time per sweep
        (it skips test predictions between stored draws). `gbart` with sparse off was not run long,
        so `dbarts`' extra moves and the missing sparse prior cannot be separated.
      - **Proposed (supersedes 10 short `gbart` chains; decision open).** In `deploy_gbart()`:
        `dbarts::bart2`, 200 trees, 4 chains of burn-in 3,000 + 5,000 kept, on all the pixels, 25
        stored draws per chain (100 layers, so 11 and 12F are unchanged); PSIS-LOO replaces the
        holdout metric, and R-hat is written for every fit. 9–32 pooled effective draws. About 14×
        `gbart`'s time, 6× 07's; S57 about 6.4 h on one core, about 1.8 h with a chain per core
        (`n.threads` untested). Backfill means will move (Jack pine fit 2.38 to 2.78, observed
        1.59), so results change, not only intervals. R-hat stays near 1.3–1.4 at n ≈ 3,800 and
        nothing above that was run at 200 trees, so first pilot about 10 subbasins on Fir across
        sizes, S57 and S98 included, for time, memory and R-hat. `dbarts` is not on Fir yet.
        Follow-up numbers: `Rscripts/misc/bart_production_settings_followup.R`,
        `bart_production_mbart_check.R`, `bart_mixing_followup.R`.
- [ ] **impactexploreR still shows 1995–2020** (deployed 2026-11-02; left alone by the user's
      choice when the pre-2020 outputs were deleted). Its bundle keeps the 565 pre-2020 maps and
      their `results.rds` rows. To make it 2020-only, 20 needs its `years` changed (it stops
      without the deleted pre-2020 CSVs) and the app's year slider hidden for a single year.

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

---

## Withdrawn: pre-2020 results (2026-11-02)

Moved from "Where we left off" when the 1995–2015 outputs were deleted (item 2 there). The
numbers below came from unvalidated year models and must not be reported.

2. **Pre-2020 densities (1995–2015) for the 113 reported species** (`Rscripts/species_113.txt`:
   species_todo.txt's 146 minus BAM's 8 fully withheld and the 25 wetland/grassland/open
   species). 12A, 12C, 12D and 12H read `YEARS` (staged 2026-10-22, whole closure
   checksum-synced). Each year uses its own observed stack (`{bcr}_{year}.tif`), Elly's year-y
   predictions, a q99 frozen from year y's observed landscape (as V5 truncates each year),
   year y's backfill, and the 2020 footprint (as 07 did). 12E, 12F and 12G are unchanged, so
   `code_md5` matches the 2020 run.
   - Done: `prep`: 12A `66813898` (39 tasks OOM-killed at 16G, resumed at 32G as `66820731`),
     12C `66813899`; logs and Fir's `truncation_params.rds` in `cluster_logs/pre2020_prep/`.
     1,316 species × BCR per year (6,580 in all), each with Elly's prediction tif, a weight and
     a finite q99. q99 / 2020's q99: median 1.00, 5–95% 0.61–1.47; the 284 BCR-years beyond ½×
     or 2× are irruptive or budworm species (EVGR can80 1995 42.6×, BBWA can80 2000–05 ~5×),
     likely real change in Elly's year models, not a defect.
   - Done: `smoke` (can42; AMRO, ATSP, CAJA; 2 bootstraps; 12D / 12H: 1995 `66834279` /
     `66834280`, 2000 `66834281` / `66834282`, 2005 `66834283` / `66834284`, 2010 `66834285` /
     `66834286`, 2015 `66834287` / `66834288`; logs in `cluster_logs/pre2020_smoke/`). All 15
     tasks: identity gate 9/9, complete weight > 0 footprint 100%, `nice.`; all five 12H `nice.`.
     can42 has only 9 footprint pixels, but their observed birds differ between years (CAJA
     58.7 in 1995 vs 102.3 in 2005, same 2 bootstraps), so passing the gate every year shows
     each year read its own stack. `code_md5` equals the 2020 production files'.
   - Done: 12D + 12H for all 113 species × 5 years (2026-10-22 to 10-24). Wave 1 (78 species AMCR–PUFI,
     891 tasks per year; 12H 1995 `66836901`, 2000 `66836904`, 2005 `66836907`, 2010
     `66836910`, 2015 `66836913`); wave 2 (35 species RBGR–YRWA, 425 per year; Fir's job IDs
     had wrapped: 12H `74263`, `74266`, `74302`, `74305`, `74308`). Logs in
     `cluster_logs/pre2020_full/`. All 6,580 tasks `nice.`, each species × BCR × year once, no
     OOM or timeout. Identity gate all 32 bootstraps at n/n pixels in 6,565; the other 15 are
     BAWW, BBWA and RBGR can41 × 5 years (no weight > 0 footprint pixels, impact 0, as in
     2020). Every other task kept 100% of its weight > 0 footprint pixels. All ten 12H `nice.`,
     no warnings: 255 tables + Shapley samples for 113 species in every year.
   - Done: 14B, 14C and 18 run one year each (`YEARS`, default 2020; staged 2026-10-24, 12E
     unchanged so `code_md5` holds). 14B's .sh submits one 14B + 14C and one 18 per year.
     2020 keeps its file names; other years' CSVs and maps carry `_{year}`, so 14D and
     `writing/` still read 2020. Each year's map uses its species' 2020 colour scale. Tested
     locally on the six species' 2020 tables relabelled 1995: all six CSVs identical to
     production 2020, maps identical but for the title.
   - Done (2026-10-24): 14B + 14C / 18 for 1995 `289752` / `289753`, 2000 `289754` /
     `289755`, 2005 `289756` / `289757`, 2010 `289758` / `289759`, 2015 `289760` / `289761`;
     logs in `cluster_logs/pre2020_report/`, the 30 CSVs in `sector_effects/` and 565 maps in
     `output_figures/` pulled. Every year: 113 species, sample means = tables' Shapley values
     (max diff 2.2e-9), per-sample additivity residual 0, 14C's checks passed, 113 maps with
     every pixel raster present and each grid within 0.5 bird of 14B's national change.
     Results: each year's species % changes correlate 0.96 with 2020's; median % −3.4 (1995),
     −2.7, −2.3, −3.0, −2.8, −2.6 (2020); 61–66 species with the whole interval below 0 each
     year. CAWA −15.7% (1995) to −13.3% (2020), OVEN −10.5% to −8.8%, OSFL −4.9% to −3.4%.
     Outliers, not examined: EVGR's 1995 observed population is 5.2× 2020's (irruptive; see
     the q99 note above); WITU (41–66k birds) is +1 to +21% before 2020, +93% in 2020.
3. **Rosenberg et al. (2019) Fig. 1 style figures, one file per panel** (2026-10-26). Their
   panels A–C sum ALL species of a group (increasing and declining); D is the share declining.
   `19A` (Fir) writes `sector_effects/population_samples.rds` (per-sample observed and
   no-industry population, species × year), since 14C's percentiles cannot be summed or
   differenced; tested locally on the six 2020 species (means and every 5/95 percentile equal
   14C's). `19B` (local) plots panel A, net change since 1995, observed and without industry,
   each from its own 1995 total. Done: 19A job `603212` (7 min, every check passed; log in
   `cluster_logs/rosenberg_19A/`), rds pulled, 19B → `output_figures/net_change_1995_2020.png`.
   2020 vs 1995: observed +58 M (5–95% +12 to +99 M), no industry +66 M (+12 to +113 M). Per
   sample, no industry − observed is −20 M (1995) to −12 M (2020), every interval below 0;
   the trends differ by +8 M by 2020 (+0.1 to +15 M), within 0 in every earlier year.
   Done: panel B, `19C` → `net_change_by_habitat_1995_2020.png` (Rosenberg's B is one
   trajectory per biome; S1 in the manuscript PDF shows the layout). Birds industry costs
   (no industry − observed), 1995 → 2020: conifer +54 → +45 M, deciduous +28 → +27 M,
   mixedwood +6 → +7 M, forest no type −2 M, non-forest −106 → −90 M; so the near-zero total is
   forest losses cancelled by non-forest gains. Trend differences by 2020: conifer −9 M (−14 to
   −6), non-forest +17 M (+12 to +21), the rest within 0.
   Done: panel C, `19D` → `percent_change_by_habitat_1995_2020.png`. 1995→2020 observed / no
   industry: conifer −4.7 / −5.1%, mixedwood +14.7 / +14.7%, deciduous +8.1 / +7.0%, forest no
   type −7.4 / −7.1%, non-forest +9.0 / +12.0%. Per sample (no industry − observed): conifer
   −0.4 points (−0.8 to −0.2), deciduous −1.1 (−1.6 to −0.7), non-forest +3.0 (+2.2 to +3.7),
   mixedwood and forest no type within 0.
   Done: panel D, `19E` → `species_declining_by_habitat_1995_2020.png` (declining = mean change
   < 0; solid = whole interval < 0). Observed / no industry: conifer 12 / 13 of 27 (8 / 7 sure),
   mixedwood 2 / 2 of 10 (none sure), deciduous 4 / 4 of 22, forest no type 2 / 2 of 6,
   non-forest 18 / 15 of 48 (12 / 10 sure). Decline only as observed: AMRO, EAKI, EAPH, ROPI,
   TRES; only without industry: ATTW, HOWR, WEWP. All four panels done.
   From 14C's means: the 113-species total is 3.02 bn (1995) to 3.08 bn (2020) observed, and
   the no-industry total is LOWER every year (by 20 M in 1995, 12 M in 2020: AMRO, BHCO, AMGO,
   SOSP, COGR, BARS, HOSP gain 120 M on the footprint, the rest lose 108 M), so the two lines
   nearly overlap. The +104 M observed jump 2015→2020 is 45 M PISI (irruptive; Anna and Diana
   drop the irruptives, Loose ends).
5. **How much of the 2020 footprint existed earlier (`17`, rewritten 2026-10-27).** Every
   pre-2020 year ran on the 2020 footprint. The vegetation part corrects itself (SCANFI and
   VLCE are year-specific, so land not yet cleared is backfilled forest for forest); the
   direct part does not (V5's `CanHF` is 2020's in every year's stack: −65 M every year, 113
   species, 14B's sign). 17 measures f_y = summed canopy-height deficit against the intact
   surroundings (11 km) over each sector's 2020 direct footprint, year y / 2020; SCANFI and
   VLCE treed share as checks, `canroad` for roads. Gate passed: earlier SCANFI years are
   separate maps (3–6% of heights identical to 2020) and see dated clearing that is neither
   fire nor harvest (oil sands mines treed 0.83 → 0.61 from 2005; Eastmain-1 water 0.29 →
   0.48 by 2010; Romaine is not in HP's 2020 dams). f_1995 (height): built 0.98, crop 0.94,
   pasture 0.93, rail 0.96, reservoirs 0.99, roads ≥ 1 (`canroad` 2001 vs 2020: 1.03),
   mines 0.59, oil_gas 0.59. Direct part × f_y per sector × BCR: 1995 gap −20.3 → −18.0 M;
   1995→2020 trend difference +8.1 → +5.7 M (+4.5 to +7.3 M over the three indicators and a
   2015 reference, since SCANFI height steps down 2015→2020 in built, roads and rail);
   observed change +58.3 → +60.6 M. Means only: the +0.1 to +15 M interval was not
   recomputed. Not covered: footprint gone by 2020 (abandoned farmland), roads before 2001.
   The 1995–2015 mosaics are now local too (8.3 GB).
