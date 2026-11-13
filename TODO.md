# TODO (last updated 2026-11-13)

> **Where we are.** 2020 is done for all 146 species on the C2d backfill; 113 are reported (14D).
> CSVs are in `sector_effects/`, maps in `output_figures/change_maps/`. Only 2020 is reported
> (pre-2020 outputs deleted 2026-11-02). The project is in write-up: log small new problems under
> Loose ends and do not reopen finished stages for them. `writing/` is private (never commit or
> push it); `rds_files/` and `sector_effects*/` CSVs stay local.
>
> Full records of finished work: `git show 30bf904:TODO.md` (to 2026-11-12: the pre-2020 runs and
> their deletion, Figure 1, 17, the full VATH / PR1 / PR2 / BC / BART-mixing investigations),
> `f4b0c36` (to 2026-10-20), `c91a863` (to 2026-10-13), `e619087` (to 2026-10-04), `3fcef6d`
> (to 2026-09-28). The pre-cleanup CLAUDE.md: `git show 30bf904:CLAUDE.md`.

## Next

- [ ] **Fir cleanup (2026-11-13), for the user to run** (Globus delete was blocked here): old-named
      2020 maps (copies in local `change_maps/`), the withdrawn 19A log, `refit_vath/`, the 3,370
      pre-2020 07 logs and `v4_bootstraps/` (196 GB):
      `cd /home/mannfred/scratch/impact_assessment && rm -f output_figures/change_grid_*.png Rscripts/slurm-603212.out logs/Y{1995,2000,2005,2010,2015}_S*.log && rm -rf refit_vath ../v4_bootstraps && echo FIR_CLEAN_DONE`
- [ ] **Supplementary results**: paste `writing/supp_results_2026-10-31.docx` (113 species: map,
      then BCR sector table from `results_bcr.csv`; built by `writing/supp_results_docx.R`, ~28
      min) into the Google Doc "Supp. Mat. for: Declining boreal birds" after the Supplementary
      Methods, replacing the 138-species tables. The Drive connector cannot edit the Doc's body.
- [ ] Stage 18 to Fir before its next run (it now writes `output_figures/change_maps/` with the
      year in every name; Fir has the old version).

## Results (2020, C2d backfill, BAM release filter applied)

From 14C. Change = observed − no-industry (negative = birds lost to industry); % of the
no-industry population.

| Species | Observed | No industry | Change (5–95%) | % | Largest sectors (%) |
|---|---|---|---|---|---|
| CAWA | 4.47 M | 5.15 M | −683k (−782k, −564k) | −13.3 | roads −5.7, pasture −3.0, crop −2.7 |
| OVEN | 37.65 M | 41.27 M | −3.62 M (−4.14 M, −2.75 M) | −8.8 | crop −3.3, pasture −3.0, built −1.3 |
| OSFL | 10.74 M | 11.11 M | −372k (−537k, −217k) | −3.3 | roads −2.1 |
| LEYE | 5.79 M | 5.83 M | −41k (−61k, −25k) | −0.7 | roads −0.6 |
| GRSP | 0.96 M | 0.90 M | +59k (−293k, +408k) | +6.6 | roads +7.1, pasture +2.6, crop −2.5 |
| BOBO | 6.76 M | 3.63 M | +3.14 M (+2.96 M, +3.33 M) | +86.4 | roads +30.3, crop +26.6, pasture +20.2 |

- Of the 138 species before 14D's filter: 88 decline (72 with the whole interval below 0), 50
  increase (38 above 0); median −1.7%. Largest declines EATO −62.2%, FISP −43.3%, GCFL −30.7%;
  largest increases are open-country or human-associated (BARS +156.9%, ROPI +151.7%). Roads are
  the largest sector for 98 species. EATO's and FISP's no-industry populations (2.7× and 1.8×
  today's, mostly can13) were not examined.
- Withheld by BAM (14B drops them): CAWA can40; GRSP can10, can61; BOBO can10; LEYE can10; every
  BCR of BGGN, BWWA, DUNL, GWWA, NESP, SPGR, YBCU, YTVI.
- Direct vs vegetation (OL #9; 14B's sign, cf − obs): CAWA +138k / +545k; OVEN −1.68 M / +5.30 M.
- 10C flags 4 of 674 subbasins (50, 72, 430, 481): 22k for CAWA, −71k for BOBO, < 4k otherwise.
- Sensitivity (pre-C2d, CAWA and OVEN, not rerun; production is the headline): without the q99
  cap v(N) +1.16% / +0.05%; densities < 0.005 birds/ha set to 0: −2.7% / −0.5% nationally, but
  CAWA can11 −81% (the prairie density floor) while can13 and can61 rise.
- Production logs: `cluster_logs/07_2020_c2d/`, `c2d_6spp/`, `c2d_allspp_2026-10-14/`.

## Open before the write-up

- [ ] **Ecological check of the farmland BCRs (can11–13).** Forest species gain more birds there
      without industry than they have today (no-industry birds on the all-sector footprint per
      bird today, C2d): OVEN can11 4.9, can13 2.1; CAWA can11 2.6, can13 4.3; OSFL can13 5.9.
      Accepted for the all-species run (user, 2026-10-13), to be handled here. The user rejected a
      kNN hot-deck backfill ("a new method on top of BART that is a bit arbitrary").
- [ ] **VATH's industry-positive result is the bird model's direct part** (OL #9, not a bug).
      14C: +332k (+2.9%), mostly roads, can10 +287k and can5 +73k; 14B split (obs − cf): direct
      +409k, vegetation −77k. In can10 ~93% of the gain is the step from no footprint to 10% of
      the observed footprint, below the lowest survey decile; in can5 it builds up gradually. 58
      of 113 species have a positive direct road part. **For the write-up: report VATH with its
      split.** Figures: `output_figures/direct_part/` (`Rscripts/misc/figures/figure_footprint_scaling.R`: 0 → observed
      scaling curves for RBNU, SWTH, RCKI, OCWA, TOWA, VATH, the format going forward;
      `VATH_zero_footprint_surveys_can10.png`).
      Checked: the zero-footprint surveys (CanHF_1km = 0) do show a VATH deficit, model-free,
      against their 5 nearest matched surveys (density ratio 0.54, 90% CI 0.44–0.64; 0.63 within
      ±100 m elevation), clearly at middle and high elevation but not below ~1,000 m, where the
      model still applies its ~30% drop. Refuted: the 22 roadside BBS stops (refit without them:
      +4k). Season does not differ. Anna's predicted order of the drop (can9 > can10 > can5) holds
      for none of the six species. Untested: navigable waterways (scored into CanHF, not one of
      our 8 sectors).
- [ ] **V5 bootstrap models pair offsets with the wrong surveys (UPSTREAM; not yet reported to
      BAM, email drafted 2026-11-05).** V5 `main` a53c6ab: `03:826-832` moves non-Canadian surveys
      to the end of `visit` but not of `offsets`; `06:105`, `05.Tune:145` and `12.Validate:127`
      pair them by position. Fix, one line at `04:364`:
      `offsets <- offsets.use[match(visit$id, offsets.use$id), ]`. Scope (bootstrap 1, rows with
      another survey's offset): 13 of 19 Canadian BCRs, 56.4% of Canadian training rows (can5
      99.9%, can81 98.6, can14 98.4, can12 98.1, can82 95.9, can71 88.8, can70 84.6, can13 81.8,
      can40 81.7, can3 80.6, can72 79.1, can9 26.4, can10 25.3; 0 in can11, 41, 42, 60, 61, 80);
      also the tuning and validation AUCs. Reproducer for BAM:
      `Rscripts/misc/v5_upstream/PR1_visit_mismatch_test.R`. Refits with V5's exact construction reproduce V5's
      stored trees and predictions (to 6e-8), so the arms differ only in offsets. Matched − V5,
      observed population, 8 paired bootstraps (`Rscripts/misc/v5_upstream/PR1_offset_refit_observed.{R,sh}`;
      CSVs and logs in `cluster_logs/pr1/`):
      can81 OVEN +0.97%, SWTH +0.84%, REVI +2.35%; can5 VATH +4.9% (interval spans 0), RCKI
      −5.7% (spans 0), TOWA +11.0% (8/8 up). In can5 the offsets are near noise and 69–74% of
      cells move > 10%. VATH can10 (5 seeds): direct part +317k → +396k, so the VATH drop is not
      an artefact of the misalignment. Not tested: the effect on obs − cf. Fir keeps the fits
      (`offset_test/`) and the V5 training table (`data/raw_data/04_NM5.0_data_stratify.Rdata`,
      also local) to extend it.
- [ ] **Ten V5 training covariates are attached to the wrong surveys (UPSTREAM; not yet reported,
      email drafted 2026-11-05).** The eight ERA layers, CCNL_1km and MODISLCC_1km (every
      temporally matched GEE layer with `Extraction == "point"`). `03.ExtractCovariates.R`
      extracts by year (`:676`), stacks the year blocks (`:769`), then `cbind()`s them by position
      (`:770`) onto `loc.gee`, which is not in year order; radius layers join by id and are fine.
      `03_AK:564-566` has the same code (untested); committed 03 does not parse (the j loop's `}`
      is missing after line 766). Evidence (VATH can10): these ten correlate −0.20 to 0.08 with
      the stack at the survey, the other 40 at 0.30–1.00; undoing the sort recovers r 0.96–1.00.
      An ERA value belongs to a survey a median 1,251 km away. Prediction uses true values. The
      ten carry a median 10.7% of relative influence (5–95% 1.6–24.3) over 867 Canadian models;
      PUFI can61 95%, RUBL can12 51%. Effect on predictions not measured (needs a refit).
      Reproducer for BAM: `Rscripts/misc/v5_upstream/PR2_covariate_mismatch.R`.
- [ ] **BC valley bottoms: is the backfill too thin where training land is high?** BC's direct
      footprint sits low (median 848 m vs training 1,175 m; 41% below its subbasin's training 5th
      percentile). Backfilled vegetation vs intact land in the same subbasin and 100 m band, both
      at footprint 0 (birds/km²): OCWA +8.6%, RBNU +6.9%, RCKI +10.6%, SWTH +3.6%, TOWA −9.5%,
      VATH −23.3% (widest in valley bottoms). In birds, intact minus bf is +255k VATH and +275k
      TOWA, the two industry-positive species. A what-if, not a correction: intact valley land is
      a selected sample. Trap: compare with intact land predicted at footprint 0, never observed
      density. **Next:** compare backfilled canopy height, closure and conifer share with intact
      land at the same elevation (needs the BART mosaics, so on Fir). Scripts
      `Rscripts/misc/bc_valley_bottoms/`; tables `output_figures/bc_backfill/`;
      `data/derived_data/bc_d0_*` and `data/raw_data/{dem_bc_1km.tif,bc_boundary.gpkg}`.
- [ ] **Limitations sentence on the 2020-only decision.** BAM's minutes: the year models
      "underestimate the magnitude of both positive and negative trends". If that is attenuated
      habitat responses (BRT shrinkage, 1 km covariate error), our counterfactual differences are
      attenuated too, so the industry effects are conservative; attenuation from drivers outside
      the covariates cancels in a one-year obs − cf contrast.

## Loose ends

- [ ] **BART does not mix at 08B's settings (decision open).** 8 production fits: median pixel
      R-hat 1.92 (50 trees), 1.69 (200); one chain's 100 stored draws are worth ~1.3–1.6 effective
      draws. Footprint-shaped holdout vs 08B's random one: RMSE ×1.5–1.7, 90% interval covers
      81% / 63% (observed / cascaded inputs). Effect: BART is 1.4% of the variance of national
      change for the median species (NOCA 56%); chain-to-chain error, missing from the intervals,
      would widen SDs 2–13% for the median species. No `gbart` tweak fixes mixing; long `gbart`
      runs memorise the training land. **Proposed:** `deploy_gbart()` → `dbarts::bart2`, 200
      trees, no sparse prior, 4 chains of burn-in 3,000 + 5,000 kept, all pixels, 25 stored draws
      per chain (100 layers, so 11 and 12F unchanged), PSIS-LOO for the holdout metric, R-hat
      written per fit. ~6× 07's time; backfill means will move, so results change. First pilot
      ~10 subbasins on Fir across sizes (S57, S98 included); `dbarts` is not on Fir yet. `mbart2`
      mixes fine. Scripts `Rscripts/misc/bart_mixing/`; CSVs
      `data/derived_data/rds_files/bart_production_settings/`.
- [ ] **impactexploreR still shows 1995–2020** (user's choice). To make it 2020-only, change 20's
      `years` (it stops without the deleted pre-2020 CSVs) and hide the app's year slider.
- [ ] **Anna's extreme densities do not reach our tables** (her land-change stacks never passed
      V5's truncation; we apply densmax and the frozen q99 to obs, d0 and bf). Open, reporting
      only: Anna and Diana drop CLSW can14, LISP can14, WAVI can5 and the irruptives PISI, EVGR,
      PIGR, WWCR, RECR, all `include = TRUE` in 14D. To drop pairs, add them to 14B's
      `withheld_models`.
- [ ] **BART and kNN differ where training land is thin or far** (BOBO can13 subbasin 98; GRSP
      can11 subbasins 55, 60). Same direction, no coding defect: report as method sensitivity
      with 10C's flags.
- [ ] 12A's memory grows with every BCR-year in one task (39 of 113 OOM-killed at 16G after 58–69
      BCR-years). If it runs many years again, give it 32G or cap GDAL's cache
      (`terra::setGDALconfig("GDAL_CACHEMAX", "1024")`, untested).
- [ ] Optional: run 12A on Fir for CAWA into a scratch dir and diff against the G:-built outputs.
- [ ] Optional, cheaper 12D: 16 distinct BART draws per bootstrap instead of 100 with replacement
      (~64 → ~25 core-hours for CAWA + OVEN; all-sector mean moves 0.08–0.18% of the impact).
- [ ] `15C_singletons_plot.R` map panel is broken (`predictions_coalitions/` deleted);
      regenerate via `save_arrays_ids` if wanted.

## Species lists

- `Rscripts/species_113.txt`: the 113 reported species (14D's `include`);
  `run_species_batches.sh`'s default. Only they have results on either end.
- `Rscripts/species_todo.txt`: 146 species, 14D's candidate list (BOBO, CAWA, GRSP, LEYE, OSFL,
  OVEN commented out: they ran first). Left out: DEJU, DOWO, EUST, FOSP (no 2020
  `07_predictions` on Fir), ROSA (no Canadian models).
- 12H needs ~0.027 GB per per-BCR file, so its 32G holds a wave of ~1,000.

## Staging discipline: read before any cluster run

Stage the entry point's whole **`source()` closure**, data files included, and checksum-sync all
of it (a 0-byte transfer proves equality). Unedited files have been stale on Fir four times.

- Never restage an `.R` file while a job running it is active (Rscript rereads its file at the
  end; Globus rewrites it in place).
- `.sh` files must be LF: `tr -dc '\r' < f | wc -c` must print 0.
- Globus refuses sources outside the local endpoint's root: stage from inside the repo.
- Before staging a 12F change, run `Rscript Rscripts/misc/harnesses/verify_12f_can60_local.R` (~4 min, on
  Fir's CAWA can60 inputs in `cluster_logs/localtest/` and `cluster_logs/smoke7/`). Its baseline
  predates C2f, so compare the new 12F against the previous 12F on the same inputs: tables and
  arrays `identical()` unless predictions changed on purpose (Fir vs Windows: ~3e-16 relative).
- Closures: 06 = `06` `.R`+`.sh`, `CAfire/CAfire_{year}_masked.tif`,
  `isric_soil/isric_soil_covariates_masked.tif`.
  07/11 = `07` `.R`+`.sh`, `08A`, `08B_*`, `09_*`, `11` `.R`+`.sh`.
  12D/12H = `12D` `.R`+`.sh`, `12E`, `12F`, `12G` `.R`+`.cpp`, `12H` `.R`+`.sh`.
  12A = `12A` `.R`+`.sh`, `12B`. 12C = `12C` `.R`+`.sh`.
  14B = `14B` `.R`+`.sh`, `14C`, `18` `.R`+`.sh`, `12E`, `data/raw_data/species_names.csv`,
  `data/raw_data/ModelReleaseDecisions.xlsx`, `data/derived_data/rds_files/extrapolation_flags.csv`.

## Globus from this laptop

- Write filters as `--opt=value` (`"--include=slurm-<id>_*.out" "--exclude=*"`): a bare `*` is
  glob-expanded on Windows. Filters do not keep subdirectories out of a recursive pull.
- From Git Bash, a Fir path as its own argument (`globus ls`, `globus rename`) needs
  `MSYS_NO_PATHCONV=1`.
- Stuck ACTIVE with `GC_NOT_CONNECTED`: start
  `C:\Program Files (x86)\Globus Connect Personal\bin\globus_connect_personal.exe`.
- A log still being written fails checksum verification: pass `--no-verify-checksum`.
- Globus Connect Personal cannot read Claude's scratchpad under `AppData`: stage from the project.
- Re-login: `globus session update computecanada.ca --no-local-server` in a PowerShell window.
- Listing a directory of thousands of files (Fir `logs/`) can time out; retry.
- Fir's `sacct -S/-E` times are Pacific (Globus listings are UTC). Logs land in `Rscripts/` on
  Fir. Job IDs wrapped in 2026-10 (ceiling 67,043,328): pull logs by the IDs sbatch printed. A
  filter's `[3-9]` is not a numeric range.
- The laptop's C: drive is 236 GB and has filled before: check `Get-PSDrive C` before a large pull.

## Workstream codes

Code comments and logs name finished work by these codes (records in the `git show` snapshots
above): **A** CAfire fix (Fix A, Fix B; OL #5); **B** V5 packaging conformance (gates G1–G4);
**C1** the 12D rework (identity gate, C++ tree walk, per-bootstrap reduction); **C2a** per-sample
Shapley sums; **C2b** 100 BART draws with replacement; **C2c** 10C as area of applicability;
**C2d** modelled blanks, no water class, drawn land-cover classes (the production backfill);
**C2e** footprint covariates at 0 (OL #9); **C2f** two masks per sector (OL #10); **C3** no-q99
sensitivity; **C4** low-density floor (`DENS_MIN`); **C5** pixel-level change maps.
**PR1** / **PR2**: the two upstream V5 defects above (offsets, covariates).

V5 renumbered its scripts in `f082866` (`10.Package` → `10.Truncate` + `11.Package`, and so on);
`06_bootstraps` and `07_predictions`, the only V5 outputs we read, were not renamed.
