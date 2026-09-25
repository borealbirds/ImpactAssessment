# TODO — consolidated plan (last updated 2026-09-24)

Open work is the cluster backfill re-run (**A**) and the converge steps (**C**). Workstream
**B** (conform density to current V5 packaging) is finished and validated — gates G1–G4 all
passed. The diagnosis behind A is `CLAUDE.md` Open Limitation #5; refuted hypotheses live in
the `project_pipeline_history` memory.

## Critical path

```
cluster:  [A7 DONE] ─► [D DONE] ─► [smoke 7 PASSED] ─► [C1 12D+12H DONE] ─┐
local:    [B, G1–G4 DONE]                                                   ├─► [C2 means DONE; SDs open] ─► C3
```

**Next action: decide what goes into a C1 rerun** (see C2's open items). C1 passed on all 25
tasks and C2's Shapley means are in `sector_effects/`, but `shapley_sd` needs per-sample
coalition arrays from 12F. The BART-draws change (S = 16) and per-species draw seeding would
fit in the same rerun, which costs ~1 h of wall time.

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

**Current state (2026-09-24)**: `/Rscripts` on Fir matches `bbe252b`, including the workstream-D
files (`12D` `.R`+`.sh`, `12F`, `12G` `.R`+`.cpp`, `12H` `.R`+`.sh`). `density_tables/` on Fir
holds C1's merged production output (and `by_bcr/` its 25 per-BCR files).
`/Rscripts/12*` is now ten files (`12A` is local-only by design).

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

## D. 12D audit: NaN fix, speedups, per-BCR jobs (2026-09-24)

Audit prompted by scaling to hundreds of species × 6 years. Findings, all measured locally on
CAWA can11's real models and stack (harnesses in the session scratchpad; see memory
`project_gbm_nan_routing`):

- **Correctness — fixed.** terra returns missing cells as `NaN`; gbm sends only `NA_real_` down
  its missing branch. V5's observed predictions used the missing branch (NaN→NA reproduces them
  at 100% of test pixels; 12F's NaN input at 1.3%). So since Fix A, bf-side pixels with a missing
  observed-only covariate were predicted by a different rule than the obs side: CAWA can11, 3.2% of
  weight > 0 footprint pixels, ≈ +10% of the observed footprint birds as a fake impact. 12F now
  recodes, and an **identity gate** re-proves obs/bf parity against `observed_bootstraps.tif` on
  every run. Details in CLAUDE.md "12F restructure invariants".
- **Where time goes.** 98.7% of `predict.gbm` is gbm's compiled tree walk (26 ns per tree per
  pixel; CAWA can11 = 9150 trees). R overhead is 1.3%, so no language change would help.
- **12G** — a C++ tree walk, 3.6–4.4× faster and `identical()` to `predict.gbm` (every worker
  re-checks this).
- **No `M`.** Workers reduce their own bootstrap to per-coalition subbasin sums (bit-identical:
  `rowsum` is column-independent), removing the ~17 GB can11 matrix and its 2–3 copies.
- **Per-BCR tasks.** 12D is one array task per species × BCR at 8 cores / 48G / 3 h (was one
  16-core / 384G job per species running BCRs in series); 12H merges in the old BCR order.
  Alliance bills memory as core-equivalents (~3.9 GB each on Fir), so 384G on 16 cores was
  charged like ~98 cores.
- Measured and rejected: merging identical trees (8969 of 9150 unique → 1.02×); fewer
  coalitions (the cost is the superset field, which even a with/without-industry trend needs).

- [x] **Regression check (local, 2026-09-24)**, real CAWA can60 inputs Globus'd from Fir into
      `cluster_logs/localtest/`, 2 boots. (1) HEAD's 12F run locally is `identical()` to the
      cluster's smoke-6 tables and arrays (cid 2, 9, 256), so the local environment reproduces Fir.
      (2) The rewrite with the NaN recode + identity gate patched out is `identical()` to HEAD on all
      255 tables and 9 arrays: 12G and the per-worker reduction change nothing. (3) The per-BCR
      path (12D save → 12H `combine_bcr_results`) is `identical()` to the single-job path.
      (4) Identity gate: 2469/2469 match V5 for both bootstraps (1041 with a missing covariate).
      NaN-fix effect on can60 is small (full-coalition impact −0.016%, largest singleton −2%); can11
      is where it bites. Sampling 320 s → 75 s on one core (4.3×).
- [x] **Staged** 2026-09-24 (checksum sync, LF verified): `12D` `.R`+`.sh`, `12E` (0 bytes — equal),
      `12F`, `12G` `.R`+`.cpp`, `12H` `.R`+`.sh`. Rcpp is in Fir's R library.
- [x] **Smoke 7 PASSED.** 12D task 1 PASSED (job `61340018`, CAWA can60, 2 boots, 8 cores): `sourceCpp`
      compiled on Fir, every worker's fast-vs-`predict.gbm` check passed, identity gate 2469/2469
      for both bootstraps, sampling 17 s (smoke 6: 55 s), `nice.`. Its per-BCR file, pulled back,
      is `identical()` to the local rewrite on all 255 tables and 9 arrays, and its `code_md5`
      equals HEAD. 12H (`61340019`) correctly refused: `TEST_BCR=can60` also matches OVEN's can60,
      so the table has 2 tasks and only task 1 had run. Fixed: CLAUDE.md smoke is now
      `--array=1-2`; 12H's resubmit hint now carries `TEST_BCR`/`TEST_N_BOOT` (without them
      "task 2" would have been full-bootstrap CAWA can11), and 12H refuses a merge whose
      `TEST_N_BOOT` differs from the files'. All three behaviours checked locally.
      **Task 2 PASSED** (job `61342817_2`, OVEN can60): 22226 / 22226 weight > 0 complete, identity
      gate 2508/2508 for both bootstraps (1082 predicted pixels carry a missing covariate), sampling
      24 s, `nice.`. **12H PASSED** (job `61342818`): both files from one `code_md5`, `n_boot` 2,
      `wrote 255 coalition tables` + `wrote 9 array files` for each species, `nice.`. The smoke
      left real can60 files in `density_tables/`, `arrays/` and `by_bcr/` — C1's wipe covers all three.
- [x] **C1 attempt 1 (job `61344103`, 2026-09-24) — cancelled after the first 11 tasks exposed two
      defects.** Timings were excellent where tasks finished: CAWA can10 32 boots in 7.5 min
      (sampling 5 min on 8 cores; smoke 5's old code took 17 min for 2 boots), can40 1.5 min,
      can60 2 min, can71 1 min. Every weight > 0 superset was ≥ 99.998% complete.
      (1) **OOM at 64G** — can12, can13, can80, the last 80 s into sampling on only 56k complete
      rows, so not data size. Workers inherit the parent's GC trigger (inflated to ~10–15 GB by
      the draw reads), and 255 `sc[kr, ]` copies per bootstrap fed each one garbage up to it.
      Fixed: `coal_rowsum()` (C++, no copy; `identical()` to the old path on 200 random cases),
      weight/NA audit folded into the per-draw loop, `gc(full = FALSE)` per draw, one fork per
      bootstrap. A killed worker is now reported as such (it returned `NULL`, which the old check
      turned into an opaque `vapply` error).
      (2) **Identity gate, can14 bootstrap 32: 3/2986 off V5** (max rel diff 1%). Cause: 12F
      factor-converted categoricals with bootstrap 1's `var.levels` for all 32, but levels differ
      per bootstrap. Proven locally on all 837,600 can14 cells: bootstrap-1 levels → 155 cells
      off V5; each model's own levels → 0. Fixed with `as_model_factors()`; the complete-case gate
      now tests the raw backfilled class. The gate did its job — without it this would have gone
      into the tables unseen (any bootstrap 2–31 was never checked).
      Regression: new 12F on local can60 is `identical()` to Fir's smoke 7 (255 tables, 9 arrays),
      and at 4 boots `identical()` to HEAD run locally. Against Fir's C1 can60, 4 of 400 cid-2
      array cells differ by 1 ulp in bootstraps 3–4 — for HEAD and the new code alike, so it is
      Fir-vs-Windows floating point, not the change; tables were identical.
      **Full tally of attempt 1** (21 of 25 tasks ran before it was stopped): OOM in 9 — CAWA
      can11/12/13/61/80/81, OVEN can10/11/12; gate failure in 4, all bootstrap 32 — CAWA can14,
      OVEN can14/41/70. All four gate failures proven locally over every cell of the BCR: bootstrap-1
      levels put 155 / 219 / 12,230 / 3,752 cells off V5 (max rel diff 1.9% / 271% / 99.8% / 27%);
      each model's own levels put 0 off. OVEN can70's bootstrap 1 is the ODD one out (it lacks
      MODISLCC_1km classes 3/16 that all 31 others know). Old per-BCR outputs saved in
      `cluster_logs/c1_old/` for a Fir-vs-Fir check after the rerun: CAWA can10/40/60/71 and OVEN
      can13/60 must come back `identical()`; OVEN can40 must NOT — its 10 lost weight > 0 pixels
      were lost only to MODISLCC classes unknown to bootstrap 1, and the new gate keeps them
      (likewise OVEN can41's 6, which never finished).
- [ ] **Decision (yours): BART draws per bootstrap.** 100 scenarios drawn with replacement give
      ~64 distinct draws per bootstrap. Smoke-6 arrays (CAWA can60, 2 boots): for the full
      coalition the bootstrap spread is 5× the BART spread, so scenario count barely matters;
      for a small single sector BART dominates, but 10–20 distinct draws per bootstrap still leave
      Monte Carlo error < 1% of the impact. 16 distinct draws (no replacement) would cut gbm calls
      ~4× but changes numbers within MC noise. Confirm on C1's 32-boot arrays before adopting.
      **Confirmed on C1's national arrays (2026-09-24), all 9 array coalitions × 2 species.** BART
      is 3–28% of the per-sample variance (CAWA) and 10–34% (OVEN); the MC error of the mean is set
      by the 32 bootstraps (σ_B/√32 ≈ 17% of the reported SD, ~1–6% of the impact) plus a term
      from the 100 stored draws that every design shares. Moving to S distinct draws raises that
      MC error by at most 0.9% (S = 16) or 2.3% (S = 8), worst case OVEN mines. Empirically,
      subsampling 16 of the 100 scenarios per bootstrap moves the all-sector mean by 0.08% (CAWA)
      and 0.18% (OVEN) of the impact. Cost: sampling was 4.25 of C1's 5.17 task-hours, so S = 16
      should cut a C1 to ~2 task-hours (~25 billed core-equivalent-hours at 48G, from ~64).
      Side finding: `chosen_k` is seeded per species × BCR, so a subbasin straddling two BCRs is
      given independent draws on each side, though draw j is ONE joint posterior sample across the
      whole subbasin. Means are unaffected; the BART spread of straddling subbasins is slightly
      understated. Seeding on species × bootstrap only would give every BCR the same draw indices.
      Harness: session scratchpad `draws_analysis.R`.
- [x] **Profiled C1 and resized 12D to 8 cores / 48G / 3 h** (was 64G / 12 h). `sacct` over all 25
      tasks: 5.17 h of task time in total; longest CAWA can11 42 min, OVEN can11 38, OVEN can61 33;
      16 of 25 under 15 min. Peak MaxRSS 28.8 GiB (OVEN can61); can11/can61 22–29 GiB, mid-sized
      BCRs 16–23, small 4–17. 48G = 1.7× the peak (40G's 1.4× was judged too thin: MaxRSS is
      sampled every 30 s, and other species' models may carry more backfilled covariates). 3 h =
      4.3× the longest task, and ≤ 3 h jobs start sooner on Alliance clusters. Billed cost of a C1
      (Alliance bills max(cores, mem ÷ ~3.9 GB) × elapsed): 85 core-equivalent-hours at 64G → ~64 at 48G.
      **Savings vs the smoke-5 code** (one 16-core / 384G job per species, BCRs in series): per
      bootstrap 13–14× less compute, measured on CAWA can10 (17 → 1.2 core-min) and can11
      (130 → 9.8); prep per BCR 5–7 → 1–4 min. Projected old C1: ~19 h longest job, ~510 CPU
      core-hours, ~3,100 billed core-equivalent-hours (~2,100 if unmeasured BCRs sped up only the
      8× the separately measured parts explain); new: 42 min, 41 core-hours, ~64 billed.
      The C++ walk (3.6–4.4×) and distinct-draw + weight-0 skipping (1.87×) account for ~8× of the
      13–14×; the rest is not isolated.
- [ ] Out of 12D's scope but on the multi-year path: `07` (backfill) must re-run per year at the
      same 128G-per-core ratio; pre-2020 years need historical footprint layers, not the 2020 mask.

## C. Converge

- [x] **C1 PASSED (attempt 2: 12D `61367890`, 12H `61367891`, 2026-09-24, Fir at `bbe252b`).**
      All 25 tasks `nice.`, 0 errors / OOM / dead workers at 64G, including the 9 that OOM-ed in
      attempt 1. Identity gate: all 32 bootstraps reproduce V5 on every task (1,258–3,122
      observed-design pixels each; the odd-class over-sample was non-empty in 17 of 25, 138 pixels
      in can14). Weight > 0 superset 100% complete everywhere except can14 (130,502 / 130,504, both
      species), whose 2 lost pixels lie in no subbasin zone, so they could never enter a table.
      The bf-field NA audit never fired, so drawing `complete_mask` from draw 1 is safe as
      guarded. 12H: 25 per-BCR files, one `code_md5`, 255 tables + 9 array files per species.
      **Fir-vs-Fir against attempt 1** (`cluster_logs/c1_old` vs `c1_new`): CAWA can10/40/60/71
      and OVEN can13/60 `identical()` in every table and array. OVEN can40 differs exactly where
      predicted: only the 128 coalitions containing `roads` (cid 129–256), in 3 subbasins
      (192, 193, 263), with `obs_total` untouched. Those are the 10 road pixels bootstrap 1's levels
      had dropped. Effect on the all-sector BCR impact: −1277.78 → −1275.79 (0.16%).
      Speedup: CAWA can10 sampling 4.7 min vs 5 min in attempt 1, so the memory fixes cost nothing.
      Old tables were stale on four counts: pre-weighting, pre-gate-change (A),
      pre-truncation-conformance (B), pre-NaN-fix (D).
- [x] **C2 run (2026-09-24): Shapley MEANS are good; `shapley_sd` is not fit to report.**
      Pulled C1's 510 tables + 18 arrays + `extrapolation_flags.csv` (local gpkg checksum-equal to
      Fir's; merged tables `identical()` to the 7 per-BCR files held locally). 14B then ran in 7 min.
      **Fixed in 14B:** 12F marks "no kept pixels" by obs_on mean `NaN` with sd `NA` (sd of NaNs is
      NA), and 14B zeroed only `is.nan()`, so every subbasin with an empty coalition carried NA
      Shapley SDs and every national SD was NA. 14B now zeroes on the NaN-mean marker, stops if bf
      is non-zero there, and stops on any other NA. Means were `identical()` before and after.
      National (CAWA can40 withheld): CAWA v(N) = +574,377 birds on 4.47 M observed (12.9%): roads
      33.6%, pasture 22.4%, crop 21.0%, built 17.6%, rail 3.4%, mines 1.0%, dams 0.5%,
      oil_gas 0.4%. OVEN v(N) = +3,527,058 on 37.65 M (9.4%): crop 32.9%, pasture 31.7%,
      built 20.6%, roads 10.8%, rail 3.0%, mines 0.6%, oil_gas 0.6%, dams −0.2%. Additivity
      residual −6.2 / +2.4 birds is rounding only (every subbasin is within 0.5).
      **Open, for review before anything is reported:**
      (1) **`shapley_sd` is wrong in both directions.** It assumes subbasins are independent, but
      within a BCR they share the same 32 bird models: against the joint (bootstrap × scenario)
      arrays, 14B-style SDs of v(S) are 0.37–0.89× the true ones. And it treats v(S ∪ j) and v(S)
      as independent, though they come from the same samples and differ only on j's exclusive
      pixels, so a small sector inherits the big coalitions' noise: CAWA dams ±4,911 vs roads
      ±5,482, on means of 2,978 vs 192,957. Correct SDs need Shapley computed per sample, i.e. 12F
      saving BCR-level `[n_boot × n_scen]` sums for all 255 coalitions, not just 9 (~20 MB per BCR
      uncompressed), which means a C1 rerun.
      (2) **Extrapolation flags are uninformative:** 10C flags 667 / 667 subbasins (`ks_max > 0.5`
      OR Mahalanobis exceedance > 0.3; min `ks_max` is 0.52; median exceedance 0.93).
      (3) **Very large and negative BCR impacts deserve ecological scrutiny:** CAWA can11 +227%
      and can13 +301% of observed, OVEN can11 +275% and can13 +89% (the most converted BCRs).
      Negative: OVEN can12 −828k (−8.5%), can81 −201k (−6.3%), can80 −68k.
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
