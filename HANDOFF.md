# Handoff — 2026-09-23

Read `CLAUDE.md` (pipeline, Open Limitations) and `TODO.md` (task list) first. This file is only
the live state they do not capture. Delete it once it is stale.

## Where things stand

- **A6 is done for 2020.** `07` backfilled 674/674 subbasins; `11` built 19/19 BCR mosaics
  (job `60806856`, 2026-09-21 → 09-22, every `.out` ends `done.`).
- **Next action: the `12D` smoke, which is also A7** (the only test of Fix A). The user was about
  to submit it; check whether it has run before doing anything else:
  ```bash
  cd /home/mannfred/scratch/impact_assessment/Rscripts && sbatch --array=1 --time=01:00:00 --mem=192G --export=ALL,TEST_BCR=can60,TEST_N_BOOT=2 12D_repredict_all_coalitions.sh
  ```
  Pass condition: `12F`'s line `complete superset pixels: N / M (X%)` is far above the old 1–3%.
  Read it with:
  ```bash
  grep -h 'complete superset pixels\|nice\.\|Error' /home/mannfred/scratch/impact_assessment/Rscripts/slurm-*_1.out 2>/dev/null | tail -n 5 | awk '{print "|SMOKE| "$0}'
  ```
  If it fails, stop and diagnose — do not launch C1 (2 × 384 G × 24 h).
- **After a passing smoke, wipe then C1.** The smoke writes can60-only tables into
  `density_tables/` and `arrays/`; `12D:131` skips empty coalitions without writing, so leftovers
  could survive into C2:
  ```bash
  rm -f /home/mannfred/scratch/impact_assessment/data/derived_data/density_tables/*.rds /home/mannfred/scratch/impact_assessment/data/derived_data/density_tables/arrays/*.rds ; cd /home/mannfred/scratch/impact_assessment/Rscripts && sbatch 12D_repredict_all_coalitions.sh
  ```
- The user is working on **year 2020 only**; other years come later (see TODO for what that needs).
  Test species are **CAWA and OVEN** (not OSFL).

## Changed this session, not yet on the cluster

- `12A_observed.R`: BCR code now parsed from the `06_bootstraps` filename instead of loading each
  `b.list` (up to 654 MB each) for `attr(b.list[[1]], "bcr")`.
- `12F_predict_species_all_coalitions.R`: new guard `stop()`s if a `b.list`'s `bcr` attribute
  disagrees with its filename (12F loads `b.list` anyway, so it costs nothing).
- **Stage `12F` to the cluster before C1** (it is in `12D`'s `source()` closure; Globus, one file,
  `--sync-level checksum`, never `--batch`). A smoke run with the old `12F` is harmless — the guard
  is additive. `12A` is still local-only, so it needs no staging yet.

## Next project after C1: move `12A` to the cluster (for 100+ species)

Evidence already gathered (details in `TODO.md` → Loose ends): Elly's
`/home/mannfred/projects/def-ecknight/NationalModels/output/` has 151 species in both
`06_bootstraps` and `07_predictions`; the 2020 predictions and the bootstraps pair 1:1 (2,868
each, 1,686 Canadian); 4,504/4,505 files common with G: match in size (32 species compared);
`q.out` lists exactly the same 151 species. `12A` never reads V5's truncated products, so
whether Elly re-ran the truncation revision there is irrelevant.

Remaining step: run `12A` on the cluster for CAWA into a scratch output dir and diff against the
existing G:-built `observed_bootstraps.tif` + `truncation_params.rds`. If equal, switch `cc <- TRUE`,
write a `.sh`, stage `12B_v5_truncate.R`, and remove the Globus step from `CLAUDE.md`. For 100+
species, `species_vec` must come from `list.dirs(06_bootstraps)` (commented out at `12A:123`).

## Working conventions (easy to get wrong)

- **No SSH from this machine.** The user runs cluster commands and pastes output. Their terminal
  drops newlines and escapes underscores, so give **single-line** commands whose output lines
  start with a unique marker (`|SMOKE| `, `|DONE| `…). Watch for dropped spaces too — a paste of
  `--mem=192G--export=...` produced `Invalid --mem specification`.
- **Globus can list the cluster read-only without SSH**, including Elly's directory:
  `globus ls -l "8dec4129-9ab4-451d-a45f-5b4b8471f7a3:/home/mannfred/projects/def-ecknight/NationalModels/output/..."`.
  A recursive listing of `07_predictions` takes ~8 min.
- **Do not bulk-list or stat G:** (`find`, `ls -l` over big folders). Google Drive for desktop
  starts downloading tens of GB of tifs. G: may also be unmounted; ask before relying on it.
- **Completion oracles**: `sacct COMPLETED 0:0` proves nothing for `07` or `11`. `07`: the
  `_confusion.rds` + mtime. `11`: the `.out` ending `masked and written to …` / `done.`.
  `12D`: last line `nice.`.
- Commit directly to `main` and push (user-authorized). Other uncommitted files in the tree
  (`.gitignore`, `CLAUDE.md`, `15C`, `writing/`) are the user's — do not sweep them into commits.
