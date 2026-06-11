# Design: 12C restructure — compute the backfilled field once, reduce 255 coalitions cheaply

**Created:** 2026-06-02
**Motivation:** scale 12B/12C from 2 species to ~60. Profiling of the clean 510-task run
(`collect_seff.sh` → `analyze_seff.R`) shows peak RSS ~306 G, peak wall ~8.5 h, zero OOM/TIMEOUT,
and — critically — **memory is ~flat across coalition size** (size-1 already ~215 G). Cost is the
coalition-INVARIANT species×BCR setup, not the coalition footprint. See
`memory/project_12b_seff_profile.md`.

## Key insight: the per-pixel backfilled prediction is coalition-INDEPENDENT

In `12C_predict_species_bcr.R` the joint-sampling seed (line ~334) is:

```r
set.seed((sum(utf8ToInt(paste0(species, bcr_code))) + i * 1000L + k) %% .Machine$integer.max)
chosen <- sample(n_draws, 1)
```

It depends only on **species, BCR, bootstrap i, scenario k** — not the coalition. The predicted
design-matrix row for pixel p (disturbance→0, categorical→backfill@p, biotic-cont→draw `chosen`@p,
else observed@p) also has no coalition dependence. **The coalition only chooses which pixels are
masked in; never the backfilled value at a pixel.** So for fixed (species, BCR, i, k) the bf
density field over pixels is identical across all 255 coalitions. Today every coalition job
recomputes it (incl. 3200 gbm predicts/BCR) from scratch → the redundancy the flat-memory
profile reveals.

## Restructure (per species × BCR)

Expensive work ONCE over the **superset** = all-8-sectors mask (any sector footprint ∧ CanHF≥1,
i.e. cid=256):

1. Read `b.list`, `stack_obs`, `stack_bf` mosaic, obs bootstraps — once.
2. Extract obs covariates, BART draws (`expm1`, non-finite→NA, `pmax(.,0)`), categorical backfill
   at superset pixels — once.
3. Build `complete_mask` over superset — once (coalition-independent per pixel).
4. Joint BRT×BART sampling over superset pixels → bf field `M[n_super × (32·100)]` and obs field
   `O[n_super × 32]`. **3200 gbm predicts/BCR happen once, not ×255.** Keep the `qsp` per-pixel cap
   (`pmin(pred, qsp)`) here — it feeds the density table.

Then loop 255 coalitions as pure reductions (no raster I/O, no gbm):

5. Precompute 8 boolean sector-membership vectors over superset (once). Coalition pixel set = OR
   of its sectors' vectors; `keep = coalition ∧ complete_mask`.
6. `bf_on_coalition  = rowsum(M[keep,], zones[keep])` → [n_sub × (32·100)] → `combine_stats`.
   `obs_on_coalition = rowsum(O[keep,], zones[keep])` (scenario-invariant; replicate across k).
   `obs_total` is coalition-independent → compute ONCE (full-raster zonal per bootstrap).
   Multiply densities by 100 (birds/ha → birds/km²) exactly where the current `agg()` does.

## Bit-identical guarantee + verification

Because the seed is coalition-free, the once-computed field equals exactly what each coalition job
produces today. **Verification harness:** rerun a singleton (e.g. cid=129 roads) and the all-8
(cid=256) old-way vs new-way; assert the 10-column density tables match to floating point.

## Job-axis options (decide before coding)

- **A — 1 job/species** (array task = species): loops all BCRs, does superset-compute + 255-coalition
  reduce per BCR, writes 255 `.rds` at end. 60 jobs for 60 species. No stitch step; matches current
  output layout. Memory = one BCR's superset field (~today's size-8 peak ~267 G).
- **B — 1 job/species×BCR**: finer parallelism, shorter/smaller jobs, easy per-BCR retry; needs a
  trivial reduce step to stitch per-coalition tables across BCRs. ~720 jobs for 60 species.

## Correctness constraints (all preserved)

- complete-case mask identical between obs and bf (single `keep` drives both).
- BART draws log1p-scaled → `expm1` before use; non-finite → NA, never 0.
- categorical `var.levels` indexed by `match()` against `var.names`, never by name.
- `qsp` cap stays in the gbm step (feeds table). q99/q0 caps only ever touched the inspection
  rasters, never the density table.
- INTERLEAVE=BAND grouped draw reads (now once/BCR).

## Side simplification

`predictions_coalitions/.../backfilled_{mean,sd}.tif` are optional-inspection outputs no downstream
script reads (12D/12F/14B consume only `.rds` tables). Emit them only for the 9 `save_arrays`
coalitions (full + 8 singletons), identical to today when produced — drops 255×→9× raster writes.

## Projected savings

Mosaic reads ×255→×1/BCR; gbm predicts ÷255; memory unchanged; jobs 15,300→60 (A) or ~720 (B).
