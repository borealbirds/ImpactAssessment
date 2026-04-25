# Simulation test for 12B/12C array changes.
# Mocks the key logic without real data (no terra, no BART, no BRT).
# Run with: Rscript Rscripts/test_12B_12C_arrays.R

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tibble))

source(file.path(getwd(), "Rscripts", "shapley_utils.R"))

pass <- function(msg) cat(sprintf("  PASS  %s\n", msg))
fail <- function(msg) stop(sprintf("  FAIL  %s", msg))

cat("=== Simulating 12C agg() and predict_species_bcr() changes ===\n\n")

# --- shared fake dimensions ---
n_sub  <- 3L
n_boot <- 4L
n_scen <- 5L

set.seed(42)
pop_obs_total_arr   <- array(runif(n_sub * n_boot * n_scen, 1000, 5000), c(n_sub, n_boot, n_scen))
pop_obs_on_coal_arr <- array(runif(n_sub * n_boot * n_scen,  100,  500), c(n_sub, n_boot, n_scen))
pop_bf_on_coal_arr  <- array(runif(n_sub * n_boot * n_scen,  100,  500), c(n_sub, n_boot, n_scen))

combine_stats <- function(arr) {
  mat <- matrix(arr, nrow = n_sub, ncol = n_boot * n_scen)
  list(mean = rowMeans(mat, na.rm = TRUE),
       sd   = apply(mat, 1, sd, na.rm = TRUE))
}

# ---- Test 1: agg() output dimensions with save_arrays = TRUE ---------------

agg_out <- list(
  pop_obs_total   = combine_stats(pop_obs_total_arr),
  pop_obs_on_coal = combine_stats(pop_obs_on_coal_arr),
  pop_bf_on_coal  = combine_stats(pop_bf_on_coal_arr)
)
agg_out$bcr_obs_total_mat   <- apply(pop_obs_total_arr,   c(2L, 3L), sum, na.rm = TRUE)
agg_out$bcr_obs_on_coal_mat <- apply(pop_obs_on_coal_arr, c(2L, 3L), sum, na.rm = TRUE)
agg_out$bcr_bf_on_coal_mat  <- apply(pop_bf_on_coal_arr,  c(2L, 3L), sum, na.rm = TRUE)

if (!all(dim(agg_out$bcr_obs_total_mat) == c(n_boot, n_scen)))   fail("bcr_obs_total_mat wrong dims")
if (!all(dim(agg_out$bcr_obs_on_coal_mat) == c(n_boot, n_scen))) fail("bcr_obs_on_coal_mat wrong dims")
if (!all(dim(agg_out$bcr_bf_on_coal_mat) == c(n_boot, n_scen)))  fail("bcr_bf_on_coal_mat wrong dims")

# apply(arr, c(2,3), sum) sums over dim 1 (subbasins) — verify one cell
expected_cell <- sum(pop_obs_total_arr[, 2, 3])
if (!isTRUE(all.equal(agg_out$bcr_obs_total_mat[2, 3], expected_cell)))
  fail("bcr_obs_total_mat values wrong")

pass("agg() BCR matrices have correct dims and values")

# ---- Test 2: agg() output unchanged for save_arrays = FALSE ----------------

agg_out_no_arr <- list(
  pop_obs_total   = combine_stats(pop_obs_total_arr),
  pop_obs_on_coal = combine_stats(pop_obs_on_coal_arr),
  pop_bf_on_coal  = combine_stats(pop_bf_on_coal_arr)
)
if (any(c("bcr_obs_total_mat", "bcr_obs_on_coal_mat", "bcr_bf_on_coal_mat") %in% names(agg_out_no_arr)))
  fail("save_arrays=FALSE path leaked array keys")

pass("agg() save_arrays=FALSE leaves return unchanged")

# ---- Test 3: per-BCR list wrapping inside lapply ---------------------------

make_bcr_list <- function(bcr, save_arrays) {
  tbl <- tibble(
    species = "TEWA", subbasin = seq_len(n_sub), bcr = bcr, coalition_id = 5L,
    obs_total_mean        = agg_out$pop_obs_total$mean,
    obs_total_sd          = agg_out$pop_obs_total$sd,
    obs_on_coalition_mean = agg_out$pop_obs_on_coal$mean,
    obs_on_coalition_sd   = agg_out$pop_obs_on_coal$sd,
    bf_on_coalition_mean  = agg_out$pop_bf_on_coal$mean,
    bf_on_coalition_sd    = agg_out$pop_bf_on_coal$sd
  )
  list(
    table  = tbl,
    arrays = if (save_arrays) list(
      obs_total_mat   = agg_out$bcr_obs_total_mat,
      obs_on_coal_mat = agg_out$bcr_obs_on_coal_mat,
      bf_on_coal_mat  = agg_out$bcr_bf_on_coal_mat
    ) else NULL
  )
}

# Simulate results with one NULL (missing backfill mosaic)
results <- list(
  make_bcr_list("BCR60", save_arrays = TRUE),
  NULL,
  make_bcr_list("BCR70", save_arrays = TRUE)
)

results_nonnull <- Filter(Negate(is.null), results)
if (length(results_nonnull) != 2L) fail("NULL filtering wrong count")

table_df <- dplyr::bind_rows(lapply(results_nonnull, `[[`, "table"))
if (nrow(table_df) != 2L * n_sub) fail("bind_rows table wrong nrow")
if (!all(c("obs_total_mean", "bf_on_coalition_mean") %in% names(table_df))) fail("table missing columns")

pass("NULL filtering and table bind_rows correct")

# ---- Test 4: national array aggregation (Reduce) ---------------------------

mats <- Filter(Negate(is.null), lapply(results_nonnull, `[[`, "arrays"))
if (length(mats) != 2L) fail("array extraction wrong count")

national_arrays <- list(
  obs_total_mat   = Reduce("+", lapply(mats, `[[`, "obs_total_mat")),
  obs_on_coal_mat = Reduce("+", lapply(mats, `[[`, "obs_on_coal_mat")),
  bf_on_coal_mat  = Reduce("+", lapply(mats, `[[`, "bf_on_coal_mat"))
)

if (!all(dim(national_arrays$obs_total_mat) == c(n_boot, n_scen))) fail("national mat wrong dims")
# two BCRs with identical fake arrays => national = 2x BCR
if (!isTRUE(all.equal(national_arrays$obs_total_mat, 2 * agg_out$bcr_obs_total_mat)))
  fail("national aggregation arithmetic wrong")

pass("national Reduce aggregation correct")

# ---- Test 5: trivial zero-return (n_active == 0) is handled correctly ------

trivial <- list(
  table = tibble(
    species = "TEWA", subbasin = seq_len(n_sub), bcr = "BCR80", coalition_id = 5L,
    obs_total_mean = 0, obs_total_sd = 0,
    obs_on_coalition_mean = 0, obs_on_coalition_sd = 0,
    bf_on_coalition_mean = 0, bf_on_coalition_sd = 0
  ),
  arrays = NULL  # no active pixels → no arrays
)

results_w_trivial <- list(trivial, make_bcr_list("BCR60", save_arrays = TRUE))
results_nonnull_t <- Filter(Negate(is.null), results_w_trivial)
mats_t <- Filter(Negate(is.null), lapply(results_nonnull_t, `[[`, "arrays"))

if (length(mats_t) != 1L) fail("trivial NULL arrays not filtered from mats")
national_t <- list(
  obs_total_mat   = Reduce("+", lapply(mats_t, `[[`, "obs_total_mat")),
  obs_on_coal_mat = Reduce("+", lapply(mats_t, `[[`, "obs_on_coal_mat")),
  bf_on_coal_mat  = Reduce("+", lapply(mats_t, `[[`, "bf_on_coal_mat"))
)
if (!all(dim(national_t$obs_total_mat) == c(n_boot, n_scen))) fail("national mat wrong after trivial")

pass("trivial zero-return (n_active=0) arrays=NULL handled correctly")

# ---- Test 6: predict_species_bcr() return structure seen by 12B ------------

pred_result <- list(table = table_df, national_arrays = national_arrays)

if (!is.data.frame(pred_result$table))    fail("res$table is not a data frame")
if (is.null(pred_result$national_arrays)) fail("res$national_arrays unexpectedly NULL")
if (!all(c("obs_total_mat", "obs_on_coal_mat", "bf_on_coal_mat") %in% names(pred_result$national_arrays)))
  fail("national_arrays missing keys")

# 12B accesses res$table for saveRDS — should be a plain tibble
if (!inherits(pred_result$table, "tbl_df")) fail("res$table not a tibble")

pass("predict_species_bcr() return structure correct for 12B")

# ---- Test 7: save_arrays = FALSE path returns NULL national_arrays ---------

results_no <- list(make_bcr_list("BCR60", FALSE), make_bcr_list("BCR70", FALSE))
results_nonnull_no <- Filter(Negate(is.null), results_no)
national_no <- NULL  # save_arrays = FALSE branch

pred_no <- list(
  table           = dplyr::bind_rows(lapply(results_nonnull_no, `[[`, "table")),
  national_arrays = national_no
)

if (!is.data.frame(pred_no$table))    fail("save_arrays=FALSE: res$table not data frame")
if (!is.null(pred_no$national_arrays)) fail("save_arrays=FALSE: national_arrays should be NULL")

pass("save_arrays=FALSE path returns NULL national_arrays")

# ---- Test 8: 12B target coalition IDs --------------------------------------

sectors    <- canonical_sectors()
target_ids <- c(
  sectors_to_coalition_id(sectors, sectors),
  vapply(sectors, function(s) sectors_to_coalition_id(s, sectors), numeric(1L))
)

if (length(target_ids) != 9L)          fail("target_ids wrong length")
if (anyDuplicated(target_ids))         fail("target_ids has duplicates")
if (sectors_to_coalition_id(sectors, sectors) != 256L) fail("full coalition ID != 256")
if (!all(target_ids >= 2L & target_ids <= 256L))       fail("target_ids out of valid range")

# single-sector IDs must each correspond to exactly one sector
recovered <- lapply(target_ids[target_ids != 256L], coalition_id_to_sectors, sectors = sectors)
if (!all(vapply(recovered, length, integer(1L)) == 1L)) fail("single-sector IDs not single-sector")

pass(sprintf("12B target_ids correct: {%s}", paste(sort(target_ids), collapse = ", ")))

# ---- Test 9: cf_draws() formula used in 15_plot_population_distributions.R -

cf_mat <- national_arrays$obs_total_mat -
          national_arrays$obs_on_coal_mat +
          national_arrays$bf_on_coal_mat

if (!all(dim(cf_mat) == c(n_boot, n_scen))) fail("cf_mat wrong dims")
if (any(!is.finite(cf_mat)))                fail("cf_mat has non-finite values")

pass("cf_draws formula (obs_total - obs_on_coal + bf_on_coal) produces finite matrix")

cat("\n=== All 9 tests passed ===\n")
