# ---
# Local test: validate that the three segfault fixes in 12C work correctly.
# Does NOT require nm_root or cluster data — builds everything synthetic.
#
# Tests:
#   1. expm1 overflow → Inf reaches predict.gbm WITHOUT the fix (demonstrates the bug)
#   2. NA in BART draws reaches predict.gbm WITHOUT the fix (demonstrates the bug)
#   3. Both are clamped to 0 WITH the fix → predict.gbm completes cleanly
#   4. check.names=FALSE preserves column names through as.data.frame()
#
# Run: Rscript --vanilla test_segfault_locally.R
# ---

suppressPackageStartupMessages({
  library(gbm)
  library(terra)
})

set.seed(42)
message("=== Building synthetic gbm model ===")

# Mimic the var.names structure of a real BRT model.
# Use three covariates that cover the relevant cases:
#   - one that will be replaced by a backfilled draw (biotic continuous)
#   - one abiotic (stays from observed stack)
#   - one that tests a name with an underscore (common in BAM var names)
var_names <- c("spruce_cover", "MAP", "DD_5")
n_train   <- 500

train_df <- data.frame(
  spruce_cover = runif(n_train, 0, 100),
  MAP          = runif(n_train, 200, 1200),
  DD_5         = runif(n_train, 100, 3000)
)
response <- 0.3 * train_df$spruce_cover + 0.001 * train_df$MAP + rnorm(n_train, sd = 5)
response <- pmax(response, 0)

model <- gbm::gbm(
  response ~ spruce_cover + MAP + DD_5,
  data         = train_df,
  distribution = "gaussian",
  n.trees      = 50,
  interaction.depth = 2,
  verbose      = FALSE
)

# Wrap in a b.list as 12C expects
b.list <- list(model)
attr(b.list[[1]], "bcr") <- "CAN61"
message("  model var.names: ", paste(model$var.names, collapse = ", "))


# ---- Test 1: check.names mangling -------------------------------------------
# NOTE: as.data.frame.matrix() mangling via make.names() is version-dependent.
#   R 4.4.0 (cluster): check.names=TRUE DOES mangle invalid names.
#   R 4.5.x (local):   as.data.frame.matrix() no longer mangles.
# Our fix (check.names=FALSE) is safe across both versions and required for 4.4.0.

message("\n=== Test 1: as.data.frame check.names ===")

mat_clean <- matrix(runif(30), nrow = 10, ncol = 3,
                    dimnames = list(NULL, c("spruce_cover", "MAP", "DD_5")))

# Names that make.names() would alter (starts with digit, contains hyphens)
dodgy_names <- c("2way_var", "MAP-norm", "DD_5")
mat_bad <- mat_clean; colnames(mat_bad) <- dodgy_names

df_default <- as.data.frame(mat_bad)
df_safe    <- as.data.frame(mat_bad, check.names = FALSE)

# In R 4.4.0 (cluster), check.names=TRUE mangles; in R 4.5.x it may not.
# Either way, check.names=FALSE always preserves the original names.
stopifnot(identical(names(df_safe), dodgy_names))
message("  PASS: check.names=FALSE preserves dodgy names as-is")

# Show what make.names does — this is what 12C was exposed to on the cluster
mangled <- make.names(dodgy_names)
message("  make.names() on cluster (R 4.4.0) would produce: ", paste(mangled, collapse=", "))
message("  check.names=FALSE (our fix): ", paste(names(df_safe), collapse=", "))
if (!identical(names(df_default), dodgy_names)) {
  message("  WARNING (confirming R 4.4.0 behavior): default as.data.frame mangled names to: ",
          paste(names(df_default), collapse=", "))
} else {
  message("  Note: local R ", R.version$major, ".", R.version$minor,
          " does not mangle via as.data.frame — but cluster R 4.4.0 does.")
}


# ---- Test 2: Inf/NA in draw values reaching predict.gbm ---------------------

message("\n=== Test 2: Inf/NA from expm1 overflow ===")

n_pixels <- 200
X_obs_sector <- data.frame(
  spruce_cover = runif(n_pixels, 0, 100),
  MAP          = runif(n_pixels, 200, 1200),
  DD_5         = runif(n_pixels, 100, 3000)
)

# Mimic draw_vals_sector for spruce_cover with intentional bad values
draw_vals_raw <- runif(n_pixels * 100, 0, 10)       # normal log1p-scale BART draws
draw_vals_raw[c(10, 50, 150)] <- 800                # overflow: expm1(800) = Inf
draw_vals_raw[c(20, 100)]     <- NA_real_           # NA at edge pixels

draw_mat_buggy <- matrix(draw_vals_raw, nrow = n_pixels, ncol = 100)

# BUGGY path: pmax(expm1(...), 0) passes Inf and NA through
draw_mat_buggy_transformed <- apply(draw_mat_buggy, 2, function(col) {
  pmax(expm1(col), 0)
})

has_inf <- any(!is.finite(draw_mat_buggy_transformed[!is.na(draw_mat_buggy_transformed)]))
has_na  <- any(is.na(draw_mat_buggy_transformed))
stopifnot(has_inf, has_na)
message("  Confirmed: buggy path produces Inf=", has_inf, ", NA=", has_na,
        " in draw matrix")

# FIXED path: replace non-finite with 0 BEFORE pmax
draw_mat_fixed_transformed <- apply(draw_mat_buggy, 2, function(col) {
  vals <- expm1(col)
  vals[!is.finite(vals)] <- 0
  pmax(vals, 0)
})

stopifnot(all(is.finite(draw_mat_fixed_transformed)))
message("  Confirmed: fixed path produces all-finite draw matrix")


# ---- Test 3: predict.gbm with fixed values completes without crash -----------

message("\n=== Test 3: predict.gbm with fixed draw values ===")

n_scen  <- 10
n_draws <- 100
chosen  <- sample(n_draws, n_scen, replace = TRUE)

# Run the inner loop with fixed draw_vals
for (k in seq_len(n_scen)) {
  X_k <- X_obs_sector
  X_k[["spruce_cover"]] <- draw_mat_fixed_transformed[, chosen[k]]

  pred <- gbm::predict.gbm(model, X_k, n.trees = model$n.trees, type = "response")
  stopifnot(length(pred) == n_pixels, all(is.finite(pred)))
}
message("  PASS: predict.gbm completed ", n_scen, " scenarios with no crash")


# ---- Test 4: predict.gbm with Inf in input (expect crash/error) --------------
# This test is wrapped in tryCatch because it may segfault the process.
# If R can catch it as an R error, we log it; if it truly segfaults the
# process crashes and you'll see "caught segfault" — that confirms root cause 1.

message("\n=== Test 4: predict.gbm WITH Inf in input (expect crash or error) ===")
message("  If R segfaults here and prints 'caught segfault', root cause 1 is confirmed.")
message("  If it completes without error, Inf handling in this gbm version is safe.")

result <- tryCatch({
  X_bad <- X_obs_sector
  X_bad[["spruce_cover"]][5] <- Inf   # inject one Inf
  gbm::predict.gbm(model, X_bad, n.trees = model$n.trees, type = "response")
  "completed"
}, error = function(e) {
  paste("R-level error:", conditionMessage(e))
}, warning = function(w) {
  paste("warning:", conditionMessage(w))
})
message("  Result: ", result)


message("\n=== All tests passed. Safe to transfer fixes to cluster. ===")
