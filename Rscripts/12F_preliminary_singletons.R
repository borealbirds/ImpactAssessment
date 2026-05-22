# ---
# title: Preliminary single-sector standalone impacts (interim, pre-Shapley)
# author: Mannfred Boehm
# ---
# Computes each sector's STANDALONE impact from the 8 single-sector coalition
# density tables (run AFTER 12D_combine.R fills the obs columns).
#
#   v(S) = cf(S) - obs ; per subbasin = bf_on_coalition - obs_on_coalition
#   sector standalone impact = sum over all subbasins (all BCRs) of v({sector})
#
# The SUM of the 8 standalone impacts is an UPPER BOUND on the all-8 industry
# impact (overlapping footprints -> jointly removing recovers fewer birds than
# the sum of removing each alone; submodularity). Coalition 64 (6 sectors,
# excludes rail+roads) is reported as a lower reference floor. The true total
# (= grand coalition 256, still pending) lies between them; the full Shapley
# decomposition will resolve the overlap.
#
# Uncertainty here is APPROXIMATE: subbasins treated independent, and within a
# subbasin Var(bf - obs) ~ Var(bf) + Var(obs) (ignores the positive bf/obs
# bootstrap covariance, so this OVERESTIMATES SD -> conservative). Exact joint
# SD needs the *_arrays.rds, which were not saved in bf-only mode.
#
# Run on a Fir login node:
#   module load StdEnv/2023 r/4.4.0
#   cd /home/mannfred/scratch/impact_assessment/Rscripts && Rscript 12F_preliminary_singletons.R
# ---

ia_dir <- "/home/mannfred/scratch/impact_assessment"
dt_dir <- file.path(ia_dir, "data", "derived_data", "density_tables")
source(file.path(ia_dir, "Rscripts", "12E_shapley_utils.R"))

sectors     <- canonical_sectors()
species_vec <- c("CAWA", "OVEN")
year        <- 2020
floor_id    <- 64L  # largest completed coalition: 6 sectors, excludes rail+roads

singleton_ids <- vapply(sectors,
                         function(s) sectors_to_coalition_id(s, sectors),
                         numeric(1L))

read_dt <- function(species, cid) {
  f   <- file.path(dt_dir, sprintf("%s_%d_coalition_%d.rds",         species, year, cid))
  fbf <- file.path(dt_dir, sprintf("%s_%d_coalition_%d_bf_only.rds", species, year, cid))
  if (file.exists(f))   return(list(dt = readRDS(f),   bf_only = FALSE, path = f))
  if (file.exists(fbf)) return(list(dt = readRDS(fbf), bf_only = TRUE,  path = fbf))
  list(dt = NULL, bf_only = NA, path = f)
}

summ <- function(dt) {
  imp <- dt$bf_on_coalition_mean - dt$obs_on_coalition_mean          # v(S), per subbasin
  v   <- dt$bf_on_coalition_sd^2 + dt$obs_on_coalition_sd^2           # approx, conservative
  # true gap = bf has a value (subbasin has footprint for this sector) but obs is NA.
  # subbasins with NO footprint have bf NA too — those are no-ops, not gaps.
  list(n_sub       = nrow(dt),
       n_bcr       = length(unique(dt$bcr)),
       impact_mean = sum(imp, na.rm = TRUE),
       impact_var  = sum(v,   na.rm = TRUE),
       obs_total   = sum(dt$obs_total_mean, na.rm = TRUE),
       n_na        = sum(!is.na(dt$bf_on_coalition_mean) & is.na(dt$obs_on_coalition_mean)))
}

fmt <- function(x) formatC(x, format = "f", big.mark = ",", digits = 0)
pct <- function(num, den) if (is.finite(den) && den != 0) sprintf("%+.1f%%", 100 * num / den) else "NA"

pooled <- list(sum8_m = 0, sum8_v = 0, floor_m = 0, floor_v = 0, obs = 0, ok = TRUE)

for (sp in species_vec) {
  cat(strrep("=", 78), "\n", sp, " ", year, "  (units: birds; impact = birds lost to footprint,\n",
      "  i.e. recovered if that sector were removed)\n", strrep("=", 78), "\n", sep = "")

  rows  <- list()
  miss  <- character(0)
  sum_m <- 0; sum_v <- 0; obs_ref <- NA_real_

  for (s in sectors) {
    cid <- singleton_ids[[s]]
    r   <- read_dt(sp, cid)
    if (is.null(r$dt)) { miss <- c(miss, sprintf("%s (cid %d): NO TABLE", s, cid)); next }
    if (r$bf_only)     { miss <- c(miss, sprintf("%s (cid %d): still _bf_only (12D not done)", s, cid)) }
    a <- summ(r$dt)
    if (a$n_na > 0)    { miss <- c(miss, sprintf("%s (cid %d): %d subbasins with NA obs", s, cid, a$n_na)) }
    rows[[s]] <- a
    sum_m <- sum_m + a$impact_mean
    sum_v <- sum_v + a$impact_var
    if (is.na(obs_ref)) obs_ref <- a$obs_total
  }

  ord <- names(sort(vapply(rows, `[[`, numeric(1L), "impact_mean"), decreasing = TRUE))
  cat(sprintf("%-32s %14s %12s %9s %6s\n", "sector", "impact", "~SD", "%obs", "BCRs"))
  cat(strrep("-", 78), "\n")
  for (s in ord) {
    a <- rows[[s]]
    cat(sprintf("%-32s %14s %12s %9s %6d\n",
                s, fmt(a$impact_mean), fmt(sqrt(a$impact_var)),
                pct(a$impact_mean, a$obs_total), a$n_bcr))
  }
  cat(strrep("-", 78), "\n")
  cat(sprintf("%-32s %14s %12s %9s\n", "SUM of 8 (UPPER-BOUND headline)",
              fmt(sum_m), fmt(sqrt(sum_v)), pct(sum_m, obs_ref)))

  fl <- read_dt(sp, floor_id)
  if (!is.null(fl$dt)) {
    fa <- summ(fl$dt)
    cat(sprintf("%-32s %14s %12s %9s   [6 sectors, no rail/roads]\n",
                "coalition 64 (LOWER floor)",
                fmt(fa$impact_mean), fmt(sqrt(fa$impact_var)), pct(fa$impact_mean, fa$obs_total)))
    pooled$floor_m <- pooled$floor_m + fa$impact_mean
    pooled$floor_v <- pooled$floor_v + fa$impact_var
  } else {
    cat("coalition 64 (LOWER floor)        : NO TABLE\n"); pooled$ok <- FALSE
  }
  cat(sprintf("%-32s %14s\n", "observed national population", fmt(obs_ref)))

  if (length(miss)) {
    cat("\n  ! data gaps (totals computed with na.rm, treat as incomplete):\n")
    for (m in miss) cat("    -", m, "\n")
    pooled$ok <- FALSE
  }
  cat("\n")

  pooled$sum8_m <- pooled$sum8_m + sum_m
  pooled$sum8_v <- pooled$sum8_v + sum_v
  pooled$obs    <- pooled$obs    + obs_ref
}

cat(strrep("=", 78), "\n", "POOLED  (CAWA + OVEN; species independent)\n", strrep("=", 78), "\n", sep = "")
cat(sprintf("%-34s %14s  %s\n", "SUM-of-8 upper-bound total impact",
            fmt(pooled$sum8_m), sprintf("~SD %s  (%s of obs)",
            fmt(sqrt(pooled$sum8_v)), pct(pooled$sum8_m, pooled$obs))))
cat(sprintf("%-34s %14s  %s\n", "coalition-64 lower floor",
            fmt(pooled$floor_m), sprintf("~SD %s  (%s of obs)",
            fmt(sqrt(pooled$floor_v)), pct(pooled$floor_m, pooled$obs))))
cat(sprintf("%-34s %14s\n", "observed national population", fmt(pooled$obs)))
cat("\ntrue all-8 industry impact (coalition 256, pending) lies between the\n",
    "floor and the upper bound; full Shapley resolves the footprint overlap.\n", sep = "")
if (!pooled$ok) cat("\nNOTE: data gaps above — rerun after 12D_combine fully completes.\n")
