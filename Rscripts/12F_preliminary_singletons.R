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
  # "active" = subbasin where the sector has REAL footprint contributing birds.
  # After the 12C keep_mask fix, subbasins whose every footprint pixel is
  # incomplete-case give obs_on_coalition_mean = 0 (and bf likewise), so
  # `!is.na(imp)` overcounts them. Require obs > 0 to count as active.
  active     <- !is.na(dt$obs_on_coalition_mean) & dt$obs_on_coalition_mean > 0
  bcr_imp    <- tapply(imp, dt$bcr, sum, na.rm = TRUE)                # per-BCR totals
  bcr_active <- tapply(active, dt$bcr, any)
  n_sub_act  <- sum(active)
  n_bcr_act  <- sum(bcr_active, na.rm = TRUE)
  bcr_obs    <- tapply(dt$obs_total_mean, dt$bcr, sum, na.rm = TRUE)
  sum_imp    <- sum(imp, na.rm = TRUE)
  obs_all    <- sum(dt$obs_total_mean, na.rm = TRUE)
  obs_actsub <- sum(dt$obs_total_mean[active], na.rm = TRUE)
  obs_actbcr <- sum(bcr_obs[bcr_active], na.rm = TRUE)
  obs_fp     <- sum(dt$obs_on_coalition_mean, na.rm = TRUE)            # birds on footprint pixels
  pct        <- function(num, den) if (is.finite(den) && den > 0) 100 * num / den else NA_real_
  list(n_sub             = nrow(dt),
       n_bcr             = length(unique(dt$bcr)),
       n_sub_active      = n_sub_act,
       n_bcr_active      = n_bcr_act,
       impact_mean       = sum_imp,
       impact_var        = sum(v, na.rm = TRUE),
       pct_footprint     = pct(sum_imp, obs_fp),
       pct_sub_active    = pct(sum_imp, obs_actsub),
       pct_bcr_active    = pct(sum_imp, obs_actbcr),
       denom_footprint   = obs_fp,
       denom_sub_active  = obs_actsub,
       denom_bcr_active  = obs_actbcr,
       obs_total         = obs_all,
       n_na              = sum(!is.na(dt$bf_on_coalition_mean) & is.na(dt$obs_on_coalition_mean)))
}

fmt <- function(x) formatC(x, format = "f", big.mark = ",", digits = 0)
pct <- function(num, den) if (is.finite(den) && den != 0) sprintf("%+.1f%%", 100 * num / den) else "NA"

pooled <- list(sum8_m = 0, sum8_v = 0, floor_m = 0, floor_v = 0, obs = 0, ok = TRUE)

# tidy-CSV accumulator: one row per (species, scope) -------------------------
csv_rows <- list()
add_row  <- function(sp, scope, label, cid, bf_only, a) {
  impact_sd <- sqrt(a$impact_var)
  sd_pct    <- function(den) if (is.finite(den) && den > 0)
                               100 * impact_sd / den else NA_real_
  pct_natl  <- if (is.finite(a$obs_total) && a$obs_total > 0)
                  100 * a$impact_mean / a$obs_total else NA_real_
  pct_natl_sd <- sd_pct(a$obs_total)
  # columns ordered so spatial scope increases left-to-right:
  # footprint -> subbasin-active -> BCR-active -> national
  csv_rows[[length(csv_rows) + 1L]] <<- data.frame(
    species           = sp,
    year              = year,
    scope             = scope,                       # "sector" | "coalition_64" | "sum_of_8"
    label             = label,
    coalition_id      = cid,
    bf_only           = bf_only,
    impact_mean       = a$impact_mean,
    impact_sd         = impact_sd,
    pct_footprint     = a$pct_footprint,
    pct_footprint_sd  = sd_pct(a$denom_footprint),
    pct_sub_active    = a$pct_sub_active,
    pct_sub_active_sd = sd_pct(a$denom_sub_active),
    pct_bcr_active    = a$pct_bcr_active,
    pct_bcr_active_sd = sd_pct(a$denom_bcr_active),
    pct_obs_natl      = pct_natl,
    pct_obs_natl_sd   = pct_natl_sd,
    n_sub_active      = a$n_sub_active,
    n_sub             = a$n_sub,
    n_bcr_active      = a$n_bcr_active,
    n_bcr             = a$n_bcr,
    obs_total_natl    = a$obs_total,
    stringsAsFactors  = FALSE
  )
}

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
    add_row(sp, "sector", s, cid, r$bf_only, a)
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

  # synthetic sum-of-8 row for the CSV (no n_sub/n_bcr counts at this scope)
  add_row(sp, "sum_of_8", "sum_of_8", NA_integer_, NA,
          list(impact_mean = sum_m, impact_var = sum_v,
               pct_footprint = NA_real_, pct_sub_active = NA_real_, pct_bcr_active = NA_real_,
               denom_footprint = NA_real_, denom_sub_active = NA_real_, denom_bcr_active = NA_real_,
               n_sub = NA_integer_, n_sub_active = NA_integer_,
               n_bcr = NA_integer_, n_bcr_active = NA_integer_,
               obs_total = obs_ref))

  fl <- read_dt(sp, floor_id)
  if (!is.null(fl$dt)) {
    fa <- summ(fl$dt)
    cat(sprintf("%-32s %14s %12s %9s   [6 sectors, no rail/roads]\n",
                "coalition 64 (LOWER floor)",
                fmt(fa$impact_mean), fmt(sqrt(fa$impact_var)), pct(fa$impact_mean, fa$obs_total)))
    add_row(sp, "coalition_64", "coalition_64", floor_id, fl$bf_only, fa)
    pooled$floor_m <- pooled$floor_m + fa$impact_mean
    pooled$floor_v <- pooled$floor_v + fa$impact_var
  } else {
    fa <- NULL
    cat("coalition 64 (LOWER floor)        : NO TABLE\n"); pooled$ok <- FALSE
  }
  cat(sprintf("%-32s %14s\n", "observed national population", fmt(obs_ref)))

  # ---- impact as % of observed, at three spatial scales ---------------------
  # footprint = total_impact / obs on footprint pixels only
  #             (most intense: "% of birds on the footprint that are displaced")
  # sub act   = total_impact / obs over subbasins WITH footprint
  # BCR act   = total_impact / obs over BCRs WITH footprint
  fpct <- function(x) if (is.na(x)) "NA" else sprintf("%+7.2f%%", x)
  cat("\n", sprintf("%-32s %10s %10s %10s %5s %5s\n",
                    "  impact (% of obs at scale)", "footprint",
                    "sub act", "BCR act", "nsub*", "nbcr*"), sep = "")
  cat(strrep("-", 82), "\n")
  for (s in ord) {
    a <- rows[[s]]
    cat(sprintf("%-32s %10s %10s %10s %5d %5d\n",
                s,
                fpct(a$pct_footprint),
                fpct(a$pct_sub_active), fpct(a$pct_bcr_active),
                a$n_sub_active, a$n_bcr_active))
  }
  if (!is.null(fa)) {
    cat(strrep("-", 82), "\n")
    cat(sprintf("%-32s %10s %10s %10s %5d %5d\n",
                "coalition 64 (LOWER floor)",
                fpct(fa$pct_footprint),
                fpct(fa$pct_sub_active), fpct(fa$pct_bcr_active),
                fa$n_sub_active, fa$n_bcr_active))
  }
  cat("  * nsub/nbcr = # subbasins/BCRs with footprint (denominators for 'act')\n")

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

# write tidy CSV ------------------------------------------------------------
out_csv <- file.path(ia_dir, "logs", "12F_singletons_summary.csv")
dir.create(dirname(out_csv), showWarnings = FALSE, recursive = TRUE)
csv_df <- do.call(rbind, csv_rows)
write.csv(csv_df, out_csv, row.names = FALSE)
cat("\nwrote", out_csv, "  (", nrow(csv_df), "rows )\n")
