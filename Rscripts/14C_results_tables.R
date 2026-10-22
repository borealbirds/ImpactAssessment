# ---
# title: Species results tables
# author: Mannfred Boehm
# ---
# Two tables, one row per species, from 12H's per-sample subbasin Shapley values:
#   results_summary.csv  — observed and no-industry populations, the change industry made
#                          (-v(N)), and that change as % of the no-industry population, each
#                          with its 5th/95th percentiles
#   results_sectors.csv  — each sector's change (-Shapley value) as % of the no-industry
#                          population, with its 5th/95th percentiles; the sectors sum to the total %
# and results_bcr.csv, both of those for each species x BCR, % of the BCR's no-industry population.
# Reported as change = observed - counterfactual, so negative = birds lost to industry, as the
# maps show it (18). The % is relative to the no-industry population: dividing by the observed
# one instead makes "-16.9%" read as a 16.9% decline when it is 14.5% (CAWA), and GRSP's near
# doubling (+94%) read as +48%.
#
# Percentiles are over the 3200 (bootstrap, scenario) samples. A % divides each sample by ITS OWN
# no-industry total (its bootstrap's observed total + its v(N)), never by the mean: both come from
# the same bird model, so dividing quantiles by the mean misstates the interval.
#
# Run after 14B. It keeps only the BCRs in 14B's shapley_bcr.csv (14B's release filter) and
# stops unless its means equal 14B's national and BCR values. Sensitivity runs (DROP_Q99, DENS_MIN)
# read and write their own directories, as in 14B.
# ---

# inside a Slurm job this is the cluster; anywhere else, the local RProject
cc    <- nzchar(Sys.getenv("SLURM_JOB_ID"))
local <- !cc

if (cc)            { ia_dir <- "/home/mannfred/scratch/impact_assessment" }
if (!cc && local)  { ia_dir <- getwd() }
if (!cc && !local) { ia_dir <- file.path("G:/Shared drives/BAM_NationalModels5", "data", "Extras",
                                          "sandbox_data", "impactassessment_sandbox") }

source(file.path(ia_dir, "Rscripts", "12E_shapley_utils.R"))

sens   <- sensitivity_run()
dt_dir <- file.path(ia_dir, "data/derived_data", paste0("density_tables", sens$suffix))
se_dir <- file.path(ia_dir, "data/derived_data", paste0("sector_effects", sens$suffix))
year   <- 2020

# grouped for reading: agriculture, settlement, transport, energy and extraction
sectors <- c("crop", "pasture", "built", "roads", "rail", "oil_gas", "mines",
             "dam_and_associated_reservoir")
stopifnot(setequal(sectors, canonical_sectors()))

nat <- read.csv(file.path(se_dir, "shapley_national.csv"))
bcr <- read.csv(file.path(se_dir, "shapley_bcr.csv"))
nat <- nat[nat$year == year, ]
bcr <- bcr[bcr$year == year, ]

q05 <- function(x) if (all(is.na(x))) NA_real_ else unname(quantile(x, 0.05, na.rm = TRUE))
q95 <- function(x) if (all(is.na(x))) NA_real_ else unname(quantile(x, 0.95, na.rm = TRUE))

# one row of results from phi [sector x sample] and each sample's observed total. A BCR whose
# weight is 0 on every subbasin (BAWW can41) has no birds either way, so its % is NA
results_row <- function(phi, obs) {
  tot <- colSums(phi)
  cf  <- obs + tot                                                       # no-industry total
  chg <- -tot                                                            # observed - cf
  pct <- function(x) ifelse(cf > 0, 100 * x / cf, NA_real_)
  pct_mean <- function(x) if (mean(cf) > 0) round(100 * mean(x) / mean(cf), 3) else NA_real_
  r <- data.frame(
    observed         = round(mean(obs)),
    without_industry = round(mean(cf)),
    change           = round(mean(chg)),
    change_q05       = round(q05(chg)),
    change_q95       = round(q95(chg)),
    change_pct       = pct_mean(chg),
    change_pct_q05   = round(q05(pct(chg)), 3),
    change_pct_q95   = round(q95(pct(chg)), 3))
  for (s in sectors) {
    r[[s]]                 <- pct_mean(-phi[s, ])
    r[[paste0(s, "_q05")]] <- round(q05(pct(-phi[s, ])), 3)
    r[[paste0(s, "_q95")]] <- round(q95(pct(-phi[s, ])), 3)
  }
  r
}
sector_cols <- c(rbind(sectors, paste0(sectors, "_q05"), paste0(sectors, "_q95")))

# BCR 3, the subunits of BCR 4 (can40-42), 5, 6 (can60-61), 7, 8, then 9 to 14
bcr_key <- function(b) { n <- as.integer(sub("^can", "", b)); ifelse(n >= 40, n / 10, n) }

summary_rows <- list()
sector_rows  <- list()
bcr_rows     <- list()
for (sp in sort(unique(nat$species))) {
  ss <- readRDS(file.path(dt_dir, paste0(sp, "_", year, "_shapley_samples.rds")))
  stopifnot(identical(dimnames(ss$phi)[[2]], canonical_sectors()))
  keep_bcrs <- unique(bcr$bcr[bcr$species == sp])
  keep      <- ss$bcr %in% keep_bcrs

  phi  <- colSums(ss$phi[keep, , , drop = FALSE])                       # [sector x sample]
  boot <- rep(seq_len(ss$n_boot), each = ss$n_scen)                      # as in 14B
  obs  <- colSums(ss$obs_total[keep, , drop = FALSE])[boot]              # per sample

  # 14B rounds its national totals to whole birds, so compare to within a bird (a relative 1e-6
  # failed EATO's 198,987.4 vs 198,987)
  chk <- nat[nat$species == sp, ]
  if (abs(mean(colSums(phi)) - chk$total_HF_impact[1]) > 1 ||
      abs(mean(obs) - chk$obs_population[1]) > 1 ||
      any(abs(rowMeans(phi)[sectors] - chk$shapley_mean[match(sectors, chk$sector)]) > 1))
    stop(sp, ": the Shapley samples do not match 14B's shapley_national.csv - rerun 14B")

  r <- results_row(phi, obs)
  summary_rows[[sp]] <- cbind(data.frame(species = sp, n_bcrs = length(keep_bcrs)),
                              r[setdiff(names(r), sector_cols)])
  sector_rows[[sp]]  <- cbind(data.frame(species = sp), r[sector_cols])

  # the same per BCR, for the per-species BCR tables
  for (b in keep_bcrs[order(bcr_key(keep_bcrs))]) {
    kb    <- ss$bcr == b
    phi_b <- colSums(ss$phi[kb, , , drop = FALSE])
    obs_b <- colSums(ss$obs_total[kb, , drop = FALSE])[boot]
    chk_b <- bcr[bcr$species == sp & bcr$bcr == b, ]
    if (abs(mean(colSums(phi_b)) - chk_b$total_HF_impact[1]) > 1 ||
        abs(mean(obs_b) - chk_b$obs_population[1]) > 1)
      stop(sp, " ", b, ": the Shapley samples do not match 14B's shapley_bcr.csv - rerun 14B")
    bcr_rows[[paste(sp, b)]] <- cbind(data.frame(species = sp, bcr = b), results_row(phi_b, obs_b))
  }
  message(sp, ": ", length(keep_bcrs), " BCRs, change = ", r$change, " (",
          round(r$change_pct, 1), "% of the no-industry population)")
}

ord <- order(vapply(summary_rows, `[[`, numeric(1), "change_pct"))
results_summary <- do.call(rbind, summary_rows[ord])
results_sectors <- do.call(rbind, sector_rows[ord])
results_bcr     <- do.call(rbind, bcr_rows)
results_bcr     <- results_bcr[order(match(results_bcr$species, names(summary_rows)[ord])), ]
write.csv(results_summary, file.path(se_dir, "results_summary.csv"), row.names = FALSE)
write.csv(results_sectors, file.path(se_dir, "results_sectors.csv"), row.names = FALSE)
write.csv(results_bcr,     file.path(se_dir, "results_bcr.csv"),     row.names = FALSE)
message("wrote results_summary.csv, results_sectors.csv and results_bcr.csv to ", se_dir)
