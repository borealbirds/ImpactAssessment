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
# stops unless its means equal 14B's national values. Sensitivity runs (DROP_Q99, DENS_MIN)
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

q05 <- function(x) unname(quantile(x, 0.05))
q95 <- function(x) unname(quantile(x, 0.95))

summary_rows <- list()
sector_rows  <- list()
for (sp in sort(unique(nat$species))) {
  ss <- readRDS(file.path(dt_dir, paste0(sp, "_", year, "_shapley_samples.rds")))
  stopifnot(identical(dimnames(ss$phi)[[2]], canonical_sectors()))
  keep_bcrs <- unique(bcr$bcr[bcr$species == sp])
  keep      <- ss$bcr %in% keep_bcrs

  phi  <- colSums(ss$phi[keep, , , drop = FALSE])                       # [sector x sample]
  boot <- rep(seq_len(ss$n_boot), each = ss$n_scen)                      # as in 14B
  obs  <- colSums(ss$obs_total[keep, , drop = FALSE])[boot]              # per sample
  tot  <- colSums(phi)
  cf   <- obs + tot                                                      # no-industry total
  chg  <- -tot                                                           # observed - cf
  pct  <- 100 * chg / cf

  chk <- nat[nat$species == sp, ]
  if (abs(mean(tot) - chk$total_HF_impact[1]) > 1 ||
      abs(mean(obs) / chk$obs_population[1] - 1) > 1e-6 ||
      any(abs(rowMeans(phi)[sectors] - chk$shapley_mean[match(sectors, chk$sector)]) > 1))
    stop(sp, ": the Shapley samples do not match 14B's shapley_national.csv - rerun 14B")

  summary_rows[[sp]] <- data.frame(
    species          = sp,
    n_bcrs           = length(keep_bcrs),
    observed         = round(mean(obs)),
    without_industry = round(mean(cf)),
    change           = round(mean(chg)),
    change_q05       = round(q05(chg)),
    change_q95       = round(q95(chg)),
    change_pct       = round(100 * mean(chg) / mean(cf), 3),
    change_pct_q05   = round(q05(pct), 3),
    change_pct_q95   = round(q95(pct), 3))

  r <- data.frame(species = sp)
  for (s in sectors) {
    p <- -100 * phi[s, ] / cf
    r[[s]]                   <- round(-100 * mean(phi[s, ]) / mean(cf), 3)
    r[[paste0(s, "_q05")]]   <- round(q05(p), 3)
    r[[paste0(s, "_q95")]]   <- round(q95(p), 3)
  }
  sector_rows[[sp]] <- r
  message(sp, ": ", length(keep_bcrs), " BCRs, change = ", round(mean(chg)), " (",
          round(100 * mean(chg) / mean(cf), 1), "% of the no-industry population)")
}

ord <- order(vapply(summary_rows, `[[`, numeric(1), "change_pct"))
results_summary <- do.call(rbind, summary_rows[ord])
results_sectors <- do.call(rbind, sector_rows[ord])
write.csv(results_summary, file.path(se_dir, "results_summary.csv"), row.names = FALSE)
write.csv(results_sectors, file.path(se_dir, "results_sectors.csv"), row.names = FALSE)
message("wrote results_summary.csv and results_sectors.csv to ", se_dir)
