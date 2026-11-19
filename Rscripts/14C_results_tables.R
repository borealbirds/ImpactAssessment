# ---
# title: Species results tables
# author: Mannfred Boehm
# ---

# from 12H's per-sample subbasin Shapley values, over the BCRs 14B kept (its release filter):
#   results_summary.csv  per species: observed and no-industry populations, the change industry
#                        made, and that change as % of the no-industry population
#   results_sectors.csv  per species: each sector's change as % of the no-industry population
#   results_bcr.csv      both, per species x BCR
# each with 5th/95th percentiles over the 3200 (bootstrap, scenario) samples. change = observed -
# counterfactual (negative = birds lost to industry), the opposite of 14B's sign. a % divides
# each sample by its OWN no-industry total, never by the mean (both come from the same bird
# model). % of the no-industry population, not the observed one, which would overstate a
# decline. stops unless its means equal 14B's to within a bird.
# runs after 14B (14B's .sh runs both); YEARS as in 14B

ia_dir <- if (nzchar(Sys.getenv("SLURM_JOB_ID"))) "/home/mannfred/scratch/impact_assessment" else getwd()
dt_dir <- file.path(ia_dir, "data/derived_data/density_tables")
se_dir <- file.path(ia_dir, "data/derived_data/sector_effects")

source(file.path(ia_dir, "Rscripts", "12E_shapley_utils.R"))

year   <- as.integer(Sys.getenv("YEARS", "2020"))
se_csv <- function(name) file.path(se_dir, paste0(name, if (year == 2020L) "" else paste0("_", year), ".csv"))

# grouped for reading: agriculture, settlement, transport, energy and extraction
sectors <- c("crop", "pasture", "built", "roads", "rail", "oil_gas", "mines", "dam_and_associated_reservoir")
stopifnot(setequal(sectors, canonical_sectors()))

nat <- read.csv(se_csv("shapley_national"))
bcr <- read.csv(se_csv("shapley_bcr"))
nat <- nat[nat$year == year, ]
bcr <- bcr[bcr$year == year, ]

q05 <- function(x) if (all(is.na(x))) NA_real_ else unname(quantile(x, 0.05, na.rm = TRUE))
q95 <- function(x) if (all(is.na(x))) NA_real_ else unname(quantile(x, 0.95, na.rm = TRUE))

# one row of results from phi [sector x sample] and each sample's observed total. a BCR whose
# weight is 0 everywhere (BAWW can41) has no birds either way, so its % is NA
results_row <- function(phi, obs) {
  tot <- colSums(phi)
  cf  <- obs + tot   # no-industry total
  chg <- -tot        # observed - cf
  pct      <- function(x) ifelse(cf > 0, 100 * x / cf, NA_real_)
  pct_mean <- function(x) if (mean(cf) > 0) round(100 * mean(x) / mean(cf), 3) else NA_real_
  r <- data.frame(observed         = round(mean(obs)),
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

# BCR 3, the parts of BCR 4 (can40-42), 5, 6 (can60-61), 7, 8, then 9 to 14
bcr_key <- function(b) { n <- as.integer(sub("^can", "", b)); ifelse(n >= 40, n / 10, n) }

summary_rows <- sector_rows <- bcr_rows <- list()
for (sp in sort(unique(nat$species))) {
  ss <- readRDS(file.path(dt_dir, paste0(sp, "_", year, "_shapley_samples.rds")))
  stopifnot(identical(dimnames(ss$phi)[[2]], canonical_sectors()))
  keep_bcrs <- unique(bcr$bcr[bcr$species == sp])
  boot <- rep(seq_len(ss$n_boot), each = ss$n_scen)   # each sample's bootstrap, as in 14B

  # phi [sector x sample] and the observed total per sample, over the BCRs in `b`
  over_bcrs <- function(b) {
    k <- ss$bcr %in% b
    list(phi = colSums(ss$phi[k, , , drop = FALSE]), obs = colSums(ss$obs_total[k, , drop = FALSE])[boot])
  }

  # 14B rounds its totals to whole birds, so compare to within a bird
  x   <- over_bcrs(keep_bcrs)
  chk <- nat[nat$species == sp, ]
  if (abs(mean(colSums(x$phi)) - chk$total_HF_impact[1]) > 1 || abs(mean(x$obs) - chk$obs_population[1]) > 1 ||
      any(abs(rowMeans(x$phi)[sectors] - chk$shapley_mean[match(sectors, chk$sector)]) > 1))
    stop(sp, ": the Shapley samples do not match 14B's shapley_national.csv - rerun 14B")

  r <- results_row(x$phi, x$obs)
  summary_rows[[sp]] <- cbind(data.frame(species = sp, n_bcrs = length(keep_bcrs)), r[setdiff(names(r), sector_cols)])
  sector_rows[[sp]]  <- cbind(data.frame(species = sp), r[sector_cols])

  for (b in keep_bcrs[order(bcr_key(keep_bcrs))]) {
    xb    <- over_bcrs(b)
    chk_b <- bcr[bcr$species == sp & bcr$bcr == b, ]
    if (abs(mean(colSums(xb$phi)) - chk_b$total_HF_impact[1]) > 1 || abs(mean(xb$obs) - chk_b$obs_population[1]) > 1)
      stop(sp, " ", b, ": the Shapley samples do not match 14B's shapley_bcr.csv - rerun 14B")
    bcr_rows[[paste(sp, b)]] <- cbind(data.frame(species = sp, bcr = b), results_row(xb$phi, xb$obs))
  }
  message(sp, ": change ", r$change, " (", round(r$change_pct, 1), "% of the no-industry population)")
}

# species from the largest decline (as %) to the largest increase
ord <- order(vapply(summary_rows, `[[`, numeric(1), "change_pct"))
results_bcr <- do.call(rbind, bcr_rows)
write.csv(do.call(rbind, summary_rows[ord]), se_csv("results_summary"), row.names = FALSE)
write.csv(do.call(rbind, sector_rows[ord]),  se_csv("results_sectors"), row.names = FALSE)
write.csv(results_bcr[order(match(results_bcr$species, names(summary_rows)[ord])), ], se_csv("results_bcr"), row.names = FALSE)
