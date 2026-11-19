# ---
# title: Per-sample species populations, every year
# author: Mannfred Boehm
# created: October 26, 2026
# ---

# for the Rosenberg et al. (2019) Fig. 1-style figures (19B on): each species' observed and
# no-industry population over the BCRs 14B kept, per (bootstrap, scenario) sample, for every
# year in YEARS (default 1995-2020). 12H's Shapley sample files stay on Fir; this one file holds
# what the figures need and is pulled with the CSVs:
#   sector_effects/population_samples.rds: species, years, n_boot, n_scen,
#     obs [species x year x bootstrap]  (the observed population does not vary by scenario)
#     cf  [species x year x sample]     (sample s is bootstrap ceiling(s / n_scen), as in 14B)
# totals over species and changes between years are taken sample by sample; bootstrap i is the
# same bird model in every year. species: those in 14B's tables in every year, or SPECIES among
# them. stops unless its means equal 14C's results_summary to within a bird.
# on Fir after 14B/14C for every year: sbatch 19A_population_samples.sh

ia_dir <- if (nzchar(Sys.getenv("SLURM_JOB_ID"))) "/home/mannfred/scratch/impact_assessment" else getwd()
dt_dir <- file.path(ia_dir, "data/derived_data/density_tables")
se_dir <- file.path(ia_dir, "data/derived_data/sector_effects")
years  <- as.integer(strsplit(Sys.getenv("YEARS", "1995,2000,2005,2010,2015,2020"), ",")[[1]])
se_csv <- function(name, y) file.path(se_dir, paste0(name, if (y == 2020L) "" else paste0("_", y), ".csv"))

# 14B's kept BCRs per species, per year
bcr <- lapply(years, function(y) { d <- read.csv(se_csv("shapley_bcr", y)); d[d$year == y, ] })
species <- sort(Reduce(intersect, lapply(bcr, function(d) unique(d$species))))
if (nzchar(Sys.getenv("SPECIES"))) {
  want <- strsplit(Sys.getenv("SPECIES"), ",")[[1]]
  if (length(miss <- setdiff(want, species)))
    stop("not in 14B's tables for every year: ", paste(miss, collapse = ", "))
  species <- sort(want)
}
message(length(species), " species x ", length(years), " years")

obs <- cf <- NULL
for (j in seq_along(years)) {
  y    <- years[j]
  summ <- read.csv(se_csv("results_summary", y))
  for (sp in species) {
    ss   <- readRDS(file.path(dt_dir, paste0(sp, "_", y, "_shapley_samples.rds")))
    keep <- ss$bcr %in% bcr[[j]]$bcr[bcr[[j]]$species == sp]
    if (is.null(obs)) {
      n_boot <- ss$n_boot
      n_scen <- ss$n_scen
      obs <- array(NA_real_, c(length(species), length(years), n_boot),
                   dimnames = list(species, years, NULL))
      cf  <- array(NA_real_, c(length(species), length(years), n_boot * n_scen),
                   dimnames = list(species, years, NULL))
    }
    if (ss$n_boot != n_boot || ss$n_scen != n_scen)
      stop(sp, " ", y, ": ", ss$n_boot, " x ", ss$n_scen, " samples, not ", n_boot, " x ", n_scen)

    boot <- rep(seq_len(n_boot), each = n_scen)                          # as in 14B and 14C
    o    <- colSums(ss$obs_total[keep, , drop = FALSE])                   # per bootstrap
    v    <- colSums(ss$phi[keep, , , drop = FALSE], dims = 2)             # v(N) per sample
    obs[sp, j, ] <- o
    cf[sp, j, ]  <- o[boot] + v

    chk <- summ[summ$species == sp, ]
    if (nrow(chk) != 1 || abs(mean(o) - chk$observed) > 1 ||
        abs(mean(cf[sp, j, ]) - chk$without_industry) > 1)
      stop(sp, " ", y, ": the Shapley samples do not match 14C's ", basename(se_csv("results_summary", y)))
    rm(ss)
  }
  message(y, ": observed ", round(sum(apply(obs[, j, , drop = FALSE], 1, mean)) / 1e6, 1),
          " M, no industry ", round(sum(apply(cf[, j, , drop = FALSE], 1, mean)) / 1e6, 1), " M")
}

out <- file.path(se_dir, "population_samples.rds")
saveRDS(list(species = species, years = years, n_boot = n_boot, n_scen = n_scen, obs = obs, cf = cf), out)
message("wrote ", out)
