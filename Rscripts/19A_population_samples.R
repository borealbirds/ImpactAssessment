# ---
# title: Per-sample species populations, every year
# author: Mannfred Boehm
# created: October 26, 2026
# ---
# For the Rosenberg et al. (2019) Fig. 1 style figures (19B on): each species' observed and
# no-industry population over the BCRs 14B kept (release filter), per (bootstrap, scenario)
# sample, for every year in YEARS (default 1995-2020). 12H's Shapley sample files stay on Fir
# (up to ~150 MB each); this one file holds what the figures need and is pulled with the CSVs:
#   sector_effects/population_samples.rds, a list of
#     species, years, n_boot, n_scen
#     obs [species x year x bootstrap]   observed population (it does not vary by scenario)
#     cf  [species x year x sample]      no-industry population; sample s is bootstrap
#                                        ceiling(s / n_scen), as in 14B and 14C
# Totals over species, and changes between years, are taken sample by sample, never from 14C's
# percentiles. Bootstrap i is the same bird model in every year, so a change between years
# pairs it.
#
# Species: those in 14B's tables in every year (the 113 reported species for 1995-2020), or
# SPECIES (comma-separated) among them. Stops unless its means equal 14C's results_summary to
# within a bird (14C rounds to whole birds).
#
# Run on Fir after 14B/14C for every year: sbatch 19A_population_samples.sh
# ---

# inside a Slurm job this is the cluster; anywhere else, the local RProject
cc    <- nzchar(Sys.getenv("SLURM_JOB_ID"))
local <- !cc

if (cc)            { ia_dir <- "/home/mannfred/scratch/impact_assessment" }
if (!cc && local)  { ia_dir <- getwd() }
if (!cc && !local) { ia_dir <- file.path("G:/Shared drives/BAM_NationalModels5", "data", "Extras",
                                          "sandbox_data", "impactassessment_sandbox") }

dt_dir <- file.path(ia_dir, "data/derived_data/density_tables")
se_dir <- file.path(ia_dir, "data/derived_data/sector_effects")
years  <- as.integer(strsplit(Sys.getenv("YEARS", "1995,2000,2005,2010,2015,2020"), ",")[[1]])
if (anyNA(years)) stop("YEARS must be comma-separated years; got '", Sys.getenv("YEARS"), "'")
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
