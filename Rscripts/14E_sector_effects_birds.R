# ---
# title: Sector effects in birds, one row per species
# author: Mannfred Boehm
# ---

# 14C's results_sectors.csv in birds instead of %: each sector's national change (-Shapley
# value, so negative = birds lost to the sector) with its 5th/95th percentiles, in 14C's row
# order. a pivot of 14B's shapley_national.csv, so it runs locally after 14C. negating swaps the
# tails, q05(-v) = -q95(v), exactly, as R's default (type 7) quantiles are symmetric.
# YEARS as in 14C -> data/derived_data/sector_effects/results_sectors_birds.csv

se_dir <- file.path(getwd(), "data/derived_data/sector_effects")
year   <- as.integer(Sys.getenv("YEARS", "2020"))
se_csv <- function(name) file.path(se_dir, paste0(name, if (year == 2020L) "" else paste0("_", year), ".csv"))

# 14C's order: agriculture, settlement, transport, energy and extraction
sectors <- c("crop", "pasture", "built", "roads", "rail", "oil_gas", "mines", "dam_and_associated_reservoir")

nat <- read.csv(se_csv("shapley_national"))
nat <- nat[nat$year == year, ]
pct <- read.csv(se_csv("results_sectors"))   # 14C, for the row order
stopifnot(setequal(nat$species, pct$species), setequal(nat$sector, sectors), !anyDuplicated(nat[c("species", "sector")]))

out <- data.frame(species = pct$species)
for (s in sectors) {
  ns <- nat[nat$sector == s, ][match(out$species, nat$species[nat$sector == s]), ]
  out[[s]]                 <- round(-ns$shapley_mean)
  out[[paste0(s, "_q05")]] <- round(-ns$shapley_q95)
  out[[paste0(s, "_q95")]] <- round(-ns$shapley_q05)
}

# the sectors add up to the total change (to within the rounding of 14B's CSVs), and every
# sign agrees with 14C's %
tot <- nat[!duplicated(nat$species), ]
tot <- -tot$total_HF_impact[match(out$species, tot$species)]
sum_unrounded <- vapply(out$species, function(sp) -sum(nat$shapley_mean[nat$species == sp]), numeric(1))
bad <- out$species[abs(sum_unrounded - tot) > 1]
if (length(bad)) stop("sectors do not add up to the total change for: ", toString(bad))
for (s in sectors) {
  clash <- sign(out[[s]]) * sign(pct[[s]]) < 0
  if (any(clash, na.rm = TRUE)) stop(s, ": sign differs from results_sectors.csv for ", toString(out$species[which(clash)]))
}

write.csv(out, se_csv("results_sectors_birds"), row.names = FALSE)
