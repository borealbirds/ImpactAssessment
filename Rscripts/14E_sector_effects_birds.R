# ---
# title: Sector effects in birds, one row per species
# author: Mannfred Boehm
# ---
# results_sectors_birds.csv: each sector's national change in birds, with its 5th/95th
# percentiles, in the layout of 14C's results_sectors.csv (which gives the same as % of the
# no-industry population), and in its row order. Change = -Shapley value, as in 14C, so
# negative = birds lost to the sector.
#
# A pivot of 14B's shapley_national.csv, so it runs locally after 14C (14C's own input,
# density_tables/, stays on Fir). Negating a sample swaps its tails: q05(-v) = -q95(v). That is
# exact, because R's default (type 7) quantiles are symmetric. One year, YEARS (default 2020), as in
# 14C; another year's table carries _{year}. Production only, so no sensitivity runs.
# ---

ia_dir <- getwd()
se_dir <- file.path(ia_dir, "data/derived_data/sector_effects")
year   <- suppressWarnings(as.integer(trimws(Sys.getenv("YEARS", "2020"))))
if (is.na(year))
  stop("YEARS must name one year, as for 14C; got '", Sys.getenv("YEARS"), "'")
se_csv <- function(name) file.path(se_dir, paste0(name, if (year == 2020L) "" else paste0("_", year), ".csv"))

# 14C's order: agriculture, settlement, transport, energy and extraction
sectors <- c("crop", "pasture", "built", "roads", "rail", "oil_gas", "mines",
             "dam_and_associated_reservoir")

nat <- read.csv(se_csv("shapley_national"))
nat <- nat[nat$year == year, ]
pct <- read.csv(se_csv("results_sectors"))                              # 14C, for the row order
stopifnot(setequal(nat$species, pct$species), setequal(nat$sector, sectors),
          !anyDuplicated(nat[c("species", "sector")]))

out <- data.frame(species = pct$species)
for (s in sectors) {
  ns <- nat[nat$sector == s, ][match(out$species, nat$species[nat$sector == s]), ]
  out[[s]]                 <- round(-ns$shapley_mean)
  out[[paste0(s, "_q05")]] <- round(-ns$shapley_q95)
  out[[paste0(s, "_q95")]] <- round(-ns$shapley_q05)
}

# the sectors add up to the total change (14B rounds the total to whole birds and the Shapley
# values to 2 decimals), and every sign agrees with 14C's %
tot <- nat[!duplicated(nat$species), ]
tot <- -tot$total_HF_impact[match(out$species, tot$species)]
sum_unrounded <- vapply(out$species, function(sp)
  -sum(nat$shapley_mean[nat$species == sp]), numeric(1))
bad <- out$species[abs(sum_unrounded - tot) > 1]
if (length(bad)) stop("sectors do not add up to the total change for: ", toString(bad))
for (s in sectors) {
  clash <- sign(out[[s]]) * sign(pct[[s]]) < 0
  if (any(clash, na.rm = TRUE))
    stop(s, ": sign differs from results_sectors.csv for ", toString(out$species[which(clash)]))
}

write.csv(out, se_csv("results_sectors_birds"), row.names = FALSE)
message("wrote ", basename(se_csv("results_sectors_birds")), ": ", nrow(out), " species x ",
        length(sectors), " sectors, ", year)
