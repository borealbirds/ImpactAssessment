# ---
# title: Net change in abundance, 1995-2020
# author: Mannfred Boehm
# created: October 26, 2026
# ---
# Our version of Rosenberg et al. (2019, Science 366:120) Fig. 1A, one figure per panel (19B on):
# the summed population of the reported species (species_113.txt), observed and without
# industry, as the net change since the first year, each relative to its OWN first-year total.
# As in Rosenberg's, the total sums every species, increasing and declining, so the most abundant
# species dominate it.
#
# Uncertainty from 19A's per-sample populations: each sample's species total minus the SAME
# sample's first-year total (bootstrap i is one bird model in every year), then the 5th/95th
# percentiles, as in 14C. The observed population varies only by bootstrap (32 samples), the
# no-industry one by bootstrap and BART draw (3200). 12F seeds the BART draw by species, so the
# band treats each species' backfill draw as independent of the others'.
#
# Every year's counterfactual removes the 2020 footprint (07 ran 1995-2015 with
# HF_MASK_YEAR=2020), so the two lines differ by how that footprint's effect changed over time,
# not by the footprint's growth.
#
# Run locally after pulling sector_effects/population_samples.rds (19A, on Fir):
#   Rscript Rscripts/19B_net_change_trajectory.R
# Writes output_figures/net_change_1995_2020.png and its values to
# sector_effects/net_change_1995_2020.csv.
# ---

suppressPackageStartupMessages(library(ggplot2))

ia_dir <- getwd()   # the local RProject
se_dir <- file.path(ia_dir, "data", "derived_data", "sector_effects")
ps     <- readRDS(file.path(se_dir, "population_samples.rds"))

# shared by every panel (19B on); validated with the dataviz skill's checker
scen_cols <- c("Observed" = "#2a78d6", "Without industry" = "#eb6834")

species <- trimws(readLines(file.path(ia_dir, "Rscripts", "species_113.txt")))
species <- species[nzchar(species) & !startsWith(species, "#")]
if (length(miss <- setdiff(species, ps$species)))
  stop("population_samples.rds lacks ", paste(miss, collapse = ", "), " - rerun 19A")

years <- ps$years
y0    <- which.min(years)
stub  <- paste0("net_change_", years[y0], "_", max(years))

# summed over species: [year x bootstrap] and [year x sample], then each sample's change
obs   <- colSums(ps$obs[species, , , drop = FALSE])
cf    <- colSums(ps$cf[species, , , drop = FALSE])
d_obs <- sweep(obs, 2, obs[y0, ])
d_cf  <- sweep(cf, 2, cf[y0, ])

band <- function(m, scenario) data.frame(
  year     = years,
  scenario = scenario,
  total    = round(rowMeans(if (scenario == "Observed") obs else cf)),
  change   = round(rowMeans(m)),
  q05      = round(apply(m, 1, quantile, 0.05)),
  q95      = round(apply(m, 1, quantile, 0.95)))
d <- rbind(band(d_obs, "Observed"), band(d_cf, "Without industry"))
d <- d[order(d$year), ]
write.csv(d, file.path(se_dir, paste0(stub, ".csv")), row.names = FALSE)
print(d, row.names = FALSE)

p <- ggplot(d, aes(year, change / 1e6, colour = scenario, fill = scenario)) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_ribbon(aes(ymin = q05 / 1e6, ymax = q95 / 1e6), alpha = 0.18, colour = NA) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2) +
  scale_colour_manual(values = scen_cols) +
  scale_fill_manual(values = scen_cols) +
  scale_x_continuous(breaks = years) +
  labs(x = NULL, y = paste0("Net change in abundance since ", years[y0], "\n(millions of birds)"),
       colour = NULL, fill = NULL,
       caption = paste0(length(species), " species. Shading: 5th to 95th percentiles.")) +
  theme_classic(base_size = 11) +
  theme(legend.position      = "inside",
        legend.position.inside = c(0.02, 0.98),
        legend.justification = c(0, 1),
        axis.text            = element_text(colour = "grey25"),
        axis.line            = element_line(colour = "grey60", linewidth = 0.3),
        axis.ticks           = element_line(colour = "grey60", linewidth = 0.3),
        plot.caption         = element_text(colour = "grey40", size = 8, hjust = 0))

out <- file.path(ia_dir, "output_figures", paste0(stub, ".png"))
ggsave(out, p, width = 6, height = 4, dpi = 300)
message("wrote ", out)
