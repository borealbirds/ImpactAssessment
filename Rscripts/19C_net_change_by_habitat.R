# ---
# title: Net change in abundance by habitat group, 1995-2020
# author: Mannfred Boehm
# created: October 26, 2026
# ---
# Our version of Rosenberg et al. (2019, Science 366:120) Fig. 1B: 19B's net change since the
# first year, observed and without industry, for each of 14D's habitat groups. One panel per
# group, all on one y-axis so the groups' sizes compare, as Rosenberg's biomes share one axis.
# Each group sums all its species, increasing and declining.
#
# Uncertainty as in 19B: each sample's group total minus the SAME sample's first-year total, then
# the 5th/95th percentiles. Colours as 19B.
#
# Run locally after 14D, with 19A's population_samples.rds pulled:
#   Rscript Rscripts/19C_net_change_by_habitat.R
# Writes output_figures/net_change_by_habitat_1995_2020.png and its values to
# sector_effects/net_change_by_habitat_1995_2020.csv.
# ---

suppressPackageStartupMessages(library(ggplot2))

ia_dir <- getwd()   # the local RProject
se_dir <- file.path(ia_dir, "data", "derived_data", "sector_effects")
ps     <- readRDS(file.path(se_dir, "population_samples.rds"))
master <- read.csv(file.path(ia_dir, "data", "raw_data", "anna_data", "master_species_list.csv"))

scen_cols <- c("Observed" = "#2a78d6", "Without industry" = "#eb6834")   # as 19B

# 14D's groups of the reported species, forest types first
groups <- c("conifer" = "Conifer", "mixedwood" = "Mixedwood", "deciduous" = "Deciduous",
            "forest no type" = "Forest, no type", "non-forest" = "Non-forest")
master <- master[master$include, ]
if (length(odd <- setdiff(master$group, names(groups))))
  stop("14D group(s) not in this script: ", paste(odd, collapse = ", "))
if (length(miss <- setdiff(master$species, ps$species)))
  stop("population_samples.rds lacks ", paste(miss, collapse = ", "), " - rerun 19A")

years <- ps$years
y0    <- which.min(years)
stub  <- paste0("net_change_by_habitat_", years[y0], "_", max(years))

# m: a group's total [year x bootstrap] or [year x sample]
band <- function(m, scenario) {
  d <- sweep(m, 2, m[y0, ])
  data.frame(year     = years,
             scenario = scenario,
             total    = round(rowMeans(m)),
             change   = round(rowMeans(d)),
             q05      = round(apply(d, 1, quantile, 0.05)),
             q95      = round(apply(d, 1, quantile, 0.95)))
}

rows <- list()
for (g in names(groups)) {
  s <- master$species[master$group == g]
  r <- rbind(band(colSums(ps$obs[s, , , drop = FALSE]), "Observed"),
             band(colSums(ps$cf[s, , , drop = FALSE]),  "Without industry"))
  rows[[g]] <- cbind(data.frame(group = g, n_species = length(s)), r)
}
d <- do.call(rbind, rows)
write.csv(d[order(match(d$group, names(groups)), d$year), ],
          file.path(se_dir, paste0(stub, ".csv")), row.names = FALSE)

n_sp    <- tapply(d$n_species, d$group, `[`, 1)[names(groups)]
d$panel <- factor(paste0(groups[d$group], " (", n_sp[d$group], " species)"),
                  levels = paste0(groups, " (", n_sp, " species)"))

p <- ggplot(d, aes(year, change / 1e6, colour = scenario, fill = scenario)) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_ribbon(aes(ymin = q05 / 1e6, ymax = q95 / 1e6), alpha = 0.18, colour = NA) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.5) +
  facet_wrap(~ panel, ncol = 3) +
  scale_colour_manual(values = scen_cols) +
  scale_fill_manual(values = scen_cols) +
  scale_x_continuous(breaks = years, labels = function(x) ifelse(x %% 10 == 0, x, "")) +
  labs(x = NULL, y = paste0("Net change in abundance since ", years[y0], "\n(millions of birds)"),
       colour = NULL, fill = NULL,
       caption = paste0(nrow(master), " species. Shading: 5th to 95th percentiles.")) +
  theme_classic(base_size = 10) +
  theme(legend.position        = "inside",
        legend.position.inside = c(0.835, 0.25),
        legend.justification   = c(0.5, 0.5),
        strip.background       = element_blank(),
        strip.text             = element_text(hjust = 0, face = "bold"),
        panel.spacing          = unit(1, "lines"),
        axis.text              = element_text(colour = "grey25"),
        axis.line              = element_line(colour = "grey60", linewidth = 0.3),
        axis.ticks             = element_line(colour = "grey60", linewidth = 0.3),
        plot.caption           = element_text(colour = "grey40", size = 8, hjust = 0))

out <- file.path(ia_dir, "output_figures", paste0(stub, ".png"))
ggsave(out, p, width = 8, height = 5.5, dpi = 300)
message("wrote ", out)
