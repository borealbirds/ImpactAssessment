# ---
# title: Share of species declining by habitat group, 1995-2020
# author: Mannfred Boehm
# created: October 26, 2026
# ---
# Our version of Rosenberg et al. (2019, Science 366:120) Fig. 1D: the share of each of 14D's
# habitat groups' species that declined from the first year to the last, observed and without
# industry, as paired bars (Rosenberg's pies, one per biome, cannot pair two scenarios). Groups
# in 19C's order.
#
# A species declines when its mean change is negative, as Rosenberg counts negative trend
# trajectories. The solid part of a bar is the species whose whole 5th-95th percentile interval of
# change is below 0, the pale rest those declining on the mean only. A species' change is taken
# per sample (bootstrap i is one bird model in every year): 32 observed samples, 3200 without
# industry. Colours as 19B.
#
# Run locally after 14D, with 19A's population_samples.rds pulled:
#   Rscript Rscripts/19E_species_declining_by_habitat.R
# Writes output_figures/species_declining_by_habitat_1995_2020.png and its values to
# sector_effects/species_declining_by_habitat_1995_2020.csv.
# ---

suppressPackageStartupMessages(library(ggplot2))

ia_dir <- getwd()   # the local RProject
se_dir <- file.path(ia_dir, "data", "derived_data", "sector_effects")
ps     <- readRDS(file.path(se_dir, "population_samples.rds"))
master <- read.csv(file.path(ia_dir, "data", "raw_data", "anna_data", "master_species_list.csv"))

scen_cols <- c("Observed" = "#2a78d6", "Without industry" = "#eb6834")   # as 19B

# as 19C
groups <- c("conifer" = "Conifer", "mixedwood" = "Mixedwood", "deciduous" = "Deciduous",
            "forest no type" = "Forest, no type", "non-forest" = "Non-forest")
master <- master[master$include, ]
if (length(odd <- setdiff(master$group, names(groups))))
  stop("14D group(s) not in this script: ", paste(odd, collapse = ", "))
if (length(miss <- setdiff(master$species, ps$species)))
  stop("population_samples.rds lacks ", paste(miss, collapse = ", "), " - rerun 19A")

years <- ps$years
y0    <- which.min(years)
y1    <- which.max(years)
stub  <- paste0("species_declining_by_habitat_", years[y0], "_", years[y1])

# each species' change, per sample: [species x bootstrap] and [species x sample]
sp      <- master$species
chg_obs <- ps$obs[sp, y1, ] - ps$obs[sp, y0, ]
chg_cf  <- ps$cf[sp, y1, ]  - ps$cf[sp, y0, ]
declining <- function(chg) data.frame(species   = rownames(chg),
                                      mean_down = rowMeans(chg) < 0,
                                      all_down  = rowMeans(chg) < 0 & apply(chg, 1, quantile, 0.95) < 0)
dec <- rbind(cbind(scenario = "Observed",         declining(chg_obs)),
             cbind(scenario = "Without industry", declining(chg_cf)))
dec$group <- master$group[match(dec$species, master$species)]

d <- aggregate(cbind(n_declining = mean_down, n_interval_below_0 = all_down) ~ group + scenario, dec, sum)
d$n_species <- as.vector(table(master$group)[d$group])
d$pct_declining        <- round(100 * d$n_declining / d$n_species, 1)
d$pct_interval_below_0 <- round(100 * d$n_interval_below_0 / d$n_species, 1)
d <- d[order(match(d$group, names(groups)), d$scenario), ]
write.csv(d, file.path(se_dir, paste0(stub, ".csv")), row.names = FALSE)
print(d, row.names = FALSE)

flip <- merge(dec[dec$scenario == "Observed", c("species", "group", "mean_down")],
              dec[dec$scenario == "Without industry", c("species", "mean_down")],
              by = "species", suffixes = c("_obs", "_cf"))
flip <- flip[flip$mean_down_obs != flip$mean_down_cf, ]
if (nrow(flip))
  message("declining in one scenario only: ",
          paste0(flip$species, " (", flip$group, ", ",
                 ifelse(flip$mean_down_obs, "observed", "without industry"), ")", collapse = ", "))

# first group on top, and within a group the observed bar above the no-industry one
lab        <- paste0(groups[d$group], " (", d$n_species, ")")
d$label    <- factor(lab, levels = rev(unique(paste0(groups, " (", table(master$group)[names(groups)], ")"))))
d$scenario <- factor(d$scenario, levels = rev(names(scen_cols)))
dodge      <- position_dodge(width = 0.8)
certainty  <- c("Mean change below 0" = 0.35, "Whole interval below 0" = 1)

p <- ggplot(d, aes(y = label, fill = scenario)) +
  geom_vline(xintercept = 50, colour = "grey75", linewidth = 0.3, linetype = "dashed") +
  geom_col(aes(x = pct_declining, alpha = names(certainty)[1]), position = dodge, width = 0.72) +
  geom_col(aes(x = pct_interval_below_0, alpha = names(certainty)[2]), position = dodge, width = 0.72) +
  scale_fill_manual(values = scen_cols, breaks = names(scen_cols)) +
  scale_alpha_manual(values = certainty, breaks = names(certainty)) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 25), expand = c(0, 0)) +
  guides(fill = guide_legend(order = 1),
         alpha = guide_legend(order = 2, override.aes = list(fill = "grey35"))) +
  labs(x = paste0("Species declining, ", years[y0], " to ", years[y1], " (%)"), y = NULL,
       fill = NULL, alpha = NULL,
       caption = paste0(nrow(master), " species; number per group in brackets. ",
                        "Interval: 5th to 95th percentiles of the change.")) +
  theme_classic(base_size = 10) +
  theme(legend.position      = "top",
        legend.justification = "left",
        legend.box           = "vertical",
        legend.box.just      = "left",
        legend.spacing.y     = unit(0, "pt"),
        axis.text            = element_text(colour = "grey25"),
        axis.text.y          = element_text(size = 10),
        axis.line            = element_line(colour = "grey60", linewidth = 0.3),
        axis.ticks           = element_line(colour = "grey60", linewidth = 0.3),
        axis.line.y          = element_blank(),
        axis.ticks.y         = element_blank(),
        plot.margin          = margin(5.5, 12, 5.5, 5.5),
        plot.caption         = element_text(colour = "grey40", size = 8, hjust = 0))

out <- file.path(ia_dir, "output_figures", paste0(stub, ".png"))
ggsave(out, p, width = 6, height = 4.3, dpi = 300)
message("wrote ", out)
