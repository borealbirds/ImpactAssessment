# ---
# title: Percent change in abundance by habitat group, 1995-2020
# author: Mannfred Boehm
# created: October 26, 2026
# ---
# Our version of Rosenberg et al. (2019, Science 366:120) Fig. 1C: each of 14D's habitat groups'
# change in abundance from the first year to the last, as % of its first-year total, observed and
# without industry (each relative to its OWN first-year total, as in 19B and 19C). Groups in
# 19C's order.
#
# As in 14C: the bar is the mean change over the mean first-year total; the whiskers are the
# 5th/95th percentiles of each sample's change over THAT sample's first-year total (bootstrap i
# is one bird model in every year). Colours as 19B.
#
# Run locally after 14D, with 19A's population_samples.rds pulled:
#   Rscript Rscripts/19D_percent_change_by_habitat.R
# Writes output_figures/percent_change_by_habitat_1995_2020.png and its values to
# sector_effects/percent_change_by_habitat_1995_2020.csv.
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
stub  <- paste0("percent_change_by_habitat_", years[y0], "_", years[y1])

# m: a group's total [year x bootstrap] or [year x sample]
pct_row <- function(m, scenario) {
  chg <- m[y1, ] - m[y0, ]
  pct <- 100 * chg / m[y0, ]
  data.frame(scenario    = scenario,
             first_total = round(mean(m[y0, ])),
             last_total  = round(mean(m[y1, ])),
             pct         = round(100 * mean(chg) / mean(m[y0, ]), 2),
             q05         = round(unname(quantile(pct, 0.05)), 2),
             q95         = round(unname(quantile(pct, 0.95)), 2))
}

rows <- list()
for (g in names(groups)) {
  s <- master$species[master$group == g]
  r <- rbind(pct_row(colSums(ps$obs[s, , , drop = FALSE]), "Observed"),
             pct_row(colSums(ps$cf[s, , , drop = FALSE]),  "Without industry"))
  rows[[g]] <- cbind(data.frame(group = g, n_species = length(s)), r)
}
d <- do.call(rbind, rows)
write.csv(d, file.path(se_dir, paste0(stub, ".csv")), row.names = FALSE)
print(d, row.names = FALSE)

# first group on top, and within a group the observed bar above the no-industry one
lab        <- paste0(groups[d$group], " (", d$n_species, ")")
d$label    <- factor(lab, levels = rev(unique(lab)))
d$scenario <- factor(d$scenario, levels = rev(names(scen_cols)))
dodge      <- position_dodge(width = 0.8)

p <- ggplot(d, aes(x = pct, y = label, fill = scenario)) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_col(position = dodge, width = 0.72) +
  geom_errorbar(aes(xmin = q05, xmax = q95), position = dodge, width = 0.25,
                colour = "grey25", linewidth = 0.35) +
  scale_fill_manual(values = scen_cols, breaks = names(scen_cols)) +
  labs(x = paste0("Change in abundance, ", years[y0], " to ", years[y1], " (%)"), y = NULL,
       fill = NULL,
       caption = paste0(nrow(master), " species; number per group in brackets. ",
                        "Whiskers: 5th to 95th percentiles.")) +
  theme_classic(base_size = 10) +
  theme(legend.position = "top",
        legend.justification = "left",
        axis.text            = element_text(colour = "grey25"),
        axis.text.y          = element_text(size = 10),
        axis.line            = element_line(colour = "grey60", linewidth = 0.3),
        axis.ticks           = element_line(colour = "grey60", linewidth = 0.3),
        axis.line.y          = element_blank(),
        axis.ticks.y         = element_blank(),
        plot.caption         = element_text(colour = "grey40", size = 8, hjust = 0))

out <- file.path(ia_dir, "output_figures", paste0(stub, ".png"))
ggsave(out, p, width = 6, height = 4, dpi = 300)
message("wrote ", out)
