# ---
# title: Percent change in abundance by habitat group, 1995-2020
# author: Mannfred Boehm
# created: October 26, 2026
# ---

# our version of Rosenberg et al. (2019) Fig. 1C: each of 14D's habitat groups' change from the
# first year to the last, as % of its first-year total, observed and without industry. groups top
# to bottom in the order 19C's observed lines end; colour = group, observed solid and filled,
# without industry dashed and open (19B's legend serves). as 14C: the point is the mean change
# over the mean first-year total, the interval the 5th/95th percentiles of each sample's change
# over THAT sample's first-year total.
# runs locally after 14D, with 19A's population_samples.rds pulled
#   -> output_figures/figure1/fig1C_percent_change_by_habitat_1995_2020.png
#   -> sector_effects/percent_change_by_habitat_1995_2020.csv

library(ggplot2)

ia_dir <- getwd()
se_dir <- file.path(ia_dir, "data", "derived_data", "sector_effects")
ps     <- readRDS(file.path(se_dir, "population_samples.rds"))
master <- read.csv(file.path(ia_dir, "data", "raw_data", "anna_data", "master_species_list.csv"))

scen_lines  <- c("Observed" = "solid", "Without industry" = "dashed")   # as 19C
scen_shapes <- c("Observed" = 16,      "Without industry" = 21)

# as 19C
groups   <- c("conifer" = "Conifer", "mixedwood" = "Mixedwood", "deciduous" = "Deciduous",
              "forest no type" = "Forest, no type", "non-forest" = "Non-forest")
grp_cols <- c("conifer" = "#009E73", "mixedwood" = "#CC79A7", "deciduous" = "#E69F00",
              "forest no type" = "#56B4E9", "non-forest" = "#000000")
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

# groups top to bottom as 19C's observed lines end (largest gain first), and within a group
# observed above without industry
obs     <- d[d$scenario == "Observed", ]
groups  <- groups[obs$group[order(obs$first_total - obs$last_total)]]
n_sp    <- table(master$group)[names(groups)]
d$row   <- match(d$group, rev(names(groups))) + ifelse(d$scenario == "Observed", 0.15, -0.15)
d$group <- factor(d$group, levels = names(groups))

p <- ggplot(d, aes(x = pct, y = row, colour = group, linetype = scenario, shape = scenario)) +
  geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_linerange(aes(xmin = q05, xmax = q95), linewidth = 0.6) +
  geom_point(size = 2.2, fill = "white", stroke = 0.8) +
  scale_colour_manual(values = grp_cols, guide = "none") +
  scale_linetype_manual(values = scen_lines, guide = "none") +
  scale_shape_manual(values = scen_shapes, guide = "none") +
  scale_y_continuous(breaks = rev(seq_along(groups)), labels = paste0(groups, " (", n_sp, " species)"),
                     expand = expansion(add = 0.4)) +
  labs(x = "Change in abundance (%)", y = NULL) +
  theme_classic(base_size = 11) +
  theme(axis.text    = element_text(colour = "grey25"),
        axis.text.y  = element_text(size = 11),   # the group names, which 19E's rows rely on
        axis.line    = element_line(colour = "grey60", linewidth = 0.3),
        axis.ticks   = element_line(colour = "grey60", linewidth = 0.3),
        axis.line.y  = element_blank(),
        axis.ticks.y = element_blank())

out <- file.path(ia_dir, "output_figures", "figure1", paste0("fig1C_", stub, ".png"))
ggsave(out, p, width = 9.85, height = 4, dpi = 300)   # with 19E's 2.15 in, as wide as 19B + 19C above it
message("wrote ", out)
