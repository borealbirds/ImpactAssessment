# ---
# title: Net change in abundance by habitat group, 1995-2020
# author: Mannfred Boehm
# created: October 26, 2026
# ---

# our version of Rosenberg et al. (2019) Fig. 1B: the net change since the first year, observed
# and without industry, for each of 14D's habitat groups on one plot (colour = group, line type =
# scenario as in 19B). no legend: 19B names the line types, 19G the colours. panel C of the
# figure1.pptx variant; its y title serves 19F (panel D) beside it. a labelled copy names each
# group at the tip of its observed line.
# bands: each sample's group total minus the SAME sample's first-year total, 5th/95th percentiles.
# runs locally after 14D, with 19A's population_samples.rds pulled
#   -> output_figures/figure1/fig1BC_net_change_by_habitat{,_labelled}_1995_2020.png (+ _pct)
#   -> sector_effects/net_change_by_habitat_1995_2020.csv

library(ggplot2)

ia_dir <- getwd()
se_dir <- file.path(ia_dir, "data", "derived_data", "sector_effects")
ps     <- readRDS(file.path(se_dir, "population_samples.rds"))
master <- read.csv(file.path(ia_dir, "data", "raw_data", "anna_data", "master_species_list.csv"))

scen_lines <- c("Observed" = "solid", "Without industry" = "dashed")   # as 19B

# 14D's groups of the reported species, forest types first
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
stub  <- paste0("net_change_by_habitat_", years[y0], "_", max(years))

# m: a group's total [year x bootstrap] or [year x sample]. pct: the % copy (slide 2), change as %
# of the first-year total, as 19D: mean change over mean first-year total; band from each
# sample's change over its own first-year total
band <- function(m, scenario) {
  d <- sweep(m, 2, m[y0, ])
  r <- 100 * sweep(d, 2, m[y0, ], "/")
  data.frame(year     = years,
             scenario = scenario,
             total    = round(rowMeans(m)),
             change   = round(rowMeans(d)),
             q05      = round(apply(d, 1, quantile, 0.05)),
             q95      = round(apply(d, 1, quantile, 0.95)),
             pct      = round(100 * rowMeans(d) / mean(m[y0, ]), 3),
             pct_q05  = round(apply(r, 1, quantile, 0.05), 3),
             pct_q95  = round(apply(r, 1, quantile, 0.95), 3))
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

d$group <- factor(d$group, levels = names(groups))

# ten bands overlap near 0, hence a paler fill than 19B's
panel <- function(y, lo, hi, ylab, step) {
  d$y <- y; d$lo <- lo; d$hi <- hi
  ggplot(d, aes(year, y, colour = group, fill = group, linetype = scenario,
                group = interaction(group, scenario))) +
    geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.1, colour = NA) +
    geom_line(linewidth = 0.6) +
    geom_point(size = 1.5) +
    scale_colour_manual(values = grp_cols, guide = "none") +
    scale_fill_manual(values = grp_cols, guide = "none") +
    scale_linetype_manual(values = scen_lines, guide = "none") +
    scale_x_continuous(breaks = years) +
    scale_y_continuous(breaks = scales::breaks_width(step)) +
    labs(x = "Year", y = ylab) +
    theme_classic(base_size = 11) +
    theme(axis.text  = element_text(colour = "grey25"),
          axis.line  = element_line(colour = "grey60", linewidth = 0.3),
          axis.ticks = element_line(colour = "grey60", linewidth = 0.3))
}

out <- file.path(ia_dir, "output_figures", "figure1", paste0("fig1BC_", stub, ".png"))
p_pct <- panel(d$pct, d$pct_q05, d$pct_q95,
               paste0("Change since ", years[y0], " (% of ", years[y0], " population)"), 5)
ggsave(sub("\\.png$", "_pct.png", out), p_pct, width = 6, height = 4, dpi = 300)
message("wrote ", sub("\\.png$", "_pct.png", out))

p <- panel(d$change / 1e6, d$q05 / 1e6, d$q95 / 1e6,
           paste0("Change since ", years[y0], " (millions of birds)"), 50)
ggsave(out, p, width = 6, height = 4, dpi = 300)   # as 19B
message("wrote ", out)

# the variant's copy, labelled as 19F
tips  <- d[d$year == max(years) & d$scenario == "Observed", ]
tips$y <- tips$change / 1e6
p_lab <- p +
  ggrepel::geom_text_repel(data = tips, aes(label = groups[as.character(group)]), size = 3.2, hjust = 0,
                           direction = "y", nudge_x = 0.8, xlim = c(max(years) + 0.6, Inf),
                           segment.size = 0.2, segment.colour = "grey60", min.segment.length = 0,
                           box.padding = 0.12, show.legend = FALSE) +
  expand_limits(x = max(years) + 5.5)   # room for the names
out <- sub("by_habitat_", "by_habitat_labelled_", out)
ggsave(out, p_lab, width = 6, height = 4, dpi = 300)
message("wrote ", out)
