# ---
# title: Net change in abundance of the species that changed most, 1995-2020
# author: Mannfred Boehm
# created: October 27, 2026
# ---
# A variant of our Rosenberg et al. (2019, Science 366:120) Fig. 1, for its panel D: 19C's net
# change since the first year, observed and without industry, for each SPECIES in the top decile
# of change (14D's reported species ranked by the absolute change in their mean observed
# population, first year to last, gains and losses alike). Colour = habitat group as 19C, line
# type = scenario as 19B (whose legend serves), each species' code at the tip of its observed
# line. No legend: 19G (the variant's panel B) names the colours. It sits beside 19C (panel C),
# whose y title serves both, on its own y range (`ylim`, zoomed in on these lines): 19F's own
# title is drawn transparent, so that its panel starts where 19G's, above it, does. Its x range
# is 19G's too.
#
# Values as 19C's, per species: each sample's total minus the SAME sample's first-year total,
# then the 5th/95th percentiles. The CSV keeps them; the plot draws no bands (22 would overlap).
#
# Run locally after 14D, with 19A's population_samples.rds pulled:
#   Rscript Rscripts/19F_net_change_top_species.R
# Writes output_figures/figure1/fig1D_net_change_top_species_1995_2020.png and its values to
# sector_effects/net_change_top_species_1995_2020.csv.
# ---

suppressPackageStartupMessages(library(ggplot2))

ia_dir <- getwd()   # the local RProject
se_dir <- file.path(ia_dir, "data", "derived_data", "sector_effects")
ps     <- readRDS(file.path(se_dir, "population_samples.rds"))
master <- read.csv(file.path(ia_dir, "data", "raw_data", "anna_data", "master_species_list.csv"))

top_q <- 0.9          # species with |change| at or above this quantile
ylim  <- c(-120, 70)  # y axis (millions of birds); widens if a line leaves it

scen_lines <- c("Observed" = "solid", "Without industry" = "dashed")   # as 19B

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
stub  <- paste0("net_change_top_species_", years[y0], "_", years[y1])

# the species whose mean observed population changed most, either way
chg  <- rowMeans(ps$obs[master$species, y1, ]) - rowMeans(ps$obs[master$species, y0, ])
top  <- master[abs(chg) >= quantile(abs(chg), top_q), c("species", "group")]
top$change_M <- round(chg[top$species] / 1e6, 2)
top  <- top[order(-top$change_M), ]
print(top, row.names = FALSE)

# m: a species' population [year x bootstrap] or [year x sample], as 19C's band() (pct included)
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
for (i in seq_len(nrow(top))) {
  s <- top$species[i]
  r <- rbind(band(ps$obs[s, , ], "Observed"), band(ps$cf[s, , ], "Without industry"))
  rows[[s]] <- cbind(data.frame(species = s, group = top$group[i]), r)
}
d <- do.call(rbind, rows)
write.csv(d, file.path(se_dir, paste0(stub, ".csv")), row.names = FALSE)

panel <- function(y, ylab, step, ylim = 0) {
  d$y  <- y
  tips <- d[d$year == years[y1] & d$scenario == "Observed", ]
  ggplot(d, aes(year, y, colour = group, linetype = scenario,
                group = interaction(species, scenario))) +
    geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
    geom_line(linewidth = 0.5) +
    ggrepel::geom_text_repel(data = tips, aes(label = species), size = 3.2, hjust = 0,
                             direction = "y", nudge_x = 0.8, xlim = c(years[y1] + 0.6, Inf),
                             segment.size = 0.2, segment.colour = "grey60", min.segment.length = 0,
                             box.padding = 0.12, show.legend = FALSE) +
    scale_colour_manual(values = grp_cols, guide = "none") +
    scale_linetype_manual(values = scen_lines, guide = "none") +
    scale_x_continuous(breaks = years) +
    expand_limits(x = years[y1] + 5.5, y = ylim) +   # 19G's x range
    scale_y_continuous(breaks = scales::breaks_width(step)) +
    labs(x = "Year", y = ylab) +   # as 19C
    theme_classic(base_size = 11) +
    theme(axis.title.y = element_text(colour = "transparent"),
          axis.text  = element_text(colour = "grey25"),
          axis.line  = element_line(colour = "grey60", linewidth = 0.3),
          axis.ticks = element_line(colour = "grey60", linewidth = 0.3))
}

out <- file.path(ia_dir, "output_figures", "figure1", paste0("fig1D_", stub, ".png"))
p   <- panel(d$change / 1e6, paste0("Change since ", years[y0], " (millions of birds)"), 25, ylim)   # as 19G
ggsave(out, p, width = 6, height = 4, dpi = 300)   # as 19B and 19C
message("wrote ", out)

# the % copy (slide 2), the same species: each line's change as % of its own first-year total
out <- sub("\\.png$", "_pct.png", out)
p   <- panel(d$pct, paste0("Change since ", years[y0], " (% of ", years[y0], " population)"), 20)
ggsave(out, p, width = 6, height = 4, dpi = 300)
message("wrote ", out)
