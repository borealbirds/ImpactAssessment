# ---
# title: Industry's effect on each habitat group, 1995-2020
# author: Mannfred Boehm
# created: October 27, 2026
# ---

# panel B of the figure1.pptx variant: each year's observed minus no-industry population (14C's
# sign: negative = birds lost to industry) for each of 14D's habitat groups and for all reported
# species together. 19B shows the total as the gap between its lines, 19C and 19F show change
# since the first year (which subtracts it); here it is the level. colours as 19C, lines dashed
# as 19B's no-industry line, each line named at its tip, on 19F's x range.
# bands: each sample's group total observed (its bootstrap's) minus the SAME sample's
# no-industry total, 5th/95th percentiles, as 14C. every year removes the 2020 footprint (07
# ran 1995-2015 with HF_MASK_YEAR=2020).
# runs locally after 14D, with 19A's population_samples.rds pulled
#   -> output_figures/figure1/fig1B_industry_effect_by_habitat_1995_2020{,_pct}.png
#   -> sector_effects/industry_effect_by_habitat_1995_2020.csv

library(ggplot2)

ia_dir <- getwd()
se_dir <- file.path(ia_dir, "data", "derived_data", "sector_effects")
ps     <- readRDS(file.path(se_dir, "population_samples.rds"))
master <- read.csv(file.path(ia_dir, "data", "raw_data", "anna_data", "master_species_list.csv"))

# as 19C, plus all groups together
groups   <- c("conifer" = "Conifer", "mixedwood" = "Mixedwood", "deciduous" = "Deciduous",
              "forest no type" = "Forest, no type", "non-forest" = "Non-forest", "all" = "All species")
grp_cols <- c("conifer" = "#009E73", "mixedwood" = "#CC79A7", "deciduous" = "#E69F00",
              "forest no type" = "#56B4E9", "non-forest" = "#000000", "all" = "grey55")
master <- master[master$include, ]
if (length(odd <- setdiff(master$group, names(groups))))
  stop("14D group(s) not in this script: ", paste(odd, collapse = ", "))
if (length(miss <- setdiff(master$species, ps$species)))
  stop("population_samples.rds lacks ", paste(miss, collapse = ", "), " - rerun 19A")

years <- ps$years
stub  <- paste0("industry_effect_by_habitat_", min(years), "_", max(years))
boot  <- rep(seq_len(ps$n_boot), each = ps$n_scen)   # sample s's bootstrap, as 19A

rows <- list()
for (g in names(groups)) {
  s   <- if (g == "all") master$species else master$species[master$group == g]
  cfs <- colSums(ps$cf[s, , , drop = FALSE])                  # [year x sample]
  eff <- colSums(ps$obs[s, , , drop = FALSE])[, boot] - cfs
  pct <- 100 * eff / cfs   # the % copy (slide 2): of the no-industry population, as 14C
  rows[[g]] <- data.frame(group     = g,
                          n_species = length(s),
                          year      = years,
                          effect    = round(rowMeans(eff)),
                          q05       = round(apply(eff, 1, quantile, 0.05)),
                          q95       = round(apply(eff, 1, quantile, 0.95)),
                          pct       = round(100 * rowMeans(eff) / rowMeans(cfs), 3),
                          pct_q05   = round(apply(pct, 1, quantile, 0.05), 3),
                          pct_q95   = round(apply(pct, 1, quantile, 0.95), 3))
}
d <- do.call(rbind, rows)
write.csv(d, file.path(se_dir, paste0(stub, ".csv")), row.names = FALSE)
print(d[d$year == max(years), ], row.names = FALSE)

d$group <- factor(d$group, levels = names(groups))

panel <- function(y, lo, hi, ylab, step) {
  d$y <- y; d$lo <- lo; d$hi <- hi
  tips <- d[d$year == max(years), ]
  ggplot(d, aes(year, y, colour = group, fill = group)) +
    annotate("segment", x = -Inf, xend = max(years), y = 0, yend = 0,   # clear of the names near 0
             colour = "grey70", linewidth = 0.3) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.1, colour = NA) +   # as 19C
    geom_line(linewidth = 0.6, linetype = "dashed") +
    geom_point(size = 1.5) +
    ggrepel::geom_text_repel(data = tips, aes(label = groups[as.character(group)]), size = 3.2, hjust = 0,
                             direction = "y", nudge_x = 0.8, xlim = c(max(years) + 0.6, Inf),
                             segment.size = 0.2, segment.colour = "grey60", min.segment.length = 0,
                             box.padding = 0.12, show.legend = FALSE) +
    scale_colour_manual(values = grp_cols, guide = "none") +
    scale_fill_manual(values = grp_cols, guide = "none") +
    scale_x_continuous(breaks = years) +
    expand_limits(x = max(years) + 5.5) +   # room for the names; 19F's x range too
    scale_y_continuous(breaks = scales::breaks_width(step)) +
    labs(x = "Year", y = ylab) +
    theme_classic(base_size = 11) +
    theme(axis.text  = element_text(colour = "grey25"),
          axis.line  = element_line(colour = "grey60", linewidth = 0.3),
          axis.ticks = element_line(colour = "grey60", linewidth = 0.3))
}

out <- file.path(ia_dir, "output_figures", "figure1", paste0("fig1B_", stub, ".png"))
p   <- panel(d$effect / 1e6, d$q05 / 1e6, d$q95 / 1e6, "Effect of industry (millions of birds)", 25)
ggsave(out, p, width = 6, height = 4, dpi = 300)   # as 19B, 19C and 19F
message("wrote ", out)

out <- sub("\\.png$", "_pct.png", out)
p   <- panel(d$pct, d$pct_q05, d$pct_q95, "Effect of industry (% of no-industry population)", 5)
ggsave(out, p, width = 6, height = 4, dpi = 300)
message("wrote ", out)
