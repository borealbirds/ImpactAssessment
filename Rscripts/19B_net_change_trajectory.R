# ---
# title: Population of the reported species, 1995-2020
# author: Mannfred Boehm
# created: October 26, 2026
# ---

# our version of Rosenberg et al. (2019, Science 366:120) Fig. 1A: the summed population of the
# reported species (species_113.txt), observed and without industry, each year. not the net
# change since the first year, which would start both lines at 0 and hide the industry effect
# (the gap between them, 19G's "All species"). the y axis does not start at 0 (the gap is < 1%
# of the total). its line-type legend serves all four panels (19C-19E draw none).
# bands: each year's 5th/95th percentiles of the total over 19A's samples (observed: 32
# bootstraps; without industry: 3200, BART draws treated as independent across species).
# every year's counterfactual removes the 2020 footprint (07 ran 1995-2015 with
# HF_MASK_YEAR=2020): the gap is that footprint's effect under each year's conditions.
# runs locally after pulling 19A's population_samples.rds
#   -> output_figures/figure1/fig1A_net_change_1995_2020{,_pct}.png, sector_effects/net_change_1995_2020.csv

library(ggplot2)

ia_dir <- getwd()
se_dir <- file.path(ia_dir, "data", "derived_data", "sector_effects")
ps     <- readRDS(file.path(se_dir, "population_samples.rds"))

# two reds a lightness step apart (L* 54 and 31), told apart by line type as well
scen_cols  <- c("Observed" = "#D55E00", "Without industry" = "#8F1D0A")
scen_lines <- c("Observed" = "solid",   "Without industry" = "dashed")

species <- trimws(readLines(file.path(ia_dir, "Rscripts", "species_113.txt")))
species <- species[nzchar(species) & !startsWith(species, "#")]
if (length(miss <- setdiff(species, ps$species)))
  stop("population_samples.rds lacks ", paste(miss, collapse = ", "), " - rerun 19A")

years <- ps$years
y0    <- which.min(years)
stub  <- paste0("net_change_", years[y0], "_", max(years))

# summed over species: [year x bootstrap] and [year x sample]
obs  <- colSums(ps$obs[species, , , drop = FALSE])
cf   <- colSums(ps$cf[species, , , drop = FALSE])
boot <- rep(seq_len(ps$n_boot), each = ps$n_scen)   # sample s's bootstrap, as 19A
base <- obs[y0, ]                                   # each bootstrap's first-year observed total

# pct: the % copy (slide 2), both scenarios over the first year's OBSERVED total so the gap stays.
# As 14C, the line is the mean total over the mean first-year total and the band each sample over
# its own bootstrap's (b: each column's bootstrap)
band <- function(m, scenario, b) {
  r <- 100 * sweep(m, 2, base[b], "/")
  data.frame(
    year     = years,
    scenario = scenario,
    total    = round(rowMeans(m)),
    q05      = round(apply(m, 1, quantile, 0.05)),
    q95      = round(apply(m, 1, quantile, 0.95)),
    change   = round(rowMeans(m) - mean(m[y0, ])),   # since the first year
    pct      = round(100 * rowMeans(m) / mean(base), 3),
    pct_q05  = round(apply(r, 1, quantile, 0.05), 3),
    pct_q95  = round(apply(r, 1, quantile, 0.95), 3))
}
d <- rbind(band(obs, "Observed", seq_len(ps$n_boot)), band(cf, "Without industry", boot))
d <- d[order(d$year), ]
write.csv(d, file.path(se_dir, paste0(stub, ".csv")), row.names = FALSE)
print(d, row.names = FALSE)

panel <- function(y, lo, hi, ylab, step, acc) {
  d$y <- y; d$lo <- lo; d$hi <- hi
  ggplot(d, aes(year, y, colour = scenario, fill = scenario, linetype = scenario)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.18, colour = NA, show.legend = FALSE) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 2, show.legend = FALSE) +
    scale_colour_manual(values = scen_cols) +
    scale_fill_manual(values = scen_cols) +
    scale_linetype_manual(values = scen_lines) +
    guides(colour = guide_legend(override.aes = list(linewidth = 1.1))) +   # it serves all four panels
    scale_x_continuous(breaks = years) +
    scale_y_continuous(breaks = scales::breaks_width(step), labels = scales::label_number(accuracy = acc)) +
    expand_limits(y = max(d$hi) + 0.3 * diff(range(d$lo, d$hi))) +   # room for the legend
    labs(x = "Year", y = ylab, colour = NULL, fill = NULL, linetype = NULL) +
    theme_classic(base_size = 11) +
    theme(legend.position        = "inside",
          legend.position.inside = c(0.02, 0.98),
          legend.justification   = c(0, 1),
          legend.key.width       = unit(3.5, "lines"),
          legend.text            = element_text(size = 12),
          axis.text              = element_text(colour = "grey25"),
          axis.line              = element_line(colour = "grey60", linewidth = 0.3),
          axis.ticks             = element_line(colour = "grey60", linewidth = 0.3))
}

out <- file.path(ia_dir, "output_figures", "figure1", paste0("fig1A_", stub, ".png"))
p   <- panel(d$total / 1e9, d$q05 / 1e9, d$q95 / 1e9, "Population (billions of birds)", 0.05, 0.01)
ggsave(out, p, width = 6, height = 4, dpi = 300)
message("wrote ", out)

out <- sub("\\.png$", "_pct.png", out)
p   <- panel(d$pct, d$pct_q05, d$pct_q95, paste0("Population (% of observed ", years[y0], ")"), 1, 1)
ggsave(out, p, width = 6, height = 4, dpi = 300)
message("wrote ", out)
