# ---
# title: Share of species declining by habitat group, 1995-2020
# author: Mannfred Boehm
# created: October 26, 2026
# ---

# our version of Rosenberg et al. (2019) Fig. 1D: the share of each habitat group's species that
# declined from the first year to the last, observed and without industry, as pies on 19D's rows
# (same order and height, no row labels: 19D's name them); observed left of no-industry, rims
# solid and dashed as 19C's line types. a species declines when its mean change is negative (as
# Rosenberg counts negative trends); the solid wedge is the species whose whole 5th-95th
# percentile interval is below 0, the pale one those declining on the mean only. a species'
# change is taken per sample: 32 observed samples, 3200 without industry.
# runs locally after 14D, with 19A's population_samples.rds pulled
#   -> output_figures/figure1/fig1D_species_declining_by_habitat_1995_2020.png
#   -> sector_effects/species_declining_by_habitat_1995_2020.csv

library(ggplot2)

ia_dir <- getwd()
se_dir <- file.path(ia_dir, "data", "derived_data", "sector_effects")
ps     <- readRDS(file.path(se_dir, "population_samples.rds"))
master <- read.csv(file.path(ia_dir, "data", "raw_data", "anna_data", "master_species_list.csv"))

scen_lines <- c("Observed" = "solid", "Without industry" = "dashed")   # as 19C

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

# groups top to bottom as 19D (19C's observed lines' end, largest gain first)
groups <- groups[names(sort(tapply(rowMeans(chg_obs), master$group, sum), decreasing = TRUE))]

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

# pies of radius r at (column, row), first group on top as in 19D; wedges clockwise from
# 12 o'clock: whole interval below 0, then mean only, then the rest left blank
r         <- 0.38
d$x       <- ifelse(d$scenario == "Observed", 1, 2.4)
d$y       <- match(d$group, rev(names(groups)))
certainty <- c("Whole interval below 0" = 1, "Mean change below 0" = 0.35)
arc       <- function(from, to) 2 * pi * seq(from, to, length.out = 2 + ceiling(100 * (to - from)))

wedges <- do.call(rbind, lapply(seq_len(nrow(d)), function(i) {
  cut <- c(0, d$n_interval_below_0[i], d$n_declining[i]) / d$n_species[i]
  do.call(rbind, lapply(which(diff(cut) > 0), function(k) {
    a <- arc(cut[k], cut[k + 1])
    data.frame(id = paste(i, k), group = d$group[i], tier = names(certainty)[k],
               x = d$x[i] + c(0, r * sin(a)), y = d$y[i] + c(0, r * cos(a)))
  }))
}))
rims <- do.call(rbind, lapply(seq_len(nrow(d)), function(i) {
  a <- arc(0, 1)
  data.frame(id = i, group = d$group[i], scenario = d$scenario[i],
             x = d$x[i] + r * sin(a), y = d$y[i] + r * cos(a))
}))

# 19D's theme, y range and bottom axis (its text and ticks kept but invisible, for their height),
# so its rows and this panel's fall at the same heights; no legend: the caption explains the shading
p <- ggplot() +
  geom_polygon(data = wedges, aes(x, y, group = id, fill = group, alpha = tier)) +
  geom_path(data = rims, aes(x, y, group = id, colour = group, linetype = scenario), linewidth = 0.5) +
  geom_text(data = d, aes(x + r + 0.08, y, label = paste0(n_declining, "/", n_species)),
            hjust = 0, size = 3, colour = "grey25") +
  scale_fill_manual(values = grp_cols, guide = "none") +
  scale_colour_manual(values = grp_cols, guide = "none") +
  scale_linetype_manual(values = scen_lines, guide = "none") +
  scale_alpha_manual(values = certainty, guide = "none") +
  scale_x_continuous(breaks = 1, labels = "0") +
  scale_y_continuous(breaks = NULL) +
  coord_fixed(xlim = c(1 - r - 0.1, 2.4 + r + 0.6), ylim = c(0.45, length(groups) + 0.55),
              expand = FALSE) +
  labs(x = "Species declining", y = NULL) +
  theme_classic(base_size = 11) +
  theme(axis.text    = element_text(colour = "grey25"),
        axis.text.x  = element_text(colour = NA),
        axis.ticks.x = element_line(colour = NA),
        axis.line    = element_blank(),
        axis.ticks.y = element_blank())

out <- file.path(ia_dir, "output_figures", "figure1", paste0("fig1D_", stub, ".png"))
ggsave(out, p, width = 2.15, height = 4, dpi = 300, bg = "white")
message("wrote ", out)
