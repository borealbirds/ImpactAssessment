# ---
# title: Population distribution plots, observed vs counterfactual
# author: Mannfred Boehm
# ---

# densities of the national population over the (bootstrap, scenario) samples, from 12H's
# arrays (density_tables/arrays/{species}_{year}_coalition_{id}_arrays.rds):
#   {species}_{year}_all_sectors.png     observed vs all sectors removed
#   {species}_{year}_single_sectors.png  observed vs each sector removed alone (8 panels)
# counterfactual population = obs_total - obs_on_coal + bf_on_coal, per sample
#   -> data/derived_data/figures/population_distributions/

library(tidyverse)

ia_dir  <- getwd()
arr_dir <- file.path(ia_dir, "data", "derived_data", "density_tables", "arrays")
fig_dir <- file.path(ia_dir, "data", "derived_data", "figures", "population_distributions")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

source(file.path(ia_dir, "Rscripts", "12E_shapley_utils.R"))

sectors    <- canonical_sectors()
full_id    <- sectors_to_coalition_id(sectors, sectors)
single_ids <- setNames(vapply(sectors, function(s) sectors_to_coalition_id(s, sectors), numeric(1L)), sectors)

arr_files <- list.files(arr_dir, pattern = "_arrays\\.rds$", full.names = TRUE)
if (length(arr_files) == 0) stop("no array files in ", arr_dir, " (12H writes them)")
parsed    <- do.call(rbind, regmatches(basename(arr_files), regexec("^(.+)_([0-9]{4})_coalition_([0-9]+)_arrays\\.rds$", basename(arr_files))))
arr_index <- data.frame(path = arr_files, species = parsed[, 2], year = parsed[, 3], coalition_id = as.integer(parsed[, 4]))

# national population per sample: observed (the same in every scenario of a bootstrap) and counterfactual
obs_draws <- function(arr) as.vector(arr$obs_total_mat)
cf_draws  <- function(arr) as.vector(arr$obs_total_mat - arr$obs_on_coal_mat + arr$bf_on_coal_mat)

pop_plot <- function(d, cf_label, title, linewidth) {
  ggplot(d, aes(x = population, fill = scenario, colour = scenario)) +
    geom_density(alpha = 0.35, linewidth = linewidth) +
    scale_fill_manual(values   = setNames(c("#4477AA", "#EE6677"), c("Observed", cf_label))) +
    scale_colour_manual(values = setNames(c("#4477AA", "#EE6677"), c("Observed", cf_label))) +
    scale_x_continuous(labels = scales::comma) +
    labs(title = title, x = "National population (birds)", y = "Density", fill = NULL, colour = NULL) +
    theme_classic(base_size = 11) +
    theme(strip.background = element_blank(), strip.text = element_text(face = "bold"), legend.position = "bottom")
}



# ---------------------------------------------------
# one pair of figures per species x year

for (sy in split(arr_index, paste(arr_index$species, arr_index$year))) {

  sp <- sy$species[1]
  yr <- sy$year[1]

  # observed vs all sectors removed
  if (full_id %in% sy$coalition_id) {
    arr <- readRDS(sy$path[sy$coalition_id == full_id])
    d   <- data.frame(population = c(obs_draws(arr), cf_draws(arr)),
                      scenario   = rep(c("Observed", "All sectors removed"), each = length(obs_draws(arr))))
    ggsave(file.path(fig_dir, paste0(sp, "_", yr, "_all_sectors.png")),
           pop_plot(d, "All sectors removed", paste(sp, yr, "— national population"), linewidth = 0.7),
           width = 7, height = 4.5, dpi = 180)
  }

  # observed vs each sector removed alone
  singles <- sy[sy$coalition_id %in% single_ids, ]
  if (nrow(singles) == 0) next
  obs_vec <- obs_draws(readRDS(singles$path[1]))
  d <- bind_rows(lapply(seq_len(nrow(singles)), function(i) {
    data.frame(population = c(obs_vec, cf_draws(readRDS(singles$path[i]))),
               scenario   = rep(c("Observed", "Sector removed"), each = length(obs_vec)),
               sector     = names(single_ids)[single_ids == singles$coalition_id[i]])
  }))
  d$sector <- factor(d$sector, levels = sort(unique(d$sector)))
  ggsave(file.path(fig_dir, paste0(sp, "_", yr, "_single_sectors.png")),
         pop_plot(d, "Sector removed", paste(sp, yr, "— national population by sector"), linewidth = 0.6) +
           facet_wrap(~ sector, scales = "free", ncol = 2),
         width = 10, height = 10, dpi = 180)
}
