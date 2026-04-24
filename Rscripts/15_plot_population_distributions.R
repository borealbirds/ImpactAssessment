# ---
# title: Population distribution plots — observed vs counterfactual
# author: Mannfred Boehm
# ---
# Reads the national bootstrap x scenario arrays saved by 12B for target
# coalitions and plots empirical density distributions of population size.
#
# Two figure types per species:
#   1. Observed vs all-sectors-removed (one panel)
#   2. Observed vs each single-sector-removed (8 panels, one per sector)
#
# Inputs (density_tables/arrays/):
#   {species}_{year}_coalition_{id}_arrays.rds
#   Each file contains a named list:
#     $obs_total_mat   — n_boot x n_scen matrix, observed national population
#     $obs_on_coal_mat — n_boot x n_scen matrix, observed birds on coalition pixels
#     $bf_on_coal_mat  — n_boot x n_scen matrix, backfilled birds on coalition pixels
#   Counterfactual population = obs_total - obs_on_coal + bf_on_coal (per cell)
#
# Outputs (data/derived_data/figures/population_distributions/):
#   {species}_{year}_all_sectors.png
#   {species}_{year}_single_sectors.png
# ---

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

# ---- Execution context -------------------------------------------------------

cc    <- FALSE
local <- TRUE

if (cc)            { ia_dir <- "/home/mannfred/scratch/impact_assessment" }
if (!cc && local)  { ia_dir <- getwd() }
if (!cc && !local) { ia_dir <- file.path("G:/Shared drives/BAM_NationalModels5", "data", "Extras",
                                          "sandbox_data", "impactassessment_sandbox") }

# ---- Source utilities --------------------------------------------------------

source(file.path(ia_dir, "Rscripts", "shapley_utils.R"))

# ---- Paths -------------------------------------------------------------------

arr_dir <- file.path(ia_dir, "data", "derived_data", "density_tables", "arrays")
fig_dir <- file.path(ia_dir, "data", "derived_data", "figures", "population_distributions")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Identify target coalition IDs ------------------------------------------

sectors    <- canonical_sectors()
full_id    <- sectors_to_coalition_id(sectors, sectors)
single_ids <- setNames(
  vapply(sectors, function(s) sectors_to_coalition_id(s, sectors), integer(1L)),
  sectors
)

# ---- Discover available array files -----------------------------------------

arr_files <- list.files(arr_dir, pattern = "_arrays\\.rds$", full.names = TRUE)
if (length(arr_files) == 0) stop("No array files found in ", arr_dir,
                                 "\nRe-run 12B on the cluster first.")

parsed <- regmatches(basename(arr_files),
  regexec("^(.+)_([0-9]{4})_coalition_([0-9]+)_arrays\\.rds$", basename(arr_files)))
parsed <- do.call(rbind, lapply(parsed, function(x) x[2:4]))
arr_index <- data.frame(
  path         = arr_files,
  species      = parsed[, 1],
  year         = parsed[, 2],
  coalition_id = as.integer(parsed[, 3]),
  stringsAsFactors = FALSE
)

species_years <- unique(arr_index[, c("species", "year")])
message("Found array files for ", nrow(species_years), " species x year combinations")

# ---- Helper: extract flat vector of national population draws ---------------

# For a given arrays RDS, compute the counterfactual population as a flat
# vector across all (boot x scenario) draws.
cf_draws <- function(arr) {
  as.vector(arr$obs_total_mat - arr$obs_on_coal_mat + arr$bf_on_coal_mat)
}

obs_draws <- function(arr) {
  # obs_total does not vary across BART scenarios — use column means to get
  # the 32-element bootstrap distribution, then replicate for comparability
  as.vector(arr$obs_total_mat)
}

# ---- Plot theme -------------------------------------------------------------

theme_pop <- theme_classic(base_size = 11) +
  theme(
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold"),
    legend.position  = "bottom"
  )

obs_colour <- "#4477AA"
cf_colour  <- "#EE6677"

# ---- Main loop --------------------------------------------------------------

for (sy in seq_len(nrow(species_years))) {

  sp <- species_years$species[sy]
  yr <- species_years$year[sy]
  message("\n=== ", sp, " ", yr, " ===")

  idx <- arr_index$species == sp & arr_index$year == yr
  avail <- arr_index[idx, ]

  # ---- Figure 1: observed vs all sectors removed ---------------------------

  full_row <- avail[avail$coalition_id == full_id, ]

  if (nrow(full_row) == 1) {

    arr_full <- readRDS(full_row$path)

    df_full <- data.frame(
      population = c(obs_draws(arr_full), cf_draws(arr_full)),
      scenario   = rep(c("Observed", "All sectors removed"),
                       each = length(obs_draws(arr_full)))
    )

    p_full <- ggplot(df_full, aes(x = population, fill = scenario, colour = scenario)) +
      geom_density(alpha = 0.35, linewidth = 0.7) +
      scale_fill_manual(values   = c("Observed" = obs_colour, "All sectors removed" = cf_colour)) +
      scale_colour_manual(values = c("Observed" = obs_colour, "All sectors removed" = cf_colour)) +
      scale_x_continuous(labels = scales::comma) +
      labs(
        title    = paste(sp, yr, "— national population"),
        subtitle = "Counterfactual: all industrial sectors removed",
        x        = "National population (birds)",
        y        = "Density",
        fill     = NULL, colour = NULL
      ) +
      theme_pop

    out_path <- file.path(fig_dir, paste0(sp, "_", yr, "_all_sectors.png"))
    ggsave(out_path, p_full, width = 7, height = 4.5, dpi = 180)
    message("  wrote ", basename(out_path))

  } else {
    message("  WARNING: no full-coalition array for ", sp, " ", yr, " — skipping figure 1")
  }

  # ---- Figure 2: observed vs each single sector removed -------------------

  single_rows <- avail[avail$coalition_id %in% single_ids, ]

  if (nrow(single_rows) == 0) {
    message("  WARNING: no single-sector arrays for ", sp, " ", yr, " — skipping figure 2")
    next
  }

  # read observed baseline from any available single-sector file
  arr_ref    <- readRDS(single_rows$path[1])
  obs_vec    <- obs_draws(arr_ref)

  panel_list <- vector("list", nrow(single_rows))

  for (i in seq_len(nrow(single_rows))) {
    cid     <- single_rows$coalition_id[i]
    sec     <- names(single_ids)[single_ids == cid]
    arr_sec <- readRDS(single_rows$path[i])

    panel_list[[i]] <- data.frame(
      population = c(obs_vec, cf_draws(arr_sec)),
      scenario   = rep(c("Observed", "Sector removed"), each = length(obs_vec)),
      sector     = sec,
      stringsAsFactors = FALSE
    )
  }

  df_single <- bind_rows(panel_list)
  df_single$sector <- factor(df_single$sector, levels = sort(unique(df_single$sector)))

  p_single <- ggplot(df_single, aes(x = population, fill = scenario, colour = scenario)) +
    geom_density(alpha = 0.35, linewidth = 0.6) +
    facet_wrap(~ sector, scales = "free", ncol = 2) +
    scale_fill_manual(values   = c("Observed" = obs_colour, "Sector removed" = cf_colour)) +
    scale_colour_manual(values = c("Observed" = obs_colour, "Sector removed" = cf_colour)) +
    scale_x_continuous(labels = scales::comma) +
    labs(
      title    = paste(sp, yr, "— national population by sector"),
      subtitle = "Each panel: counterfactual with that sector's footprint removed",
      x        = "National population (birds)",
      y        = "Density",
      fill     = NULL, colour = NULL
    ) +
    theme_pop

  out_path <- file.path(fig_dir, paste0(sp, "_", yr, "_single_sectors.png"))
  ggsave(out_path, p_single, width = 10, height = 10, dpi = 180)
  message("  wrote ", basename(out_path))

}

message("\nDone. Figures written to ", fig_dir)
