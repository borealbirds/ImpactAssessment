# ---
# title: Map the 16 subbasins sampled in 16_lowHF_isnot_a_proxy.R
# Reproduces subbasin_index from the same seed without re-running the RF pipeline.
# ---

library(terra)
library(tidyterra)
library(tidyverse)

root   <- "G:/Shared drives/BAM_NationalModels5"
ia_dir <- file.path(root, "data", "Extras", "sandbox_data", "impactassessment_sandbox")

# --- Reproduce subbasin_index ------------------------------------------------
# Use one layer of the covariate stack as CRS/resolution reference (lazy-loaded)
stack_ref <- terra::rast(file.path(ia_dir, "covariates_mosaiced_2020.tif"))[[1]]

all_subbasins_subset <- terra::vect(file.path(ia_dir, "hydrobasins_masked_merged_subset.gpkg"))
all_subbasins_subset <- terra::project(all_subbasins_subset, stack_ref)

lowhf_mask_raw  <- terra::rast(file.path(ia_dir, "hirshpearson", "CanHF_1km_lessthan1.tif"))
highhf_mask_raw <- terra::rast(file.path(ia_dir, "hirshpearson", "CanHF_1km_morethan1.tif"))
lowhf_mask_raw  <- terra::project(lowhf_mask_raw,  stack_ref, method = "near")
highhf_mask_raw <- terra::project(highhf_mask_raw, stack_ref, method = "near")

min_n    <- 100
n_target <- 64
set.seed(42)
candidate_order <- sample(seq_len(nrow(all_subbasins_subset)))
subbasin_index  <- integer(0)

for (idx in candidate_order) {
  sub_poly <- all_subbasins_subset[idx]
  lo_sub <- tryCatch(terra::crop(lowhf_mask_raw,  sub_poly) |> terra::mask(sub_poly), error = function(e) NULL)
  hi_sub <- tryCatch(terra::crop(highhf_mask_raw, sub_poly) |> terra::mask(sub_poly), error = function(e) NULL)
  n_lo <- if (!is.null(lo_sub)) terra::global(lo_sub, "notNA")[[1]] else 0
  n_hi <- if (!is.null(hi_sub)) terra::global(hi_sub, "notNA")[[1]] else 0
  if (n_lo >= min_n && n_hi >= min_n) subbasin_index <- c(subbasin_index, idx)
  if (length(subbasin_index) == n_target) break
}
cat(sprintf("Reproduced %d subbasins\n", length(subbasin_index)))

subbasin_s      <- all_subbasins_subset[subbasin_index]
subbasin_labels <- seq_along(subbasin_index)

sub_centroids       <- terra::centroids(subbasin_s)
sub_centroids$label <- subbasin_labels

# --- Map ---------------------------------------------------------------------
bam_boundary <- terra::vect(file.path(root, "Regions", "BAM_BCR_NationalModel_UnBuffered.shp"))
bam_boundary <- bam_boundary[bam_boundary$subUnit != 23, ]

p_map <- ggplot() +
  geom_spatvector(data = bam_boundary, fill = NA, colour = "grey") +
  geom_spatvector(data = subbasin_s, fill = NA, color = "#D55E00", linewidth = 0.5) +
  geom_spatvector_text(
    data = sub_centroids,
    aes(label = label),
    nudge_x = 20000,
    nudge_y = 100000,
    size = 3) +
  coord_sf(crs = crs(subbasin_s)) +
  theme_minimal()

out_fig_dir <- file.path(getwd(), "output_figures", "lowHF_isnotaproxy")
dir.create(out_fig_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(out_fig_dir, "sampled_subbasins_map.png"),
       p_map, width = 1500, height = 1000, units = "px", dpi = 150)
cat("Map saved.\n")

# --- Update rf_summary.csv with label column ---------------------------------
summary_path <- file.path(getwd(), "data/derived_data/rf_summary.csv")
summary_df   <- read.csv(summary_path)

label_key  <- tibble(subbasin = all_subbasins_subset$first_HYBAS_ID[subbasin_index],
                     label    = subbasin_labels)
summary_df <- dplyr::left_join(summary_df, label_key, by = "subbasin") |>
              dplyr::relocate(label, .before = subbasin)

write.csv(summary_df, summary_path, row.names = FALSE)
cat("rf_summary.csv updated with label column.\n")
