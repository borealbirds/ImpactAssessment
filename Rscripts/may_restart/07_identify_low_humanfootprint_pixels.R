# ---
# title: Impact Assessment: generate low and high human footprint layers, and merge data-sparse subbasins 
# author: Mannfred Boehm
# created: September 16, 2025
# ---

library(terra)
library(tidyverse)

root <- "G:/Shared drives/BAM_NationalModels5"
ia_dir <- file.path(root, "data", "Extras", "sandbox_data", "impactassessment_sandbox")

# identify and store HF layers
# hf_layers <- intersect(c("CanHF_1km", "CanHF_5x5"), names(stack_y))
# hf_stack <- if (length(hf_layers)) stack_y[[hf_layers]] else NULL


# set low HF to <1 
# from Hirsh-Pearson: we found that 82% of Canada’s land areas had a 
# HF < 1 and therefore were considered intact
lowhf_mask <- 
  terra::lapp(stack_y[["CanHF_1km"]],
              \(v) ifelse(!is.finite(v), NA, ifelse(v < 1, 1, NA))) |> 
  as.factor()

#  per-subbasin lowhf counts 
# problem: many subbasins don't have any pixels with HF < 1
counts_df <- terra::extract(
  lowhf_mask,
  all_subbasins,
  fun   = function(x) sum(!is.na(x)),  # count non-NA cells
  ID    = TRUE)
mean(counts_df$CanHF_1km); range(counts_df$CanHF_1km)

# all_subbasins$sub_count <- counts_df$CanHF_1km[match(seq_len(nrow(all_subbasins)), counts_df$ID)]

# tag each point with its basin
ij <- terra::intersect(as.points(lowhf_mask, na.rm = TRUE), all_subbasins)

# paste counts back onto the SpatVector (by HYBAS_ID)
ij$sub_count <- all_subbasins$sub_count[ match(ij$HYBAS_ID, all_subbasins$HYBAS_ID) ]

plot_df <- cbind(terra::crds(ij, df = TRUE), sub_count = ij$sub_count)

# plot (exported at 1000 x 751)
ggplot() +

# points colored by subbasin density
geom_point(data = slice_sample(plot_df, prop=0.01), aes(x = x, y = y, colour = sub_count), size = 0.05) +
scale_color_gradient(low="#56B4E9", high="#CC79A7") +

# BCR and basin outlines
geom_spatvector(data = bam_boundary, fill = NA, colour = "grey") +
geom_spatvector(data = all_subbasins, fill = NA, color = "grey25", linewidth = 0.3) +

coord_sf(crs = crs(all_subbasins)) +
theme_minimal()

