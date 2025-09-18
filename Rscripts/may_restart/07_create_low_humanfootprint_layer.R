# ---
# title: Impact Assessment: generate low and high human footprint layers, and merge data-sparse subbasins 
# author: Mannfred Boehm
# created: September 16, 2025
# ---

library(terra)
library(tidyterra)
library(tidyverse)

root <- "G:/Shared drives/BAM_NationalModels5"
ia_dir <- file.path(root, "data", "Extras", "sandbox_data", "impactassessment_sandbox")


# ---------------------------------------------------------------
# OBJECTIVE 1: create and save a low HF raster for model training

# import any arbitrary covariate stack to access CanHF_1km 
# (HF layer is the same across years)
CanHF_1km <- rast(file.path(ia_dir, "covariates_mosaiced_2020.tif"))[["CanHF_1km"]]

# set low HF to <1 
# from Hirsh-Pearson: "we found that 82% of Canada’s land areas had a 
# HF < 1 and therefore were considered intact"
lowhf_mask <- 
  terra::lapp(CanHF_1km, \(x) ifelse(!is.finite(x), NA, ifelse(x < 1, 1, NA))) |> 
  as.factor()

terra::writeRaster(lowhf_mask, file.path(ia_dir, "CanHF_1km_lessthan1.tif"))


# ---------------------------------------------------------------
# OBJECTIVE 2: identify subbasins with zero or few low HF pixels,


# import subbains multi-polygon
all_subbasins <- terra::vect(file.path(ia_dir, "hydrobasins_subset_2020.gpkg"))

# per-subbasin lowhf counts 
# note: many subbasins don't have any pixels with HF < 1
counts_df <- terra::extract(
  lowhf_mask,
  all_subbasins,
  fun   = function(x) sum(!is.na(x)),  # count non-NA cells
  ID    = TRUE)

sum(counts_df$CanHF_1km) #1,917,452
quantile(counts_df$CanHF_1km)

# ---------------------------------------------------------------
# OBJECTIVE 3: visualize subbasins with their low HF pixel density

# assigns a low HF pixel count to each subbasin
all_subbasins$sub_count <- counts_df$CanHF_1km[match(seq_len(nrow(all_subbasins)), counts_df$ID)]

# join each of the 1.9 million points to its subbasin ID
ij <- terra::intersect(as.points(lowhf_mask, na.rm = TRUE), all_subbasins)

# paste subbasin-level lowHF counts back onto each point
# every point is aware of if it's subbasin pixel density
ij$sub_count <- all_subbasins$sub_count[ match(ij$HYBAS_ID, all_subbasins$HYBAS_ID) ]

# generate dataframe: columns are lat, long (for each point) 
# and subbasin pixel density for that point
plot_df <- cbind(terra::crds(ij, df = TRUE), sub_count = ij$sub_count)

# import BCR boundaries for plotting
bam_boundary <- terra::vect(file.path(root, "Regions", "BAM_BCR_NationalModel_UnBuffered.shp"))


# plot (exported at 1000 x 751)
ggplot() +

# points colored by subbasin density
geom_point(data = slice_sample(plot_df, prop=0.04), aes(x = x, y = y, colour = sub_count), size = 0.01) +
scale_color_gradient(low="#56B4E9", high="#CC79A7") +

# BCR and basin outlines
geom_spatvector(data = bam_boundary, fill = NA, colour = "grey") +
geom_spatvector(data = all_subbasins, fill = NA, color = "grey25", linewidth = 0.3) +

coord_sf(crs = crs(all_subbasins)) +
theme_minimal()


# ---------------------------------------------------------------
# OBJECTIVE 4: merge low pixel density subbasins with 
# neighbouring subbasins, but only if they're in the same BCR


# find which BCR(s) a subbasin centroid is nearest to (including inside of)
hits <- terra::nearest(centroids(all_subbasins), bam_boundary)

# assign BCR (always returns one)
all_subbasins$BCR <- bam_boundary$subUnit[hits$to_id]




# sample size threshold is 25 percentile
threshold <- quantile(counts_df$CanHF_1km)[2]


