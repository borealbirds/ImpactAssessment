# ---
# title: Impact Assessment: subset hydrobasins to those with industry footprints
# author: Mannfred Boehm
# created: August 19, 2025
# ---

library(terra)
library(tidyterra)
library(tidyverse)

root <- "G:/Shared drives/BAM_NationalModels5"
ia_dir <- file.path(root, "data", "Extras", "sandbox_data", "impactassessment_sandbox")


# -----------------------------------------------------------
# PART I
# create polygons for mines and oil&gas

# use %d as a placeholder to insert year_y (2020)
year <- 2020
mines_filepath <- file.path(ia_dir, sprintf("mincan_mines_%d_masked.tif", year))

# import mining layer for year_y
mines_y <- terra::rast(mines_filepath)

# convert mines patches to a multi-polygon (one polygon per mine)
# this is necessary for terra::relate() to find the associated subbasin
mines_polygons_y <- terra::as.polygons(mines_y, dissolve = TRUE, na.rm = TRUE)
names(mines_polygons_y) <- "patch_id"

# import oil & gas layer (Hirsh-Pearson)
oilgas <- terra::rast(file.path(ia_dir, "hirshpearson_oil_gas_masked.tif"))
oilgas_polygons <- terra::as.polygons(oilgas, dissolve = TRUE, na.rm = TRUE)
names(oilgas_polygons) <- "patch_id"

# merge mines and oilgas footprints
combined_poly <- rbind(mines_polygons_y[mines_polygons_y$patch_id == 1],
                       oilgas_polygons[oilgas_polygons$patch_id == 1])
              
terra::writeVector(combined_poly, file.path(ia_dir, "combined_industry_footprint.gpkg"), overwrite=TRUE)



# -----------------------------------------------------------
# PART II 
# subset subbasins to those with industry footprints

# import industry footprint (can skip PART I)
combined_poly <- terra::vect(file.path(ia_dir, "combined_industry_footprint.gpkg"))

# import subbasins multi-polygon (merged so that every subbasin has sufficient low HF pixels)
all_subbasins_merged <- terra::vect(file.path(ia_dir, "hydrobasins_masked_merged.gpkg"))

# identify subbasins that intersect mines or oilgas
hits1 <- terra::relate(all_subbasins_merged, combined_poly, relation = "intersects")  
hits2 <- unique(which(hits1, arr.ind = TRUE)[, 1])                          
all_subbasins_subset <- all_subbasins_merged[hits2]                                       

terra::writeVector(all_subbasins_subset, file.path(ia_dir, "hydrobasins_masked_merged_subset.gpkg"), overwrite=TRUE)




# -----------------------------------------------------------
# PART III
# count low HF training pixels per subbasin

# restrict analysis to BCRs that we have covariate data for
bam_template <- terra::rast(file.path(root, "PredictionRasters", "Biomass", "SCANFI", "1km", "SCANFIBalsamFir_1km_2020.tif"))
bam_boundary <- terra::vect(file.path(root, "Regions", "BAM_BCR_NationalModel_UnBuffered.shp"))
bam_boundary <- bam_boundary[bam_boundary$subUnit != 23, ]

# import Hirsh-Pearson HF layer
CanHF_1km <- 
  terra::rast(file.path(root, "CovariateRasters", "Disturbance", "cum_threat2020.02.18.tif")) |> 
  terra::project(x = _, y = bam_template) |>
  terra::resample(x = _, y = bam_template) |> 
  terra::crop(x = _, y = bam_template) |> 
  terra::mask(x = _, mask = bam_template)

names(CanHF_1km) <- "CanHF_1km"

# set low HF to <1 
# from Hirsh-Pearson: "we found that 82% of Canada’s land areas had a 
# HF < 1 and therefore were considered intact"
lowhf_mask <- 
  terra::lapp(CanHF_1km, \(x) ifelse(!is.finite(x), NA, ifelse(x < 1, 1, NA))) |> 
  as.factor()

counts_df <- 
  terra::extract(lowhf_mask, all_subbasins_subset, fun=function(x) sum(!is.na(x)), ID=TRUE) |> 
  dplyr::as_tibble() |>
  dplyr::arrange(CanHF_1km)

quantile(counts_df$CanHF_1km) 
# 0%     25%     50%     75%    100% 
# 905.0  2518.5  4318.0  8654.0 26261.0 

# assign a low HF pixel count to each subbasin
all_subbasins_subset$sub_count <- counts_df$CanHF_1km[match(seq_len(nrow(all_subbasins_subset)), counts_df$ID)]

# join each of the 5 million points to its subbasin ID
ij <- terra::intersect(as.points(lowhf_mask, na.rm = TRUE), all_subbasins_subset)

# paste subbasin-level lowHF counts back onto each point
# every point is aware of if it's subbasin pixel density
# match by group_key (present post-merge), not HYBAS_ID
ij$sub_count <- all_subbasins_subset$sub_count[
  match(ij$group_key, all_subbasins_subset$group_key)
]
# generate dataframe: rows are low HF points, 
# columns are lat, long, subbasin pixel density for that point
plot_df <- cbind(terra::crds(ij, df = TRUE), sub_count = ij$sub_count)

# plot (exported at 1000 x 751)
ggplot() +
  
  # points colored by subbasin density
  geom_point(data = slice_sample(plot_df, prop=0.04), aes(x = x, y = y, colour = sub_count), size = 0.01) +
  scale_color_gradient(low="#56B4E9", high="#CC79A7") +
  
  # BCR and basin outlines
  geom_spatvector(data = bam_boundary, fill = NA, colour = "grey") +
  geom_spatvector(data = all_subbasins_subset, fill = NA, color = "grey25", linewidth = 0.3) +
  
  coord_sf(crs = crs(all_subbasins_subset)) +
  theme_minimal()



# -------------------------------------------------------------------
# PART IV: get subbasin stats
# compute area of each polygon (km^2)
areas_km2 <- terra::expanse(all_subbasins_subset, unit = "km")

# summary stats
subbasin_stats <- tibble(
  total_subbasins = length(areas_km2),
  total_area = sum(areas_km2),
  mean_area  = mean(areas_km2),
  sd_area    = sd(areas_km2),
  median_area = median(areas_km2),
  min_area   = min(areas_km2),
  max_area   = max(areas_km2)
)

# subbasin_stats
# # A tibble: 1 × 7
# total_subbasins total_area mean_area sd_area median_area min_area max_area
# <int>      <dbl>     <dbl>   <dbl>       <dbl>    <dbl>    <dbl>
#   1             295   4068351.    13791.  17513.       9408.    1355.  179656.
