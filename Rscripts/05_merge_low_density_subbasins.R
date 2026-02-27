# ---
# title: Impact Assessment: create low/high HF layers and merge data-sparse subbasins 
# author: Mannfred Boehm
# created: September 16, 2025
# ---

library(terra)
library(tidyterra)
library(tidyverse)

root <- "G:/Shared drives/BAM_NationalModels5"
ia_dir <- file.path(root, "data", "Extras", "sandbox_data", "impactassessment_sandbox")



# restrict analysis to BCRs that we have covariate data for
bam_template <- terra::rast(file.path(root, "PredictionRasters", "Biomass", "SCANFI", "1km", "SCANFIBalsamFir_1km_2020.tif"))

bam_boundary <- terra::vect(file.path(root, "Regions", "BAM_BCR_NationalModel_UnBuffered.shp"))
bam_boundary <- bam_boundary[bam_boundary$subUnit != 23, ]



# ---------------------------------------------------------------
# identify subbasins with zero or few low HF pixels,

# import low HF raster
lowhf_mask <- terra::rast(file.path(ia_dir, "hirshpearson", "CanHF_1km_lessthan1.tif"))

# import subbasins multi-polygon
all_subbasins <- terra::vect(file.path(ia_dir, "hydrobasins_masked.gpkg"))


# per-subbasin lowhf counts 
# note: many subbasins don't have any pixels with HF < 1
# only count non-NA cells
counts_df <- 
  terra::extract(lowhf_mask, all_subbasins, fun=function(x) sum(!is.na(x)), ID=TRUE) |> 
  dplyr::as_tibble() |>
  dplyr::arrange(CanHF_1km)
  

sum(counts_df$CanHF_1km) #5,755,889 low HF pixels
quantile(counts_df$CanHF_1km) 
# 0%   25%   50%   75%  100% 
# 0   892  2949  6742 29787 



# ---------------------------------------------------------------
# visualize low HF pixel density per subbasin

# assign a low HF pixel count to each subbasin
all_subbasins$sub_count <- counts_df$CanHF_1km[match(seq_len(nrow(all_subbasins)), counts_df$ID)]

# join each of the 5 million points to its subbasin ID
ij <- terra::intersect(as.points(lowhf_mask, na.rm = TRUE), all_subbasins)

# paste subbasin-level lowHF counts back onto each point
# every point is aware of if it's subbasin pixel density
ij$sub_count <- all_subbasins$sub_count[ match(ij$HYBAS_ID, all_subbasins$HYBAS_ID) ]

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
geom_spatvector(data = all_subbasins, fill = NA, color = "grey25", linewidth = 0.3) +

coord_sf(crs = crs(all_subbasins_subset)) +
theme_minimal()



# ---------------------------------------------------------------
# merge low pixel density subbasins with 
# nearest subbasins, but only if they're in the same BCR

threshold <- quantile(counts_df$CanHF_1km, probs = 0.25, na.rm = TRUE)

merge_by_nearest_subbasin <- function(all_subbasins, bam_boundary, counts_df, threshold) {
  # counts aligned to row order (ID from terra::extract == row index of all_subbasins)
  cnt <- { v <- integer(nrow(all_subbasins)); v[counts_df$ID] <- counts_df$CanHF_1km; v }
  
  # map each subbasin to its BCR (max overlap)
  bcr_by_row <- {
    ij <- terra::intersect(all_subbasins["HYBAS_ID"], bam_boundary["subUnit"])
    ij$area <- terra::expanse(ij)
    m <- dplyr::as_tibble(ij) |>
      dplyr::group_by(HYBAS_ID) |>
      dplyr::slice_max(area, n = 1, with_ties = FALSE) |>
      dplyr::ungroup() |>
      dplyr::select(HYBAS_ID, subUnit)
    m$subUnit[match(all_subbasins$HYBAS_ID, m$HYBAS_ID)]
  }
  
  # dissolve helper: returns geometry-only; fix invalid only if needed
  dissolve_members <- function(idx) {
    g <- terra::aggregate(all_subbasins[idx])[, 0]
    iv <- tryCatch(terra::is.valid(g), error = function(e) FALSE)
    if (!all(iv, na.rm = TRUE)) {
      g <- tryCatch(terra::buffer(g, 0), error = function(e) g)
    }
    g
  }
  
  # split processing by BCR (treat NA as its own key)
  split_key <- ifelse(is.na(bcr_by_row), "__NA__", as.character(bcr_by_row))
  split_idx <- split(seq_len(nrow(all_subbasins)), f = split_key, drop = FALSE)
  
  group_key <- integer(nrow(all_subbasins))
  next_key  <- 1L
  
  for (key in names(split_idx)) {
    idx <- split_idx[[key]]
    if (!length(idx)) next
    
    # state within this BCR
    members   <- lapply(idx, function(i) i)         # list of *global* row-index vectors
    sum_cnt   <- cnt[idx]                           # counts per local cluster
    geom_list <- lapply(idx, function(i) all_subbasins[i, 0])  # geometry-only per cluster
    
    # anchor picker: largest deficit first
    pick_anchor <- function() {
      under <- which(sum_cnt < threshold)
      if (!length(under)) return(NA_integer_)
      under[order(threshold - sum_cnt[under], decreasing = TRUE)][1]
    }
    
    while (TRUE) {
      i <- pick_anchor()
      if (is.na(i)) break
      
      repeat {
        if (sum_cnt[i] >= threshold) break
        pool <- if (length(sum_cnt) > 1) setdiff(seq_along(sum_cnt), i) else integer(0)
        if (!length(pool)) break
        
        # distance from current cluster i to all other clusters in this BCR
        gi <- geom_list[[i]]
        gp <- do.call(rbind, geom_list[pool])
        d  <- as.numeric(terra::distance(gi, gp))
        nb <- pool[which.min(d)]
        
        # merge i and nb -> update only the merged cluster, drop nb
        members[[i]] <- c(members[[i]], members[[nb]])
        sum_cnt[i]   <- sum_cnt[i] + sum_cnt[nb]
        geom_list[[i]] <- dissolve_members(members[[i]])
        
        members   <- members[-nb]
        sum_cnt   <- sum_cnt[-nb]
        geom_list <- geom_list[-nb]
        
        if (nb < i) i <- i - 1L
      }
    }
    
    # assign group IDs for this BCR
    for (g in seq_along(members)) {
      group_key[members[[g]]] <- next_key
      next_key <- next_key + 1L
    }
  }
  
  tibble::tibble(HYBAS_ID = all_subbasins$HYBAS_ID, group_key = group_key)
}

# merge low density subbasins until there are no subbasin densities below `threshold`
# group_key = the merged unit ID
# subbasins with the same `group_key` will be merged together
grouping <- merge_by_nearest_subbasin(all_subbasins, bam_boundary, counts_df, threshold)

# dissolve polygons by group 
merged_units <- 
  all_subbasins |>
  dplyr::left_join(grouping, by = "HYBAS_ID") |>
  terra::aggregate(by = "group_key", fun = "first")


# check low HF density by merged group
counts_df2 <- 
  counts_df |>
  dplyr::transmute(ID = grouping$group_key[ID], CanHF_1km = CanHF_1km) |>
  dplyr::group_by(ID) |>
  dplyr::summarise(CanHF_1km = sum(CanHF_1km, na.rm = TRUE), .groups = "drop") |>
  dplyr::arrange(CanHF_1km)

# per-subbasin lowhf counts 
sum(counts_df2$CanHF_1km) #5,719,097 (same as counts_df)


quantile(counts_df2$CanHF_1km) 
# 0%     25%     50%     75%    100% 
# 902.0  2468.5  4351.0  8420.5 30102.0 

# ---------------------------------------------------------------
# visualize merged subbasins

# dissolve to merged units
all_subbasins_merged <- terra::aggregate(
  dplyr::left_join(all_subbasins, grouping, by = "HYBAS_ID"),
  by  = "group_key",
  fun = "first"
)

# assign a low HF pixel count to each subbasin
# map by group_key, not by row index
all_subbasins_merged$sub_count <- counts_df2$CanHF_1km[
  match(all_subbasins_merged$group_key, counts_df2$ID)
]


# join each of the 5 million points to its subbasin ID
ij2 <- terra::intersect(as.points(lowhf_mask, na.rm = TRUE), all_subbasins_merged)

# paste subbasin-level lowHF counts back onto each point
# every point is aware of if it's subbasin pixel density
# paste group-level lowHF counts back onto each point — match by group_key
ij2$sub_count <- all_subbasins_merged$sub_count[
  match(ij2$group_key, all_subbasins_merged$group_key)
]

# generate dataframe: rows are low HF points, 
# columns are lat, long, subbasin pixel density for that point
plot_df2 <- cbind(terra::crds(ij2, df = TRUE), sub_count = ij2$sub_count)


# visualize new subbasin boundaries (exported at 1000 x 751)
ggplot() +
  
  # points colored by subbasin density
  geom_point(data = slice_sample(plot_df2, prop=0.04), aes(x = x, y = y, colour = sub_count), size = 0.01) +
  scale_color_gradient(low="#56B4E9", high="#CC79A7") +
  
  # BCR and basin outlines
  geom_spatvector(data = bam_boundary, fill = NA, colour = "grey") +
  geom_spatvector(data = all_subbasins_merged, fill = NA, color = "grey25", linewidth = 0.3) +
  
  coord_sf(crs = crs(all_subbasins_merged)) +
  theme_minimal()


# save 
terra::writeVector(all_subbasins_merged, file.path(ia_dir, "hydrobasins_masked_merged.gpkg"), overwrite=TRUE)


#  get subbasin stats
# compute area of each polygon (km^2)
areas_km2 <- terra::expanse(all_subbasins_merged, unit = "km")

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
#   1             887   8613815.     9711.  11704.       6557.     902.  179656.



# -----------------------------------------------------------
# subset subbasins to just those with high human footprint (ie HF >= 1)
# because we don't need to train or backfill areas without disturbance

# import raster with high human footprint
highhf_mask <- terra::rast(file.path(ia_dir, "hirshpearson", "CanHF_1km_morethan1.tif"))

# import subbasins multi-polygon (merged so that every subbasin has sufficient low HF pixels)
all_subbasins_merged <- terra::vect(file.path(ia_dir, "hydrobasins_masked_merged.gpkg"))

# identify subbasins that intersect with high human footprint
# first, count non-NA high-HF cells per subbasin
counts_df3 <- terra::extract(
  highhf_mask, all_subbasins_merged,
  fun = function(x) sum(!is.na(x)), ID = TRUE) |> 
  dplyr::as_tibble()

# keep only those subbasins with at least 1 high-HF pixel
keep_ids <- counts_df3$ID[counts_df3$CanHF_1km > 0]  
all_subbasins_subset <- all_subbasins_merged[keep_ids]

# save 
terra::writeVector(all_subbasins_subset, file.path(ia_dir, "hydrobasins_masked_merged_subset.gpkg"), overwrite=TRUE)


# -------------------------------------------------------------------
# visualize high human footprint pixels per subbasin

sum(counts_df3$CanHF_1km) #1,787,681 high HF pixels
quantile(counts_df3$CanHF_1km)
# 0%      25%      50%      75%     100% 
# 0.0      1.5    220.0   1304.0 109057.0 

# assign a high HF pixel count to each subbasin
all_subbasins_subset$sub_count <- counts_df3$CanHF_1km[match(seq_len(nrow(all_subbasins_subset)), counts_df3$ID)]

# join each of the 1 million points to its subbasin ID
ij <- terra::intersect(as.points(highhf_mask, na.rm = TRUE), all_subbasins_subset)

# paste subbasin-level high HF counts back onto each point
# every point is aware of if it's subbasin pixel density
ij$sub_count <- all_subbasins_subset$sub_count[ match(ij$first_HYBAS_ID, all_subbasins_subset$first_HYBAS_ID) ]

# generate dataframe: rows are high HF points, 
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


# -----------------------------------------------------------
# visualise low HF pixels for the subsetted subbasin set


# import subbsetted subbasin boundaries
all_subbasins_subset <- terra::vect(file.path(ia_dir, "hydrobasins_masked_merged_subset.gpkg"))

# import low HF raster
lowhf_mask <- 
  terra::rast(file.path(ia_dir, "hirshpearson", "CanHF_1km_lessthan1.tif")) |> 
  terra::crop(x = _, y = all_subbasins_subset) |> 
  terra::mask(x = _, mask = all_subbasins_subset)


# rasterize subbasins onto the lowHF grid (use a stable ID; here 'group_key')
zones <- rasterize(all_subbasins_subset, lowhf_mask, field = "group_key")

# zonal sum: since lowhf_mask is 1s, the sum = number of low-HF pixels
z <- 
  terra::zonal(lowhf_mask, zones, fun = "sum", na.rm = TRUE) |>
  dplyr::as_tibble() |>
  dplyr::rename(sub_count = CanHF_1km)

# make points from the low-HF mask (1s only); sample ~4%
pts <- terra::as.points(lowhf_mask, na.rm = TRUE)
pts_s <- pts[sample(1:nrow(pts), size = round(nrow(pts) * 0.04)), ]

# 2) Attach subbasin group_key to each point (lookup from zones raster)
pts_s$group_key <- terra::extract(zones, pts_s, ID = FALSE)[, 1]

# 3) Join subbasin counts
pts_df <- cbind(terra::crds(pts_s), data.frame(pts_s)) |>
  as_tibble() |>
  left_join(z, by = "group_key") |>
  mutate(sub_count = ifelse(is.na(sub_count), 0L, sub_count))

# visualize new subbasin boundaries (exported at 1000 x 751)
ggplot() +
  
  # points colored by subbasin density
  geom_point(data = pts_df, aes(x = x, y = y, colour = sub_count), size = 0.01) +
  scale_color_gradient(low="#56B4E9", high="#CC79A7") +
  
  # BCR and basin outlines
  geom_spatvector(data = bam_boundary, fill = NA, colour = "grey") +
  geom_spatvector(data = all_subbasins_subset, fill = NA, color = "grey25", linewidth = 0.3) +
  
  coord_sf(crs = crs(all_subbasins_subset)) +
  theme_minimal()

