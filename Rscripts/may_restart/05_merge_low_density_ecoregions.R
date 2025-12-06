# ---
# title: Impact Assessment: create low/high HF layers and merge data-sparse ecoregions
# author: Mannfred Boehm
# created: December 5, 2025
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
# identify ecoregions with zero or few low HF pixels,

# import low HF raster
lowhf_mask <- terra::rast(file.path(ia_dir, "hirshpearson", "CanHF_1km_lessthan1.tif"))

# import ecoregions multi-polygon
ecoregions <- terra::vect(file.path(ia_dir, "ecoregions_masked.gpkg"))


# per-ecoregion lowhf counts 
# note: some ecoregions don't have any pixels with HF < 1
# only count non-NA cells
counts_df <- 
  terra::extract(lowhf_mask, ecoregions, fun=function(x) sum(!is.na(x)), ID=TRUE) |> 
  dplyr::as_tibble() |>
  dplyr::arrange(CanHF_1km)


sum(counts_df$CanHF_1km) #5,756,093 low HF pixels
quantile(counts_df$CanHF_1km, probs=seq(0,1,0.1)) 
#0%      10%      20%      30%      40%      50%      60%      70% 
#  0.0    165.0   1915.0   3405.2   6742.6  12993.5  21676.4  34097.4 
#80%      90%     100% 
#54258.2  90028.9 173519.0 



# ---------------------------------------------------------------
# visualize low HF pixel density per ecoregion

# assign a low HF pixel count to each subbasin
ecoregions$REGION_COUNT <- counts_df$CanHF_1km[match(seq_len(nrow(ecoregions)), counts_df$ID)]

# join each of the 5 million points to its ecoregion ID
ij <- terra::intersect(as.points(lowhf_mask, na.rm = TRUE), ecoregions)

# paste subbasin-level lowHF counts back onto each point
# every point is aware of if it's subbasin pixel density
ij$REGION_COUNT <- ecoregions$REGION_COUNT[ match(ij$REGION_ID, ecoregions$REGION_ID) ]
saveRDS(ij, file="ij2.rds")
ij<-readRDS(file="ij2.rds")

# generate dataframe: rows are low HF points, 
# columns are lat, long, subbasin pixel density for that point
plot_df <- cbind(terra::crds(ij, df = TRUE), REGION_COUNT = ij$REGION_COUNT)

# plot (exported at 1000 x 751)
ggplot() +
  
  # points colored by subbasin density
  geom_point(data = slice_sample(plot_df, prop=0.01), aes(x = x, y = y, colour = REGION_COUNT), size = 0.01) +
  scale_color_gradient(low="#56B4E9", high="#CC79A7") +
  
  # BCR and basin outlines
  geom_spatvector(data = bam_boundary, fill = NA, colour = "#999999", linewidth = 0.3) +
  geom_spatvector(data = ecoregions, fill = NA, color = "#000000", linewidth = 0.2) +
  
  coord_sf(crs = crs(ecoregions)) +
  theme_minimal()



# ---------------------------------------------------------------
# merge low pixel density ecoregions with 
# nearest ecoregions, but only if they're in the same BCR

threshold <- quantile(counts_df$CanHF_1km, probs = 0.20, na.rm = TRUE)

merge_by_nearest_ecoregion <- function(ecoregions, bam_boundary, counts_df, threshold) {
  
  t_start <- Sys.time()
  cat("\n=== start merge_by_nearest_ecoregion ===\n")
  cat("start time:", format(t_start), "\n\n")
  
  cat("Checking for nested ecoregions (interior-point in polygon)...\n")
  
  # Get a guaranteed interior point per polygon
  intpts <- terra::centroids(ecoregions, inside = TRUE)
  
  # Extract polygon containing each interior point
  hits <- terra::extract(ecoregions["REGION_ID"], intpts)
  is_nested <- which(hits$REGION_ID != own_ids)
  print(names())  # for debugging
  
  # Find the ID column (row index of point)
  id_col <- intersect(names(hits), c("ID", "id", "id.x", "id.y"))
  if (length(id_col) != 1)
    stop("Could not identify centroid/point ID column in extract() output.")
  
  # REGION_ID column
  reg_col <- intersect(names(hits), c("REGION_ID", "REGION_ID_1", "REGION_ID_2"))
  if (length(reg_col) != 1)
    stop("Could not identify REGION_ID column in extract() output.")
  
  point_row  <- hits[[id_col]]      # which interior point (== polygon index)
  container  <- hits[[reg_col]]     # REGION_ID of containing polygon
  
  own_ids <- ecoregions$REGION_ID
  
  contained_map <- integer(length(own_ids))
  contained_map[point_row] <- container
  
  # A polygon is nested if its interior point lies inside a *different* polygon
  is_nested <- which(contained_map != 0 & contained_map != own_ids)
  
  if (length(is_nested) == 0) {
    cat("No nested ecoregions detected.\n")
  } else {
    cat("Found", length(is_nested), "nested ecoregion(s).\n")
  }
  
  group_key <- integer(nrow(ecoregions))
  next_key  <- 1L
  
  for (key in names(split_idx)) {
    
    idx <- split_idx[[key]]
    if (!length(idx)) next
    
    # remove nested polygons from the active cluster set
    idx <- idx[cnt[idx] > 0]
    if (!length(idx)) {
      cat("  - no positive-count clusters in this BCR; skipping.\n")
      next
    }
    
    cat("=== processing BCR:", key, " | ecoregions:", length(idx), "===\n")
    
    # state within this BCR
    members   <- lapply(idx, function(i) i)         # list of *global* row-index vectors
    sum_cnt   <- cnt[idx]                           # counts per local cluster
    geom_list <- lapply(idx, function(i) ecoregions[i, 0])  # geometry-only per cluster
    
    # anchor picker: largest deficit first
    pick_anchor <- function() {
      under <- which(sum_cnt < threshold)
      if (!length(under)) return(NA_integer_)
      under[order(threshold - sum_cnt[under], decreasing = TRUE)][1]
    }
    
    # merging loop
    iter <- 1L
    repeat {
      i <- pick_anchor()
      if (is.na(i)) break
      
      cat("  · Iteration", iter,
          "| under-threshold cluster index:", i,
          "| current size:", sum_cnt[i], "\n")
      iter <- iter + 1L
      
      repeat {
        if (sum_cnt[i] >= threshold) {
          cat("    Cluster now meets threshold (", sum_cnt[i], "). Moving on.\n")
          break
        }
        
        pool <- if (length(sum_cnt) > 1) setdiff(seq_along(sum_cnt), i) else integer(0)
        if (!length(pool)) {
          cat("    no neighbors left to merge.\n")
          break
        }
        
        # distances
        gi <- geom_list[[i]]
        gp <- do.call(rbind, geom_list[pool])
        d  <- as.numeric(terra::distance(gi, gp))
        nb <- pool[which.min(d)]
        
        cat("    -> merging with neighbor", nb,
            "| neighbor size:", sum_cnt[nb],
            "| distance:", round(min(d), 2), "\n")
        
        # merge clusters
        members[[i]] <- c(members[[i]], members[[nb]])
        sum_cnt[i]   <- sum_cnt[i] + sum_cnt[nb]
        
        geom_list[[i]] <- dissolve_members(members[[i]])
        
        # drop neighbor
        members   <- members[-nb]
        sum_cnt   <- sum_cnt[-nb]
        geom_list <- geom_list[-nb]
        
        if (nb < i) i <- i - 1L
      }
    }
    
    # assign group keys
    for (g in seq_along(members)) {
      group_key[members[[g]]] <- next_key
      next_key <- next_key + 1L
    }
    
    cat("Finished BCR", key, "\n\n")
  }
  
  t_end <- Sys.time()
  cat("=== finished merge_by_nearest_ecoregion ===\n")
  cat("tnd time:", format(t_end), "\n")
  cat("total runtime:", round(as.numeric(t_end - t_start), 1), "seconds\n\n")
  
  tibble::tibble(REGION_ID = ecoregions$REGION_ID, group_key = group_key)
}

# merge low density ecoregions until there are no ecoregions densities below `threshold`
# group_key = the merged unit ID
# ecoregions with the same `group_key` will be merged together
# polygon merging is expensive this takes ~90 minutes
grouping <- merge_by_nearest_ecoregion(ecoregions, bam_boundary, counts_df, threshold)

# dissolve polygons by group 
merged_units <- 
  ecoregions |>
  dplyr::left_join(grouping, by = "REGION_ID") |>
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

