# ---
# title: Impact Assessment: merge data-sparse subbasins within BCRs
# author: Mannfred Boehm
# created: September 16, 2025
# ---

# subbasins with few low HF (training) pixels are merged with their nearest neighbour in the
# same BCR until each has at least the 25th percentile of low HF pixels; subbasins with no
# high HF pixel are then dropped (nothing to backfill). the result, 674 subbasins, is the set
# every later script indexes. the existing file was built in February 2026 from the low HF
# mask as it was before 03's 2026-09-28 fix (the numbers below are from that run), so a rerun
# will not reproduce it, and the final writeVector() refuses to overwrite it.

library(terra)
library(tidyterra)
library(tidyverse)

root    <- "G:/Shared drives/BAM_NationalModels5"
raw_dir <- file.path("data", "raw_data")

bam_boundary <- terra::vect(file.path(root, "Regions", "BAM_BCR_NationalModel_UnBuffered.shp"))
bam_boundary <- bam_boundary[bam_boundary$subUnit != 23, ]   # no covariate data for BCR 23



# ---------------------------------------------------------------
# count low HF pixels per subbasin

lowhf_mask    <- terra::rast(file.path(raw_dir, "hirshpearson", "CanHF_1km_lessthan1.tif"))
all_subbasins <- terra::vect(file.path(raw_dir, "hydrobasins_masked.gpkg"))

# many subbasins have no pixels with HF < 1
counts_df <-
  terra::extract(lowhf_mask, all_subbasins, fun = function(x) sum(!is.na(x)), ID = TRUE) |>
  dplyr::as_tibble() |>
  dplyr::arrange(CanHF_1km)

sum(counts_df$CanHF_1km) #5,755,889 low HF pixels
quantile(counts_df$CanHF_1km)
# 0%   25%   50%   75%  100%
# 0   892  2949  6742 29787



# ---------------------------------------------------------------
# map low HF pixels, coloured by their subbasin's count (exported at 1000 x 751)

all_subbasins$sub_count <- counts_df$CanHF_1km[match(seq_len(nrow(all_subbasins)), counts_df$ID)]

# the ~5 million low HF pixels as points, each carrying its subbasin's count
ij <- terra::intersect(terra::as.points(lowhf_mask, na.rm = TRUE), all_subbasins)
plot_df <- cbind(terra::crds(ij, df = TRUE), sub_count = ij$sub_count)

ggplot() +
  geom_point(data = slice_sample(plot_df, prop = 0.04), aes(x = x, y = y, colour = sub_count), size = 0.01) +
  scale_color_gradient(low = "#56B4E9", high = "#CC79A7") +
  geom_spatvector(data = bam_boundary, fill = NA, colour = "grey") +
  geom_spatvector(data = all_subbasins, fill = NA, color = "grey25", linewidth = 0.3) +
  coord_sf(crs = crs(all_subbasins)) +
  theme_minimal()



# ---------------------------------------------------------------
# merge low pixel density subbasins with their nearest subbasin, within the same BCR

threshold <- quantile(counts_df$CanHF_1km, probs = 0.25, na.rm = TRUE)

# returns each subbasin's group_key: subbasins sharing a key are merged
merge_by_nearest_subbasin <- function(all_subbasins, bam_boundary, counts_df, threshold) {

  # low HF counts in row order (extract's ID is the row index)
  cnt <- integer(nrow(all_subbasins))
  cnt[counts_df$ID] <- counts_df$CanHF_1km

  # each subbasin's BCR: the one it overlaps most
  ij <- terra::intersect(all_subbasins["HYBAS_ID"], bam_boundary["subUnit"])
  ij$area <- terra::expanse(ij)
  m <- dplyr::as_tibble(ij) |>
    dplyr::group_by(HYBAS_ID) |>
    dplyr::slice_max(area, n = 1, with_ties = FALSE) |>
    dplyr::ungroup()
  bcr_by_row <- m$subUnit[match(all_subbasins$HYBAS_ID, m$HYBAS_ID)]

  # geometry of a set of subbasins dissolved into one, repaired only if invalid
  dissolve_members <- function(idx) {
    g  <- terra::aggregate(all_subbasins[idx])[, 0]
    iv <- tryCatch(terra::is.valid(g), error = function(e) FALSE)
    if (!all(iv, na.rm = TRUE)) g <- tryCatch(terra::buffer(g, 0), error = function(e) g)
    g
  }

  # one pass per BCR (subbasins without a BCR form their own group)
  split_key <- ifelse(is.na(bcr_by_row), "__NA__", as.character(bcr_by_row))
  split_idx <- split(seq_len(nrow(all_subbasins)), f = split_key, drop = FALSE)

  group_key <- integer(nrow(all_subbasins))
  next_key  <- 1L

  for (key in names(split_idx)) {
    idx <- split_idx[[key]]
    if (!length(idx)) next

    # clusters in this BCR: member row indices, low HF count, geometry
    members   <- as.list(idx)
    sum_cnt   <- cnt[idx]
    geom_list <- lapply(idx, function(i) all_subbasins[i, 0])

    # the cluster furthest below the threshold goes first
    pick_anchor <- function() {
      under <- which(sum_cnt < threshold)
      if (!length(under)) return(NA_integer_)
      under[order(threshold - sum_cnt[under], decreasing = TRUE)][1]
    }

    repeat {
      i <- pick_anchor()
      if (is.na(i) || length(sum_cnt) == 1) break   # done, or one cluster left under the threshold

      # merge cluster i with its nearest cluster (nb) until it reaches the threshold
      while (sum_cnt[i] < threshold && length(sum_cnt) > 1) {
        pool <- setdiff(seq_along(sum_cnt), i)
        d    <- as.numeric(terra::distance(geom_list[[i]], do.call(rbind, geom_list[pool])))
        nb   <- pool[which.min(d)]

        members[[i]]   <- c(members[[i]], members[[nb]])
        sum_cnt[i]     <- sum_cnt[i] + sum_cnt[nb]
        geom_list[[i]] <- dissolve_members(members[[i]])

        members   <- members[-nb]
        sum_cnt   <- sum_cnt[-nb]
        geom_list <- geom_list[-nb]
        if (nb < i) i <- i - 1L
      }
    }

    for (g in seq_along(members)) {
      group_key[members[[g]]] <- next_key
      next_key <- next_key + 1L
    }
  }

  tibble::tibble(HYBAS_ID = all_subbasins$HYBAS_ID, group_key = group_key)
}

grouping <- merge_by_nearest_subbasin(all_subbasins, bam_boundary, counts_df, threshold)

# low HF count per merged subbasin
counts_df2 <-
  counts_df |>
  dplyr::transmute(ID = grouping$group_key[ID], CanHF_1km = CanHF_1km) |>
  dplyr::group_by(ID) |>
  dplyr::summarise(CanHF_1km = sum(CanHF_1km, na.rm = TRUE), .groups = "drop") |>
  dplyr::arrange(CanHF_1km)

sum(counts_df2$CanHF_1km) # equals sum(counts_df$CanHF_1km)
quantile(counts_df2$CanHF_1km)
# 0%     25%     50%     75%    100%
# 902.0  2468.5  4351.0  8420.5 30102.0



# ---------------------------------------------------------------
# map merged subbasins (exported at 1000 x 751)

all_subbasins_merged <- terra::aggregate(
  dplyr::left_join(all_subbasins, grouping, by = "HYBAS_ID"),
  by  = "group_key",
  fun = "first"
)

all_subbasins_merged$sub_count <- counts_df2$CanHF_1km[match(all_subbasins_merged$group_key, counts_df2$ID)]

ij2 <- terra::intersect(terra::as.points(lowhf_mask, na.rm = TRUE), all_subbasins_merged)
plot_df2 <- cbind(terra::crds(ij2, df = TRUE), sub_count = ij2$sub_count)

ggplot() +
  geom_point(data = slice_sample(plot_df2, prop = 0.04), aes(x = x, y = y, colour = sub_count), size = 0.01) +
  scale_color_gradient(low = "#56B4E9", high = "#CC79A7") +
  geom_spatvector(data = bam_boundary, fill = NA, colour = "grey") +
  geom_spatvector(data = all_subbasins_merged, fill = NA, color = "grey25", linewidth = 0.3) +
  coord_sf(crs = crs(all_subbasins_merged)) +
  theme_minimal()

terra::writeVector(all_subbasins_merged, file.path(raw_dir, "hydrobasins_masked_merged.gpkg"), overwrite = TRUE)

areas_km2 <- terra::expanse(all_subbasins_merged, unit = "km")
tibble(total_subbasins = length(areas_km2),
       total_area      = sum(areas_km2),
       mean_area       = mean(areas_km2),
       sd_area         = sd(areas_km2),
       median_area     = median(areas_km2),
       min_area        = min(areas_km2),
       max_area        = max(areas_km2))
# total_subbasins total_area mean_area sd_area median_area min_area max_area
#             887   8613815.     9711.  11704.       6557.     902.  179656.



# -----------------------------------------------------------
# keep only subbasins with high human footprint (HF >= 1): the rest have nothing to backfill

highhf_mask <- terra::rast(file.path(raw_dir, "hirshpearson", "CanHF_1km_morethan1.tif"))
all_subbasins_merged <- terra::vect(file.path(raw_dir, "hydrobasins_masked_merged.gpkg"))

counts_df3 <-
  terra::extract(highhf_mask, all_subbasins_merged, fun = function(x) sum(!is.na(x)), ID = TRUE) |>
  dplyr::as_tibble()

has_highhf <- counts_df3$CanHF_1km > 0
all_subbasins_subset <- all_subbasins_merged[counts_df3$ID[has_highhf]]

# overwrite = FALSE: this file defines the production subbasins (see header)
terra::writeVector(all_subbasins_subset, file.path(raw_dir, "hydrobasins_masked_merged_subset.gpkg"), overwrite = FALSE)

sum(counts_df3$CanHF_1km) #1,787,681 high HF pixels
quantile(counts_df3$CanHF_1km)
# 0%      25%      50%      75%     100%
# 0.0      1.5    220.0   1304.0 109057.0



# -------------------------------------------------------------------
# map high HF pixels, coloured by their subbasin's count (exported at 1000 x 751)

all_subbasins_subset$sub_count <- counts_df3$CanHF_1km[has_highhf]

ij <- terra::intersect(terra::as.points(highhf_mask, na.rm = TRUE), all_subbasins_subset)
plot_df <- cbind(terra::crds(ij, df = TRUE), sub_count = ij$sub_count)

ggplot() +
  geom_point(data = slice_sample(plot_df, prop = 0.04), aes(x = x, y = y, colour = sub_count), size = 0.01) +
  scale_color_gradient(low = "#56B4E9", high = "#CC79A7") +
  geom_spatvector(data = bam_boundary, fill = NA, colour = "grey") +
  geom_spatvector(data = all_subbasins_subset, fill = NA, color = "grey25", linewidth = 0.3) +
  coord_sf(crs = crs(all_subbasins_subset)) +
  theme_minimal()



# -----------------------------------------------------------
# map low HF pixels of the kept subbasins, coloured by their subbasin's count (exported at 1000 x 751)

all_subbasins_subset <- terra::vect(file.path(raw_dir, "hydrobasins_masked_merged_subset.gpkg"))

lowhf_mask <-
  terra::rast(file.path(raw_dir, "hirshpearson", "CanHF_1km_lessthan1.tif")) |>
  terra::crop(all_subbasins_subset) |>
  terra::mask(all_subbasins_subset)

# the mask is all 1s, so each subbasin's sum is its low HF count
zones <- terra::rasterize(all_subbasins_subset, lowhf_mask, field = "group_key")
z <-
  terra::zonal(lowhf_mask, zones, fun = "sum", na.rm = TRUE) |>
  dplyr::as_tibble() |>
  dplyr::rename(sub_count = CanHF_1km)

# a 4% sample of the low HF pixels, each tagged with its subbasin's count
pts   <- terra::as.points(lowhf_mask, na.rm = TRUE)
pts_s <- pts[sample(1:nrow(pts), size = round(nrow(pts) * 0.04)), ]
pts_s$group_key <- terra::extract(zones, pts_s, ID = FALSE)[, 1]

pts_df <-
  cbind(terra::crds(pts_s), data.frame(pts_s)) |>
  as_tibble() |>
  left_join(z, by = "group_key") |>
  mutate(sub_count = ifelse(is.na(sub_count), 0L, sub_count))

ggplot() +
  geom_point(data = pts_df, aes(x = x, y = y, colour = sub_count), size = 0.01) +
  scale_color_gradient(low = "#56B4E9", high = "#CC79A7") +
  geom_spatvector(data = bam_boundary, fill = NA, colour = "grey") +
  geom_spatvector(data = all_subbasins_subset, fill = NA, color = "grey25", linewidth = 0.3) +
  coord_sf(crs = crs(all_subbasins_subset)) +
  theme_minimal()
