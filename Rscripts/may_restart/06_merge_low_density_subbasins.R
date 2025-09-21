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
# OBJECTIVE 1: create low HF raster to count low HF pixels per subbasin



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



# ---------------------------------------------------------------
# OBJECTIVE 2: identify subbasins with zero or few low HF pixels,


# import subbasins multi-polygon
all_subbasins <- terra::vect(file.path(ia_dir, "hydrobasins_masked.gpkg"))


# per-subbasin lowhf counts 
# note: many subbasins don't have any pixels with HF < 1
# only count non-NA cells
counts_df <- 
  terra::extract(lowhf_mask, all_subbasins, fun=function(x) sum(!is.na(x)), ID=TRUE) |> 
  dplyr::as_tibble() |>
  dplyr::arrange(CanHF_1km)
  

sum(counts_df$CanHF_1km) #5,755,891 low HF pixels
quantile(counts_df$CanHF_1km) 
# quantile(counts_df$CanHF_1km) 
# 0%   25%   50%   75%  100% 
# 0   899  2964  6740 29787 



# ---------------------------------------------------------------
# OBJECTIVE 3: visualize low HF pixel density per subbasin

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

coord_sf(crs = crs(all_subbasins)) +
theme_minimal()



# ---------------------------------------------------------------
# OBJECTIVE 4: merge low pixel density subbasins with 
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
sum(counts_df2$CanHF_1km) #5755891 (same as counts_df)

# range now 30102 (from 902-30102), but 25 percentile is now 2468 (from 899)
quantile(counts_df2$CanHF_1km) 


# ---------------------------------------------------------------
# OBJECTIVE 5: visualize merged subbasins

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
