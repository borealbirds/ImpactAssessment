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

bam_boundary <- terra::vect(file.path(root, "Regions", "BAM_BCR_NationalModel_UnBuffered.shp"))


# ---------------------------------------------------------------
# OBJECTIVE 1: create and save a low HF raster for model training

# import an arbitrary bam template for re-projecting HF layer
bam_template <- terra::rast(file.path(root, "PredictionRasters", "Biomass", "SCANFI", "1km", "SCANFIBalsamFir_1km_2020.tif"))


# import Hirsh-Pearson HF layer
CanHF_1km <- 
  terra::rast(file.path(root, "CovariateRasters", "Disturbance", "cum_threat2020.02.18.tif")) |> 
  terra::project(x = _, y = bam_template) |>
  terra::resample(x = _, y = bam_template) |> 
  terra::crop(x = _, y = bam_template) |> 
  terra::mask(x = _, mask = bam_template)
  
names(CanHF_1km) <- "CanHF_1km"


# ---------------------------------------------------------------
# OBJECTIVE 2: identify subbasins with zero or few low HF pixels,


# import subbasins multi-polygon
all_subbasins <- terra::vect(file.path(ia_dir, "hydrobasins_masked.gpkg"))


# set low HF to <1 
# from Hirsh-Pearson: "we found that 82% of Canada’s land areas had a 
# HF < 1 and therefore were considered intact"
lowhf_mask <- 
  terra::lapp(CanHF_1km, \(x) ifelse(!is.finite(x), NA, ifelse(x < 1, 1, NA))) |> 
  as.factor()

# per-subbasin lowhf counts 
# note: many subbasins don't have any pixels with HF < 1
counts_df <- terra::extract(
  lowhf_mask,
  all_subbasins,
  fun   = function(x) sum(!is.na(x)),  # count non-NA cells
  ID    = TRUE)

sum(counts_df$CanHF_1km) #5,552,391 low HF pixels
quantile(counts_df$CanHF_1km) 



# ---------------------------------------------------------------
# OBJECTIVE 3: visualize low HF pixel density per subbasin

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

# sample size threshold is 25 percentile (1085 pixels)
threshold <- quantile(counts_df$CanHF_1km)[2]


# 1) neighbor matrix (touching polygons)
nb_mat <- terra::relate(all_subbasins, all_subbasins, relation = "touches")

# row-wise: indices of touching neighbors
nb <- apply(nb_mat, 1, function(r) which(r))

merge_id <- seq_len(nrow(all_subbasins))
taken    <- rep(FALSE, nrow(all_subbasins))
ord      <- order(all_subbasins$sub_count)

for (i in ord) {
  if (taken[i] || all_subbasins$sub_count[i] >= threshold) next
  
  b   <- all_subbasins$BCR[i]
  grp <- i
  tot <- all_subbasins$sub_count[i]
  
  fr <- nb[[i]]
  fr <- fr[ all_subbasins$BCR[fr] == b & !taken[fr] & !(fr %in% grp) ]
  
  while (tot < threshold && length(fr) > 0) {
    j <- fr[ which.min(all_subbasins$sub_count[fr]) ]
    grp <- c(grp, j)
    tot <- tot + all_subbasins$sub_count[j]
    taken[j] <- TRUE
    
    nbj <- nb[[j]]
    nbj <- nbj[ all_subbasins$BCR[nbj] == b & !taken[nbj] & !(nbj %in% grp) ]
    fr  <- unique(c(setdiff(fr, j), nbj))
  }
  
  merge_id[grp] <- i
  taken[grp]    <- TRUE
}

all_subbasins$merge_id <- merge_id

# 2) dissolve by merge_id; sum numeric fields (drop BCR before summing)
tmp <- all_subbasins; tmp$BCR <- NULL
merged_subs <- aggregate(tmp, by = "merge_id", fun = sum)

# 3) reattach BCR from seed members (all members share a BCR)
bcr_map <- unique(data.frame(merge_id = merge_id, BCR = all_subbasins$BCR))
merged_subs <- left_join(merged_subs, bcr_map, by = "merge_id")

# quick check
table(merged_subs$sub_count >= threshold)

# 4) per-subbasin lowhf counts 
# note: many subbasins don't have any pixels with HF < 1
counts_df2 <- terra::extract(
  lowhf_mask,
  merged_subs,
  fun   = function(x) sum(!is.na(x)),  # count non-NA cells
  ID    = TRUE)

sum(counts_df2$CanHF_1km) #1,917,452 (same as counts_df)
quantile(counts_df2$CanHF_1km) # range remains 0-29466, but 25 percentile is now 2404
min(counts_df2$CanHF_1km)


# visualize new subbasin boundaries (exported at 1000 x 751)
ggplot() +
  
  # points colored by subbasin density
  geom_point(data = slice_sample(plot_df, prop=0.04), aes(x = x, y = y, colour = sub_count), size = 0.01) +
  scale_color_gradient(low="#56B4E9", high="#CC79A7") +
  
  # BCR and basin outlines
  geom_spatvector(data = bam_boundary, fill = NA, colour = "grey") +
  geom_spatvector(data = merged_subs, fill = NA, color = "grey25", linewidth = 0.3) +
  
  coord_sf(crs = crs(merged_subs)) +
  theme_minimal()

# SHOULD BE DONE BY THIS POINT






# attach counts and flag under-threshold
merged_subs$lowhf_n   <- counts_df2$CanHF_1km[match(seq_len(nrow(merged_subs)), counts_df2$ID)]
merged_subs$under_thr <- merged_subs$lowhf_n < threshold

# quick count
table(merged_subs$under_thr)

# map
ggplot() +
  geom_spatvector(data = merged_subs, aes(fill = under_thr), color = "grey30", linewidth = 0.2) +
  geom_spatvector(data = bam_boundary, fill = NA, color = "black", linewidth = 0.6) +
  
  scale_fill_manual(values = c(`TRUE` = "#F8766D", `FALSE` = "white"),
                    name = paste0("< ", threshold, " low-HF px")) +
  coord_sf(crs = crs(merged_subs)) +
  theme_minimal()




terra::writeRaster(lowhf_mask, file.path(ia_dir, "CanHF_1km_lessthan1.tif"))

