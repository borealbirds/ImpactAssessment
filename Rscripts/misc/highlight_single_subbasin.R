# plotting and highlighting specific subbasins
library(terra)
library(tidyterra)

# inspect subbasin indexes
# centroids for labeling
sub_centroids <- terra::centroids(all_subbasins_subset)

ggplot() +
  geom_spatvector(data = all_subbasins_subset,
                  fill = NA, color = "grey25", linewidth = 0.3) +
  
  # text labels at centroids
  geom_spatvector_text(
    data = sub_centroids,
    aes(label = seq_along(all_subbasins_subset)), 
    size = 1.5
  ) +
  
  coord_sf(crs = crs(all_subbasins_subset)) +
  theme_minimal()


# highlight specific subbasins
# 301 = great slave lake, 240 = west coast, 57 = prairies, 491= labrador, 100=NS
ggplot() +
  
  geom_spatvector(data = all_subbasins_subset, fill = NA, color = "grey25", linewidth = 0.3) +
  geom_spatvector(data = all_subbasins_subset[c(301,240,57,491,100)], fill = NA, color = "#CC79A7", linewidth = 1.2) +
  coord_sf(crs = crs(all_subbasins_subset)) +
  theme_minimal()



# visualize backfilling for a specific subbasin
root <- "G:/Shared drives/BAM_NationalModels5"
ia_dir <- file.path(root, "data", "Extras", "sandbox_data", "impactassessment_sandbox")

year <- 2020
stack_y <- terra::rast(file.path(ia_dir, sprintf("covariates_mosaiced_%d.tif", year)))

all_subbasins_subset <- terra::vect(file.path(ia_dir, "hydrobasins_masked_merged_subset.gpkg"))
all_subbasins_subset <- terra::project(x=all_subbasins_subset, y=stack_y)

subbasin_index <- 301
subbasin_s <- all_subbasins_subset[subbasin_index]

# crop covariate stack to subbasin
cov_s <- 
  stack_y |>
  terra::crop(x = _, y = subbasin_s) |>
  terra::mask(x = _, mask = subbasin_s)

test <- rast(file.path(ia_dir, "bart_models/2020/subbasin_301/subbasin_301_backfill.tif"))


my_colours <- c(
  "#bae4b3",  # ID 1: bryoid
  "#bae4b3",  # ID 2: herb
  "#E69F00",  # ID 3: rock
  "#bae4b3",  # ID 4: shrub
  "#7FFF00",  # ID 5: broadleaf (treed, darker shade)
  "#006644",  # ID 6: conifer (treed, darker shade)
  "#9ACD32",  # ID 7: mixed (treed, darker shade)
  "#0072B2"   # ID 8: water
)

plot(cov_s$SCANFI_1km, col=my_colours)
plot(test$SCANFI_1km, col=my_colours, add=TRUE)
