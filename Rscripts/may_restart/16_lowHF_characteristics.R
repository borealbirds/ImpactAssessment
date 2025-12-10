# ---
# title: Impact Assessment: low HF is not a proxy for pre-industrialization
# author: Mannfred Boehm
# created: December 9, 2025
# ---

library(terra)
library(tidyterra)
library(tidyverse)

root <- "G:/Shared drives/BAM_NationalModels5"
ia_dir <- file.path(root, "data", "Extras", "sandbox_data", "impactassessment_sandbox")
year <- 2020

# ---------------------------------------------------
# import pre-mosaiced covariate stack for year_y
stack_y <- terra::rast(file.path(ia_dir, sprintf("covariates_mosaiced_%d.tif", year)))

# import subbasin boundaries and project to current stack
all_subbasins_subset <- terra::vect(file.path(ia_dir, "hydrobasins_masked_merged_subset.gpkg"))
all_subbasins_subset <- terra::project(x=all_subbasins_subset, y=stack_y)

# subset to a random sample of subbasins
set.seed(123)
subbasin_index <- sample(1:length(all_subbasins_subset), size = 20)
subbasin_s <- all_subbasins_subset[subbasin_index]

# crop covariate stack to subbasin
cov_s <- 
  stack_y |>
  terra::crop(x = _, y = subbasin_s) |>
  terra::mask(x = _, mask = subbasin_s)


# ---------------------------------------------------
# visualize randomly selected subbasins

# get BCR boundaries to orient our subbasins in space
bam_boundary <- terra::vect(file.path(root, "Regions", "BAM_BCR_NationalModel_UnBuffered.shp"))
bam_boundary <- bam_boundary[bam_boundary$subUnit != 23, ]

# plot (exported at 1000 x 751)
ggplot() +
  
  # BCR and basin outlines
  geom_spatvector(data = bam_boundary, fill = NA, colour = "grey") +
  geom_spatvector(data = subbasin_s, fill = NA, color = "#D55E00", linewidth = 0.3) +
  
  coord_sf(crs = crs(subbasin_s)) +
  theme_minimal()



# ---------------------------------------------------
# analyse low HF traits in randomly selected subbasins

# convert stacks from random subbasins to dataframe
cov_df <- terra::as.data.frame(cov_s)
saveRDS(cov_df, file=file.path(getwd(), "data/derived_data/rds_files/covariate_dataframe_random_subbasins_2020.rds"))


# import low hf layer and project to current stack
lowhf_mask <- terra::rast(file.path(ia_dir, "hirshpearson", "CanHF_1km_lessthan1.tif"))
lowhf_mask <- terra::project(x=lowhf_mask, y=cov_s, method = "near")
cov_lowhf  <- terra::mask(cov_s, lowhf_mask) 
cov_lowhf_df <- terra::as.data.frame(cov_lowhf)
saveRDS(cov_lowhf_df, file=file.path(getwd(), "data/derived_data/rds_files/lowhf_covariate_df_random_subbasins_2020.rds"))

# import high hf layer and project to current stack
highhf_mask <- terra::rast(file.path(ia_dir, "hirshpearson", "CanHF_1km_morethan1.tif"))
highhf_mask <- terra::project(x=highhf_mask, y=cov_s, method = "near")
cov_highhf  <- terra::mask(cov_s, highhf_mask) 
cov_highhf_df <- terra::as.data.frame(cov_highhf)

saveRDS(cov_highhf, file=file.path(getwd(), "data/derived_data/rds_files/highhf_covariate_df_random_subbasins_2020.rds"))


