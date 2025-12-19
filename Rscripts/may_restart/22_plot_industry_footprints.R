# ---
# title: Impact Assessment: plot industry footprint
# author: Mannfred Boehm
# created: October 2, 2025
# ---


library(terra)
library(tidyterra)
library(tidyverse)

# set root path
root <- "G:/Shared drives/BAM_NationalModels5"
ia_dir <- file.path(root, "data", "Extras", "sandbox_data", "impactassessment_sandbox")

# BCR boundaries
bam_boundary <- terra::vect(file.path(root, "Regions", "BAM_BCR_NationalModel_Unbuffered.shp")) 

# combined industry polygon
combined_poly <- terra::vect(file.path(ia_dir, "combined_industry_footprint.gpkg"))

# oil and gas pixels
oilgas_no_urban <- terra::rast(file.path(ia_dir, "hirshpearson_oil_gas_masked.tif"))

# mining pixels for year 2020
mines_no_urban <- terra::rast(file.path(ia_dir, paste0("mincan_", "mines_", "2020_", "masked.tif")))


# plotting
# 1) Build a categorical raster: 1 = Oil&Gas, 2 = Mines (2020), 3 = Overlap
ogb <- terra::lapp(oilgas_no_urban, \(v) ifelse(is.finite(v) & v > 0, 1L, 0L))
mib <- terra::lapp(mines_no_urban,   \(v) ifelse(is.finite(v) & v > 0, 1L, 0L))
ind_cat <- ogb + 2L * mib
ind_cat[ind_cat == 0] <- NA
ind_fac <- as.factor(ind_cat)
levels(ind_fac) <- data.frame(
  value = c(1, 2, 3),
  label = c("Oil & Gas", "Mines (2020)", "Overlap")
)

# 2) Plot: pixels + outlines
ggplot() +
  geom_spatraster(data = ind_fac) +
  scale_fill_manual(
    values = c(
      "Oil & Gas"   = "#E69F00",
      "Mines (2020)"= "#CC79A7",
      "Overlap"     = "#009E73"
    ),
    na.value = NA,
    name = "Industry pixels"
  ) +
  #geom_spatvector(data = combined_poly, fill = NA, color = "grey80", linewidth = 0.4) +
  geom_spatvector(data = bam_boundary,  fill = NA, color = "grey25", linewidth = 0.3) +
  coord_sf(crs = crs(bam_boundary)) +
  theme_minimal()
