# ---
# title: Impact Assessment: comparing MinCan to Hirsh-Pearson 
# author: Mannfred Boehm
# created: September 15, 2025
# ---

library(terra)
library(tidyverse)

root <- "G:/Shared drives/BAM_NationalModels5"
ia_dir <- file.path(root, "data", "Extras", "sandbox_data", "impactassessment_sandbox")


hp_mines <- terra::rast(file.path(ia_dir, "hirshpearson_mines.tif"))


