# ---
# title: Impact Assessment: comparing MinCan to Hirsh-Pearson 
# author: Mannfred Boehm
# created: September 15, 2025
# ---

# will use MinCan 2020 ****IF**** Hirsh-Pearson is a subset of MinCan
library(terra)
library(tidyverse)

root <- "G:/Shared drives/BAM_NationalModels5"
ia_dir <- file.path(root, "data", "Extras", "sandbox_data", "impactassessment_sandbox")

mincan <- terra::rast(file.path(ia_dir, "mincan_mines_2020_masked.tif"))

hp_mines <- 
  terra::rast(file.path(ia_dir, "hirshpearson_mines.tif")) |> 
  terra::project(x=_, y=mincan) |> 
  terra::crop(x=_, y=mincan) |> 
  terra::mask(x=_,)

terra::writeRaster(hp_mines, file.path(ia_dir, "hirshpearson_mines_masked.tif"))


# find pixels where hp_mines has mines but mincan does not
diff <- hp_mines_aligned == 1 & mincan != 1

# check if any violations exist
is_subset <- global(diff, "max", na.rm = TRUE)[[1]] == 0

if (is_subset) {
  message("hp_mines is a subset of mincan")
} else {
  message("hp_mines has mines outside mincan")
}



