library(terra)
library(tidyverse)

# set working directory
root <- "G:/Shared drives/BAM_NationalModels5"

# import boundary shapefile
bcr_boundary <- terra::vect(file.path(root, "Regions", "BAM_BCR_NationalModel_Unbuffered.shp"))

# index the raster files ending with "_2020.tif" in the 'stacks' folder.
test_files <- list.files(
  path = file.path(root, "gis", "stacks"),
  pattern = "_2020\\.tif$",
  full.names = TRUE
)

# define an output directory for inspecting plots
output_directory <- file.path(
  root, "data", "Extras", "sandbox_data",
  "impactassessment_sandbox", "inspect_categorical_covariates"
)

# create a function that processes a single raster stack file
inspect_layers <- function(test_file, boundary, output_dir) {
  
  # read in a BCR's raster stack
  rast_stack <- terra::rast(test_file)
  
  # define layers to plot 
  layers_to_plot <- c("MODISLCC_5x5", "MODISLCC_1km", "SCANFI_1km", "VLCE_1km", "NLCD_1km", "ABoVE_1km")
  
  # extract a base name from the file for naming the PNGs
  base_name <- tools::file_path_sans_ext(basename(test_file))
  
  # loop through each layer and plot if it exists
  for (layer_name in layers_to_plot) {
    if (layer_name %in% names(rast_stack)) {
      
      # convert layer to a factor (categorical)
      layer <- terra::as.factor(rast_stack[[layer_name]])
      
      # define the output file name
      output_file <- file.path(output_dir, paste0(base_name, "_", layer_name, ".png"))
      
      # open the PNG device
      png(filename = output_file, width = 990, height = 586)
      
      # plot the layer as a factor and add the boundary
      plot(layer, main = paste(base_name, layer_name))
      lines(boundary)
      
      # close the device to save the file
      dev.off()
      
    } else {
      warning(paste("layer", layer_name, "not found in", test_file))
    } # close ifelse
  } # close loop
}# close function

# generate plots for inspection
lapply(X = test_files, FUN = inspect_layers, boundary = bcr_boundary, output_dir = output_directory)
