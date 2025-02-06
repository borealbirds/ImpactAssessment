# ---
# title: National Models 5.0 - kriging interpolation for backfilling human footprints
# author: Mannfred Boehm
# created: January 15, 2025
# ---



root1 <- "G:/Shared drives/BAM_NationalModels5/"
root2 <- here()

# ------------------------------
# import BCR shapefile: project, crop, and mask
bcr <- 
  terra::vect(file.path(root1, "Regions", "BAM_BCR_NationalModel.shp")) |> 
  (\(x) terra::subset(x=x, subset=x$subUnit == 61))() |> # trailing () is needed to execute \(x) with the piped object
  terra::project(y="epsg:5072") |> 
  (\(x) terra::crop(x=x, y=terra::ext(x)))() |> #create bounding box defined by BCR61 polygon
  (\(x) terra::mask(x=x, mask=x))() #assign NA to pixels in bounding box but outside the BCR61 polygon


# check reproject/crop/mask
test1 <- terra::vect(file.path(root, "Regions", "BAM_BCR_NationalModel.shp"))
test2 <- terra::subset(test1, test1$subUnit == 61)
test3 <- bcr

plot(test1, col = "lightgray", main = "Original BCRs") #og shapefile
plot(test2, col = "blue", add = TRUE, lwd = 2) #reprojected BCR61
plot(bcr, col = "green", add = TRUE, lwd = 2) #cropped and masked



# ------------------------------
# import Canada Human Footprint (CHF) layer 
chf <- terra::rast(file.path(root2, "data", "raw_data", "cum_threat2020.02.18.tif"))

# crop+mask before reprojecting to save computing time
bcr_buffered <- 
  terra::buffer(bcr, width = 1000) |> 
  terra::project(y=crs(chf))


chf61 <-
  chf |> 
  (\(x) terra::crop(x=x, y=ext(bcr_buffered)))() |> #using the extent of BCR61 plus a buffer to account for different projections
  (\(x) terra::mask(x=x, mask=bcr_buffered))() |> #assign NA to pixels in bounding box but outside the BCR61 polygon
  terra::project(y="epsg:5072") |> 
  (\(x) terra::crop(x=x, y=ext(bcr)))() |>  #now that we're in right projection, do full crop and mask
  (\(x) terra::mask(x=x, mask=bcr))() #assign NA to pixels in bounding box but outside the BCR61 polygon

# explicitly align extents
ext(chf61) <- ext(bcr) 

# mask again to ensure exact polygon boundaries
# now ext(bcr) and ext(chf61) are the same
chf61 <- terra::mask(chf61, bcr)

# saveRDS(chf61, file=here("data/derived_data/rds_files/chf61.rds"))


# ------------------------------
# backfill





