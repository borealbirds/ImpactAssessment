# Elevation on the national mosaic's 1 km grid, over BC and its surroundings.
# AWS Terrain Tiles via elevatr at zoom 7 (~0.7 km at 55 N), area-averaged to 1 km.
suppressMessages({library(terra); library(elevatr); library(sf)})
out <- file.path("data", "raw_data")
st <- rnaturalearth::ne_states(country = "Canada", returnclass = "sf")
bc <- st[st$name == "British Columbia", ]
st_write(bc[, "name"], file.path(out, "bc_boundary.gpkg"), delete_dsn = TRUE, quiet = TRUE)

bb <- st_bbox(bc)
print(bb)
loc <- data.frame(x = c(bb["xmin"] - 2, bb["xmax"] + 2), y = c(bb["ymin"] - 0.5, bb["ymax"] + 1.5))
e <- get_elev_raster(locations = loc, prj = "EPSG:4326", z = 7, src = "aws", clip = "bbox", verbose = FALSE)
e <- rast(e)
print(e)

mos  <- rast("data/raw_data/covariates_mosaiced/covariates_mosaiced_2020.tif")[[1]]
bcv  <- project(vect(bc), crs(mos))
tmpl <- crop(mos, buffer(bcv, 300000), snap = "out")
dem  <- project(e, tmpl, method = "average", threads = FALSE)
names(dem) <- "elev_m"
writeRaster(dem, file.path(out, "dem_bc_1km.tif"), overwrite = TRUE, gdal = "COMPRESS=DEFLATE")
print(dem)
print(global(dem, c("min", "max", "mean"), na.rm = TRUE))
cat("DEM_DONE\n")
