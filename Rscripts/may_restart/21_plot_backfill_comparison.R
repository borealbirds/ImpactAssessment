library(terra)
library(tidyterra) 
library(ggplot2)
library(maptiles)

# set root path
root <- "G:/Shared drives/BAM_NationalModels5"
ia_dir <- file.path(root, "data", "Extras", "sandbox_data", "impactassessment_sandbox")

# from "create_pop_summary_table.R"
# cawa_diff_max
# A tibble: 1 × 5
# subbasin backfilled observed abs_diff rel_diff
# <dbl>      <dbl>    <dbl>    <dbl>    <dbl>
#   1 7060245400      3869.    2000.    1869.     93.5

# plot observed raster layer---------------------------
# import all merged subbasins
all_subbasins <- vect(file.path(ia_dir, "hydrobasins_masked_merged.gpkg"))

# Athabasca Oil Sands is in BCR60, 61, and 80
aoi_poly <- as.polygons(ext(c(-114.487, -112.278, 57.289, 57.803)), crs = "EPSG:4326")
aoi_poly <- project(aoi_poly, crs(all_subbasins))
crs(ath_poly)


# find subbasins overlapping with oil sands
aoi_subbasins <- {
  hits <- terra::relate(all_subbasins, aoi_poly, relation = "intersects")
  all_subbasins[rowSums(hits) > 0, ]
}

# OR find subbasins overlapping with max pop change
aoi_subbasins <-  all_subbasins_subset[
  all_subbasins_subset$first_HYBAS_ID == cawa_diff_max$subbasin, ]

# choose some layer to plot
bcr60 <- 
  rast(file.path(root, "gis", "stacks", "can60_2020.tif"))[["SCANFIprcD_1km"]] |> 
  crop(x=_, y=aoi_subbasins) |> 
  mask(x=_, mask=aoi_subbasins)

bcr61 <- 
  rast(file.path(root, "gis", "stacks", "can61_2020.tif"))[["SCANFI_1km"]] |> 
  crop(x=_, y=aoi_subbasins) |> 
  mask(x=_, mask=aoi_subbasins)

bcr80 <- 
  rast(file.path(root, "gis", "stacks", "can80_2020.tif"))[["SCANFI_1km"]] |> 
  crop(x=_, y=aoi_subbasins) |> 
  mask(x=_, mask=aoi_subbasins)

bcr11 <-
  rast(file.path(root, "gis", "stacks", "can11_2020.tif"))[["SCANFI_1km"]] |> 
  crop(x=_, y=aoi_subbasins) |> 
  mask(x=_, mask=aoi_subbasins)


# athabasca
r_mos <- mosaic(mosaic(bcr60, bcr61), bcr80) 

# peace river
r_mos <- mosaic(mosaic(bcr60, bcr61), bcr11)

# prairies
r_mos <- mosaic(mosaic(bcr80, bcr61), bcr11) 

# import industry footprint polygon
combined_poly <- vect(file.path(ia_dir, "combined_industry_footprint.gpkg"))
combined_poly <- terra::intersect(
  terra::project(combined_poly, crs(aoi_subbasins)),
  terra::aggregate(aoi_subbasins)  # dissolve subbasins to avoid dup splits
)

ggplot() +
  geom_spatraster(data = as.factor(r_mos)) +
  scale_fill_viridis_c(na.value = NA, name = "SCANFIprcD_1km") +
  #geom_spatvector(data = sb_ath, fill = NA, linewidth = 0.5) +
  geom_spatvector(data = combined_poly, fill = NA, linewidth = 0.5, colour="white") +
  coord_sf() + theme_minimal()



# compare to backfilled layers---------------------------
# find subbasin indexes for sb_ath
all_subbasins_subset <- vect(file.path(ia_dir, "hydrobasins_masked_merged_subset.gpkg"))
m   <- which(terra::relate(all_subbasins_subset, aoi_subbasins, "intersects"), arr.ind = TRUE)
idx_in_order_of_sb_ath <- m[order(m[,2]), 1] # [1] 154 207 274 275 276 277 279 280 284

sub_ids <- all_subbasins_subset$first_HYBAS_ID[idx_in_order_of_sb_ath]

mosaic_backfilled_onevar <- function(sub_ids, year, layer_nm,
                                     subbasins_path, ia_dir, template = NULL) {
  sub_sv <- terra::vect(subbasins_path)
  idx_chr <- match(sub_ids, sub_sv$first_HYBAS_ID)
  
  ptab <- tibble::tibble(
    HYBAS_ID = as.character(sub_ids),
    subbasin_index = idx_chr,
    path = file.path(
      ia_dir, "xgboost_models", paste0("year=", year),
      paste0("subbasin=", idx_chr),
      paste0("backfilled_stack_subbasin-", sprintf("%03d", idx_chr), ".tif"))
  ) |>
    dplyr::filter(!is.na(subbasin_index), file.exists(path))
  
  if (nrow(ptab) == 0) return(NULL)
  
  rlist <- lapply(ptab$path, \(p) terra::rast(p)[[layer_nm]])
  
  if (!is.null(template)) {
    rlist <- lapply(rlist, \(r) if (!terra::compareGeom(r, template, stopOnError = FALSE))
      terra::resample(r, template, method = "near") else r)
  } else {
    ref <- rlist[[1]]
    rlist <- lapply(rlist, \(r) if (!terra::compareGeom(r, ref, stopOnError = FALSE))
      terra::resample(r, ref, method = "near") else r)
  }
  
  out <- Reduce(terra::cover, rlist)
  names(out) <- layer_nm
  out
}

layer_nm <- "SCANFI_1km"

# observed (already built)
# r_mos <- mosaic(r60, r61)      # from the previous step
obs_ath <- terra::mask(terra::crop(r_mos[[layer_nm]], aoi_subbasins), aoi_subbasins)

# backfilled mosaic
r_bf <- mosaic_backfilled_onevar(
  sub_ids        = sub_ids,
  year           = 2020,
  layer_nm       = layer_nm,
  subbasins_path = file.path(ia_dir, "hydrobasins_masked_merged_subset.gpkg"),
  ia_dir         = ia_dir,
  template       = r_mos[[layer_nm]]        # ensures same grid as observed
)

bf_ath <- terra::mask(terra::crop(r_bf, aoi_subbasins), aoi_subbasins)
bf_ath <- terra::cover(bf_ath, obs_ath)

ggplot() +
  geom_spatraster(data = bf_ath) +
  scale_fill_viridis_c(na.value = NA, name = "SCANFIprcD_1km") +
  # geom_spatvector(data = sb_ath, fill = NA, linewidth = 0.5) +
  geom_spatvector(data = combined_poly, fill = NA, linewidth = 0.5, colour="white") +
  coord_sf() + theme_minimal()




# satellite basemap over the Athabasca subbasins --------------
sat <- maptiles::get_tiles(
  x        = sb_ath,                # your AOI (SpatVector)
  provider = "Esri.WorldImagery",      # high-res imagery
  zoom     = 10,                       # bump to 11–12 if you need more detail
  crop     = TRUE
)

sat <- project(sat, r_mos, method = "near") |> 
      crop(x=_, y=sb_ath) |> 
      mask(x=_, mask=sb_ath) 

ggplot() +
  geom_spatraster_rgb(data = sat) +
  #geom_spatvector(data = sb_ath, fill = NA, linewidth = 0.5) +
  geom_spatvector(data = combined_poly, fill = NA, linewidth = 0.5, colour="white") +
  coord_sf() + theme_minimal()



# CAWA population visualization-----------------------------------
dir_ca <- file.path(ia_dir, "density_predictions", "CAWA", "2020")

# precedence mosaic for CAWA observed/backfilled tiles
mosaic_with_precedence <- function(dir_ca, scenario = c("observed","backfilled"),
                                   template, precedence_codes = NULL) {
  scenario <- match.arg(scenario)
  
  # pick files for scenario
  paths <- list.files(dir_ca, pattern = "CAWA.*2020.*\\.tif$", full.names = TRUE, ignore.case = TRUE)
  tag   <- if (scenario == "observed") "observ" else "backfill"
  paths <- paths[grepl(tag, basename(paths), ignore.case = TRUE)]
  stopifnot(length(paths) > 0)
  
  # extract tile code like "can60" -> 60
  codes <- as.integer(sub(".*can(\\d{2}).*", "\\1", basename(paths)))
  ord <- if (is.null(precedence_codes)) order(codes) else order(match(codes, precedence_codes), codes)
  paths <- paths[ord]
  
  # align each tile to the template grid (CRS + resolution + origin)
  align_to_template <- function(r, tmpl) {
    if (as.character(crs(r)) != as.character(crs(tmpl))) r <- project(r, tmpl, method = "bilinear")
    if (!compareGeom(r, tmpl, stopOnError = FALSE))       r <- resample(r, tmpl, method = "bilinear")
    r
  }
  rs <- lapply(paths, \(p) align_to_template(rast(p), template))
  
  # precedence mosaic: first covers later
  out <- Reduce(cover, rs)
  names(out) <- paste0("CAWA_", scenario)
  out
}

cawa <- 
  terra::rast(file.path(root, "output", "08_mosaics_can", "CAWA_2020.tif")) |> 
  crop(x=_, y=sb_ath) |> 
  mask(x=_, mask=sb_ath)

cawa_mean <- mean(cawa, na.rm = TRUE)

obs <- mosaic_with_precedence(dir_ca, "observed",   template = cawa_mean, precedence_codes = c(60, 61))
bf  <- mosaic_with_precedence(dir_ca, "backfilled", template = cawa_mean, precedence_codes = c(60, 61))

obs_ath <- terra::mask(terra::crop(obs, sb_ath), sb_ath) |> terra::trim()
bf_ath  <- terra::mask(terra::crop(bf,  sb_ath), sb_ath) |> terra::trim()


ggplot() +
  geom_spatraster(data = bf_ath) +
  scale_fill_viridis_c(na.value = NA, name = paste(spp, "density")) +
  geom_spatvector(data = combined_poly, fill = NA, linewidth = 0.5, colour="white") +
  coord_sf() + theme_minimal()

  
