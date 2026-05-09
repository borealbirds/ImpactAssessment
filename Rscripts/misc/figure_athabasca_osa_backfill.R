# ---
# title: 2x3 figure: Athabasca Oil Sands — context | satellite | CanHF / observed | backfilled | uncertainty
# author: Mannfred Boehm
# ---
#

library(terra)
library(tidyterra)
library(ggplot2)
library(maptiles)
library(patchwork)

cc    <- FALSE
local <- TRUE

if (cc)            { ia_dir <- "/home/mannfred/scratch/impact_assessment" }
if (!cc && local)  { ia_dir <- getwd() }
if (!cc && !local) { ia_dir <- file.path("G:/Shared drives/BAM_NationalModels5", "data", "Extras", "sandbox_data", "impactassessment_sandbox") }

year <- 2020

SUBBASIN_IDS <- c(356, 360, 376)

# --- AOI: 100 x 100 km centred on 57°09'15.2"N 111°32'05.0"W ---
center_wgs84 <- matrix(c(-111.534722, 57.154222), nrow = 1,
                       dimnames = list(NULL, c("x", "y")))
center_pt    <- project(vect(center_wgs84, crs = "EPSG:4326"), "EPSG:5072")
cx   <- crds(center_pt)[, "x"]
cy   <- crds(center_pt)[, "y"] - 8000  # 8 km south
half <- 50000  # metres

aoi      <- ext(cx - half, cx + half, cy - half, cy + half)
aoi_poly <- as.polygons(aoi, crs = "EPSG:5072")

# --- Load covariate stack ---
cov_stack <- rast(file.path(ia_dir, "data", "raw_data", "covariates_mosaiced",
                             sprintf("covariates_mosaiced_%d.tif", year)))

# --- Load and merge backfill draws across subbasins ---
# draws stored in log1p space; expm1() back-transforms before averaging
backfill_stacks <- lapply(SUBBASIN_IDS, function(id) {
  rast(file.path(ia_dir, "data", "derived_data", "bart_models",
                 as.character(year),
                 sprintf("subbasin_%d", id),
                 sprintf("subbasin_%d_backfill.tif", id)))
})

draw_layers <- grep("^SCANFIheight_1km_draw_", names(backfill_stacks[[1]]), value = TRUE)

draw_merged <- Reduce(
  terra::merge,
  lapply(backfill_stacks, function(s) s[[draw_layers]])
)

biomass_bf <- app(draw_merged, function(x) mean(expm1(x))) |>
  project(y = cov_stack, method = "bilinear") |>
  crop(y = aoi)
biomass_sd <- app(draw_merged, function(x) sd(expm1(x))) |>
  project(y = cov_stack, method = "bilinear") |>
  crop(y = aoi)

# --- Backfill-eligible pixels: highhf_mask intersected with union of all subbasins ---
highhf_mask   <- rast(file.path(ia_dir, "data", "raw_data", "hirshpearson", "CanHF_1km_morethan1.tif"))
all_subbasins <- vect(file.path(ia_dir, "data", "raw_data", "hydrobasins_masked_merged_subset.gpkg"))

subbasins_union <- aggregate(project(all_subbasins[SUBBASIN_IDS], cov_stack))

highhf_proj <- project(highhf_mask, cov_stack, method = "near")
highhf_si   <- crop(mask(highhf_proj, subbasins_union), aoi)  # 1 on eligible pixels, NA elsewhere

# Mask SD to backfill-eligible pixels only
biomass_sd <- mask(biomass_sd, highhf_si)

# --- Top-left: Regional context satellite (~2000 km × 2000 km, southern extent ~49°N) ---
south_pt     <- project(
  vect(matrix(c(-111.65, 49.0), nrow = 1, dimnames = list(NULL, c("x", "y"))), crs = "EPSG:4326"),
  "EPSG:5072"
)
sy_ctx       <- crds(south_pt)[, "y"]
aoi_ctx      <- ext(cx - 1e6, cx + 1e6, sy_ctx, sy_ctx + 2e6)
aoi_ctx_poly <- as.polygons(aoi_ctx, crs = "EPSG:5072")

mine_pt <- data.frame(x = cx, y = cy)

sat_ctx <- get_tiles(
  x        = project(aoi_ctx_poly, "EPSG:4326"),
  provider = "Esri.WorldImagery",
  zoom     = 4,
  crop     = TRUE
)
sat_ctx <- project(sat_ctx, "EPSG:5072", method = "near") |> crop(y = aoi_ctx)

sb_x0_ctx  <- cx - 900000
sb_x1_ctx  <- sb_x0_ctx + 500000
sb_y_ctx   <- sy_ctx + 275000
tick_h_ctx <- 25000

p_ctx <- ggplot() +
  geom_spatraster_rgb(data = sat_ctx) +
  geom_point(data = mine_pt, aes(x = x, y = y),
             shape = 0, fill = "red", colour = "white", size = 3.5, stroke = 0.9) +
  annotate("segment",
           x = sb_x0_ctx, xend = sb_x1_ctx, y = sb_y_ctx, yend = sb_y_ctx,
           colour = "white", linewidth = 0.8) +
  annotate("segment",
           x = sb_x0_ctx, xend = sb_x0_ctx, y = sb_y_ctx - tick_h_ctx, yend = sb_y_ctx + tick_h_ctx,
           colour = "white", linewidth = 0.8) +
  annotate("segment",
           x = sb_x1_ctx, xend = sb_x1_ctx, y = sb_y_ctx - tick_h_ctx, yend = sb_y_ctx + tick_h_ctx,
           colour = "white", linewidth = 0.8) +
  annotate("text",
           x = (sb_x0_ctx + sb_x1_ctx) / 2, y = sb_y_ctx - 56000,
           label = "500 km", colour = "white", size = 6.4, vjust = 1) +
  labs(title = "Athabasca Oil Sands") +
  coord_sf(crs = "EPSG:5072") +
  theme_void() +
  theme(plot.title = element_text(size = 20, hjust = 0.5))

sat <- get_tiles(
  x        = project(aoi_poly, "EPSG:4326"),
  provider = "Esri.WorldImagery",
  zoom     = 11,
  crop     = TRUE
)
sat <- project(sat, "EPSG:5072", method = "near") |> crop(y = aoi)

canhf <- cov_stack[["CanHF_1km"]] |> crop(y = aoi)

biomass_obs       <- cov_stack[["SCANFIheight_1km"]] |> crop(y = aoi)
biomass_composite <- cover(biomass_bf, biomass_obs)  # backfilled on eligible pixels, observed elsewhere

val_range_biomass <- range(c(values(biomass_obs, na.rm = TRUE), values(biomass_composite, na.rm = TRUE), values(biomass_sd, na.rm = TRUE)))

# --- Build panels ---
# Scale bar geometry in EPSG:5072 coordinates (25 km = 25000 m)
sb_x0  <- cx - 47000
sb_x1  <- sb_x0 + 25000
sb_y   <- cy - 36250  # 13.75% from panel bottom, matching context panel scale bar height
tick_h <- 1500

p_sat <- ggplot() +
  geom_spatraster_rgb(data = sat) +
  annotate("segment",
           x = sb_x0, xend = sb_x1, y = sb_y, yend = sb_y,
           colour = "white", linewidth = 0.8) +
  annotate("segment",
           x = sb_x0, xend = sb_x0, y = sb_y - tick_h, yend = sb_y + tick_h,
           colour = "white", linewidth = 0.8) +
  annotate("segment",
           x = sb_x1, xend = sb_x1, y = sb_y - tick_h, yend = sb_y + tick_h,
           colour = "white", linewidth = 0.8) +
  annotate("text",
           x = (sb_x0 + sb_x1) / 2, y = sb_y - 2500,
           label = "25 km", colour = "white", size = 6.4, vjust = 1) +
  labs(title = "Athabasca Oil Sands") +
  coord_sf(crs = "EPSG:5072") +
  theme_void() +
  theme(plot.title = element_text(size = 20, hjust = 0.5))

p_hf <- ggplot() +
  geom_spatraster(data = canhf) +
  scale_fill_viridis_c(
    name      = "Human\nfootprint",
    na.value  = NA,
    option    = "A",
    direction = -1
  ) +
  labs(title = "Human Footprint Score, 2020") +
  coord_sf(crs = "EPSG:5072") +
  theme_void() +
  theme(plot.title      = element_text(size = 20, hjust = 0.5),
        legend.position = "right",
        legend.title    = element_text(size = 14),
        legend.text     = element_text(size = 12))

p_obs <- ggplot() +
  geom_spatraster(data = biomass_obs) +
  scale_fill_viridis_c(
    limits   = val_range_biomass,
    name     = "Tree height (m)",
    na.value = NA,
    option   = "D",
    guide    = "none"
  ) +
  labs(title = "Observed tree height, 2020") +
  coord_sf(crs = "EPSG:5072") +
  theme_void() +
  theme(plot.title = element_text(size = 20, hjust = 0.5))

p_bf <- ggplot() +
  geom_spatraster(data = biomass_composite) +
  scale_fill_viridis_c(
    limits   = val_range_biomass,
    name     = "Tree height (m)",
    na.value = NA,
    option   = "D",
    guide    = "none"
  ) +
  labs(title = "Backfilled tree height, 2020") +
  coord_sf(crs = "EPSG:5072") +
  theme_void() +
  theme(plot.title = element_text(size = 20, hjust = 0.5))

p_uncert <- ggplot() +
  geom_spatraster(data = biomass_sd) +
  scale_fill_viridis_c(
    limits   = val_range_biomass,
    name     = "Tree height (m)",
    na.value = NA,
    option   = "D"
  ) +
  labs(title = "Backfill uncertainty") +
  coord_sf(crs = "EPSG:5072") +
  theme_void() +
  theme(plot.title      = element_text(size = 20, hjust = 0.5),
        legend.position = "right",
        legend.title    = element_text(size = 14),
        legend.text     = element_text(size = 12))

# --- Assemble 2x3 and save ---
fig <- (p_ctx | p_sat | p_hf) / (p_obs | p_bf | p_uncert)

ggsave(
  filename = file.path(ia_dir, "output_figures", "figure_athabasca_osa_backfill.png"),
  plot     = fig,
  width    = 15,
  height   = 10,
  dpi      = 300
)
