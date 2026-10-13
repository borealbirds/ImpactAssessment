# ---
# title: Impact Assessment: map the population change due to industrial footprint per subbasin
# author: Mannfred Boehm
# ---
# Reads 14B's subbasin output and maps, per subbasin, what industry did to each species:
# change = observed (2020) - counterfactual (no industry) = -total_HF_impact, so red is a
# decline due to industry and blue an increase. Subbasins that straddle BCRs are summed
# across them by HYBAS_ID. Grey: no subbasin with footprint, outside the species' modelled
# BCRs, or withheld by BAM (14B drops them, e.g. CAWA can40).
#   (1) absolute change, birds per km2 of subbasin, one figure per species on its own scale
#   (2) percent change relative to the no-industry population, CAWA and OVEN side by side
# Display CRS is EPSG:3978 (display only).
# A sensitivity run is mapped by setting its switch, as for 14B (12E sensitivity_run()):
#   DENS_MIN=0.005 Rscript Rscripts/18_map_population_change.R
# reads sector_effects_dmin0.005/ and writes the same figures with a _dmin0.005 suffix. The
# absolute maps keep the production run's colour scale, so the two can be compared directly.
# ---

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(ggplot2)
})


# set paths ------------------------------------------------------

cc    <- FALSE
local <- TRUE

if (cc)            { ia_dir <- "/home/mannfred/scratch/impact_assessment" }
if (!cc && local)  { ia_dir <- getwd() }
if (!cc && !local) { ia_dir <- file.path("G:/Shared drives/BAM_NationalModels5", "data", "Extras", "sandbox_data", "impactassessment_sandbox") }

source(file.path(ia_dir, "Rscripts", "12E_shapley_utils.R"))
sens <- sensitivity_run()

shap_path  <- file.path(ia_dir, "data", "derived_data", paste0("sector_effects", sens$suffix), "shapley_subbasin.csv")
prod_path  <- file.path(ia_dir, "data", "derived_data", "sector_effects", "shapley_subbasin.csv")
basin_path <- file.path(ia_dir, "data", "raw_data", "hydrobasins_masked_merged_subset.gpkg")
bcr_path   <- file.path(ia_dir, "data", "raw_data", "Regions", "BAM_BCR_NationalModel_Unbuffered.shp")
fig_dir    <- file.path(ia_dir, "output_figures")
dir.create(fig_dir, showWarnings = FALSE)


# subbasin change ------------------------------------------------

# one row per species x bcr x subbasin (total_HF_impact repeats across the 8 sector rows),
# then sum straddling subbasins across BCRs
subbasin_change <- function(path) {
  read.csv(path) |>
    distinct(species, bcr, HYBAS_ID, obs_population, total_HF_impact) |>
    group_by(species, HYBAS_ID) |>
    summarise(obs = sum(obs_population), impact = sum(total_HF_impact), .groups = "drop") |>
    mutate(cf     = obs + impact,
           change = -impact,                                 # obs - cf: what industry did
           pct    = ifelse(cf > 0, 100 * change / cf, NA_real_))
}
sub <- subbasin_change(shap_path)

bcr_can <- st_read(bcr_path, quiet = TRUE) |>
  filter(country == "can")
bcr <- st_transform(bcr_can, 3978)

# Keep only the Canadian part of each subbasin: the birds are counted inside the Canadian BCRs
# only (12C's in-BCR weight), so the area that turns them into birds/km2 must be Canadian too.
canada <- bcr_can |> st_transform(5072) |> st_make_valid() |> st_union()
basins <- st_read(basin_path, quiet = TRUE) |>
  select(HYBAS_ID = first_HYBAS_ID) |>
  st_make_valid() |>
  st_intersection(canada) |>
  mutate(area_km2 = as.numeric(st_area(geom)) / 1e6)

m <- basins |>
  inner_join(sub, by = "HYBAS_ID") |>
  mutate(dens = change / area_km2) |>
  st_transform(3978)

lab <- c(CAWA = "Canada Warbler", OVEN = "Ovenbird", OSFL = "Olive-sided Flycatcher",
         GRSP = "Grasshopper Sparrow", LEYE = "Lesser Yellowlegs", BOBO = "Bobolink")
miss <- setdiff(unique(sub$species), names(lab))
if (length(miss) > 0L) stop("add common names to `lab` for: ", paste(miss, collapse = ", "))
tot <- sub |> group_by(species) |> summarise(o = sum(obs), c = sum(cf))
# the % is 14C's change_pct: what industry did (obs - cf) as % of the no-industry population
headline <- setNames(sprintf("%s\n%.2fM birds in 2020 vs %.2fM without industry (%+.1f%%)",
                             lab[tot$species], tot$o / 1e6, tot$c / 1e6,
                             100 * (tot$o - tot$c) / tot$c), tot$species)


# map helpers ----------------------------------------------------

red <- "#c0392b"; blue <- "#2c6fbb"; mid <- "#f0efec"; land <- "#dcdcd8"

# symmetric scale capped at +-lim; the end labels read <= / >= because values beyond are clamped
capped_labels <- function(lim, fmt) {
  function(x) ifelse(abs(x) >= lim,
                     paste0(ifelse(x < 0, "\u2264", "\u2265"), sprintf(fmt, x)),
                     sprintf(fmt, x))
}

# round x up to 1, 2, 2.5 or 5 x 10^k
nice_ceiling <- function(x) {
  k <- 10^floor(log10(x))
  k * c(1, 2, 2.5, 5, 10)[which(x / k <= c(1, 2, 2.5, 5, 10))[1]]
}

change_map <- function(data, fill_var, lim, legend, labels_fun) {
  ggplot() +
    geom_sf(data = bcr, fill = land, colour = "white", linewidth = 0.3) +
    geom_sf(data = data, aes(fill = pmax(pmin({{ fill_var }}, lim), -lim)), colour = NA) +
    geom_sf(data = bcr, fill = NA, colour = "white", linewidth = 0.3) +
    scale_fill_gradient2(low = red, mid = mid, high = blue, midpoint = 0,
                         limits = c(-lim, lim), breaks = c(-1, -0.5, 0, 0.5, 1) * lim,
                         labels = labels_fun,
                         na.value = land, name = legend) +
    coord_sf(datum = NA) +
    theme_void(base_size = 11) +
    theme(legend.position = "bottom", legend.key.width = unit(2.2, "cm"),
          legend.key.height = unit(0.35, "cm"), legend.title.position = "top",
          legend.title = element_text(hjust = 0.5),
          plot.title = element_text(size = 11, hjust = 0.5, margin = margin(b = 4)),
          strip.text = element_text(size = 11, margin = margin(b = 4)),
          plot.background = element_rect(fill = "white", colour = NA),
          plot.margin = margin(12, 12, 12, 12))
}


# (1) absolute change, one figure per species ---------------------

# cap each species at the 98th percentile of |change| in the PRODUCTION run, so a few small,
# dense subbasins don't wash out the rest and a sensitivity run shares production's scale
lims <- subbasin_change(prod_path) |>
  inner_join(st_drop_geometry(basins), by = "HYBAS_ID") |>
  group_by(species) |>
  summarise(lim = nice_ceiling(quantile(abs(change / area_km2), 0.98, na.rm = TRUE)))
lims <- setNames(lims$lim, lims$species)

for (sp in sort(unique(m$species))) {
  d <- filter(m, species == sp)
  lim <- lims[[sp]]
  fmt <- if (lim <= 1) "%+.2f" else if (lim <= 10) "%+.1f" else "%+.0f"
  cat(sprintf("%s: change density range %.2f to %.2f birds/km2, scale capped at +-%g\n",
              sp, min(d$dens), max(d$dens), lim))
  p <- change_map(d, dens, lim, "Population change (birds per km\u00b2)", capped_labels(lim, fmt)) +
    labs(title = headline[[sp]])
  ggsave(file.path(fig_dir, sprintf("change_map_abs_%s%s.png", sp, sens$suffix)), p,
         width = 7.5, height = 6.6, dpi = 200)
}


# (2) percent change, both species --------------------------------

p <- change_map(filter(m, species %in% c("CAWA", "OVEN")), pct, 50, "Percent population change",
                capped_labels(50, "%+.0f%%")) +
  facet_wrap(~species, ncol = 2, labeller = as_labeller(headline))
ggsave(file.path(fig_dir, paste0("change_map_pct_CAWA_OVEN", sens$suffix, ".png")), p,
       width = 12, height = 6.2, dpi = 200)


# (3) pixel-level change on a 10 km grid --------------------------
# 12F writes, per species x BCR, the mean observed and no-industry density at every footprint
# pixel the tables sum (birds/km2, weighted; density_tables{suffix}/pixel_change/). Off the
# footprint nothing changes. Each 1 km pixel's change (obs - cf, in birds) goes to the 10 km cell
# that holds it (every BCR grid has origin 0,0 at 1 km, so cells nest exactly), and each cell's
# sum is divided by its area inside the species' modelled BCRs. Colour then follows the
# footprint, not subbasin outlines. Cells are drawn whole; everything outside Canada is painted
# over, and BCRs the species is not modelled in (or BAM withholds) are drawn grey on top.

cell_m  <- 10000
px_dir  <- function(suffix) file.path(ia_dir, "data", "derived_data", paste0("density_tables", suffix), "pixel_change")
bcr_tab <- read.csv(sub("shapley_subbasin", "shapley_bcr", shap_path)) |> distinct(species, bcr)  # after the release filter
can_3978 <- st_transform(canada, 3978)

grid_change <- function(sp, suffix) {
  bcrs  <- bcr_tab$bcr[bcr_tab$species == sp]
  files <- file.path(px_dir(suffix), paste0(sp, "_2020_", bcrs, ".tif"))
  # a BCR whose all-sector coalition has no footprint pixel writes no raster: it has no change
  cat(sprintf("%s%s: pixel rasters for %d of %d modelled BCRs\n", sp, suffix, sum(file.exists(files)), length(files)))
  parts <- lapply(files[file.exists(files)], function(f) {
    r  <- terra::rast(f)
    v  <- terra::values(r, mat = TRUE)
    ok <- which(!is.na(v[, 1]))
    xy <- terra::xyFromCell(r, ok)
    data.frame(gx = floor(xy[, 1] / cell_m), gy = floor(xy[, 2] / cell_m), change = v[ok, 1] - v[ok, 2])
  })
  if (length(parts) == 0L) stop("no pixel_change rasters for ", sp, " in ", px_dir(suffix))
  bind_rows(parts) |> group_by(gx, gy) |> summarise(change = sum(change), .groups = "drop")
}

grid_density <- function(sp, cells) {
  modelled <- bcr_can |>
    filter(paste0(country, subUnit) %in% bcr_tab$bcr[bcr_tab$species == sp]) |>
    st_transform(5072) |> st_make_valid() |> st_union()
  e    <- terra::ext(terra::vect(modelled))
  tmpl <- terra::rast(terra::ext(floor(e$xmin / cell_m) * cell_m, ceiling(e$xmax / cell_m) * cell_m,
                                 floor(e$ymin / cell_m) * cell_m, ceiling(e$ymax / cell_m) * cell_m),
                      resolution = cell_m, crs = "EPSG:5072")
  # km2 modelled = 1 km pixels touching the modelled BCRs, counted as 12C's in-BCR weight counts
  # them (touches = TRUE), so every pixel 12F predicts lies in a cell with area. (cover = TRUE
  # estimates from sub-cell points and gave 0 to cells a BCR edge barely clips: BOBO lost 4 birds.)
  fine <- terra::rast(terra::ext(tmpl), resolution = 1000, crs = "EPSG:5072")
  area <- terra::rasterize(terra::vect(modelled), fine, touches = TRUE, background = 0) |>
    terra::aggregate(cell_m / 1000, fun = "sum")
  area <- terra::ifel(area == 0, NA, area)
  idx  <- terra::cellFromXY(tmpl, cbind((cells$gx + 0.5) * cell_m, (cells$gy + 0.5) * cell_m))
  keep <- !is.na(idx) & !is.na(terra::values(area, mat = FALSE)[idx])
  if (sum(abs(cells$change[!keep])) > 1e-6 * sum(abs(cells$change)))
    stop(sp, ": ", round(sum(abs(cells$change[!keep]))), " birds of change fall outside the modelled BCRs")
  chg <- terra::rast(tmpl); terra::values(chg) <- 0
  chg[idx[keep]] <- cells$change[keep]
  dens <- chg / area
  names(dens) <- "dens"
  list(dens = dens, total = sum(cells$change[keep]))
}

grid_map <- function(sp, g, lim, legend, labels_fun) {
  modelled <- bcr$subUnit %in% sub("^can", "", bcr_tab$bcr[bcr_tab$species == sp])
  # equal-valued neighbours (mostly the unchanged 0s) are merged, and each cell is outlined in
  # its own colour, or the squares' anti-aliased edges show as seams
  cells    <- st_as_sf(terra::as.polygons(g$dens, aggregate = TRUE, na.rm = TRUE)) |> st_transform(3978)
  bb       <- st_bbox(bcr)
  outside  <- st_difference(st_as_sfc(bb) |> st_buffer(100000), st_union(st_make_valid(can_3978)))
  ggplot() +
    geom_sf(data = bcr, fill = land, colour = NA) +
    geom_sf(data = cells, aes(fill = pmax(pmin(dens, lim), -lim), colour = pmax(pmin(dens, lim), -lim)),
            linewidth = 0.1) +
    geom_sf(data = bcr[!modelled, ], fill = land, colour = NA) +
    geom_sf(data = outside, fill = "white", colour = NA) +
    geom_sf(data = bcr, fill = NA, colour = "white", linewidth = 0.3) +
    scale_fill_gradient2(low = red, mid = mid, high = blue, midpoint = 0,
                         limits = c(-lim, lim), breaks = c(-1, -0.5, 0, 0.5, 1) * lim,
                         labels = labels_fun, na.value = land, name = legend,
                         aesthetics = c("fill", "colour")) +
    coord_sf(datum = NA, xlim = bb[c("xmin", "xmax")], ylim = bb[c("ymin", "ymax")]) +
    theme_void(base_size = 11) +
    theme(legend.position = "bottom", legend.key.width = unit(2.2, "cm"),
          legend.key.height = unit(0.35, "cm"), legend.title.position = "top",
          legend.title = element_text(hjust = 0.5),
          plot.title = element_text(size = 11, hjust = 0.5, margin = margin(b = 4)),
          plot.background = element_rect(fill = "white", colour = NA),
          plot.margin = margin(12, 12, 12, 12))
}

nat_tab <- read.csv(sub("shapley_subbasin", "shapley_national", shap_path)) |> distinct(species, total_HF_impact)

for (sp in sort(unique(bcr_tab$species))) {
  g <- grid_density(sp, grid_change(sp, sens$suffix))
  # the cells must hold the whole national change in 14B. Not `tot`: 14B rounds every subbasin
  # to a whole bird (LEYE's ~740 rows drift 4 birds); the national row is rounded once (+-0.5),
  # and the float32 rasters add far less than that
  nat <- -nat_tab$total_HF_impact[nat_tab$species == sp]
  cat(sprintf("%s: grid change %.1f birds vs 14B national %.0f\n", sp, g$total, nat))
  if (abs(g$total - nat) > 1) stop(sp, ": the pixel rasters do not add up to 14B's national change")
  # scale capped at the 98th percentile of |change| over cells that changed, in the PRODUCTION run
  gp  <- if (nchar(sens$suffix) == 0L) g else grid_density(sp, grid_change(sp, ""))
  dv  <- terra::values(gp$dens, mat = FALSE)
  lim <- nice_ceiling(quantile(abs(dv[!is.na(dv) & dv != 0]), 0.98))
  fmt <- if (lim <= 1) "%+.2f" else if (lim <= 10) "%+.1f" else "%+.0f"
  p <- grid_map(sp, g, lim, "Population change (birds per km\u00b2, 10 km cells)", capped_labels(lim, fmt)) +
    labs(title = headline[[sp]])
  ggsave(file.path(fig_dir, sprintf("change_grid_%s%s.png", sp, sens$suffix)), p, width = 7.5, height = 6.6, dpi = 600)
}
