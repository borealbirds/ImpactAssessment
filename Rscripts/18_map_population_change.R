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
#       (drawn only when both are mapped)
#   (3) the pixel-level change on a 10 km grid, one figure per species
# Display CRS is EPSG:3978 (display only).
#
# Species: every species in 14B's tables, or the ones named in SPECIES (comma-separated):
#   SPECIES=CAWA,OVEN Rscript Rscripts/18_map_population_change.R
# On Fir, 14B's .sh chains it after 14C; on its own: SPECIES=CAWA,OVEN sbatch 18_map_population_change.sh
# A species that fails is reported and skipped, and the script stops at the end if any did.
# Titles use data/raw_data/species_names.csv (BAMexploreR::spp_tbl, plus the 8 codes it
# lacks); a species missing from it is titled by its code.
#
# No sf: the geometry is terra's, and ggplot2 draws it with geom_polygon from terra::geom(),
# one group per polygon and one subgroup per ring, so holes are cut out as geom_sf cuts them.
#
# A sensitivity run is mapped by setting its switch, as for 14B (12E sensitivity_run()):
#   DENS_MIN=0.005 Rscript Rscripts/18_map_population_change.R
# reads sector_effects_dmin0.005/ and writes the same figures with a _dmin0.005 suffix. The
# absolute maps keep the production run's colour scale, so the two can be compared directly.
# ---

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})


# set paths ------------------------------------------------------

# inside a Slurm job this is the cluster; anywhere else, the local RProject
cc    <- nzchar(Sys.getenv("SLURM_JOB_ID"))
local <- !cc

if (cc)            { ia_dir <- "/home/mannfred/scratch/impact_assessment" }
if (!cc && local)  { ia_dir <- getwd() }
if (!cc && !local) { ia_dir <- file.path("G:/Shared drives/BAM_NationalModels5", "data", "Extras", "sandbox_data", "impactassessment_sandbox") }

source(file.path(ia_dir, "Rscripts", "12E_shapley_utils.R"))
sens <- sensitivity_run()

shap_path  <- file.path(ia_dir, "data", "derived_data", paste0("sector_effects", sens$suffix), "shapley_subbasin.csv")
prod_path  <- file.path(ia_dir, "data", "derived_data", "sector_effects", "shapley_subbasin.csv")
basin_path <- file.path(ia_dir, "data", "raw_data", "hydrobasins_masked_merged_subset.gpkg")
bcr_path   <- file.path(ia_dir, "data", "raw_data", "Regions", "BAM_BCR_NationalModel_Unbuffered.shp")
names_path <- file.path(ia_dir, "data", "raw_data", "species_names.csv")
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


# species --------------------------------------------------------

in_14b <- sort(unique(sub$species))
spp    <- if (nzchar(Sys.getenv("SPECIES"))) trimws(strsplit(Sys.getenv("SPECIES"), ",")[[1]]) else in_14b
miss   <- setdiff(spp, in_14b)
if (length(miss) > 0L) stop("not in ", shap_path, ": ", paste(miss, collapse = ", "))
sub <- filter(sub, species %in% spp)
message("mapping ", length(spp), " species: ", paste(spp, collapse = " "))

if (!file.exists(names_path)) stop("missing ", names_path)
nm    <- read.csv(names_path)
lab   <- setNames(nm$common_name, nm$species)
nolab <- setdiff(spp, names(lab))
if (length(nolab) > 0L) message("no common name in ", names_path, " for ", paste(nolab, collapse = ", "),
                                ": titled by code")
lab <- c(lab, setNames(nolab, nolab))

tot <- sub |> group_by(species) |> summarise(o = sum(obs), c = sum(cf))
# the % is 14C's change_pct: what industry did (obs - cf) as % of the no-industry population
headline <- setNames(sprintf("%s\n%.2fM birds in 2020 vs %.2fM without industry (%+.1f%%)",
                             lab[tot$species], tot$o / 1e6, tot$c / 1e6,
                             100 * (tot$o - tot$c) / tot$c), tot$species)


# geometry -------------------------------------------------------

# a SpatVector as the vertex table geom_polygon draws: group = one polygon (feature x part),
# subgroup = one ring of it (hole 0 = outer, 1.. = holes), the feature's attributes on every vertex.
# geom_polygon draws groups in level order, so the levels keep the features' order, as geom_sf
# draws them: sorted, the unchanged cells' outlines landed on top of changed cells' edges.
in_order <- function(x) factor(x, levels = unique(x))
vertices <- function(v) {
  g <- as.data.frame(terra::geom(v))
  a <- as.data.frame(v)
  if (ncol(a) > 0L) g <- cbind(a[g$geom, , drop = FALSE], g)
  rownames(g) <- NULL
  mutate(g, group = in_order(paste(geom, part)))
}

bcr_can <- terra::vect(bcr_path)
bcr_can <- bcr_can[bcr_can$country == "can", ]
bcr     <- terra::project(bcr_can, "EPSG:3978")
bcr_xy  <- vertices(bcr)

# Keep only the Canadian part of each subbasin: the birds are counted inside the Canadian BCRs
# only (12C's in-BCR weight), so the area that turns them into birds/km2 must be Canadian too.
# Planar area (transform = FALSE) in the equal-area EPSG:5072.
canada <- bcr_can |> terra::project("EPSG:5072") |> terra::makeValid() |> terra::aggregate()
basins <- terra::vect(basin_path)[, "first_HYBAS_ID"]
names(basins) <- "HYBAS_ID"
basins <- terra::crop(terra::makeValid(basins), canada)
basins$area_km2 <- terra::expanse(basins, unit = "km", transform = FALSE)
basin_area <- as.data.frame(basins)
basin_xy   <- vertices(terra::project(basins, "EPSG:3978"))

# the map's vertices for species s, with each subbasin's change per km2
basin_map <- function(s) {
  basin_xy |>
    inner_join(filter(sub, species %in% s), by = "HYBAS_ID", relationship = "many-to-many") |>
    mutate(dens = change / area_km2, group = in_order(paste(species, group)))
}

# everything outside Canada, painted white over the 10 km cells
bb       <- as.vector(terra::ext(bcr))
can_3978 <- terra::project(canada, "EPSG:3978")
outside  <- terra::as.polygons(terra::ext(bcr), crs = terra::crs(bcr)) |>
  terra::buffer(100000) |>
  terra::erase(terra::aggregate(terra::makeValid(can_3978)))
outside_xy <- vertices(outside)


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

map_theme <- function() {
  theme_void(base_size = 11) +
    theme(legend.position = "bottom", legend.key.width = unit(2.2, "cm"),
          legend.key.height = unit(0.35, "cm"), legend.title.position = "top",
          legend.title = element_text(hjust = 0.5),
          plot.title = element_text(size = 11, hjust = 0.5, margin = margin(b = 4)),
          strip.text = element_text(size = 11, margin = margin(b = 4)),
          plot.background = element_rect(fill = "white", colour = NA),
          plot.margin = margin(12, 12, 12, 12))
}

change_map <- function(data, fill_var, lim, legend, labels_fun) {
  ggplot(mapping = aes(x, y, group = group, subgroup = hole)) +
    geom_polygon(data = bcr_xy, fill = land, colour = "white", linewidth = 0.3) +
    geom_polygon(data = data, aes(fill = pmax(pmin({{ fill_var }}, lim), -lim)), colour = NA) +
    geom_polygon(data = bcr_xy, fill = NA, colour = "white", linewidth = 0.3) +
    scale_fill_gradient2(low = red, mid = mid, high = blue, midpoint = 0,
                         limits = c(-lim, lim), breaks = c(-1, -0.5, 0, 0.5, 1) * lim,
                         labels = labels_fun,
                         na.value = land, name = legend) +
    coord_equal() +
    map_theme()
}


# (1) absolute change, one figure per species ---------------------

# cap each species at the 98th percentile of |change| in the PRODUCTION run, so a few small,
# dense subbasins don't wash out the rest and a sensitivity run shares production's scale
lims <- subbasin_change(prod_path) |>
  filter(species %in% spp) |>
  inner_join(basin_area, by = "HYBAS_ID") |>
  group_by(species) |>
  summarise(lim = nice_ceiling(quantile(abs(change / area_km2), 0.98, na.rm = TRUE)))
lims <- setNames(lims$lim, lims$species)

abs_map <- function(sp) {
  d   <- basin_map(sp)
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

pct_map <- function() {
  p <- change_map(basin_map(c("CAWA", "OVEN")), pct, 50, "Percent population change",
                  capped_labels(50, "%+.0f%%")) +
    facet_wrap(~species, ncol = 2, labeller = as_labeller(headline))
  ggsave(file.path(fig_dir, paste0("change_map_pct_CAWA_OVEN", sens$suffix, ".png")), p,
         width = 12, height = 6.2, dpi = 200)
}


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
  modelled <- bcr_can[paste0(bcr_can$country, bcr_can$subUnit) %in% bcr_tab$bcr[bcr_tab$species == sp], ] |>
    terra::project("EPSG:5072") |> terra::makeValid() |> terra::aggregate()
  e    <- terra::ext(modelled)
  tmpl <- terra::rast(terra::ext(floor(e$xmin / cell_m) * cell_m, ceiling(e$xmax / cell_m) * cell_m,
                                 floor(e$ymin / cell_m) * cell_m, ceiling(e$ymax / cell_m) * cell_m),
                      resolution = cell_m, crs = "EPSG:5072")
  # km2 modelled = 1 km pixels touching the modelled BCRs, counted as 12C's in-BCR weight counts
  # them (touches = TRUE), so every pixel 12F predicts lies in a cell with area. (cover = TRUE
  # estimates from sub-cell points and gave 0 to cells a BCR edge barely clips: BOBO lost 4 birds.)
  fine <- terra::rast(terra::ext(tmpl), resolution = 1000, crs = "EPSG:5072")
  area <- terra::rasterize(modelled, fine, touches = TRUE, background = 0) |>
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
  modelled <- bcr_xy$subUnit %in% sub("^can", "", bcr_tab$bcr[bcr_tab$species == sp])
  # equal-valued neighbours (mostly the unchanged 0s) are merged, and each cell is outlined in
  # its own colour, or the squares' anti-aliased edges show as seams
  cells <- terra::as.polygons(g$dens, aggregate = TRUE, na.rm = TRUE) |>
    terra::project("EPSG:3978") |>
    vertices()
  ggplot(mapping = aes(x, y, group = group, subgroup = hole)) +
    geom_polygon(data = bcr_xy, fill = land, colour = NA) +
    geom_polygon(data = cells, aes(fill = pmax(pmin(dens, lim), -lim), colour = pmax(pmin(dens, lim), -lim)),
                 linewidth = 0.1) +
    geom_polygon(data = bcr_xy[!modelled, ], fill = land, colour = NA) +
    geom_polygon(data = outside_xy, fill = "white", colour = NA) +
    geom_polygon(data = bcr_xy, fill = NA, colour = "white", linewidth = 0.3) +
    scale_fill_gradient2(low = red, mid = mid, high = blue, midpoint = 0,
                         limits = c(-lim, lim), breaks = c(-1, -0.5, 0, 0.5, 1) * lim,
                         labels = labels_fun, na.value = land, name = legend,
                         aesthetics = c("fill", "colour")) +
    coord_equal(xlim = bb[1:2], ylim = bb[3:4]) +
    map_theme()
}

nat_tab <- read.csv(sub("shapley_subbasin", "shapley_national", shap_path)) |> distinct(species, total_HF_impact)

change_grid <- function(sp) {
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


# run ------------------------------------------------------------

# one species' failure is reported and skipped, so a batch maps every species it can
failed <- character(0)
for (sp in spp) {
  t0 <- Sys.time()
  ok <- tryCatch({ abs_map(sp); change_grid(sp); TRUE },
                 error = function(e) { message(sp, ": FAILED: ", conditionMessage(e)); FALSE })
  if (!ok) failed <- c(failed, sp)
  message(sprintf("%s: %.0f s", sp, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}
if (all(c("CAWA", "OVEN") %in% spp)) {
  ok <- tryCatch({ pct_map(); TRUE },
                 error = function(e) { message("CAWA/OVEN % map: FAILED: ", conditionMessage(e)); FALSE })
  if (!ok) failed <- c(failed, "CAWA/OVEN % map")
}

message("figures in ", fig_dir)
if (length(failed) > 0L) stop(length(failed), " of ", length(spp), " failed: ", paste(failed, collapse = ", "))
