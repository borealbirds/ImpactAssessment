# ---
# title: Impact Assessment: visualize scale-dependent sector impacts (CAWA, OVEN)
# author: Mannfred Boehm
# ---
# Reads logs/12F_singletons_summary.csv and renders the same sector impacts at
# three increasing geographic scales:
#   (1) footprint   (pct_footprint)   — % of obs on the sector's footprint pixels
#   (2) watershed   (pct_sub_active)  — % of obs on subbasins where the sector is active
#   (3) BCR         (pct_bcr_active)  — % of obs on BCRs where the sector is active
# Sectors on x, % impact on y, species dodged. Error bars are joint BRT×BART SD
# propagated to each scale via SD/% = impact_sd / |impact_mean|.
# ---

suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
})


# set paths ------------------------------------------------------

cc    <- FALSE
local <- TRUE

if (cc)            { ia_dir <- "/home/mannfred/scratch/impact_assessment" }
if (!cc && local)  { ia_dir <- getwd() }
if (!cc && !local) { ia_dir <- file.path("G:/Shared drives/BAM_NationalModels5", "data", "Extras", "sandbox_data", "impactassessment_sandbox") }

csv_path <- file.path(ia_dir, "logs", "12F_singletons_summary.csv")
fig_dir  <- file.path(ia_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)


# read + tidy ------------------------------------------------------

raw <- read.csv(csv_path, stringsAsFactors = FALSE)

# Relative SD (CoV in birds) — apply to each scaled %-impact to get its SD.
sectors <- raw |>
  filter(scope == "sector") |>
  mutate(
    sector = ifelse(label == "dam_and_associated_reservoir", "dams", label),
    rel_sd = ifelse(is.finite(impact_mean) & impact_mean != 0,
                    abs(impact_sd / impact_mean), NA_real_)
  )

# sector order: rank by CAWA BCR-scale impact (descending), apply to both species
sector_order <- sectors |>
  filter(species == "CAWA") |>
  arrange(desc(pct_bcr_active)) |>
  pull(sector)

sectors <- sectors |>
  mutate(
    sector  = factor(sector,  levels = sector_order),
    species = factor(species, levels = c("CAWA", "OVEN"))
  )


# pivot to long, propagate SD per scale ------------------------------

scale_labels <- c(
  pct_footprint  = "Footprint (% of obs on footprint)",
  pct_sub_active = "Watershed (% of obs on active subbasins)",
  pct_bcr_active = "BCR (% of obs on active BCRs)"
)

long <- sectors |>
  select(species, sector, rel_sd,
         pct_footprint, pct_sub_active, pct_bcr_active) |>
  pivot_longer(c(pct_footprint, pct_sub_active, pct_bcr_active),
               names_to = "scale", values_to = "value") |>
  mutate(
    sd    = abs(value) * rel_sd,
    scale = factor(scale, levels = names(scale_labels), labels = scale_labels)
  )


# plot ------------------------------------------------------

dodge <- position_dodge(width = 0.7)

p <- ggplot(long, aes(x = sector, y = value, color = species)) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_errorbar(aes(ymin = value - 2 * sd, ymax = value + 2 * sd),
                position = dodge, width = 0.35,
                linewidth = 0.8, na.rm = TRUE) +
  geom_point(position = dodge, size = 4, na.rm = TRUE) +
  facet_wrap(~ scale, ncol = 1, scales = "free_y") +
  scale_color_manual(values = c(CAWA = "#E69F00", OVEN = "#009E73")) +
  scale_y_continuous(labels = number_format(accuracy = 0.1)) +
  labs(
    x        = NULL,
    y        = "% impact",
    color    = NULL
  ) +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid        = element_blank(),
    panel.background  = element_blank(),
    axis.ticks        = element_line(linewidth = 0.3, color = "grey40"),
    axis.line         = element_line(linewidth = 0.3, color = "grey40"),
    axis.text.x       = element_text(angle = 30, hjust = 1, size = 15),
    axis.text.y       = element_text(size = 14),
    axis.title        = element_text(size = 17),
    strip.text        = element_text(face = "bold", size = 16),
    legend.text       = element_text(size = 15),
    plot.caption      = element_text(size = 11),
    legend.position   = "top"
  )


# save ------------------------------------------------------

out_png <- file.path(fig_dir, "12F_singletons_scales.png")
ggsave(out_png, p, width = 10, height = 9, dpi = 200)
cat("wrote", out_png, "\n")

out_pdf <- file.path(fig_dir, "12F_singletons_scales.pdf")
ggsave(out_pdf, p, width = 10, height = 9)
cat("wrote", out_pdf, "\n")


# =====================================================
# Block 2: spatial impact maps for one sector (CAWA)
# Three figures showing gain/loss intensity at increasing scales:
#   (1) 4 random active subbasins, pixel-level
#   (2) can81, pixel-level
#   (3) national (multi-BCR), aggregated for tractability
# Impact = backfilled - observed (per pixel).
#   positive = gain of CAWA when industry removed   (#D55E00)
#   negative = loss of CAWA when industry removed   (#0072B2)
#   neutral  = no change                             (#999999)
# Sector defaults to built (coalition_2) — only sector with a full set of
# CAWA backfilled rasters locally available. Change `sector_label` /
# `coalition_id` to run other sectors on the cluster.
# =====================================================

suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(patchwork)
})

species      <- "CAWA"
year_spatial <- 2020
sector_label <- "mines"
coalition_id <- 9

bf_dir   <- file.path(ia_dir, "data", "derived_data", "predictions_coalitions",
                      coalition_id, species)
obs_dir  <- file.path(ia_dir, "data", "derived_data", "predictions", species)
basins_p <- file.path(ia_dir, "data", "raw_data",
                      "hydrobasins_masked_merged_subset.gpkg")

bcrs <- list.files(bf_dir)
stopifnot(length(bcrs) > 0)

# Continuous diverging palette: blue (loss) -> light grey (neutral) -> orange (gain)
pal_cont <- c("#0072B2", "#F5F5F5", "#D55E00")

make_impact <- function(bcr) {
  bf  <- rast(file.path(bf_dir,  bcr, year_spatial, "backfilled_mean.tif"))
  obs <- rast(file.path(obs_dir, bcr, year_spatial, "observed_mean.tif"))
  imp <- bf - obs
  # Non-footprint pixels: backfilled_mean is NA but observed is finite — set
  # impact to 0 (no change) so the basin interior renders as neutral.
  imp[is.na(imp) & !is.na(obs)] <- 0
  names(imp) <- "impact"
  imp
}

# Per-figure symmetric limit: 99th-pct of |nonzero| values, used for fill scale.
make_lim <- function(vals) {
  nz <- vals[abs(vals) > .Machine$double.eps & is.finite(vals)]
  if (length(nz) == 0) return(1)
  as.numeric(quantile(abs(nz), 0.99, na.rm = TRUE))
}

# SpatRaster -> long data frame, clamping to [-lim, lim] for stable color mapping.
raster_to_df <- function(r, lim, panel_label = NA_character_) {
  df <- as.data.frame(r, xy = TRUE, na.rm = FALSE)
  names(df)[3] <- "impact"
  finite <- is.finite(df$impact)
  df$impact[finite] <- pmin(pmax(df$impact[finite], -lim), lim)
  if (!is.na(panel_label)) df$panel <- panel_label
  df
}

impact_theme <- theme_minimal(base_size = 14) +
  theme(panel.grid       = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background  = element_rect(fill = "white", color = NA),
        strip.background = element_rect(fill = "white", color = NA),
        axis.text        = element_blank(),
        axis.ticks       = element_blank(),
        axis.title       = element_blank(),
        strip.text       = element_text(face = "bold", size = 12),
        legend.position  = "right",
        legend.key.height = unit(1.2, "cm"))

impact_fill <- function(lim) {
  scale_fill_gradientn(
    colors    = pal_cont,
    limits    = c(-lim, lim),
    oob       = scales::squish,
    na.value  = "white",
    name      = NULL
  )
}


# Lazy-load impact rasters only for BCRs actually requested (avoids slow
# eager loading of all BCRs when only specific subbasins / can81 are needed).
impact_cache <- new.env(parent = emptyenv())
get_impact <- function(bcr) {
  if (is.null(impact_cache[[bcr]])) impact_cache[[bcr]] <- make_impact(bcr)
  impact_cache[[bcr]]
}


# ---- Pick BCR + 4 footprints (each in a distinct subbasin) ----------
# Figure 4 (footprint-scale) drives both Figure 1 (the subbasins that contain
# those footprints) and Figure 2 (the BCR they all share). Iterate BCRs in
# random order; the first BCR with >=4 mid-sized patches falling in distinct
# subbasins is chosen.

basins   <- st_read(basins_p, quiet = TRUE)
set.seed(11)
buffer_m <- 10000  # ~10 km of context around each footprint patch

choose_bcr_and_footprints <- function() {
  for (bcr in sample(bcrs)) {
    imp <- get_impact(bcr)
    fp_mask <- abs(imp) > 0
    fp_mask[!fp_mask] <- NA
    patches_r <- terra::patches(fp_mask, directions = 8, allowGaps = FALSE)
    fr <- terra::freq(patches_r)
    if (is.null(fr) || nrow(fr) == 0) next
    fr <- fr[fr$count >= 3 & fr$count <= 80, , drop = FALSE]
    if (nrow(fr) < 4) next

    basins_proj <- st_transform(basins, crs(imp))
    chosen <- list()
    used_subs <- integer(0)
    for (pid in sample(fr$value)) {
      patch_only <- patches_r == pid
      patch_only[patch_only == 0] <- NA
      patch_ext <- ext(terra::trim(patch_only))
      cx <- (patch_ext$xmin + patch_ext$xmax) / 2
      cy <- (patch_ext$ymin + patch_ext$ymax) / 2
      pt <- st_sfc(st_point(c(cx, cy)), crs = st_crs(basins_proj))
      hit <- which(lengths(st_intersects(basins_proj, pt)) > 0)
      if (length(hit) == 0) next
      sub_id <- hit[1]
      if (sub_id %in% used_subs) next
      used_subs <- c(used_subs, sub_id)
      expanded <- ext(patch_ext$xmin - buffer_m, patch_ext$xmax + buffer_m,
                      patch_ext$ymin - buffer_m, patch_ext$ymax + buffer_m)
      npix <- sum(fr$count[fr$value == pid])
      chosen[[length(chosen) + 1]] <- list(
        sub_id  = sub_id,
        raster  = crop(imp, expanded),
        polygon = basins_proj[sub_id, ],
        npix    = npix
      )
      if (length(chosen) == 4) {
        return(list(bcr = bcr, picks = chosen))
      }
    }
  }
  stop("no BCR had 4 footprints in distinct subbasins")
}

found       <- choose_bcr_and_footprints()
target_bcr  <- found$bcr
fp_picks    <- lapply(seq_along(found$picks), function(i) {
  p <- found$picks[[i]]
  list(raster = p$raster,
       label  = paste0("footprint ", i,
                       " (subbasin ", p$sub_id, ", ", p$npix, " px)"))
})
cat("chosen BCR:", target_bcr, "  subbasins:",
    paste(vapply(found$picks, function(p) p$sub_id, integer(1)), collapse = ", "),
    "\n")


# ---- Figure 1: subbasins containing the 4 footprints ------------------------

picks <- lapply(found$picks, function(p) {
  imp     <- get_impact(target_bcr)
  poly_v  <- vect(p$polygon)
  cropped <- mask(crop(imp, poly_v), poly_v)
  list(
    raster  = cropped,
    polygon = p$polygon,
    label   = paste0("subbasin ", p$sub_id, " (", target_bcr, ")")
  )
})

panel_levels <- vapply(picks, function(p) p$label, character(1))
lim1 <- make_lim(unlist(lapply(picks, function(p) values(p$raster, na.rm = TRUE)),
                        use.names = FALSE))
cat("fig1 color limit:", signif(lim1, 3), "\n")
fig1_df  <- bind_rows(lapply(picks, function(p) raster_to_df(p$raster, lim1, p$label)))
fig1_df$panel <- factor(fig1_df$panel, levels = panel_levels)

# Polygon outlines as plain x/y paths (avoids coord_sf restriction on free scales)
fig1_outlines <- bind_rows(lapply(picks, function(p) {
  cs <- sf::st_coordinates(p$polygon)
  ring_cols <- intersect(c("L1", "L2", "L3"), colnames(cs))
  ring_id   <- do.call(paste, c(as.data.frame(cs[, ring_cols, drop = FALSE]),
                                sep = "_"))
  tibble(x = cs[, "X"], y = cs[, "Y"],
         ring = paste(p$label, ring_id, sep = "::"),
         panel = p$label)
}))
fig1_outlines$panel <- factor(fig1_outlines$panel, levels = panel_levels)

fig1 <- ggplot() +
  geom_tile(data = fig1_df, aes(x = x, y = y, fill = impact)) +
  geom_path(data = fig1_outlines, aes(x = x, y = y, group = ring),
            color = "black", linewidth = 0.4, inherit.aes = FALSE) +
  facet_wrap(~ panel, ncol = 2, scales = "free") +
  impact_fill(lim1) + impact_theme +
  theme(aspect.ratio = 1)

out1 <- file.path(fig_dir, paste0("12F_spatial_watersheds_", sector_label, ".png"))
ggsave(out1, fig1, width = 10, height = 10, dpi = 200)
cat("wrote", out1, "\n")


# ---- Figure 2: target BCR (aggregated to ~5km) with BCR boundary overlay ----

target_agg <- aggregate(get_impact(target_bcr), fact = 5, fun = "mean", na.rm = TRUE)
lim2       <- make_lim(values(target_agg, na.rm = TRUE))
cat("fig2 color limit:", signif(lim2, 3), "\n")
fig2_df <- raster_to_df(target_agg, lim2)

target_subunit <- as.integer(sub("^can", "", target_bcr))
bcr_poly  <- st_read(file.path(ia_dir, "data", "raw_data", "Regions",
                               "BAM_BCR_NationalModel_Unbuffered.shp"),
                     quiet = TRUE) |>
  dplyr::filter(country == "can", subUnit == target_subunit) |>
  st_transform(crs(target_agg))

bcr_coords <- sf::st_coordinates(bcr_poly)
ring_cols2 <- intersect(c("L1", "L2", "L3"), colnames(bcr_coords))
bcr_outline <- tibble(
  x    = bcr_coords[, "X"],
  y    = bcr_coords[, "Y"],
  ring = do.call(paste, c(as.data.frame(bcr_coords[, ring_cols2, drop = FALSE]),
                          sep = "_"))
)

fig2 <- ggplot() +
  geom_tile(data = fig2_df, aes(x = x, y = y, fill = impact)) +
  geom_path(data = bcr_outline, aes(x = x, y = y, group = ring),
            color = "black", linewidth = 0.4, inherit.aes = FALSE) +
  coord_equal() +
  impact_fill(lim2) + impact_theme

out2 <- file.path(fig_dir, paste0("12F_spatial_", target_bcr, "_", sector_label, ".png"))
ggsave(out2, fig2, width = 11, height = 8, dpi = 200)
cat("wrote", out2, "\n")


# ---- Figure 3: national (aggregated to ~10km) with all-Canadian BCR overlay ----

# Eager-load all BCRs for the national merge.
for (bcr in bcrs) get_impact(bcr)
national <- do.call(terra::merge,
                    unname(mget(bcrs, envir = impact_cache)))
# Signed-max: each 10 km cell keeps the value of largest absolute magnitude,
# so isolated point-footprint impacts (mines, dams) survive aggregation.
nat_max  <- aggregate(national, fact = 10, fun = "max", na.rm = TRUE)
nat_min  <- aggregate(national, fact = 10, fun = "min", na.rm = TRUE)
national_agg <- terra::ifel(abs(nat_max) >= abs(nat_min), nat_max, nat_min)
lim3 <- make_lim(values(national_agg, na.rm = TRUE))
cat("fig3 color limit:", signif(lim3, 3), "\n")
fig3_df <- raster_to_df(national_agg, lim3)
# Signed square-root color stretch: amplifies weak western pixels so they
# remain visible against the eastern-dominated range, while preserving sign
# and zero. Legend ticks are still labeled in original birds/km2 units.
sgn_sqrt <- function(x) sign(x) * sqrt(abs(x))
fig3_df$impact_t <- sgn_sqrt(fig3_df$impact)
lim3_t   <- sgn_sqrt(lim3)
ticks_o  <- c(-lim3, -lim3 / 4, -lim3 / 16, 0, lim3 / 16, lim3 / 4, lim3)
ticks_t  <- sgn_sqrt(ticks_o)

all_bcrs_poly <- st_read(file.path(ia_dir, "data", "raw_data", "Regions",
                                   "BAM_BCR_NationalModel_Unbuffered.shp"),
                         quiet = TRUE) |>
  dplyr::filter(country == "can") |>
  st_transform(crs(national_agg))
all_bcrs_coords <- sf::st_coordinates(all_bcrs_poly)
ring_cols3 <- intersect(c("L1", "L2", "L3"), colnames(all_bcrs_coords))
all_bcrs_outline <- tibble(
  x    = all_bcrs_coords[, "X"],
  y    = all_bcrs_coords[, "Y"],
  ring = do.call(paste, c(as.data.frame(all_bcrs_coords[, ring_cols3, drop = FALSE]),
                          sep = "_"))
)

fig3 <- ggplot() +
  geom_tile(data = fig3_df, aes(x = x, y = y, fill = impact_t)) +
  geom_path(data = all_bcrs_outline, aes(x = x, y = y, group = ring),
            color = "black", linewidth = 0.3, inherit.aes = FALSE) +
  coord_equal() +
  scale_fill_gradientn(
    colors   = pal_cont,
    limits   = c(-lim3_t, lim3_t),
    breaks   = ticks_t,
    labels   = signif(ticks_o, 2),
    oob      = scales::squish,
    na.value = "white",
    name     = NULL
  ) +
  impact_theme

out3 <- file.path(fig_dir, paste0("12F_spatial_national_", sector_label, ".png"))
ggsave(out3, fig3, width = 12, height = 9, dpi = 200)
cat("wrote", out3, "\n")


# ---- Figure 4: the 4 footprints chosen above ------------------------

fp_panel_levels <- vapply(fp_picks, function(p) p$label, character(1))
lim4 <- make_lim(unlist(lapply(fp_picks, function(p) values(p$raster, na.rm = TRUE)),
                        use.names = FALSE))
cat("fig4 color limit:", signif(lim4, 3), "\n")
fig4_df <- bind_rows(lapply(fp_picks, function(p) raster_to_df(p$raster, lim4, p$label)))
fig4_df$panel <- factor(fig4_df$panel, levels = fp_panel_levels)

fig4 <- ggplot(fig4_df, aes(x = x, y = y, fill = impact)) +
  geom_tile() +
  facet_wrap(~ panel, ncol = 2, scales = "free") +
  impact_fill(lim4) + impact_theme +
  theme(aspect.ratio = 1)

out4 <- file.path(fig_dir, paste0("12F_spatial_footprints_", sector_label, ".png"))
ggsave(out4, fig4, width = 10, height = 10, dpi = 200)
cat("wrote", out4, "\n")
