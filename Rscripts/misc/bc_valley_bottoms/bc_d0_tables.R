# ---
# title: d0 vs backfill in BC, by elevation relative to each subbasin's training land
# author: Mannfred Boehm
# ---
#
# Colleagues' concern: BC's footprint sits in valley bottoms, the low-HF land BART trains on
# sits higher, so the backfill puts high-elevation vegetation in valleys and reduces birds.
#
# d0 = observed vegetation, footprint covariates at 0. bf = backfilled vegetation, footprint
# covariates at 0. Vegetation is backfilled only on the direct footprint, so bf - d0 is compared
# there, and it is the whole of 14B's vegetation part.
#
# BC = every (subbasin, BCR) unit with at least half its kept footprint pixels inside BC, so
# that the pixel sums and 14B's per-subbasin samples cover the same land. Means come from
# pixels; 5th/95th percentiles come from the samples (14B's rules): bf - d0 per (bootstrap,
# scenario) sample, each % over that sample's own bootstrap's d0. Bands have no bf interval:
# pixel_change/ stores only bf's mean.
#
# Elevation: data/raw_data/dem_bc_1km.tif (AWS Terrain Tiles via elevatr, zoom 7, area-averaged
# to the mosaic grid). V5 has no elevation covariate; BART sees climate, mTPI and roughness.
# Rank = a direct pixel's elevation as a percentile of its own subbasin's low-HF pixels, the
# land 08A trained on there, so 0 = lower than all of them.
#
# Runs locally from the RProject root, after bc_d0_pixels.R.

library(terra)
library(dplyr)
library(tidyr)

ia_dir  <- getwd()
year    <- 2020L
species <- c("RBNU", "SWTH", "RCKI", "OCWA", "TOWA", "VATH")
d_dir   <- file.path(ia_dir, "data", "derived_data", "bc_d0_pixels")
dt_dir  <- file.path(ia_dir, "data", "derived_data", "density_tables")
out_dir <- file.path(ia_dir, "output_figures", "bc_backfill")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

res <- lapply(list.files(d_dir, "\\.rds$", full.names = TRUE), readRDS)
fr  <- bind_rows(lapply(res, function(r) mutate(r$frame, species = r$species, bcr = r$bcr)))
cat(length(res), "species x BCR,", nrow(fr), "pixel rows\n")

# ---- checks against production -------------------------------------------------------
chk <- bind_rows(lapply(res, function(r) as.data.frame(t(r$check)) |> mutate(species = r$species, bcr = r$bcr)))
cat("non-direct pixels, re-predicted d0 vs 12F's stored cf: worst pair off by",
    signif(max(abs(chk$sum_d0 - chk$sum_cf) / pmax(chk$sum_cf, 1)), 3), "(relative); max pixel gap",
    signif(max(chk$max_abs, na.rm = TRUE), 3), "birds; max share of pixels off",
    signif(max(chk$frac_px_off, na.rm = TRUE), 3), "\n")

ssub <- read.csv(file.path(ia_dir, "data", "derived_data", "sector_effects", "shapley_subbasin.csv")) |>
  filter(species %in% .env$species, year == .env$year) |>
  group_by(species, bcr, subbasin) |>
  summarise(direct_part = sum(shapley_direct_mean), veg_part = sum(shapley_vegetation_mean),
            frac_outside_aoa = first(frac_outside_aoa), aoa_flag = first(extrapolation_flag), .groups = "drop")

coal <- bind_rows(lapply(species, function(sp)
  readRDS(file.path(dt_dir, sprintf("%s_%d_coalition_256.rds", sp, year))))) |>
  select(species, bcr, subbasin, obs_on_coalition_mean, bf_on_coalition_mean)

bysub <- fr |>
  group_by(species, bcr, sid) |>
  summarise(n_fp = n(), n_direct = sum(direct), share_bc = mean(in_bc), obs_fp = sum(obs), cf_fp = sum(cf),
            d0_fp = sum(d0), obs_direct = sum(obs[direct]), d0_direct = sum(d0[direct]),
            bf_direct = sum(cf[direct]), .groups = "drop") |>
  left_join(ssub, by = c("species", "bcr", "sid" = "subbasin")) |>
  left_join(coal, by = c("species", "bcr", "sid" = "subbasin"))

rel <- function(a, b) max(abs(a - b) / pmax(abs(b), 1), na.rm = TRUE)
cat("subbasins: obs vs cid-256", signif(rel(bysub$obs_fp, bysub$obs_on_coalition_mean), 3),
    "| cf vs cid-256", signif(rel(bysub$cf_fp, bysub$bf_on_coalition_mean), 3),
    "| d0 vs obs + 14B direct", signif(rel(bysub$d0_fp, bysub$obs_fp + bysub$direct_part), 3),
    "| bf - d0 on direct vs 14B vegetation", signif(rel(bysub$bf_direct - bysub$d0_direct, bysub$veg_part), 3), "\n")

# ---- which units are BC --------------------------------------------------------------
unit <- fr |> group_by(bcr, sid) |> summarise(share_bc = mean(in_bc), .groups = "drop") |> filter(share_bc >= 0.5)
cat(nrow(unit), "(subbasin, BCR) units in BC;", length(unique(unit$sid)), "subbasins\n")
print(table(unit$bcr))
fr <- semi_join(fr, unit, by = c("bcr", "sid"))
cat("share of the units' direct pixels inside BC:", round(mean(fr$in_bc[fr$direct]), 4), "\n")

# ---- elevation of training land and direct footprint ---------------------------------
dem  <- rast(file.path(ia_dir, "data", "raw_data", "dem_bc_1km.tif"))
mos  <- rast(file.path(ia_dir, "data", "raw_data", "covariates_mosaiced", sprintf("covariates_mosaiced_%d.tif", year)))
mtpi <- crop(mos[["mTPI_1km"]], dem)
subs <- vect(file.path(ia_dir, "data", "raw_data", "hydrobasins_masked_merged_subset.gpkg"))
subs$sid <- seq_len(nrow(subs))
sids <- sort(unique(unit$sid))
zone <- rasterize(subs[sids, ], dem, field = "sid")
low  <- project(rast(file.path(ia_dir, "data", "raw_data", "hirshpearson", "CanHF_1km_lessthan1.tif")), dem, method = "near")

tr <- data.frame(sid = values(zone, mat = FALSE), low = values(low, mat = FALSE),
                 elev = values(dem, mat = FALSE), mtpi = values(mtpi, mat = FALSE))
tr <- tr[!is.na(tr$sid) & !is.na(tr$low) & tr$low == 1, ]
cat("training pixels:", nrow(tr), "| without elevation:", sum(is.na(tr$elev)), "\n")
tr <- tr[!is.na(tr$elev), ]
train <- split(tr$elev, tr$sid)

xy <- unique(fr[fr$direct, c("x", "y")])
xy$elev <- terra::extract(dem, as.matrix(xy[, c("x", "y")]))[, 1]
xy$mtpi <- terra::extract(mtpi, as.matrix(xy[, c("x", "y")]))[, 1]
dir_px  <- fr |> filter(direct) |> left_join(xy, by = c("x", "y"))
stopifnot(!anyNA(dir_px$elev))

ranks <- dir_px |> distinct(sid, elev) |> group_by(sid) |>
  mutate(rank = ecdf(train[[as.character(sid[1])]])(elev),
         below_min = elev < min(train[[as.character(sid[1])]])) |> ungroup()
dir_px <- left_join(dir_px, ranks, by = c("sid", "elev"))

rank_lab <- c("below 5th pct", "5-25th", "25-50th", "50-75th", "above 75th")
dir_px <- dir_px |> mutate(
  rank_band = cut(rank, c(-Inf, 0.05, 0.25, 0.5, 0.75, Inf), labels = rank_lab, right = FALSE),
  elev_band = cut(elev, c(-Inf, 500, 1000, 1500, Inf), labels = c("< 500 m", "500-1000 m", "1000-1500 m", "> 1500 m"), right = FALSE))

# terrain: one row per pixel, not per species
one <- dir_px |> distinct(bcr, sid, cell, .keep_all = TRUE)
tr_u <- tr |> semi_join(distinct(one, sid), by = "sid")
terrain <- function(d, t, label) data.frame(
  unit = label, direct_km2 = nrow(d), training_km2 = nrow(t),
  direct_elev_median = median(d$elev), training_elev_median = median(t$elev),
  direct_mtpi_median = median(d$mtpi, na.rm = TRUE), training_mtpi_median = median(t$mtpi, na.rm = TRUE),
  pct_direct_below_training_q05 = 100 * mean(d$rank < 0.05),
  pct_direct_below_training_q25 = 100 * mean(d$rank < 0.25),
  pct_direct_below_training_median = 100 * mean(d$rank < 0.5),
  pct_direct_below_all_training = 100 * mean(d$below_min))
ter <- bind_rows(terrain(one, tr_u, "BC"), lapply(sort(unique(one$bcr)), function(b) {
  d <- one[one$bcr == b, ]; terrain(d, tr_u[tr_u$sid %in% d$sid, ], b) }))
write.csv(ter, file.path(out_dir, "terrain_direct_vs_training.csv"), row.names = FALSE)

# ---- d0 vs bf: BC and BCR, with intervals ---------------------------------------------
summ <- function(sp, keep_bcr = NULL) {
  ss   <- readRDS(file.path(dt_dir, sprintf("%s_%d_shapley_samples.rds", sp, year)))
  u    <- if (is.null(keep_bcr)) unit else unit[unit$bcr %in% keep_bcr, ]
  rows <- which(paste(ss$bcr, ss$subbasin) %in% paste(u$bcr, u$sid))
  # the release filter: only units this species has a frame for
  have <- unique(paste(dir_px$bcr[dir_px$species == sp], dir_px$sid[dir_px$species == sp]))
  rows <- rows[paste(ss$bcr[rows], ss$subbasin[rows]) %in% have]
  boot_of <- rep(seq_len(ss$n_boot), each = ss$n_scen)
  tot   <- apply(ss$phi[rows, , , drop = FALSE], 3, sum)                # [sample]
  tot_d <- apply(ss$phi_direct[rows, , , drop = FALSE], 3, sum)         # [bootstrap]
  veg_s <- tot - tot_d[boot_of]

  d0b <- 0
  for (r in res) if (r$species == sp && (is.null(keep_bcr) || r$bcr %in% keep_bcr)) {
    f   <- r$frame[r$frame$direct, ]
    use <- paste(f$sid) %in% paste(unit$sid[unit$bcr == r$bcr])
    d0b <- d0b + colSums(r$d0_boot[use, , drop = FALSE])
  }
  bf_s <- d0b[boot_of] + veg_s
  p <- dir_px[dir_px$species == sp & (if (is.null(keep_bcr)) TRUE else dir_px$bcr %in% keep_bcr), ]
  q <- function(x, pr) unname(quantile(x, pr))
  data.frame(species = sp, unit = if (is.null(keep_bcr)) "BC" else keep_bcr,
             direct_km2 = nrow(p), obs = sum(p$obs),
             d0 = sum(p$d0), d0_q05 = q(d0b, .05), d0_q95 = q(d0b, .95),
             bf = sum(p$cf), bf_q05 = q(bf_s, .05), bf_q95 = q(bf_s, .95),
             bf_minus_d0 = sum(p$cf) - sum(p$d0), diff_q05 = q(veg_s, .05), diff_q95 = q(veg_s, .95),
             diff_mean_samples = mean(veg_s),
             pct = 100 * (sum(p$cf) - sum(p$d0)) / sum(p$d0),
             pct_q05 = q(100 * veg_s / d0b[boot_of], .05), pct_q95 = q(100 * veg_s / d0b[boot_of], .95))
}
t_bc  <- bind_rows(lapply(species, summ))
t_bcr <- bind_rows(lapply(species, function(sp) bind_rows(lapply(sort(unique(dir_px$bcr[dir_px$species == sp])), function(b) summ(sp, b)))))
cat("pixel means vs sample means of bf - d0, worst:", signif(rel(t_bc$bf_minus_d0, t_bc$diff_mean_samples), 3), "\n")
write.csv(t_bc,  file.path(out_dir, "d0_vs_bf_by_species.csv"), row.names = FALSE)
write.csv(t_bcr, file.path(out_dir, "d0_vs_bf_by_species_bcr.csv"), row.names = FALSE)

# ---- by elevation --------------------------------------------------------------------
band <- function(d, ...) d |> group_by(species, ...) |>
  summarise(direct_km2 = n(), obs = sum(obs), d0 = sum(d0), bf = sum(cf), .groups = "drop") |>
  mutate(bf_minus_d0 = bf - d0, pct = 100 * bf_minus_d0 / d0,
         d0_per_km2 = d0 / direct_km2, bf_per_km2 = bf / direct_km2)
t_rank <- band(dir_px, rank_band)
t_elev <- band(dir_px, elev_band)
write.csv(t_rank, file.path(out_dir, "d0_vs_bf_by_elevation_rank.csv"), row.names = FALSE)
write.csv(band(dir_px, bcr, rank_band), file.path(out_dir, "d0_vs_bf_by_elevation_rank_bcr.csv"), row.names = FALSE)
write.csv(t_elev, file.path(out_dir, "d0_vs_bf_by_elevation.csv"), row.names = FALSE)

# ---- per subbasin ---------------------------------------------------------------------
sub_elev <- one |> group_by(bcr, sid) |>
  summarise(direct_elev_median = median(elev), pct_below_training_q05 = 100 * mean(rank < 0.05), .groups = "drop") |>
  mutate(training_elev_median = vapply(as.character(sid), function(s) median(train[[s]]), numeric(1L)),
         elev_gap = training_elev_median - direct_elev_median)
t_sub <- bysub |> semi_join(unit, by = c("bcr", "sid")) |> filter(n_direct > 0) |>
  left_join(sub_elev, by = c("bcr", "sid")) |>
  transmute(species, bcr, subbasin = sid, n_direct, obs_direct, d0_direct, bf_direct,
            bf_minus_d0 = bf_direct - d0_direct, pct = 100 * bf_minus_d0 / d0_direct,
            direct_elev_median, training_elev_median, elev_gap, pct_below_training_q05,
            frac_outside_aoa, aoa_flag)
write.csv(t_sub, file.path(out_dir, "d0_vs_bf_by_subbasin.csv"), row.names = FALSE)

saveRDS(dir_px, file.path(d_dir, "..", "bc_d0_direct_pixels.rds"))

options(width = 170, digits = 4)
cat("\nTERRAIN\n"); print(ter)
cat("\nBC\n"); print(t_bc |> select(species, direct_km2, obs, d0, bf, bf_minus_d0, diff_q05, diff_q95, pct, pct_q05, pct_q95) |> mutate(across(where(is.numeric), ~ round(.x, 1))))
cat("\nBY BCR\n"); print(t_bcr |> select(species, unit, direct_km2, d0, bf, bf_minus_d0, diff_q05, diff_q95, pct) |> mutate(across(where(is.numeric), ~ round(.x, 1))), row.names = FALSE)
cat("\nPCT BY ELEVATION RANK\n"); print(t_rank |> select(species, rank_band, pct) |> pivot_wider(names_from = rank_band, values_from = pct) |> mutate(across(where(is.numeric), ~ round(.x, 1))))
cat("\nkm2 BY ELEVATION RANK\n"); print(one |> count(rank_band))
cat("\nBIRDS (bf - d0) BY ELEVATION RANK\n"); print(t_rank |> select(species, rank_band, bf_minus_d0) |> pivot_wider(names_from = rank_band, values_from = bf_minus_d0) |> mutate(across(where(is.numeric), ~ round(.x))))
cat("\nPCT BY ELEVATION\n"); print(t_elev |> select(species, elev_band, pct) |> pivot_wider(names_from = elev_band, values_from = pct) |> mutate(across(where(is.numeric), ~ round(.x, 1))))
cat("\nkm2 BY ELEVATION\n"); print(one |> count(elev_band))
cat("TABLES_DONE\n")
