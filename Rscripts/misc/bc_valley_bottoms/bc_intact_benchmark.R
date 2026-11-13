# ---
# title: Backfilled density against intact land at the same elevation, BC
# author: Mannfred Boehm
# ---
#
# bf > d0 cannot clear the backfill: d0 is cleared land, so nearly any forest adds forest birds.
# The question is whether bf falls short of intact land at the pixel's OWN elevation. any_elev
# (the subbasin's intact land in general, which sits higher) is context only: density up there
# reflects its climate as well as its vegetation, so it is not what imported vegetation would
# give a pixel that keeps its valley climate.
#
# For each direct pixel of bc_d0_tables.R, within its own subbasin and BCR, from bc_intact_d0.R
# (intact land predicted at footprint 0, as bf is, so vegetation is what differs):
#   same_elev = mean density of low-HF pixels in its 100 m elevation band (at least min_n of
#               them; a pixel without that is "unmatched" and left out of every column)
#   any_elev  = mean density of all the subbasin's low-HF land, weighted by each band's true area
#   same_obs  = same_elev's pixels at their observed footprint, to show what the step to 0 does
# All in weighted birds/km2, over pixels and land with weight > 0 (so not water). Intact land can
# still differ from valley bottoms in other ways than elevation.
#
# Runs locally from the RProject root, after bc_intact_d0.R.

library(dplyr)
library(tidyr)

ia_dir  <- getwd()
bin_m   <- 100
min_n   <- 5L
out_dir <- file.path(ia_dir, "output_figures", "bc_backfill")

dir_px <- readRDS(file.path(ia_dir, "data", "derived_data", "bc_d0_direct_pixels.rds")) |>
  filter(w > 0) |> mutate(bin = floor(elev / bin_m))
res    <- lapply(list.files(file.path(ia_dir, "data", "derived_data", "bc_d0_pixels", "intact"),
                            full.names = TRUE), readRDS)

bins <- bind_rows(lapply(res, function(r) {
  r$frame |> group_by(species, bcr, sid, bin) |>
    summarise(same_elev = mean(d0), same_obs = mean(obs, na.rm = TRUE), n_same = n(), .groups = "drop") |>
    mutate(key = paste(sid, bin)) |> left_join(r$n_land, by = "key") |> select(-key)
}))
stopifnot(!anyNA(bins$n_land))
any_e <- bins |> group_by(species, bcr, sid) |>
  summarise(any_elev = weighted.mean(same_elev, n_land), .groups = "drop")
same  <- bins |> filter(n_same >= min_n) |> select(species, bcr, sid, bin, same_elev, same_obs)

px <- dir_px |> left_join(same, by = c("species", "bcr", "sid", "bin")) |>
  left_join(any_e, by = c("species", "bcr", "sid")) |>
  select(species, bcr, sid, rank_band, elev, obs, d0, cf, same_elev, same_obs, any_elev)

cat("direct pixels with weight > 0:", nrow(px), "| matched:", round(100 * mean(!is.na(px$same_elev)), 1), "%\n")
print(as.data.frame(px |> group_by(rank_band) |> summarise(pct_matched = round(100 * mean(!is.na(same_elev)), 1))))

comp <- function(d, ...) d |> filter(!is.na(same_elev)) |> group_by(species, ...) |>
  summarise(km2 = n(), d0 = mean(d0), bf = mean(cf), same_elev = mean(same_elev), any_elev = mean(any_elev),
            same_obs = mean(same_obs), .groups = "drop") |>
  mutate(bf_vs_same_pct = 100 * (bf / same_elev - 1), d0_vs_same_pct = 100 * (d0 / same_elev - 1),
         any_vs_same_pct = 100 * (any_elev / same_elev - 1), obs_vs_same_pct = 100 * (same_obs / same_elev - 1))

t_bc   <- comp(px)
t_bcr  <- comp(px, bcr)
t_rank <- comp(px, rank_band)
t_rb   <- comp(px, bcr, rank_band)
unm <- px |> group_by(species) |>
  summarise(km2_unmatched = sum(is.na(same_elev)), pct_km2 = 100 * mean(is.na(same_elev)),
            bf_minus_d0_unmatched = sum((cf - d0)[is.na(same_elev)]), bf_minus_d0_all = sum(cf - d0), .groups = "drop")

# Size of the gap in birds: if every matched direct pixel held intact land's same-elevation
# density instead of bf. A what-if, not a correction: intact land is a selected sample.
size <- px |> filter(!is.na(same_elev)) |> group_by(species) |>
  summarise(same_minus_bf = sum(same_elev - cf), .groups = "drop") |>
  left_join(read.csv(file.path(ia_dir, "data", "derived_data", "sector_effects", "results_summary.csv")) |>
              select(species, any_of(c("observed", "no_industry", "change", "change_pct", "pct_change"))), by = "species")
write.csv(size,   file.path(out_dir, "bf_vs_intact_in_birds.csv"), row.names = FALSE)

write.csv(t_bc,   file.path(out_dir, "bf_vs_intact_by_species.csv"), row.names = FALSE)
write.csv(t_bcr,  file.path(out_dir, "bf_vs_intact_by_species_bcr.csv"), row.names = FALSE)
write.csv(t_rank, file.path(out_dir, "bf_vs_intact_by_elevation_rank.csv"), row.names = FALSE)
write.csv(t_rb,   file.path(out_dir, "bf_vs_intact_by_elevation_rank_bcr.csv"), row.names = FALSE)
write.csv(unm,    file.path(out_dir, "bf_unmatched_pixels.csv"), row.names = FALSE)

options(width = 170, digits = 4)
r1 <- function(d) as.data.frame(mutate(d, across(where(is.numeric), ~ round(.x, 1))))
wide <- function(d, v) r1(d |> select(species, rank_band, all_of(v)) |> pivot_wider(names_from = rank_band, values_from = all_of(v)))
cat("\nBC, birds/km2\n"); print(r1(t_bc), row.names = FALSE)
cat("\nBY BCR, birds/km2\n"); print(r1(t_bcr |> select(-same_obs, -obs_vs_same_pct) |> arrange(bcr, species)), row.names = FALSE)
cat("\nbf vs same-elevation intact land (%), by elevation rank\n"); print(wide(t_rank, "bf_vs_same_pct"), row.names = FALSE)
cat("\nd0 vs same-elevation intact land (%), by elevation rank\n"); print(wide(t_rank, "d0_vs_same_pct"), row.names = FALSE)
cat("\nany-elevation vs same-elevation intact land (%), by elevation rank\n"); print(wide(t_rank, "any_vs_same_pct"), row.names = FALSE)
for (b in c("can5", "can10", "can9")) {
  cat("\n", b, ": bf vs same-elevation (%), by rank\n", sep = "")
  print(wide(filter(t_rb, bcr == b), "bf_vs_same_pct"), row.names = FALSE)
}
cat("\nunmatched pixels\n"); print(r1(unm), row.names = FALSE)
cat("\nsame-elevation intact minus bf, birds, against the national result\n"); print(r1(size), row.names = FALSE)
cat("BENCH_DONE\n")
