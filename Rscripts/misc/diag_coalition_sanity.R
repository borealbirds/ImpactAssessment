# diag_coalition_sanity.R
# Spot-check newly-fetched coalition density tables for sensible predictions.
#
# Schema (per density_tables/{sp}_{yr}_coalition_{cid}.rds, one row per subbasin):
#   species, subbasin, bcr, coalition_id,
#   obs_total_mean/sd,            # observed birds across subbasin
#   obs_on_coalition_mean/sd,     # observed on this coalition's footprint pixels
#   bf_on_coalition_mean/sd       # backfilled (counterfactual) on those pixels
# Per-subbasin standalone impact v(S) = bf_on_coalition - obs_on_coalition
# Positive impact = "birds recovered if this coalition's footprint were removed."

suppressPackageStartupMessages({ library(dplyr); library(tidyr) })

ia_dir <- getwd()
dt_dir <- file.path(ia_dir, "data", "derived_data", "density_tables")

# canonical sector order (alphabetical, matches 12E_shapley_utils.R) ---------
sectors <- sort(c("built", "crop", "dam_and_associated_reservoir",
                  "mines", "oil_gas", "pasture", "rail", "roads"))

cid_to_sectors <- function(cid) {
  bits <- as.logical(intToBits(cid - 1L)[1:length(sectors)])
  sectors[bits]
}

# coalitions of interest -----------------------------------------------------
cids <- c(2, 5, 9, 17, 18, 64)   # built, dam, mines, oil_gas, built+oil_gas, floor

summarise_one <- function(sp, cid) {
  p <- file.path(dt_dir, sprintf("%s_2020_coalition_%d.rds", sp, cid))
  if (!file.exists(p)) return(NULL)
  dt <- readRDS(p)
  imp <- dt$bf_on_coalition_mean - dt$obs_on_coalition_mean
  v   <- dt$bf_on_coalition_sd^2 + dt$obs_on_coalition_sd^2
  active <- !is.na(dt$obs_on_coalition_mean) & dt$obs_on_coalition_mean > 0

  tibble(
    species          = sp,
    coalition_id     = cid,
    sectors          = paste(cid_to_sectors(cid), collapse = "+"),
    n_sub            = nrow(dt),
    n_sub_active     = sum(active),
    n_bcr_active     = sum(tapply(active, dt$bcr, any), na.rm = TRUE),
    obs_total        = sum(dt$obs_total_mean,      na.rm = TRUE),
    obs_on_fp        = sum(dt$obs_on_coalition_mean, na.rm = TRUE),
    bf_on_fp         = sum(dt$bf_on_coalition_mean, na.rm = TRUE),
    impact_mean      = sum(imp, na.rm = TRUE),
    impact_sd        = sqrt(sum(v, na.rm = TRUE)),
    pct_footprint    = if (sum(dt$obs_on_coalition_mean, na.rm=TRUE) > 0)
                          100 * sum(imp, na.rm=TRUE) /
                          sum(dt$obs_on_coalition_mean, na.rm=TRUE) else NA_real_,
    pct_obs_natl     = if (sum(dt$obs_total_mean, na.rm=TRUE) > 0)
                          100 * sum(imp, na.rm=TRUE) /
                          sum(dt$obs_total_mean, na.rm=TRUE) else NA_real_,
    n_neg_imp_sub    = sum(imp < 0, na.rm = TRUE),    # subbasins where obs > bf
    n_pos_imp_sub    = sum(imp > 0, na.rm = TRUE),
    n_na_obs_with_bf = sum(!is.na(dt$bf_on_coalition_mean) &
                             is.na(dt$obs_on_coalition_mean)),
    frac_zero_obs    = mean(dt$obs_on_coalition_mean == 0, na.rm = TRUE)
  )
}

grid <- expand.grid(sp = c("CAWA","OVEN"), cid = cids,
                    stringsAsFactors = FALSE)
out  <- bind_rows(Map(summarise_one, grid$sp, grid$cid))

cat("============================================================\n")
cat(" coalition-level totals (national) per (species, coalition)\n")
cat("============================================================\n")
print(out %>% select(species, coalition_id, sectors,
                     n_sub_active, n_bcr_active,
                     obs_total, obs_on_fp, bf_on_fp,
                     impact_mean, impact_sd,
                     pct_footprint, pct_obs_natl),
      n = Inf, width = 160)

cat("\n--- sanity flags ---\n")
print(out %>% select(species, coalition_id, sectors,
                     n_neg_imp_sub, n_pos_imp_sub,
                     n_na_obs_with_bf, frac_zero_obs),
      n = Inf, width = 140)

# ----------------------------------------------------------------
# monotonicity: coalition_18 = built (2) + oil_gas (17). The
# standalone impact of {built, oil_gas} on its footprint should be
# >= each singleton's impact AND <= built_impact + oil_gas_impact
# (subadditivity in expectation, because adding a sector only adds
# footprint pixels; on those new pixels impact may be + or -).
cat("\n--- monotonicity check: built (cid 2) vs oil_gas (cid 17) vs union (cid 18) ---\n")
mono <- out %>% filter(coalition_id %in% c(2, 17, 18)) %>%
  select(species, coalition_id, sectors, impact_mean, pct_obs_natl) %>%
  pivot_wider(names_from = coalition_id, values_from = c(impact_mean, pct_obs_natl)) %>%
  mutate(sum_singletons = impact_mean_2 + impact_mean_17,
         union_minus_sum = impact_mean_18 - sum_singletons,
         union_minus_max = impact_mean_18 - pmax(impact_mean_2, impact_mean_17))
print(mono, width = 160)

cat("\n--- floor (coalition_64) vs sum of 6 non-rail/roads singletons ---\n")
floor_sectors <- setdiff(sectors, c("rail", "roads"))
single_cids   <- sapply(floor_sectors, function(s)
                         { bits <- sectors == s; sum(2L^(which(bits)-1L)) + 1L })
# We only fetched a subset of singletons (built, dam, mines, oil_gas).
# crop and pasture singletons not fetched -- skip the strict 6-sum comparison.
have <- intersect(single_cids, cids)
cat("singletons fetched (cid):", paste(have, collapse=","), "\n")
cat("(no crop/pasture singletons in this fetch; sum is partial)\n")
partial_sum <- out %>% filter(coalition_id %in% have) %>%
  group_by(species) %>%
  summarise(partial_sum_singletons = sum(impact_mean), .groups = "drop")
floor_imp <- out %>% filter(coalition_id == 64) %>%
  select(species, floor_impact = impact_mean)
print(left_join(floor_imp, partial_sum, by = "species") %>%
        mutate(floor_minus_partial = floor_impact - partial_sum_singletons),
      width = 160)
