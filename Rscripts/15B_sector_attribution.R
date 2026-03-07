# ---
# title: Isolated sector attribution of bird population impacts
# author: Mannfred Boehm
# ---
# For each sector, 12A/12B was re-run with biotic covariates backfilled only at
# that sector's pixels (sector > 0 AND CanHF >= 1). The resulting backfilled
# raster differs from the observed raster ONLY at those pixels, giving an
# isolated sector counterfactual.
#
# This script reads the per-sector prediction rasters and computes % impact at
# four spatial scales:
#
#   BCR      - (sum(bf - obs) * area) / (sum(obs) * area) * 100  over BCR
#   national - same numerator summed across BCRs / national obs total
#   subbasin - per-subbasin version of BCR formula
#   footprint- per-connected-component (8-connectivity) of sector pixels
#              denominator = observed population within that footprint
#
# Uncertainty (SD) is propagated via standard error propagation.
# Only the "mean" backfill scenario is used (backfilled_mean_mean.tif).
#
# Outputs (data/derived_data/rds_files/):
#   sector_bcr.csv, sector_national.csv, sector_subbasin.csv, sector_footprint.csv
# ---

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
})

# ---- Execution context -------------------------------------------------------

cc    <- FALSE
local <- TRUE

if (cc)            { ia_dir <- "/home/mannfred/scratch/impact_assessment" }
if (!cc && local)  { ia_dir <- getwd() }
if (!cc && !local) { ia_dir <- file.path("G:/Shared drives/BAM_NationalModels5", "data", "Extras",
                                          "sandbox_data", "impactassessment_sandbox") }

# ---- Paths ------------------------------------------------------------------

hirsh_dir   <- file.path(ia_dir, "data/raw_data/hirshpearson")
obs_root    <- file.path(ia_dir, "data/derived_data/predictions")
bf_root     <- file.path(ia_dir, "data/derived_data/predictions_sectors")
basin_path  <- file.path(ia_dir, "data/raw_data/hydrobasins_masked_merged_subset.gpkg")
out_dir     <- file.path(ia_dir, "data/derived_data/sector_effects")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

canHF_path  <- file.path(hirsh_dir, "CanHF_1km_morethan1.tif")

# ---- Load sector names -------------------------------------------------------

sector_files <- list.files(hirsh_dir, pattern = "\\.tif$", full.names = FALSE)
sector_files <- sector_files[!grepl("^CanHF", sector_files)]
sector_names <- tools::file_path_sans_ext(sector_files)
message("Sectors (", length(sector_names), "): ", paste(sector_names, collapse = ", "))

# ---- Load hydrobasins --------------------------------------------------------

hydrobasins <- vect(basin_path)

# discover valid sector x species x BCR x year combinations --------------
# A combination is valid if backfilled_mean_mean.tif exists and the matching
# observed_mean.tif exists in the shared predictions/ directory.
# `combos` is a dataframe where rows are unique sector x species x BCR x year tuples

combos <- do.call(rbind, Filter(Negate(is.null), lapply(sector_names, function(sec) {
  sec_dir <- file.path(bf_root, sec)
  if (!dir.exists(sec_dir)) return(NULL)
  do.call(rbind, Filter(Negate(is.null), lapply(
    list.dirs(sec_dir, full.names = FALSE, recursive = FALSE), function(sp) {
      do.call(rbind, Filter(Negate(is.null), lapply(
        list.dirs(file.path(sec_dir, sp), full.names = FALSE, recursive = FALSE), function(bcr) {
          do.call(rbind, Filter(Negate(is.null), lapply(
            list.dirs(file.path(sec_dir, sp, bcr), full.names = FALSE, recursive = FALSE), function(yr) {
              bf_f  <- file.path(bf_root,  sec, sp, bcr, yr, "backfilled_mean_mean.tif")
              obs_f <- file.path(obs_root, sp,  bcr, yr,     "observed_mean.tif")
              if (file.exists(bf_f) && file.exists(obs_f)) {
                data.frame(sector = sec, species = sp, bcr = bcr, year = yr,
                           stringsAsFactors = FALSE)
              }
            })))
        })))
    })))
})))

if (is.null(combos) || nrow(combos) == 0) stop("No valid sector x species x BCR x year combinations found.")
message("found ", nrow(combos), "unique sector x species x BCR x year tuples")

# main processing loop ---------------------------------------------------
# outer loop: sector x BCR x year  (project sector mask once per BCR grid)
# inner loop: species

bcr_year_sector_combos <- unique(combos[, c("sector", "bcr", "year")])

bcr_rows      <- vector("list", nrow(combos))
subbasin_rows <- vector("list", nrow(combos))
footprint_rows <- vector("list", nrow(combos))
result_idx    <- 1

for (i in seq_len(nrow(bcr_year_sector_combos))) {

  cur_sector <- bcr_year_sector_combos$sector[i]
  cur_bcr    <- bcr_year_sector_combos$bcr[i]
  cur_year   <- bcr_year_sector_combos$year[i]
  sp_list    <- combos$species[combos$sector == cur_sector &
                               combos$bcr    == cur_bcr    &
                               combos$year   == cur_year]

  message("sector=", cur_sector, " BCR=", cur_bcr, " year=", cur_year,
          " (", length(sp_list), " species)")

  # use first species as spatial template (observed grid is identical across species)
  template_r <- rast(file.path(obs_root, sp_list[1], cur_bcr, cur_year, "observed_mean.tif"))
  area_r     <- cellSize(template_r, unit = "ha")

  # sector mask on BCR grid
  sector_r  <- project(rast(file.path(hirsh_dir, paste0(cur_sector, ".tif"))),
                       template_r, method = "near")
  canHF_r   <- project(rast(canHF_path), template_r, method = "near")
  sector_mask <- ifel((sector_r > 0) & (canHF_r >= 1), 1, NA)

  # connected-component patches for footprint scale
  patches_r <- patches(sector_mask, directions = 8, zeroAsNA = TRUE)

  for (sp in sp_list) {

    message("  ", sp)

    # ---- Load 12B density tables ---------------------------------------------
    # Per-subbasin bootstrap-aggregated population counts: SDs here reflect
    # the actual bootstrap distribution (32 boots x 3 scenarios), not pixel-level
    # independence assumptions.
    dt_path <- file.path(ia_dir, "data/derived_data/density_tables",
                         paste0(sp, "_", cur_year, "_", cur_sector, ".rds"))
    dt_all_hf_path <- file.path(ia_dir, "data/derived_data/density_tables",
                                paste0(sp, "_", cur_year, "_all_hf.rds"))

    if (!file.exists(dt_path) || !file.exists(dt_all_hf_path)) {
      message("  density table missing for ", sp, " — skipping")
      result_idx <- result_idx + 1
      next
    }

    dt        <- dplyr::filter(readRDS(dt_path),        bcr == cur_bcr)
    dt_all_hf <- dplyr::filter(readRDS(dt_all_hf_path), bcr == cur_bcr)

    if (nrow(dt) == 0 || nrow(dt_all_hf) == 0) {
      message("  no rows for BCR ", cur_bcr, " — skipping")
      result_idx <- result_idx + 1
      next
    }

    # ---- BCR scale -----------------------------------------------------------
    # All BCR-level quantities come directly from the bootstrap x scenario
    # distributions pre-computed in 12B. Take [1] since BCR-level values are
    # identical for every subbasin row within a BCR.

    obs_population_mean <- dt$obs_total_bcr_mean[1]  # BCR-wide; BRT bootstrap SD only
    obs_population_sd   <- dt$obs_total_bcr_sd[1]

    # all-HF counterfactual from the all_hf density table
    HF_impact_mean     <- dt_all_hf$sector_impact_bcr_mean[1]  # bootstrap + BART posterior SD
    HF_impact_sd       <- dt_all_hf$sector_impact_bcr_sd[1]
    cf_population_mean <- dt_all_hf$cf_total_bcr_mean[1]
    cf_population_sd   <- dt_all_hf$cf_total_bcr_sd[1]

    HF_percent_impact_mean <- HF_impact_mean / obs_population_mean * 100
    HF_percent_impact_sd   <- 100 * sqrt((HF_impact_sd / obs_population_mean)^2 +
                                           (HF_impact_mean * obs_population_sd / obs_population_mean^2)^2)

    # sector-specific impact from sector density table (bootstrap + BART posterior SD)
    sector_impact_mean        <- dt$sector_impact_bcr_mean[1]
    sector_impact_sd          <- dt$sector_impact_bcr_sd[1]
    cf_sector_population_mean <- dt$cf_total_bcr_mean[1]
    cf_sector_population_sd   <- dt$cf_total_bcr_sd[1]

    sector_percent_impact_mean <- sector_impact_mean / obs_population_mean * 100
    sector_percent_impact_sd   <- 100 * sqrt((sector_impact_sd / obs_population_mean)^2 +
                                               (sector_impact_mean * obs_population_sd / obs_population_mean^2)^2)

    # BCR-wide footprint populations (bootstrap SD only for obs; bootstrap + posterior for cf)
    obs_population_on_footprint_mean <- dt$obs_on_sector_bcr_mean[1]
    obs_population_on_footprint_sd   <- dt$obs_on_sector_bcr_sd[1]
    cf_population_on_footprint_mean  <- dt$bf_on_sector_bcr_mean[1]
    cf_population_on_footprint_sd    <- dt$bf_on_sector_bcr_sd[1]

    bcr_rows[[result_idx]] <- data.frame(
      species                          = sp,
      bcr                              = cur_bcr,
      year                             = cur_year,
      sector                           = cur_sector,
      obs_population_mean              = round(obs_population_mean),
      obs_population_sd                = round(obs_population_sd, 1),
      cf_population_mean               = round(cf_population_mean),
      cf_population_sd                 = round(cf_population_sd, 1),
      HF_impact_mean                   = round(HF_impact_mean),
      HF_impact_sd                     = round(HF_impact_sd, 1),
      HF_percent_impact_mean           = round(HF_percent_impact_mean, 2),
      HF_percent_impact_sd             = round(HF_percent_impact_sd, 2),
      cf_sector_population_mean        = round(cf_sector_population_mean),
      cf_sector_population_sd          = round(cf_sector_population_sd, 1),
      sector_impact_mean               = round(sector_impact_mean),
      sector_impact_sd                 = round(sector_impact_sd, 1),
      sector_percent_impact_mean       = round(sector_percent_impact_mean, 2),
      sector_percent_impact_sd         = round(sector_percent_impact_sd, 2),
      obs_population_on_footprint_mean = round(obs_population_on_footprint_mean),
      obs_population_on_footprint_sd   = round(obs_population_on_footprint_sd, 1),
      cf_population_on_footprint_mean  = round(cf_population_on_footprint_mean),
      cf_population_on_footprint_sd    = round(cf_population_on_footprint_sd, 1),
      stringsAsFactors = FALSE
    )

    # ---- Subbasin scale ------------------------------------------------------
    # Use density table directly; map subbasin row-index back to HYBAS_ID.

    sub_hybas <- hydrobasins$first_HYBAS_ID[dt$subbasin]
    sect_imp  <- dt$bf_only_mean - dt$obs_on_bf_mean
    pct_sub   <- sect_imp / dt$obs_total_mean * 100
    pct_sub[!is.finite(pct_sub) | dt$obs_total_mean < 1e-6] <- NA

    subbasin_rows[[result_idx]] <- data.frame(
      species                   = sp,
      bcr                       = cur_bcr,
      year                      = cur_year,
      sector                    = cur_sector,
      subbasin_id               = sub_hybas,
      obs_population_sub        = round(dt$obs_total_mean),
      sector_impact_sub         = round(sect_imp),
      sector_percent_impact_sub = round(pct_sub, 2),
      stringsAsFactors = FALSE
    )

    # ---- Footprint scale -----------------------------------------------------
    # Still raster-based: connected components don't align with subbasins.

    obs_mean <- rast(file.path(obs_root, sp, cur_bcr, cur_year, "observed_mean.tif"))
    bf_mean  <- rast(file.path(bf_root, cur_sector, sp, cur_bcr, cur_year, "backfilled_mean_mean.tif"))
    delta    <- bf_mean - obs_mean

    sector_impact_footprint  <- zonal(delta    * area_r, patches_r, "sum", na.rm = TRUE)
    obs_population_footprint <- zonal(obs_mean * area_r, patches_r, "sum", na.rm = TRUE)
    cf_population_footprint  <- zonal(bf_mean  * area_r, patches_r, "sum", na.rm = TRUE)
    names(sector_impact_footprint)[2]  <- "sector_impact_footprint"
    names(obs_population_footprint)[2] <- "obs_population_on_footprint"
    names(cf_population_footprint)[2]  <- "cf_population_on_footprint"

    fp_tbl <- merge(obs_population_footprint, cf_population_footprint, by = names(patches_r))
    fp_tbl <- merge(fp_tbl, sector_impact_footprint, by = names(patches_r))
    fp_tbl$sector_percent_impact_footprint <- fp_tbl$sector_impact_footprint /
                                              fp_tbl$obs_population_on_footprint * 100
    fp_tbl$sector_percent_impact_footprint[!is.finite(fp_tbl$sector_percent_impact_footprint) |
                                            fp_tbl$obs_population_on_footprint < 1e-6] <- NA

    footprint_rows[[result_idx]] <- data.frame(
      species                         = sp,
      bcr                             = cur_bcr,
      year                            = cur_year,
      sector                          = cur_sector,
      footprint_id                    = fp_tbl[[1]],
      obs_population_on_footprint     = round(fp_tbl$obs_population_on_footprint),
      cf_population_on_footprint      = round(fp_tbl$cf_population_on_footprint),
      sector_impact_footprint         = round(fp_tbl$sector_impact_footprint),
      sector_percent_impact_footprint = round(fp_tbl$sector_percent_impact_footprint, 2),
      stringsAsFactors = FALSE
    )

    result_idx <- result_idx + 1

  } # close species loop
} # close outer loop

# ---- Assemble BCR table ------------------------------------------------------

bcr_df <- bind_rows(bcr_rows)

# ---- National scale (derived from BCR) ---------------------------------------

national_df <- bcr_df |>
  group_by(species, year, sector) |>
  summarise(
    obs_population_mean        = round(sum(obs_population_mean)),
    obs_population_sd          = round(sqrt(sum(obs_population_sd^2)), 1),
    cf_population_mean         = round(sum(cf_population_mean)),
    cf_population_sd           = round(sqrt(sum(cf_population_sd^2)), 1),
    HF_impact_mean             = round(sum(HF_impact_mean)),
    HF_impact_sd               = round(sqrt(sum(HF_impact_sd^2)), 1),
    HF_percent_impact_mean     = round(HF_impact_mean / obs_population_mean * 100, 2),
    HF_percent_impact_sd       = round(100 * sqrt((HF_impact_sd / obs_population_mean)^2 +
                                         (HF_impact_mean * obs_population_sd / obs_population_mean^2)^2), 2),
    cf_sector_population_mean  = round(sum(cf_sector_population_mean)),
    cf_sector_population_sd    = round(sqrt(sum(cf_sector_population_sd^2)), 1),
    sector_impact_mean         = round(sum(sector_impact_mean)),
    sector_impact_sd           = round(sqrt(sum(sector_impact_sd^2)), 1),
    sector_percent_impact_mean = round(sector_impact_mean / obs_population_mean * 100, 2),
    sector_percent_impact_sd   = round(100 * sqrt((sector_impact_sd / obs_population_mean)^2 +
                                         (sector_impact_mean * obs_population_sd / obs_population_mean^2)^2), 2),
    .groups = "drop"
  )

# ---- Write outputs -----------------------------------------------------------

write.csv(bcr_df,             file.path(out_dir, "sector_bcr.csv"),       row.names = FALSE)
write.csv(national_df,        file.path(out_dir, "sector_national.csv"),  row.names = FALSE)
write.csv(bind_rows(subbasin_rows),  file.path(out_dir, "sector_subbasin.csv"),  row.names = FALSE)
write.csv(bind_rows(footprint_rows), file.path(out_dir, "sector_footprint.csv"), row.names = FALSE)

message("writing 4 tables to ", out_dir)
message("  sector_bcr.csv       : ", nrow(bcr_df), " rows")
message("  sector_national.csv  : ", nrow(national_df), " rows")
message("  sector_subbasin.csv  : ", nrow(bind_rows(subbasin_rows)), " rows")
message("  sector_footprint.csv : ", nrow(bind_rows(footprint_rows)), " rows")
