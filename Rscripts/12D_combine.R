# ---
# title: Impact Assessment: combine bf-only tables with observed predictions
# author: Mannfred Boehm
# ---
# Run AFTER 12A_observed.R, for species that had 12B coalition jobs complete
# before observed_bootstraps.tif files were available.
#
# For each non-empty coalition that has a *_bf_only.rds density table, reads the
# canonical observed bootstraps produced by 12A (one per species x BCR), rebuilds
# the coalition sector mask, computes obs_total and obs_on_coalition per subbasin,
# and overwrites the bf_only table with the complete density table.
#
# SLURM array over species (one task per species). Iterates all 255 non-empty
# coalitions internally — workload is zonal stats only, no GBM predictions.
# ---

suppressPackageStartupMessages({
  library(terra)
  library(tidyverse)
})

# ---- Paths -------------------------------------------------------------------

nm_root <- "/home/mannfred/projects/def-ecknight/NationalModels"

cc    <- TRUE
local <- FALSE

if (cc)            { ia_dir <- "/home/mannfred/scratch/impact_assessment" }
if (!cc && local)  { ia_dir <- getwd() }
if (!cc && !local) { ia_dir <- file.path("G:/Shared drives/BAM_NationalModels5", "data", "Extras",
                                          "sandbox_data", "impactassessment_sandbox") }

# ---- Import shared spatial data ----------------------------------------------

bam_boundary         <- terra::vect(file.path(ia_dir, "data", "raw_data", "Regions",
                                               "BAM_BCR_NationalModel_Unbuffered.shp"))
all_subbasins_subset <- terra::vect(file.path(ia_dir, "data", "raw_data",
                                               "hydrobasins_masked_merged_subset.gpkg"))

bcr_subbasins_ref <- {
  hits <- terra::relate(all_subbasins_subset, bam_boundary, relation = "intersects")
  ij   <- which(hits, arr.ind = TRUE)
  tibble(
    sub_index = ij[, 1],
    HYBAS_ID  = all_subbasins_subset$first_HYBAS_ID[ij[, 1]],
    bcr_label = paste(bam_boundary$country[ij[, 2]],
                      bam_boundary$subUnit[ij[, 2]], sep = "_"),
    bcr_code  = gsub("_", "", bcr_label)
  )
}

hirsh_dir <- file.path(ia_dir, "data", "raw_data", "hirshpearson")

source(file.path(ia_dir, "Rscripts", "12E_shapley_utils.R"))

# ---- Species from SLURM ------------------------------------------------------

species_vec <- c("CAWA", "OVEN")
# species_vec <- sort(c("BANS", "BARS", "BOBO", "CAWA", "EAWP", "EVGR", "GCTH", "GRSP", "GWWA", "LEYE", "OSFL"))
year    <- 2020
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
species <- species_vec[task_id]
message(Sys.time(), " | 12D combine for species=", species)

sectors <- canonical_sectors()
dt_dir  <- file.path(ia_dir, "data", "derived_data", "density_tables")

# ---- Loop over all non-empty coalitions --------------------------------------

n_coalitions <- 2L^length(sectors) - 1L  # 255 for 8 sectors
n_combined   <- 0L
n_skipped    <- 0L

for (coalition_id in seq_len(n_coalitions)) {

  coalition <- coalition_id_to_sectors(coalition_id, sectors)
  if (length(coalition) == 0) next

  bf_only_path <- file.path(dt_dir,
                             paste0(species, "_", year, "_coalition_", coalition_id, "_bf_only.rds"))
  if (!file.exists(bf_only_path)) {
    n_skipped <- n_skipped + 1L
    next
  }

  dt        <- readRDS(bf_only_path)
  bcr_codes <- unique(dt$bcr)
  message(Sys.time(), " | coalition=", coalition_id,
          " | bf_only table: ", nrow(dt), " rows, ", length(bcr_codes), " BCRs")

  # ---- Fill obs columns per BCR --------------------------------------------

  for (bcr_code in bcr_codes) {

    obs_boot_path <- file.path(ia_dir, "data", "derived_data", "predictions",
                               species, bcr_code, year, "observed_bootstraps.tif")
    if (!file.exists(obs_boot_path)) {
      message(Sys.time(), " | coalition=", coalition_id, " | ", bcr_code,
              " | observed_bootstraps.tif not found — skipping BCR")
      next
    }

    # BCR reference stack (spatial grid)
    stack_obs <- terra::rast(file.path(nm_root, "gis/stacks",
                                        paste0(bcr_code, "_", year, ".tif")))

    # rebuild coalition sector mask for this BCR
    canHF_r    <- terra::project(terra::rast(file.path(hirsh_dir, "CanHF_1km_morethan1.tif")),
                                  stack_obs, method = "near")
    union_mask <- terra::rast(stack_obs[[1]]); terra::values(union_mask) <- 0L
    for (sec in coalition) {
      sec_r      <- terra::project(terra::rast(file.path(hirsh_dir, paste0(sec, ".tif"))),
                                    stack_obs, method = "near")
      union_mask <- terra::ifel(sec_r > 0, 1L, union_mask)
    }
    sector_mask <- terra::ifel((union_mask == 1L) & (canHF_r >= 1), 1, NA)
    rm(canHF_r, union_mask); gc()

    # subbasin zones for this BCR
    sub_ids <- bcr_subbasins_ref |>
      dplyr::filter(bcr_code == !!bcr_code) |>
      dplyr::pull(sub_index) |>
      unique()

    subbasin_zone_r <- terra::rasterize(all_subbasins_subset[sub_ids, ], stack_obs[[1]],
                                         field = "first_HYBAS_ID")
    zone_ids  <- sort(unique(terra::values(subbasin_zone_r, na.rm = TRUE)))
    hybas_ids <- all_subbasins_subset$first_HYBAS_ID[sub_ids]
    idx       <- match(hybas_ids, zone_ids)

    # load observed bootstraps and compute per-subbasin zonal sums
    obs_boot_stack <- terra::rast(obs_boot_path)
    n_boot         <- terra::nlyr(obs_boot_stack)
    n_scen         <- 100L  # must match 12C

    obs_total_mat   <- matrix(NA_real_, nrow = length(sub_ids), ncol = n_boot)
    obs_on_coal_mat <- matrix(NA_real_, nrow = length(sub_ids), ncol = n_boot)

    for (i in seq_len(n_boot)) {
      obs_r     <- obs_boot_stack[[i]] * 100
      obs_on_fp <- terra::mask(obs_r, sector_mask)
      sub_obs_total <- terra::zonal(obs_r,     subbasin_zone_r, "sum", na.rm = TRUE)
      sub_obs_on_fp <- terra::zonal(obs_on_fp, subbasin_zone_r, "sum", na.rm = TRUE)
      obs_total_mat[, i]   <- sub_obs_total[[2]][idx]
      obs_on_coal_mat[, i] <- sub_obs_on_fp[[2]][idx]
    }
    rm(obs_boot_stack, stack_obs, sector_mask, subbasin_zone_r); gc()

    # replicate bootstrap values across n_scen to match 12C agg() SD computation:
    # in the full pipeline obs values fill an n_sub x n_boot x n_scen array where
    # each bootstrap value is identical across all scenarios.
    obs_total_rep   <- obs_total_mat[,   rep(seq_len(n_boot), n_scen), drop = FALSE]
    obs_on_coal_rep <- obs_on_coal_mat[, rep(seq_len(n_boot), n_scen), drop = FALSE]

    obs_total_mean   <- rowMeans(obs_total_rep,   na.rm = TRUE)
    obs_total_sd     <- apply(obs_total_rep,   1, sd, na.rm = TRUE)
    obs_on_coal_mean <- rowMeans(obs_on_coal_rep, na.rm = TRUE)
    obs_on_coal_sd   <- apply(obs_on_coal_rep, 1, sd, na.rm = TRUE)

    # update matching rows in table
    rows <- dt$bcr == bcr_code
    dt$obs_total_mean[rows]        <- obs_total_mean
    dt$obs_total_sd[rows]          <- obs_total_sd
    dt$obs_on_coalition_mean[rows] <- obs_on_coal_mean
    dt$obs_on_coalition_sd[rows]   <- obs_on_coal_sd

    message(Sys.time(), " | coalition=", coalition_id, " | ", bcr_code,
            " | obs columns filled (", sum(rows), " subbasins)")
  }

  # ---- Save complete table (or warn if some BCRs still missing) ------------

  incomplete <- unique(dt$bcr[is.na(dt$obs_total_mean)])

  final_path <- file.path(dt_dir,
                           paste0(species, "_", year, "_coalition_", coalition_id, ".rds"))
  saveRDS(dt, final_path)

  if (length(incomplete) == 0) {
    file.remove(bf_only_path)
    message(Sys.time(), " | coalition=", coalition_id, " | complete — bf_only file removed")
    n_combined <- n_combined + 1L
  } else {
    message(Sys.time(), " | coalition=", coalition_id,
            " | WARNING: obs still NA for BCRs: ", paste(incomplete, collapse = ", "),
            " — bf_only file kept, final table written with NAs")
  }
}

message(Sys.time(), " | species=", species,
        " | combined=", n_combined,
        " | skipped (no bf_only)=", n_skipped,
        " | done.")
