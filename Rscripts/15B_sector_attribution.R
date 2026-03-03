# ---
# title: Quantitative attribution of industrial sector effects on bird populations
# author: Mannfred Boehm
# ---
# For each species x BCR x year, partitions the total population change
# (observed vs counterfactual backfilled prediction) across industrial sectors
# using each sector's Hirsh-Pearson footprint score as a proportional weight.
# Pixels where a sector's score is < 1 are treated as unimpacted by that sector.
#
# Attribution method:
#   At each pixel, the weight for sector S is score_S / sum(score_S') where
#   the sum runs only over sectors with score >= 1. Pixels where no sector
#   reaches the threshold receive zero attribution to all sectors.
#
# Three spatial scales:
#   pixel    - mean and SD of attributed % change across active-HF pixels
#              (unweighted; measures typical pixel-level impact)
#   subbasin - mean and SD of population-weighted attributed % change across
#              subbasins within the BCR (each subbasin contributes one value)
#   BCR      - single population-weighted attributed % change across the full BCR
#
# Backfill scenario used: "mean" (backfilled_mean_mean.tif)
#
# Output: data/derived_data/rds_files/sector_attribution.csv
#   rows    = species x BCR x year
#   columns = species, bcr, year, {sector}_{pixel_mean, pixel_sd,
#             subbasin_mean, subbasin_sd, bcr_impact}
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

hirsh_dir  <- file.path(ia_dir, "data/raw_data/hirshpearson")
pred_dir   <- file.path(ia_dir, "data/derived_data/predictions")
basin_path <- file.path(ia_dir, "data/raw_data/hydrobasins_masked_merged_subset.gpkg")
out_path   <- file.path(ia_dir, "data/derived_data/rds_files/sector_attribution.csv")

# ---- Load sector rasters ----------------------------------------------------
# Exclude the aggregate CanHF masks; keep only the individual sector layers

sector_files <- list.files(hirsh_dir, pattern = "\\.tif$", full.names = TRUE)
sector_files <- sector_files[!grepl("^CanHF", basename(sector_files))]
sector_names <- tools::file_path_sans_ext(basename(sector_files))

sectors_raw        <- rast(sector_files)
names(sectors_raw) <- sector_names

message("Sectors (", length(sector_names), "): ", paste(sector_names, collapse = ", "))

# ---- Load hydrobasins -------------------------------------------------------

hydrobasins <- vect(basin_path)

# ---- Discover all valid species x BCR x year combinations ------------------
# A combination is valid if both observed_mean.tif and backfilled_mean_mean.tif exist.

combos <- do.call(rbind, Filter(Negate(is.null), lapply(
  list.dirs(pred_dir, full.names = FALSE, recursive = FALSE), function(sp) {
    do.call(rbind, Filter(Negate(is.null), lapply(
      list.dirs(file.path(pred_dir, sp), full.names = FALSE, recursive = FALSE), function(bcr) {
        do.call(rbind, Filter(Negate(is.null), lapply(
          list.dirs(file.path(pred_dir, sp, bcr), full.names = FALSE, recursive = FALSE), function(yr) {
            obs_f <- file.path(pred_dir, sp, bcr, yr, "observed_mean.tif")
            bf_f  <- file.path(pred_dir, sp, bcr, yr, "backfilled_mean_mean.tif")
            if (file.exists(obs_f) && file.exists(bf_f)) {
              data.frame(species = sp, bcr = bcr, year = yr, stringsAsFactors = FALSE)
            }
          })))
      })))
  })))

if (is.null(combos) || nrow(combos) == 0) stop("No valid species x BCR x year combinations found.")
message("Found ", nrow(combos), " species x BCR x year combinations")

# ---- Main processing loop ---------------------------------------------------
# Outer loop over BCR x year so sector weights are precomputed once per grid.
# Inner loop applies the shared weights to each species' density rasters.

bcr_year_combos <- unique(combos[, c("bcr", "year")])
results         <- vector("list", nrow(combos))
result_idx      <- 1

for (i in seq_len(nrow(bcr_year_combos))) {

  cur_bcr  <- bcr_year_combos$bcr[i]
  cur_year <- bcr_year_combos$year[i]
  sp_list  <- combos$species[combos$bcr == cur_bcr & combos$year == cur_year]

  message("BCR=", cur_bcr, " year=", cur_year, " (", length(sp_list), " species)")

  # Spatial template: the prediction grid for this BCR x year (same for all species)
  template_r <- rast(file.path(pred_dir, sp_list[1], cur_bcr, cur_year, "observed_mean.tif"))

  # ---- Precompute sector weights (shared across all species in this BCR x year) ----

  message("  Projecting sector rasters to BCR grid...")
  sectors_proj   <- project(sectors_raw, template_r, method = "bilinear")

  # Zero out sector scores below the significance threshold
  sectors_thresh <- ifel(sectors_proj >= 1, sectors_proj, 0)

  # Total active footprint score per pixel across all sectors
  total_score <- sum(sectors_thresh)

  # Proportional weights: score_S / total (0 where no sector is active)
  weight_stack <- ifel(total_score > 0,
                       sectors_thresh / total_score,
                       sectors_thresh * 0)  # sectors_thresh * 0 preserves NA structure

  # Pixel area in m^2 for population-weighted aggregation
  area_r <- cellSize(template_r, unit = "m")

  # Mask for pixels where at least one sector is active (score >= 1)
  hf_mask <- ifel(total_score > 0, 1, NA)

  # Rasterize subbasin IDs onto the BCR prediction grid (crop first for speed)
  message("  Rasterizing subbasins...")
  basins_clip <- crop(hydrobasins, ext(template_r))
  subbasin_r  <- rasterize(basins_clip, template_r, field = "HYBAS_ID")

  # ---- Inner loop: species --------------------------------------------------

  for (sp in sp_list) {
    message("    ", sp)

    obs_r <- rast(file.path(pred_dir, sp, cur_bcr, cur_year, "observed_mean.tif"))
    bf_r  <- rast(file.path(pred_dir, sp, cur_bcr, cur_year, "backfilled_mean_mean.tif"))

    # Mask non-positive backfill values to avoid division errors
    bf_r <- ifel(bf_r <= 0, NA, bf_r)

    # Total pixel-level % change: (observed - backfilled) / backfilled * 100
    # Negative values indicate industry has reduced the population below baseline.
    delta_r <- obs_r - bf_r
    pct_r   <- (delta_r / bf_r) * 100

    # Attributed % change per sector per pixel (12-layer raster)
    attr_stack <- pct_r * weight_stack

    # ---- Pixel-level: mean and SD over active-HF pixels --------------------
    # Restricted to pixels where at least one sector is active.
    attr_hf   <- mask(attr_stack, hf_mask)
    pix_stats <- global(attr_hf, fun = c("mean", "sd"), na.rm = TRUE)
    # pix_stats: data.frame [n_sectors x 2], rows named by sector, cols "mean" and "sd"

    # ---- BCR-level: single population-weighted aggregate -------------------
    # sum(delta * w_S * area) / sum(bf * area) * 100  for each sector S
    num_bcr   <- global(delta_r * weight_stack * area_r, fun = "sum", na.rm = TRUE)$sum
    denom_bcr <- as.numeric(global(bf_r * area_r, fun = "sum", na.rm = TRUE))
    bcr_impacts <- (num_bcr / denom_bcr) * 100
    # num_bcr is a length-12 vector (one sum per sector layer); denom_bcr is scalar

    # ---- Subbasin-level: per-subbasin aggregate, summarised across subbasins ----
    # For each subbasin b and sector S: sum(delta * w_S * area) / sum(bf * area) * 100
    num_zonal   <- zonal(delta_r * weight_stack * area_r, subbasin_r, fun = "sum", na.rm = TRUE)
    denom_zonal <- zonal(bf_r * area_r,                   subbasin_r, fun = "sum", na.rm = TRUE)
    # num_zonal:   data.frame [n_subbasins x (1 + n_sectors)] — zone ID + one col per sector
    # denom_zonal: data.frame [n_subbasins x 2]               — zone ID + denominator sum

    sub_pct <- sweep(as.matrix(num_zonal[, -1]), 1, denom_zonal[, 2], "/") * 100
    sub_pct[!is.finite(sub_pct)] <- NA  # subbasins with zero backfill population

    subbasin_means <- colMeans(sub_pct, na.rm = TRUE)
    subbasin_sds   <- apply(sub_pct, 2, sd, na.rm = TRUE)

    # ---- Assemble output row ------------------------------------------------

    row <- data.frame(species = sp, bcr = cur_bcr, year = cur_year, stringsAsFactors = FALSE)

    for (s_idx in seq_along(sector_names)) {
      s <- sector_names[s_idx]
      row[[paste0(s, "_pixel_mean")]]    <- pix_stats$mean[s_idx]
      row[[paste0(s, "_pixel_sd")]]      <- pix_stats$sd[s_idx]
      row[[paste0(s, "_subbasin_mean")]] <- subbasin_means[s_idx]
      row[[paste0(s, "_subbasin_sd")]]   <- subbasin_sds[s_idx]
      row[[paste0(s, "_bcr_impact")]]    <- bcr_impacts[s_idx]
    }

    results[[result_idx]] <- row
    result_idx <- result_idx + 1
  }
}

# ---- Save output ------------------------------------------------------------

output <- bind_rows(results)
write.csv(output, out_path, row.names = FALSE)
message("Saved sector attribution table (", nrow(output), " rows x ",
        ncol(output), " cols) to:\n  ", out_path)
