# ---
# title: Abiotic extrapolation diagnostics for BART backfilling
# author: Mannfred Boehm
# ---
# For each subbasin, compare the abiotic covariate distributions of low-HF
# (training) pixels vs high-HF (backfill) pixels.  Flag subbasins where the
# BART model is likely extrapolating into abiotic space not represented in
# the training data.
#
# Diagnostics per subbasin x abiotic covariate:
#   - Kolmogorov-Smirnov D statistic (univariate distributional overlap)
#
# Diagnostics per subbasin (multivariate):
#   - Fraction of high-HF pixels whose Mahalanobis distance from the low-HF
#     centroid exceeds the 95th percentile of the low-HF Mahalanobis distribution
#
# Output: data/derived_data/rds_files/extrapolation_flags.csv
# ---

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(tidyr)
})

# ---- Execution context -------------------------------------------------------

cc    <- TRUE
local <- FALSE

if (cc)            { ia_dir <- "/home/mannfred/scratch/impact_assessment" }
if (!cc && local)  { ia_dir <- getwd() }
if (!cc && !local) { ia_dir <- file.path("G:/Shared drives/BAM_NationalModels5", "data", "Extras",
                                          "sandbox_data", "impactassessment_sandbox") }

# ---- Paths -------------------------------------------------------------------

lowhf_path  <- file.path(ia_dir, "data/raw_data/hirshpearson/CanHF_1km_lessthan1.tif")
highhf_path <- file.path(ia_dir, "data/raw_data/hirshpearson/CanHF_1km_morethan1.tif")
basin_path  <- file.path(ia_dir, "data/raw_data/hydrobasins_masked_merged_subset.gpkg")
out_dir     <- file.path(ia_dir, "data/derived_data/rds_files")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

year <- 2020

# ---- Load spatial data -------------------------------------------------------

all_subbasins <- vect(basin_path)
n_sub         <- nrow(all_subbasins)

lowhf_mask  <- rast(lowhf_path)
highhf_mask <- rast(highhf_path)

stack_y <- rast(file.path(ia_dir, "data/raw_data/covariates_mosaiced",
                          paste0("covariates_mosaiced_", year, ".tif")))

# ---- Define abiotic covariates -----------------------------------------------

predictor_metadata <-
  dplyr::tibble(BAMexploreR::predictor_metadata) |>
  dplyr::filter(version == "v5") |>
  dplyr::select(predictor, predictor_class) |>
  dplyr::mutate(dplyr::across('predictor', stringr::str_replace, 'Year', 'year')) |>
  dplyr::mutate(dplyr::across('predictor', stringr::str_replace, 'Method','method'))

actually_biotic <- c("Peatland_5x5", "Peatland_1km")

abiotic_preds <-
  predictor_metadata |>
  dplyr::filter(predictor_class %in% c("Annual Climate", "Climate Normals",
                                        "Topography", "Wetland", "Disturbance",
                                        "Time", "Method")) |>
  dplyr::filter(!(predictor %in% actually_biotic)) |>
  dplyr::pull(predictor)

# subset to covariates actually in the stack
abiotic_preds <- intersect(abiotic_preds, names(stack_y))
message("Abiotic covariates (", length(abiotic_preds), "): ",
        paste(abiotic_preds, collapse = ", "))

# ---- Process each subbasin ---------------------------------------------------

results <- vector("list", n_sub)

for (s in seq_len(n_sub)) {

  sub_s <- all_subbasins[s]

  # crop covariate stack to subbasin
  cov_s <- tryCatch(
    terra::mask(terra::crop(stack_y, sub_s), sub_s),
    error = function(e) NULL
  )
  if (is.null(cov_s)) next

  # build low-HF and high-HF masks for this subbasin
  lowhf_s  <- terra::mask(terra::resample(lowhf_mask,  cov_s, method = "near"), sub_s)
  highhf_s <- terra::mask(terra::resample(highhf_mask, cov_s, method = "near"), sub_s)

  # extract abiotic values at low-HF and high-HF pixels
  avail <- intersect(abiotic_preds, names(cov_s))
  if (length(avail) == 0) next

  vals_all <- terra::values(cov_s[[avail]], mat = TRUE)
  lo_idx   <- which(terra::values(lowhf_s,  mat = FALSE) == 1)
  hi_idx   <- which(terra::values(highhf_s, mat = FALSE) == 1)

  if (length(lo_idx) < 10 || length(hi_idx) < 2) next

  lo_mat <- vals_all[lo_idx, , drop = FALSE]
  hi_mat <- vals_all[hi_idx, , drop = FALSE]

  # --- Univariate: KS D statistic per covariate ---
  ks_d <- setNames(numeric(length(avail)), avail)
  for (v in avail) {
    lo_v <- lo_mat[, v]; lo_v <- lo_v[!is.na(lo_v)]
    hi_v <- hi_mat[, v]; hi_v <- hi_v[!is.na(hi_v)]
    if (length(lo_v) > 1 && length(hi_v) > 1) {
      ks_d[v] <- suppressWarnings(ks.test(lo_v, hi_v)$statistic)
    } else {
      ks_d[v] <- NA
    }
  }

  # --- Multivariate: Mahalanobis exceedance ---
  # remove columns with zero variance or all-NA in either set
  good_cols <- vapply(avail, function(v) {
    lo_v <- lo_mat[, v]; hi_v <- hi_mat[, v]
    if (all(is.na(lo_v)) || all(is.na(hi_v))) return(FALSE)
    s <- sd(lo_v, na.rm = TRUE)
    isTRUE(s > 0)
  }, logical(1))
  good_avail <- avail[good_cols]

  mahal_exceedance <- NA_real_

  if (length(good_avail) >= 2) {
    lo_complete <- lo_mat[, good_avail, drop = FALSE]
    hi_complete <- hi_mat[, good_avail, drop = FALSE]

    # remove rows with any NA
    lo_complete <- lo_complete[complete.cases(lo_complete), , drop = FALSE]
    hi_complete <- hi_complete[complete.cases(hi_complete), , drop = FALSE]

    if (nrow(lo_complete) > ncol(lo_complete) + 1 && nrow(hi_complete) > 0) {
      mu    <- colMeans(lo_complete)
      Sigma <- cov(lo_complete)

      # regularize if near-singular
      Sigma <- Sigma + diag(1e-6, ncol(Sigma))

      tryCatch({
        Sigma_inv <- solve(Sigma)

        mahal_lo <- mahalanobis(lo_complete, mu, Sigma, inverted = FALSE)
        q95      <- quantile(mahal_lo, 0.95)

        mahal_hi <- mahalanobis(hi_complete, mu, Sigma, inverted = FALSE)
        mahal_exceedance <- mean(mahal_hi > q95)
      }, error = function(e) {
        mahal_exceedance <<- NA_real_
      })
    }
  }

  results[[s]] <- data.frame(
    subbasin        = s,
    HYBAS_ID        = all_subbasins$first_HYBAS_ID[s],
    n_lowhf         = length(lo_idx),
    n_highhf        = length(hi_idx),
    ks_max          = max(ks_d, na.rm = TRUE),
    ks_median       = median(ks_d, na.rm = TRUE),
    mahal_exceedance = mahal_exceedance,
    stringsAsFactors = FALSE
  )

  if (s %% 50 == 0) message("  processed ", s, " / ", n_sub, " subbasins")
}

# ---- Assemble and flag -------------------------------------------------------

flags_df <- bind_rows(results)

# flag subbasins where extrapolation risk is elevated:
#   KS max > 0.5 (large distributional shift in at least one covariate)
#   OR Mahalanobis exceedance > 0.3 (>30% of backfill pixels outside training 95th pct)
flags_df$flag <- with(flags_df,
  (ks_max > 0.5) | (!is.na(mahal_exceedance) & mahal_exceedance > 0.3)
)

out_path <- file.path(out_dir, "extrapolation_flags.csv")
write.csv(flags_df, out_path, row.names = FALSE)

message("Wrote ", nrow(flags_df), " subbasins to ", out_path)
message("  Flagged: ", sum(flags_df$flag, na.rm = TRUE), " / ", nrow(flags_df))
