# ---
# title: Abiotic extrapolation diagnostics for BART backfilling
# author: Mannfred Boehm
# created: May 19, 2026 (rewritten 2026-09-25)
# ---

# For each subbasin: do the high-HF (backfilled) pixels lie inside the region of abiotic
# predictor space that its low-HF (training) pixels cover? This is the area of
# applicability (AOA) of Meyer & Pebesma (2021, Methods Ecol Evol 12:1620), computed on the
# natural abiotic predictors 08A's BART models are trained on: 07's `abiotic_vars` (climate,
# phenology, topography, wetland, CAfire, soil) WITHOUT its "Disturbance" class.
#
# The Disturbance class (CanHF_1km/_5x5, canroad_1km/_5x5, CCNL_1km night lights) is human
# footprint. 07 passes it to BART, and 08A backfills with the pixels' OBSERVED values, so
# every backfilled pixel lies outside the training range on it by construction (training
# pixels are CanHF < 1, backfilled ones CanHF >= 1). Left in, it would put nearly every
# backfilled pixel outside the AOA and say nothing about the abiotic overlap this script is
# for. Whether BART should see footprint covariates at all is a separate question about
# 08A (TODO.md C2c).
#
#   DI(x)     = distance from x to its nearest training pixel / mean distance between
#               training pixels, in predictor space standardised by the training SDs.
#   threshold = upper whisker (Q75 + 1.5 IQR) of the training pixels' own DI, each measured
#               to the nearest training pixel in a DIFFERENT spatial fold (k-means blocks on
#               x, y). Spatial folds, not random ones: neighbouring pixels are near-duplicates,
#               so random folds would put every training DI near 0 and every backfilled pixel
#               outside. Blocks mimic the backfill, whose pixels sit in clusters away from
#               the training pixels.
#   frac_outside_aoa = share of the subbasin's high-HF pixels with DI > threshold.
#   flag      = frac_outside_aoa > 0.5: most of the subbasin's backfill is extrapolated.
#
# Predictors a subbasin's training pixels hold constant are dropped: BART cannot split on
# them, so their high-HF values cannot move a prediction. NAs take the training median,
# as 08A's imputation does.
#
# This replaced a KS + Mahalanobis rule that flagged 667 of 667 subbasins (C2, 2026-09-24):
#   - KS measures distribution SHIFT, not extrapolation. Low- and high-HF pixels in a
#     subbasin are spatially segregated and climate normals are smooth, so some normal
#     nearly always has near-disjoint samples (min ks_max was 0.52), however small the
#     difference in degrees.
#   - Mahalanobis distance over ~40 collinear climate covariates inverts a near-singular
#     covariance, so a tiny shift along a minor axis becomes a huge distance (median
#     exceedance 0.93).
#   - Its covariate set was not the BART models': no soil and no CAfire, but Time and
#     Method, which 07 excludes.
#
# Runs locally (every input is local, and FNN is installed here, not on Fir).
# output: data/derived_data/rds_files/extrapolation_flags.csv


suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(FNN)
})

# define paths -------------------------------------------------------

cc    <- FALSE
local <- TRUE

if (cc)            { ia_dir <- "/home/mannfred/scratch/impact_assessment" }
if (!cc && local)  { ia_dir <- getwd() }
if (!cc && !local) { ia_dir <- file.path("G:/Shared drives/BAM_NationalModels5", "data", "Extras",
                                          "sandbox_data", "impactassessment_sandbox") }

lowhf_path  <- file.path(ia_dir, "data/raw_data/hirshpearson/CanHF_1km_lessthan1.tif")
highhf_path <- file.path(ia_dir, "data/raw_data/hirshpearson/CanHF_1km_morethan1.tif")
basin_path  <- file.path(ia_dir, "data/raw_data/hydrobasins_masked_merged_subset.gpkg")
out_dir     <- file.path(ia_dir, "data/derived_data/rds_files")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

year <- 2020

n_folds      <- 10L      # spatial blocks for the training DI
max_query    <- 5000L    # training pixels whose DI sets the threshold (sampled above this)
n_pairs      <- 1e5      # sampled pairs for the mean training distance
flag_cut     <- 0.5
set.seed(20260925)

# load spatial data ---------------------------------------------------

all_subbasins <- vect(basin_path)
n_sub         <- nrow(all_subbasins)

lowhf_mask  <- rast(lowhf_path)
highhf_mask <- rast(highhf_path)

stack_y <- rast(file.path(ia_dir, "data/raw_data/covariates_mosaiced",
                          paste0("covariates_mosaiced_", year, ".tif")))

# abiotic predictors: 07's `abiotic_vars` minus its "Disturbance" class ----

predictor_metadata <-
  dplyr::tibble(BAMexploreR::predictor_metadata) |>
  dplyr::filter(version == "v5") |>
  dplyr::select(predictor, predictor_class) |>
  dplyr::mutate(dplyr::across('predictor', stringr::str_replace, 'Year', 'year')) |>
  dplyr::mutate(dplyr::across('predictor', stringr::str_replace, 'Method','method'))

soil_covs <- c("cec_0-5cm_mean_1000", "cec_100-200cm_mean_1000", "cec_15-30cm_mean_1000",
               "cec_30-60cm_mean_1000", "cec_5-15cm_mean_1000", "cec_60-100cm_mean_1000",
               "soc_0-5cm_mean_1000", "soc_100-200cm_mean_1000", "soc_15-30cm_mean_1000",
               "soc_30-60cm_mean_1000", "soc_5-15cm_mean_1000", "soc_60-100cm_mean_1000")
actually_biotic <- c("Peatland_5x5", "Peatland_1km")

abiotic_preds <-
  predictor_metadata |>
  dplyr::filter(predictor_class %in% c("Annual Climate", "Climate Normals",
                                        "Topography", "Wetland")) |>
  dplyr::filter(!(predictor %in% actually_biotic)) |>
  dplyr::pull(predictor) |>
  c("CAfire", soil_covs)

missing_preds <- setdiff(abiotic_preds, names(stack_y))
if (length(missing_preds) > 0)
  message("NOTE: not in the covariate stack, so not in 08A's training either: ",
          paste(missing_preds, collapse = ", "))
abiotic_preds <- intersect(abiotic_preds, names(stack_y))
stack_y       <- stack_y[[abiotic_preds]]
message("Abiotic predictors (", length(abiotic_preds), "): ", paste(abiotic_preds, collapse = ", "))

# AOA of one subbasin's training pixels ---------------------------------

aoa_subbasin <- function(lo, hi, lo_xy) {
  varies <- vapply(seq_len(ncol(lo)), function(j)
    sum(!is.na(lo[, j])) > 1L && isTRUE(sd(lo[, j], na.rm = TRUE) > 0), logical(1))
  lo <- lo[, varies, drop = FALSE]; hi <- hi[, varies, drop = FALSE]
  if (ncol(lo) == 0L) return(NULL)
  med <- apply(lo, 2, median, na.rm = TRUE)
  for (j in seq_len(ncol(lo))) {
    lo[is.na(lo[, j]), j] <- med[j]
    hi[is.na(hi[, j]), j] <- med[j]
  }
  mu <- colMeans(lo); s <- apply(lo, 2, sd)
  lo <- scale(lo, mu, s); hi <- scale(hi, mu, s)

  # mean distance between training pixels, from sampled pairs
  a <- sample.int(nrow(lo), n_pairs, replace = TRUE)
  b <- sample.int(nrow(lo), n_pairs, replace = TRUE)
  d_bar <- mean(sqrt(rowSums((lo[a, , drop = FALSE] - lo[b, , drop = FALSE])^2))[a != b])

  # training DI, each pixel to the nearest training pixel in another spatial block
  fold <- stats::kmeans(lo_xy, centers = min(n_folds, nrow(lo) - 1L), nstart = 3L, iter.max = 50L)$cluster
  q    <- if (nrow(lo) > max_query) sort(sample.int(nrow(lo), max_query)) else seq_len(nrow(lo))
  di_train <- numeric(length(q))
  for (f in unique(fold[q])) {
    at  <- which(fold[q] == f)
    ref <- which(fold != f)
    di_train[at] <- FNN::get.knnx(lo[ref, , drop = FALSE], lo[q[at], , drop = FALSE],
                                  k = 1)$nn.dist[, 1] / d_bar
  }
  thr <- stats::quantile(di_train, 0.75, names = FALSE) + 1.5 * stats::IQR(di_train)

  di_hi <- FNN::get.knnx(lo, hi, k = 1)$nn.dist[, 1] / d_bar
  data.frame(n_predictors     = ncol(lo),
             aoa_threshold    = thr,
             di_train_median  = stats::median(di_train),
             di_median        = stats::median(di_hi),
             frac_outside_aoa = mean(di_hi > thr))
}

# process each subbasin ----------------------------------------------

results <- vector("list", n_sub)

for (s in seq_len(n_sub)) {

  sub_s <- all_subbasins[s]

  # crop covariate stack to subbasin
  cov_s <- tryCatch(
    terra::mask(terra::crop(stack_y, sub_s), sub_s),
    error = function(e) NULL
  )
  if (is.null(cov_s)) next

  # low-HF and high-HF masks for this subbasin, resampled as 08A does
  lowhf_s  <- terra::mask(terra::resample(lowhf_mask,  cov_s, method = "near"), sub_s)
  highhf_s <- terra::mask(terra::resample(highhf_mask, cov_s, method = "near"), sub_s)

  vals_all <- terra::values(cov_s, mat = TRUE)
  lo_idx   <- which(terra::values(lowhf_s,  mat = FALSE) == 1)
  hi_idx   <- which(terra::values(highhf_s, mat = FALSE) == 1)

  if (length(lo_idx) < 10 || length(hi_idx) < 1) next

  res <- aoa_subbasin(vals_all[lo_idx, , drop = FALSE], vals_all[hi_idx, , drop = FALSE],
                      terra::xyFromCell(cov_s, lo_idx))
  if (is.null(res)) next

  results[[s]] <- data.frame(subbasin = s,
                             HYBAS_ID = all_subbasins$first_HYBAS_ID[s],
                             n_lowhf  = length(lo_idx),
                             n_highhf = length(hi_idx),
                             res)

  if (s %% 50 == 0) message("  processed ", s, " / ", n_sub, " subbasins")

} # close for loop over subbasins

# assemble results and flag areas of extrapolation ----------------------

flags_df      <- bind_rows(results)
flags_df$flag <- flags_df$frac_outside_aoa > flag_cut

out_path <- file.path(out_dir, "extrapolation_flags.csv")
write.csv(flags_df, out_path, row.names = FALSE)

message("wrote ", nrow(flags_df), " subbasins to ", out_path)
message("  flagged (frac_outside_aoa > ", flag_cut, "): ", sum(flags_df$flag, na.rm = TRUE),
        " / ", nrow(flags_df))
message("  frac_outside_aoa quantiles (0, .1, .25, .5, .75, .9, 1): ",
        paste(round(stats::quantile(flags_df$frac_outside_aoa, c(0, .1, .25, .5, .75, .9, 1)), 3),
              collapse = " "))
