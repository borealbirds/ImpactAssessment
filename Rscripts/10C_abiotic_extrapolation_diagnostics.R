# ---
# title: Abiotic extrapolation diagnostics for BART backfilling
# author: Mannfred Boehm
# created: May 19, 2026 (rewritten 2026-09-25)
# ---

# per subbasin: do the backfilled (high HF) pixels lie inside the region of abiotic predictor
# space its training (low HF) pixels cover? the area of applicability (AOA) of Meyer & Pebesma
# (2021, Methods Ecol Evol 12:1620), on 07's abiotic predictors minus the footprint class
# (CanHF, canroad, CCNL): training pixels are HF < 1 and backfilled ones HF >= 1 by
# construction, so on those every backfilled pixel would be outside.
#   DI(x)     = distance from x to its nearest training pixel / mean distance between training
#               pixels, with predictors scaled by the training SDs
#   threshold = upper whisker (Q75 + 1.5 IQR) of the training pixels' DI, each measured to the
#               nearest training pixel in another spatial fold (k-means blocks on x, y). random
#               folds would put every training DI near 0, as neighbouring pixels are near-duplicates
#   flag      = more than half of the subbasin's backfilled pixels have DI > threshold
# predictors constant on a subbasin's training pixels are dropped (BART can't split on them);
# NAs take the training median, as in 08A.
# runs locally (FNN isn't on Fir) -> data/derived_data/rds_files/extrapolation_flags.csv

library(FNN)
library(terra)

ia_dir  <- getwd()
raw_dir <- file.path(ia_dir, "data", "raw_data")
out_dir <- file.path(ia_dir, "data", "derived_data", "rds_files")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

n_folds   <- 10L     # spatial blocks for the training DI
max_query <- 5000L   # training pixels whose DI sets the threshold (sampled above this)
n_pairs   <- 1e5     # sampled pairs for the mean training distance
flag_cut  <- 0.5
set.seed(20260925)

all_subbasins <- terra::vect(file.path(raw_dir, "hydrobasins_masked_merged_subset.gpkg"))
lowhf_mask    <- terra::rast(file.path(raw_dir, "hirshpearson", "CanHF_1km_lessthan1.tif"))
highhf_mask   <- terra::rast(file.path(raw_dir, "hirshpearson", "CanHF_1km_morethan1.tif"))
stack_y       <- terra::rast(file.path(raw_dir, "covariates_mosaiced", "covariates_mosaiced_2020.tif"))

# 07's abiotic predictors, minus the footprint ("Disturbance") class
pm <- dplyr::filter(BAMexploreR::predictor_metadata, version == "v5")
soil_covs <- c("cec_0-5cm_mean_1000", "cec_100-200cm_mean_1000", "cec_15-30cm_mean_1000",
               "cec_30-60cm_mean_1000", "cec_5-15cm_mean_1000", "cec_60-100cm_mean_1000",
               "soc_0-5cm_mean_1000", "soc_100-200cm_mean_1000", "soc_15-30cm_mean_1000",
               "soc_30-60cm_mean_1000", "soc_5-15cm_mean_1000", "soc_60-100cm_mean_1000")
abiotic_preds <- c(setdiff(pm$predictor[pm$predictor_class %in% c("Annual Climate", "Climate Normals", "Topography", "Wetland")],
                           c("Peatland_5x5", "Peatland_1km")),
                   "CAfire", soil_covs)
stack_y <- stack_y[[intersect(abiotic_preds, names(stack_y))]]



# ---------------------------------------------------
# AOA of one subbasin's training pixels

aoa_subbasin <- function(lo, hi, lo_xy) {

  varies <- vapply(seq_len(ncol(lo)), function(j)
    sum(!is.na(lo[, j])) > 1L && isTRUE(sd(lo[, j], na.rm = TRUE) > 0), logical(1))
  lo <- lo[, varies, drop = FALSE]
  hi <- hi[, varies, drop = FALSE]
  if (ncol(lo) == 0L) return(NULL)

  med <- apply(lo, 2, median, na.rm = TRUE)
  for (j in seq_len(ncol(lo))) {
    lo[is.na(lo[, j]), j] <- med[j]
    hi[is.na(hi[, j]), j] <- med[j]
  }
  mu <- colMeans(lo)
  s  <- apply(lo, 2, sd)
  lo <- scale(lo, mu, s)
  hi <- scale(hi, mu, s)

  # mean distance between training pixels, from sampled pairs
  a <- sample.int(nrow(lo), n_pairs, replace = TRUE)
  b <- sample.int(nrow(lo), n_pairs, replace = TRUE)
  d_bar <- mean(sqrt(rowSums((lo[a, , drop = FALSE] - lo[b, , drop = FALSE])^2))[a != b])

  # training DI: each pixel to the nearest training pixel in another spatial block
  fold <- stats::kmeans(lo_xy, centers = min(n_folds, nrow(lo) - 1L), nstart = 3L, iter.max = 50L)$cluster
  q    <- if (nrow(lo) > max_query) sort(sample.int(nrow(lo), max_query)) else seq_len(nrow(lo))
  di_train <- numeric(length(q))
  for (f in unique(fold[q])) {
    at  <- which(fold[q] == f)
    ref <- which(fold != f)
    di_train[at] <- FNN::get.knnx(lo[ref, , drop = FALSE], lo[q[at], , drop = FALSE], k = 1)$nn.dist[, 1] / d_bar
  }
  thr <- stats::quantile(di_train, 0.75, names = FALSE) + 1.5 * stats::IQR(di_train)

  di_hi <- FNN::get.knnx(lo, hi, k = 1)$nn.dist[, 1] / d_bar
  data.frame(n_predictors     = ncol(lo),
             aoa_threshold    = thr,
             di_train_median  = stats::median(di_train),
             di_median        = stats::median(di_hi),
             frac_outside_aoa = mean(di_hi > thr))
}



# ---------------------------------------------------
# every subbasin (~55 min)

results <- vector("list", nrow(all_subbasins))

for (s in seq_len(nrow(all_subbasins))) {

  sub_s <- all_subbasins[s]
  cov_s <- terra::mask(terra::crop(stack_y, sub_s), sub_s)

  # low and high HF pixels, resampled as in 08A
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

  if (s %% 50 == 0) message("processed ", s, " / ", nrow(all_subbasins), " subbasins")
}

flags_df      <- dplyr::bind_rows(results)
flags_df$flag <- flags_df$frac_outside_aoa > flag_cut
write.csv(flags_df, file.path(out_dir, "extrapolation_flags.csv"), row.names = FALSE)

message("flagged: ", sum(flags_df$flag, na.rm = TRUE), " of ", nrow(flags_df), " subbasins")
