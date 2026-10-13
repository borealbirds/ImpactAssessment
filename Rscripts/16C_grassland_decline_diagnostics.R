# ---
# title: Is the BOBO / GRSP decline a backfill artefact or a real signal?
# author: Mannfred Boehm
# created: October 8, 2026
# ---

# Without industry, BOBO falls from 6.76 M to 3.35 M birds and GRSP from 0.96 M to 0.49 M (14C),
# almost all through the backfilled vegetation (14B: vegetation part BOBO -4.45 M, GRSP -440k;
# the footprint covariates at 0 move little). This asks whether that is BART extrapolating
# from the thin southern training sets, or what the training land itself implies for the bird.
# It runs on the subbasins carrying 80% of the loss in can11 and can13 (95% of GRSP's, 83% of
# BOBO's), with 12F's bird models, caps, weights and design rules.
#
# T1  Does the counterfactual look like its own training land, to the bird? At each direct-
#     footprint pixel (where 12F's counterfactual is the backfill), replace BART's vegetation by
#     the OBSERVED vegetation of one of the pixel's k nearest low-HF training pixels in the
#     subbasin (nearest in 10C's standardised abiotic space). This is a hot-deck imputation:
#     each pixel takes one donor's whole biotic vector, so class frequencies and the
#     covariates' joint structure survive. The pixel keeps its own abiotic covariates and the
#     footprint covariates at 0, and the same 32 bootstraps predict it. Every training pixel is
#     a donor, missing values included: SCANFIheightcv_* is NA wherever canopy height is 0
#     (91-98% of the training land in GRSP's southern can11 subbasins), V5 predicts those
#     pixels down gbm's missing branch, and a complete-case donor pool would draw only from
#     the few treed ones (the 2026-10-08 run did, and read can11 as woodier than it is).
#       knn_all  donors from every training pixel of the subbasin;
#       knn_same donors only from training pixels in the subbasin's target BCR (can13's
#                subbasins train mostly on can12 land).
#     If 12F's counterfactual (`bart`, from pixel_change/) sits well below knn_all, BART made the
#     land worse for the bird than the land it learned from.
# T2  Which backfilled covariates carry the gap? One at a time, swap BART's value (posterior mean,
#     or the modal class 08B writes) into the kNN design (delta_gap = BART - kNN) and into the
#     observed design with the footprint at 0 (delta_veg, the vegetation part). Distributions of
#     every backfilled covariate in training, observed footprint, kNN and BART.
# T3  Is the bird model extrapolating? Share of values outside each covariate's survey range in
#     each bootstrap (gbm var.levels: the survey deciles 0-100%), or classes a bootstrap never
#     saw; partial dependence of the top drivers, averaged over each group's values.
# T5  Does the backfill keep its training land's zeros? Each continuous backfilled covariate is
#     quantile-mapped onto the subbasin's finite training values: BART's ranking of the pixels
#     is kept, the training land's distribution imposed (it cannot restore missing values). BART's Gaussian model of log1p(x) never returns
#     the exact zeros a near-treeless landscape has, and a bird model that splits just above 0
#     (GRSP on SCANFIprcD_5x5: flat to 0.024%, x0.5 by 0.07%) then reads every backfilled pixel
#     as treed. qm_ratio = density under the mapped design / under BART's (8 bootstraps).
# T4  Training size and abiotic distance: within each subbasin, cf/obs inside vs outside 10C's
#     area of applicability (per pixel); across every can11/can13 subbasin, footprint cf/obs
#     against the backfill-to-training ratio and frac_outside_aoa.
#
# BF_NEW_DIR=<dir> scores a candidate backfill too (a 07 test run's subbasin_{i}_backfill.tif,
# pulled into cluster_logs/grassland_diag/<dir>/): bart_new_check predicts it as 12F will, with
# the same pixels, draws and bootstraps as bart_check (production's backfill): one draw index per
# scenario shared by every covariate, drawn land-cover classes (<cov>_class_###), and NA where a
# covariate is modelled blank (<cov>_blank == 1), which gbm sends down its missing branch.
# Outputs then go to grassland_diagnostics/<dir>/, and the distributions gain group fp_bart_new.
#
# Inputs pulled from Fir into cluster_logs/grassland_diag/ (gitignored): models/
# {sp}_{bcr}.Rdata (06_bootstraps), backfill/subbasin_{i}_backfill.tif, predictions/{sp}/
# truncation_params.rds and {bcr}/2020/weight.tif. The BCR covariate stacks are
# cluster_logs/sanity/stacks/ (byte-identical in size to Fir's gis/stacks).
# Runs locally, one case at a time (~5 GB peak with BOBO can13's models loaded).
# output: data/derived_data/rds_files/grassland_diagnostics/*.csv


suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(FNN)
  library(gbm)
})

# define paths -------------------------------------------------------

ia_dir    <- getwd()
diag_dir  <- file.path(ia_dir, "cluster_logs", "grassland_diag")
stack_dir <- file.path(ia_dir, "cluster_logs", "sanity", "stacks")
hirsh_dir <- file.path(ia_dir, "data", "raw_data", "hirshpearson")
dt_dir    <- file.path(ia_dir, "data", "derived_data", "density_tables")
out_dir   <- file.path(ia_dir, "data", "derived_data", "rds_files", "grassland_diagnostics")
bf_new_dir <- Sys.getenv("BF_NEW_DIR", "")
if (nzchar(bf_new_dir)) {
  out_dir    <- file.path(out_dir, basename(bf_new_dir))
  bf_new_dir <- file.path(diag_dir, bf_new_dir)
  if (!dir.exists(bf_new_dir)) stop("BF_NEW_DIR: no ", bf_new_dir)
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

year <- 2020

# the subbasins carrying 80% of each species' loss in can11 and can13 (14B shapley_subbasin.csv)
cases <- list(
  list(species = "GRSP", bcr = "can11", subs = c(57, 56, 63, 60, 61, 62, 55)),
  list(species = "BOBO", bcr = "can11", subs = c(57, 61, 62)),
  list(species = "BOBO", bcr = "can13", subs = c(98, 107, 92, 77, 72, 82, 103)))

# smoke test: CASE=BOBO_can13 N_FP_MAX=300 Rscript Rscripts/16C_grassland_decline_diagnostics.R
if (nzchar(Sys.getenv("CASE")))
  cases <- Filter(function(cs) paste(cs$species, cs$bcr, sep = "_") %in%
                    strsplit(Sys.getenv("CASE"), ",")[[1]], cases)

k_nn       <- 10L     # donors per footprint pixel
n_rep      <- 5L      # hot-deck replicates per pixel
n_fp_max   <- as.integer(Sys.getenv("N_FP_MAX", "5000"))   # direct-footprint pixels sampled per subbasin (T1)
n_t2       <- 2000L   # of those, pixels for T2's swaps and the BART reconstruction check
n_boot_sub <- 8L      # bootstraps for T2's swaps, the BART reconstruction check and T3's PD
n_draw_chk <- 10L     # BART draws in the reconstruction check
min_donors <- 50L     # smallest same-BCR donor pool for knn_same
cut        <- 0.2     # verdict threshold: a 20% difference counts
set.seed(20261008)

# covariate classes, as 12D defines them -------------------------------------

categorical_responses <- c("ABoVE_1km", "NLCD_1km", "MODISLCC_1km", "MODISLCC_5x5", "SCANFI_1km", "VLCE_1km")
predictor_metadata <-
  dplyr::tibble(BAMexploreR::predictor_metadata) |>
  dplyr::filter(version == "v5") |>
  dplyr::select(predictor, definition, predictor_class) |>
  dplyr::mutate(dplyr::across('predictor', stringr::str_replace, 'Year', 'year')) |>
  dplyr::mutate(dplyr::across('predictor', stringr::str_replace, 'Method', 'method'))

actually_biotic_what <- c("Peatland_5x5", "Peatland_1km")
actually_biotic_df   <- tibble::tibble(predictor = actually_biotic_what, predictor_class = c("Wetland", "Wetland"))
abiotic_vars <- predictor_metadata |>
  dplyr::filter(predictor_class %in% c("Annual Climate", "Climate Normals", "Topography", "Wetland", "Disturbance", "Time", "Method")) |>
  dplyr::filter(!(predictor %in% actually_biotic_what))
biotic_continuous_vars <- predictor_metadata |>
  dplyr::filter(!(predictor_class %in% abiotic_vars$predictor_class)) |>
  dplyr::bind_rows(actually_biotic_df) |>
  dplyr::filter(!predictor %in% categorical_responses) |>
  dplyr::pull(predictor)
disturbance_vars <- predictor_metadata |>
  dplyr::filter(predictor_class == "Disturbance") |>
  dplyr::pull(predictor)

# 10C's abiotic predictors (07's abiotic_vars minus "Disturbance"), the matching space for T1
soil_covs <- c("cec_0-5cm_mean_1000", "cec_100-200cm_mean_1000", "cec_15-30cm_mean_1000",
               "cec_30-60cm_mean_1000", "cec_5-15cm_mean_1000", "cec_60-100cm_mean_1000",
               "soc_0-5cm_mean_1000", "soc_100-200cm_mean_1000", "soc_15-30cm_mean_1000",
               "soc_30-60cm_mean_1000", "soc_5-15cm_mean_1000", "soc_60-100cm_mean_1000")
abiotic_preds <-
  predictor_metadata |>
  dplyr::filter(predictor_class %in% c("Annual Climate", "Climate Normals", "Topography", "Wetland")) |>
  dplyr::filter(!(predictor %in% actually_biotic_what)) |>
  dplyr::pull(predictor) |>
  c("CAfire", soil_covs)

# 12F's per-bootstrap factor levels and 12G's fast tree walk ---------------------

source(file.path(ia_dir, "Rscripts", "12E_shapley_utils.R"))                    # canonical_sectors()
source(file.path(ia_dir, "Rscripts", "12F_predict_species_all_coalitions.R"))   # as_model_factors()
source(file.path(ia_dir, "Rscripts", "12G_gbm_tree_walk.R"))
Rcpp::sourceCpp(file.path(ia_dir, "Rscripts", "12G_gbm_tree_walk.cpp"))

# load spatial data ---------------------------------------------------

subs_v <- vect(file.path(ia_dir, "data", "raw_data", "hydrobasins_masked_merged_subset.gpkg"))
bcr_v  <- project(vect(file.path(ia_dir, "data", "raw_data", "Regions", "BAM_BCR_NationalModel_Unbuffered.shp")),
                  crs(subs_v))
bcr_v$bcr <- paste0(bcr_v$country, bcr_v$subUnit)
bcr_v$idx <- seq_len(nrow(bcr_v))

lowhf  <- rast(file.path(hirsh_dir, "CanHF_1km_lessthan1.tif"))
mosaic <- rast(file.path(ia_dir, "data", "raw_data", "covariates_mosaiced",
                         paste0("covariates_mosaiced_", year, ".tif")))
flags  <- read.csv(file.path(ia_dir, "data", "derived_data", "rds_files", "extrapolation_flags.csv"))
abiotic_preds <- intersect(abiotic_preds, names(mosaic))
sectors <- canonical_sectors()

# helpers -------------------------------------------------------------

nan_to_na <- function(X) {
  for (v in names(X)) if (is.double(X[[v]])) X[[v]][is.nan(X[[v]])] <- NA_real_
  X
}

# 10C's area of applicability, returning the standardised matrices and every DI. Training DI
# is measured to the nearest training pixel in another spatial block; threshold = its upper
# whisker. Predictors the training pixels hold constant are dropped; NAs take the training median.
aoa_standardise <- function(lo, hi, lo_xy, n_folds = 10L, max_query = 5000L, n_pairs = 1e5) {
  varies <- vapply(seq_len(ncol(lo)), function(j)
    sum(!is.na(lo[, j])) > 1L && isTRUE(sd(lo[, j], na.rm = TRUE) > 0), logical(1))
  lo <- lo[, varies, drop = FALSE]; hi <- hi[, varies, drop = FALSE]
  med <- apply(lo, 2, median, na.rm = TRUE)
  for (j in seq_len(ncol(lo))) {
    lo[is.na(lo[, j]), j] <- med[j]
    hi[is.na(hi[, j]), j] <- med[j]
  }
  mu <- colMeans(lo); s <- apply(lo, 2, sd)
  lo <- scale(lo, mu, s); hi <- scale(hi, mu, s)

  a <- sample.int(nrow(lo), n_pairs, replace = TRUE)
  b <- sample.int(nrow(lo), n_pairs, replace = TRUE)
  d_bar <- mean(sqrt(rowSums((lo[a, , drop = FALSE] - lo[b, , drop = FALSE])^2))[a != b])

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
  list(lo = lo, hi = hi, thr = thr,
       di_hi = FNN::get.knnx(lo, hi, k = 1)$nn.dist[, 1] / d_bar)
}

# summaries of one covariate in one group of pixels
dist_cont <- function(x, var, group) {
  frac_na <- mean(!is.finite(x))
  x <- x[is.finite(x)]
  data.frame(var = var, group = group, n = length(x), frac_na = frac_na,
             mean = if (length(x)) mean(x) else NA,
             q10 = if (length(x)) quantile(x, 0.1, names = FALSE) else NA,
             q50 = if (length(x)) quantile(x, 0.5, names = FALSE) else NA,
             q90 = if (length(x)) quantile(x, 0.9, names = FALSE) else NA)
}
dist_cat <- function(x, var, group) {
  x <- x[!is.na(x)]
  if (!length(x)) return(NULL)
  tb <- table(as.character(x))
  data.frame(var = var, group = group, class = names(tb), share = as.numeric(tb) / length(x))
}

t1_rows <- list(); t2_rows <- list(); dc_rows <- list(); dk_rows <- list()
t3_rows <- list(); pd_rows <- list()

# process each species x BCR ------------------------------------------------

for (cs in cases) {

  sp <- cs$species; bcr <- cs$bcr
  message(Sys.time(), " | ", sp, " ", bcr, " | loading models")

  e <- new.env(parent = emptyenv())
  load(file.path(diag_dir, "models", paste0(sp, "_", bcr, ".Rdata")), envir = e)
  b.list <- e$b.list; rm(e); gc()
  if (!identical(attr(b.list[[1]], "bcr"), bcr)) stop(sp, " ", bcr, ": model file carries another BCR")
  n_boot     <- length(b.list)
  model_vars <- b.list[[1]]$var.names
  if (!all(vapply(b.list, function(m) identical(m$var.names, model_vars), logical(1L))))
    stop(sp, " ", bcr, ": bootstraps do not share one var.names")
  cat_vars  <- intersect(model_vars, categorical_responses)
  dist_vars <- intersect(disturbance_vars, model_vars)
  boots_sub <- seq_len(min(n_boot_sub, n_boot))

  tp  <- readRDS(file.path(diag_dir, "predictions", sp, "truncation_params.rds"))[[paste(bcr, year, sep = "_")]]
  qsp <- tp$densmax
  q99 <- tp$q99
  if (!is.finite(qsp) || !is.finite(q99)) stop(sp, " ", bcr, ": no densmax/q99")

  # capped birds/ha averaged over bootstraps, exactly as 12F's workers predict
  pm <- function(X, boots = seq_len(n_boot)) {
    acc <- 0
    for (i in boots)
      acc <- acc + pmin(pmin(gbm_predict_design(b.list[[i]],
                    gbm_design(b.list[[i]], as_model_factors(X, b.list[[i]], cat_vars))), qsp), q99)
    if (anyNA(acc)) {
      bad <- which(is.na(acc))
      message("  NOTE: ", length(bad), " of ", length(acc), " predictions are NA; first row: ",
              paste(names(X), unlist(X[bad[1], ]), sep = "=", collapse = " "))
    }
    acc / length(boots)
  }

  stack_obs <- rast(file.path(stack_dir, paste0(bcr, "_", year, ".tif")))
  w_r  <- rast(file.path(diag_dir, "predictions", sp, bcr, year, "weight.tif"))
  pc_r <- rast(file.path(dt_dir, "pixel_change", paste0(sp, "_", year, "_", bcr, ".tif")))
  if (names(w_r)[1] != "weight_v3_touches") stop(sp, " ", bcr, ": stale weight.tif")
  if (!compareGeom(w_r, stack_obs[[1]], stopOnError = FALSE) ||
      !compareGeom(pc_r, stack_obs[[1]], stopOnError = FALSE))
    stop(sp, " ", bcr, ": weight.tif or pixel_change is not on the stack grid")
  w_v      <- values(w_r, mat = FALSE)
  pc_obs_v <- values(pc_r[["obs"]], mat = FALSE)
  pc_cf_v  <- values(pc_r[["cf"]], mat = FALSE)

  # the all-sector footprint and its direct part, on the stack grid, as 12F builds them
  message(Sys.time(), " | ", sp, " ", bcr, " | footprint masks")
  hf_v    <- values(project(rast(file.path(hirsh_dir, "CanHF_1km_morethan1.tif")), stack_obs[[1]],
                            method = "near"), mat = FALSE)
  any_sec <- logical(ncell(stack_obs)); any_dir <- logical(ncell(stack_obs))
  for (sec in sectors) {
    sv  <- values(project(rast(file.path(hirsh_dir, paste0(sec, ".tif"))), stack_obs[[1]],
                          method = "near"), mat = FALSE)
    dv  <- values(project(rast(file.path(hirsh_dir, paste0(sec, "_direct.tif"))), stack_obs[[1]],
                          method = "near"), mat = FALSE)
    mem <- !is.na(sv) & sv > 0
    any_sec <- any_sec | mem
    any_dir <- any_dir | (mem & !is.na(dv) & dv == 1)
  }
  super_v  <- any_sec & !is.na(hf_v) & hf_v >= 1
  direct_v <- super_v & any_dir
  rm(hf_v, any_sec, any_dir, sv, dv, mem)

  # subbasin zones on the stack grid, as 12D/12F rasterize them
  hits    <- relate(subs_v, bcr_v[bcr_v$bcr == bcr], relation = "intersects")
  sub_ids <- which(apply(hits, 1, any))
  zone_v  <- values(rasterize(subs_v[sub_ids], stack_obs[[1]], field = "first_HYBAS_ID"), mat = FALSE)

  tab <- readRDS(file.path(dt_dir, paste0(sp, "_", year, "_coalition_256.rds")))
  tab <- tab[tab$bcr == bcr, ]

  obs_vals <- values(stack_obs[[model_vars]])
  pool <- list(obs = list(), knn = list(), bart = list(), train = list())

  for (s in cs$subs) {

    hy   <- subs_v$first_HYBAS_ID[s]
    in_s <- !is.na(zone_v) & zone_v == hy
    message(Sys.time(), " | ", sp, " ", bcr, " | subbasin ", s)

    # pixel_change must add up to the cid-256 table: pins the zone raster to 12F's
    trow <- tab[tab$subbasin == s, ]
    if (nrow(trow) != 1L) stop(sp, " ", bcr, " subbasin ", s, ": not one cid-256 row")
    for (k in c("obs", "cf")) {
      pix <- sum(if (k == "obs") pc_obs_v[in_s] else pc_cf_v[in_s], na.rm = TRUE)
      ref <- if (k == "obs") trow$obs_on_coalition_mean else trow$bf_on_coalition_mean
      if (abs(pix - ref) > 1e-4 * abs(ref) + 1e-6)
        stop(sp, " ", bcr, " subbasin ", s, ": pixel_change ", k, " sums to ", pix,
             " but the cid-256 table says ", ref, " - the zone raster differs from 12F's")
    }

    # direct-footprint pixels, where 12F's counterfactual is the backfill
    fp_all <- which(in_s & direct_v & is.finite(w_v) & w_v > 0 & is.finite(pc_obs_v) & is.finite(pc_cf_v))
    fp     <- if (length(fp_all) > n_fp_max) sort(sample(fp_all, n_fp_max)) else fp_all
    fp_xy  <- xyFromCell(stack_obs, fp)

    # training pixels exactly as 08A takes them (mosaic grid)
    sub_poly <- subs_v[s]
    cov_s    <- mask(crop(mosaic[[unique(c(model_vars, abiotic_preds))]], sub_poly), sub_poly)
    lo_s     <- mask(resample(lowhf, cov_s, method = "near"), sub_poly)
    lo_cells <- which(values(lo_s, mat = FALSE) == 1)
    if (length(lo_cells) != flags$n_lowhf[flags$subbasin == s])
      stop(sp, " ", bcr, " subbasin ", s, ": ", length(lo_cells), " training pixels, 10C counted ",
           flags$n_lowhf[flags$subbasin == s])
    mv     <- values(cov_s, mat = TRUE); mv[is.nan(mv)] <- NA_real_
    lo_xy  <- xyFromCell(cov_s, lo_cells)
    fp_m   <- cellFromXY(cov_s, fp_xy)
    lo_bcr <- bcr_v$bcr[values(rasterize(bcr_v, cov_s[[1]], field = "idx"), mat = FALSE)[lo_cells]]

    # BART's backfill at the sampled pixels, as 12F reads it ------------------------
    bf_r   <- rast(file.path(diag_dir, "backfill", paste0("subbasin_", s, "_backfill.tif")))
    bnm    <- names(bf_r)
    fp_b   <- cellFromXY(bf_r, fp_xy)
    in_grid <- !is.na(fp_m) & !is.na(fp_b)
    if (!all(in_grid)) stop(sp, " ", bcr, " subbasin ", s, ": ", sum(!in_grid),
                            " footprint pixels fall outside the subbasin's mosaic or backfill grid")
    draw_v <-unique(sub("_draw_[0-9]{3}$", "", grep("_draw_[0-9]{3}$", bnm, value = TRUE)))
    mean_v <- sub("_mean$", "", grep("_mean$", bnm, value = TRUE))
    cont_bf <- intersect(setdiff(model_vars, cat_vars),
                         union(draw_v, intersect(mean_v, biotic_continuous_vars)))
    cat_bf  <- cat_vars[cat_vars %in% bnm | paste0(cat_vars, "_mean") %in% bnm]
    bf_vars <- c(cont_bf, cat_bf)

    draws <- setNames(lapply(cont_bf, function(v) {
      ml <- paste0(v, "_mean")
      if (v %in% draw_v) {
        lyr <- intersect(paste0(v, "_draw_", sprintf("%03d", 1:100)), bnm)
        m   <- expm1(as.matrix(terra::extract(bf_r[[lyr]], fp_b)))
        m[!is.finite(m)] <- NA_real_
        if (ml %in% bnm) {
          const <- terra::extract(bf_r[[ml]], fp_b)[, 1]
          fill  <- rowSums(!is.na(m)) == 0L & is.finite(const)
          m[fill, ] <- const[fill]
        }
        pmax(m, 0)
      } else {
        const <- terra::extract(bf_r[[ml]], fp_b)[, 1]
        const[!is.finite(const)] <- NA_real_
        matrix(pmax(const, 0), ncol = 1L)
      }
    }), cont_bf)
    cats <- setNames(lapply(cat_bf, function(v)
      terra::extract(bf_r[[if (v %in% bnm) v else paste0(v, "_mean")]], fp_b)[, 1]), cat_bf)

    # 12F keeps only pixels whose backfilled covariates are complete
    ok <- stats::complete.cases(do.call(cbind, c(lapply(draws, function(m) m[, 1]), cats)))
    if (mean(ok) < 0.99)
      message("  NOTE: ", sum(!ok), " of ", length(ok), " sampled pixels have incomplete backfill; dropped")
    fp <- fp[ok]; fp_xy <- fp_xy[ok, , drop = FALSE]; fp_m <- fp_m[ok]
    draws <- lapply(draws, function(m) m[ok, , drop = FALSE]); cats <- lapply(cats, `[`, ok)
    f  <- length(fp) / length(fp_all)        # sampling fraction: birds below are scaled by 1/f
    wv <- w_v[fp] * 100                      # birds/ha -> birds/km2, weighted as 12F

    # observed, d0 and BART designs ---------------------------------------------------
    X_obs <- nan_to_na(as.data.frame(obs_vals[fp, , drop = FALSE], check.names = FALSE))
    X_d0  <- X_obs; for (v in dist_vars) X_d0[[v]] <- 0
    X_bart <- X_d0
    for (v in cont_bf) X_bart[[v]] <- rowMeans(draws[[v]], na.rm = TRUE)
    for (v in cat_bf)  X_bart[[v]] <- cats[[v]]

    p_obs <- pm(X_obs)
    pc_obs <- pc_obs_v[fp]; pc_cf <- pc_cf_v[fp]
    # identity: this path must reproduce 12F's observed side (V5's clamped bootstraps)
    bad <- abs(p_obs * wv - pc_obs) > 1e-5 * abs(pc_obs) + 1e-9
    if (any(bad))
      stop(sp, " ", bcr, " subbasin ", s, ": identity failed at ", sum(bad), " of ", length(bad),
           " pixels (max rel diff ", signif(max(abs(p_obs * wv - pc_obs)[bad] / abs(pc_obs[bad])), 3),
           ") - models, covariates or weights are not 12F's")
    p_d0 <- pm(X_d0)

    # T2 and the reconstruction check use a subset (BOBO's models predict ~9k rows/s)
    si <- if (length(fp) > n_t2) sort(sample.int(length(fp), n_t2)) else seq_along(fp)
    f2 <- length(si) / length(fp_all)

    # does the BART design reproduce 12F's counterfactual? (subsampled draws and bootstraps)
    acc <- 0
    chk_draws <- sample.int(100L, n_draw_chk)
    for (d in chk_draws) {
      X <- X_d0[si, , drop = FALSE]
      for (v in cont_bf) X[[v]] <- draws[[v]][si, if (ncol(draws[[v]]) == 1L) 1L else d]
      for (v in cat_bf)  X[[v]] <- cats[[v]][si]
      acc <- acc + pm(X, boots_sub)
    }
    p_chk <- acc / n_draw_chk
    p_obs_sub <- pm(X_obs[si, , drop = FALSE], boots_sub)   # how far the 8 bootstraps sit from all 32

    # the candidate backfill, predicted as 12F will, on bart_check's pixels, draws and bootstraps
    p_new <- NULL; new_vals <- NULL; new_na <- NA_real_
    if (nzchar(bf_new_dir)) {
      nb_r <- rast(file.path(bf_new_dir, paste0("subbasin_", s, "_backfill.tif")))
      if (!compareGeom(nb_r, bf_r, stopOnError = FALSE))
        stop(sp, " ", bcr, " subbasin ", s, ": the candidate backfill is not on production's grid")
      nnm <- names(nb_r)
      fpn <- cellFromXY(nb_r, fp_xy)
      blank_at <- function(v) {
        bl <- paste0(v, "_blank")
        if (bl %in% nnm) terra::extract(nb_r[[bl]], fpn)[, 1] %in% 1 else rep(FALSE, length(fpn))
      }
      new_draws <- setNames(lapply(cont_bf, function(v) {
        lyr <- intersect(paste0(v, "_draw_", sprintf("%03d", 1:100)), nnm)
        ml  <- paste0(v, "_mean")
        m   <- if (length(lyr)) expm1(as.matrix(terra::extract(nb_r[[lyr]], fpn)))
               else matrix(NA_real_, length(fpn), 1L)
        m[!is.finite(m)] <- NA_real_
        if (ml %in% nnm) {
          const <- terra::extract(nb_r[[ml]], fpn)[, 1]
          fill  <- rowSums(!is.na(m)) == 0L & is.finite(const)
          m[fill, ] <- const[fill]
        }
        m <- pmax(m, 0)
        m[blank_at(v), ] <- NA_real_
        m
      }), cont_bf)
      new_cls <- setNames(lapply(cat_bf, function(v) {
        lyr <- paste0(v, "_class_", sprintf("%03d", 1:100))
        if (!all(lyr %in% nnm)) stop(sp, " ", bcr, " subbasin ", s, ": no drawn classes for ", v)
        m <- as.matrix(terra::extract(nb_r[[lyr]], fpn))
        m[blank_at(v), ] <- NA_real_
        m
      }), cat_bf)
      # an NA that is not a modelled blank would be a failed backfill
      blank_any <- Reduce(`|`, lapply(bf_vars, blank_at))
      new_na <- mean(!blank_any & !stats::complete.cases(do.call(cbind, c(
        lapply(new_draws, function(m) m[, 1]), lapply(new_cls, function(m) m[, 1])))))
      acc <- 0
      for (d in chk_draws) {
        X <- X_d0[si, , drop = FALSE]
        for (v in cont_bf) X[[v]] <- new_draws[[v]][si, if (ncol(new_draws[[v]]) == 1L) 1L else d]
        for (v in cat_bf)  X[[v]] <- new_cls[[v]][si, d]
        acc <- acc + pm(X, boots_sub)
      }
      p_new <- acc / n_draw_chk
      new_vals <- c(lapply(new_draws, function(m) suppressWarnings(rowMeans(m, na.rm = TRUE))),
                    lapply(new_cls, function(m) as.vector(m)))
      rm(nb_r, new_draws, new_cls)
    }

    # T1: hot-deck kNN imputation from the subbasin's own training land ----------------
    lo_bio   <- mv[lo_cells, bf_vars, drop = FALSE]
    donor_ok <- rep(TRUE, nrow(lo_bio))      # missing values are copied, as V5 sees them
    aoa      <- aoa_standardise(mv[lo_cells, abiotic_preds, drop = FALSE],
                                mv[fp_m, abiotic_preds, drop = FALSE], lo_xy)
    hot_deck <- function(don) {
      nn  <- FNN::get.knnx(aoa$lo[don, , drop = FALSE], aoa$hi, k = min(k_nn, length(don)))$nn.index
      acc <- 0; X1 <- NULL
      for (r in seq_len(n_rep)) {
        pick <- don[nn[cbind(seq_len(nrow(nn)), sample.int(ncol(nn), nrow(nn), replace = TRUE))]]
        X <- X_d0
        for (v in bf_vars) X[[v]] <- lo_bio[pick, v]
        if (r == 1L) X1 <- X
        acc <- acc + pm(X)
      }
      list(p = acc / n_rep, X1 = X1)
    }
    don_all  <- which(donor_ok)
    don_same <- which(donor_ok & !is.na(lo_bcr) & lo_bcr == bcr)
    hd_all   <- hot_deck(don_all)
    hd_same  <- if (length(don_same) >= min_donors) hot_deck(don_same) else NULL

    # the training land's own density under the target BCR's model (footprint at 0, unweighted)
    lo_i <- if (length(lo_cells) > n_fp_max) sort(sample.int(length(lo_cells), n_fp_max)) else seq_along(lo_cells)
    X_lo <- nan_to_na(as.data.frame(mv[lo_cells[lo_i], model_vars, drop = FALSE], check.names = FALSE))
    for (v in dist_vars) X_lo[[v]] <- 0
    p_lo <- pm(X_lo)
    same_lo <- !is.na(lo_bcr[lo_i]) & lo_bcr[lo_i] == bcr

    in_aoa <- aoa$di_hi <= aoa$thr
    t1_rows[[length(t1_rows) + 1L]] <- data.frame(
      species = sp, bcr = bcr, subbasin = s, HYBAS_ID = hy,
      n_fp_direct = length(fp_all), n_sample = length(fp),
      n_train = length(lo_cells), n_donors = length(don_all), n_donors_same_bcr = length(don_same),
      frac_train_complete = mean(stats::complete.cases(lo_bio)),
      frac_train_other_bcr = mean(is.na(lo_bcr) | lo_bcr != bcr),
      obs = sum(pc_obs) / f, d0 = sum(p_d0 * wv) / f, bart = sum(pc_cf) / f,
      bart_check = sum(p_chk * wv[si]) / f2, bart_check_ref = sum(pc_cf[si]) / f2,
      obs_check = sum(p_obs_sub * wv[si]) / f2, obs_check_ref = sum(pc_obs[si]) / f2,
      bart_new_check = if (is.null(p_new)) NA_real_ else sum(p_new * wv[si]) / f2,
      new_frac_na_not_blank = new_na,
      knn_all  = sum(hd_all$p * wv) / f,
      knn_same = if (is.null(hd_same)) NA_real_ else sum(hd_same$p * wv) / f,
      dens_fp_obs = mean(p_obs) * 100,
      dens_train_all = mean(p_lo) * 100,
      dens_train_same = if (any(same_lo)) mean(p_lo[same_lo]) * 100 else NA_real_,
      aoa_threshold = aoa$thr, frac_outside_aoa_sample = mean(!in_aoa),
      obs_in_aoa  = sum(pc_obs[in_aoa]) / f,  obs_out_aoa  = sum(pc_obs[!in_aoa]) / f,
      bart_in_aoa = sum(pc_cf[in_aoa]) / f,   bart_out_aoa = sum(pc_cf[!in_aoa]) / f,
      knn_in_aoa  = sum((hd_all$p * wv)[in_aoa]) / f, knn_out_aoa = sum((hd_all$p * wv)[!in_aoa]) / f)

    # T2: one covariate at a time --------------------------------------------------------
    X_k1 <- hd_all$X1
    Xs_d0 <- X_d0[si, , drop = FALSE]; Xs_k1 <- X_k1[si, , drop = FALSE]; Xs_bart <- X_bart[si, , drop = FALSE]
    ws <- wv[si]
    p_k1 <- pm(Xs_k1, boots_sub); p_d0s <- pm(Xs_d0, boots_sub); p_bm <- pm(Xs_bart, boots_sub)
    for (v in bf_vars) {
      Xa <- Xs_d0; Xa[[v]] <- Xs_bart[[v]]
      Xb <- Xs_k1; Xb[[v]] <- Xs_bart[[v]]
      t2_rows[[length(t2_rows) + 1L]] <- data.frame(
        species = sp, bcr = bcr, subbasin = s, var = v,
        type = if (v %in% cat_vars) "categorical" else "continuous",
        delta_veg = sum((pm(Xa, boots_sub) - p_d0s) * ws) / f2,
        delta_gap = sum((pm(Xb, boots_sub) - p_k1) * ws) / f2)
    }
    t2_rows[[length(t2_rows) + 1L]] <- data.frame(
      species = sp, bcr = bcr, subbasin = s, var = "(all, jointly)", type = "total",
      delta_veg = sum((p_bm - p_d0s) * ws) / f2, delta_gap = sum((p_bm - p_k1) * ws) / f2)

    # T5: BART's backfill quantile-mapped onto the training land's values ----------------
    X_qm <- Xs_bart
    for (v in cont_bf) {
      tv <- lo_bio[is.finite(lo_bio[, v]), v]
      r  <- (rank(X_qm[[v]], ties.method = "average") - 0.5) / nrow(X_qm)
      X_qm[[v]] <- stats::quantile(tv, r, type = 1, names = FALSE)
    }
    t1_rows[[length(t1_rows)]]$qm_ratio <- sum(pm(X_qm, boots_sub) * ws) / sum(p_bm * ws)

    # distributions of each backfilled covariate ----------------------------------------
    same <- !is.na(lo_bcr) & lo_bcr == bcr
    grp  <- list(train_all = lo_bio[, , drop = FALSE], train_same = lo_bio[same, , drop = FALSE])
    for (v in bf_vars) {
      vals <- Filter(Negate(is.null), list(
        train_all = grp$train_all[, v], train_same = grp$train_same[, v],
        fp_obs = X_obs[[v]], fp_knn = X_k1[[v]], fp_bart = X_bart[[v]],
        fp_bart_new = new_vals[[v]]))
      for (g in names(vals)) {
        d <- if (v %in% cat_vars) dist_cat(vals[[g]], v, g) else dist_cont(vals[[g]], v, g)
        if (is.null(d)) next
        d <- cbind(data.frame(species = sp, bcr = bcr, subbasin = s), d)
        if (v %in% cat_vars) dk_rows[[length(dk_rows) + 1L]] <- d
        else                 dc_rows[[length(dc_rows) + 1L]] <- d
      }
    }
    pool$obs[[length(pool$obs) + 1L]]   <- X_obs[bf_vars]
    pool$knn[[length(pool$knn) + 1L]]   <- X_k1[bf_vars]
    pool$bart[[length(pool$bart) + 1L]] <- X_bart[bf_vars]
    pool$train[[length(pool$train) + 1L]] <- as.data.frame(lo_bio[donor_ok, , drop = FALSE])

    rm(cov_s, lo_s, mv, bf_r, draws, cats, X_obs, X_d0, X_bart, X_lo, X_k1, hd_all, hd_same); gc()
  }

  # T3: survey range and partial dependence, pooled over the case's subbasins -----------
  pooled <- lapply(pool, function(l) dplyr::bind_rows(l))
  t2_case <- dplyr::bind_rows(t2_rows) |> dplyr::filter(species == sp, bcr == !!bcr, type != "total")
  drivers <- t2_case |> dplyr::group_by(var) |>
    dplyr::summarise(g = sum(delta_gap), vg = sum(delta_veg), .groups = "drop")
  top <- unique(c(head(drivers$var[order(drivers$g)], 5L), head(drivers$var[order(drivers$vg)], 5L)))

  for (v in names(pooled$bart)) {
    j <- match(v, model_vars)
    share_out <- function(x) mean(vapply(seq_len(n_boot), function(i) {
      lv <- b.list[[i]]$var.levels[[j]]
      x  <- x[!is.na(x)]
      if (!length(x)) return(NA_real_)
      if (v %in% cat_vars) mean(!(as.character(x) %in% as.character(lv)))
      else mean(x < min(lv) | x > max(lv))
    }, numeric(1)))
    t3_rows[[length(t3_rows) + 1L]] <- data.frame(
      species = sp, bcr = bcr, var = v, type = if (v %in% cat_vars) "categorical" else "continuous",
      top_driver = v %in% top,
      out_train = share_out(pooled$train[[v]]), out_fp_obs = share_out(pooled$obs[[v]]),
      out_fp_knn = share_out(pooled$knn[[v]]), out_fp_bart = share_out(pooled$bart[[v]]))
  }

  # partial dependence of the top drivers: exp(f) averaged over each group's values, relative
  # to the observed footprint, so the number reads as the bird model's multiplier
  for (v in top) {
    for (i in boots_sub) {
      g <- tryCatch(plot.gbm(b.list[[i]], i.var = v, n.trees = b.list[[i]]$n.trees,
                             return.grid = TRUE, continuous.resolution = 50),
                    error = function(e) NULL)
      if (is.null(g)) next
      ev <- function(x) {
        x <- x[!is.na(x)]
        if (!length(x)) return(NA_real_)
        y <- if (v %in% cat_vars) g$y[match(as.character(x), as.character(g[[v]]))]
             else stats::approx(g[[v]], g$y, xout = x, rule = 2)$y
        mean(exp(y), na.rm = TRUE)
      }
      base <- ev(pooled$obs[[v]])
      pd_rows[[length(pd_rows) + 1L]] <- data.frame(
        species = sp, bcr = bcr, var = v, boot = i,
        pd_train = ev(pooled$train[[v]]) / base, pd_fp_knn = ev(pooled$knn[[v]]) / base,
        pd_fp_bart = ev(pooled$bart[[v]]) / base)
    }
  }

  rm(b.list, obs_vals, pooled, pool, super_v, direct_v, zone_v, w_v, pc_obs_v, pc_cf_v); gc()

  # written after every case, so an error later cannot lose a finished one
  write.csv(dplyr::bind_rows(t1_rows), file.path(out_dir, "t1_substitution.csv"), row.names = FALSE)
  write.csv(dplyr::bind_rows(t2_rows), file.path(out_dir, "t2_drivers.csv"), row.names = FALSE)
  write.csv(dplyr::bind_rows(dc_rows), file.path(out_dir, "t2_dist_continuous.csv"), row.names = FALSE)
  write.csv(dplyr::bind_rows(dk_rows), file.path(out_dir, "t2_dist_categorical.csv"), row.names = FALSE)
  write.csv(dplyr::bind_rows(t3_rows), file.path(out_dir, "t3_survey_range.csv"), row.names = FALSE)
  write.csv(dplyr::bind_rows(pd_rows), file.path(out_dir, "t3_partial_dependence.csv"), row.names = FALSE)
}

t1 <- dplyr::bind_rows(t1_rows)
t2 <- dplyr::bind_rows(t2_rows)
t3 <- dplyr::bind_rows(t3_rows)
pd <- dplyr::bind_rows(pd_rows)

# T4 across subbasins: footprint cf/obs against training size and AOA (cid-256 tables) -------

t4 <- dplyr::bind_rows(lapply(unique(vapply(cases, `[[`, "", "species")), function(sp)
  readRDS(file.path(dt_dir, paste0(sp, "_", year, "_coalition_256.rds"))))) |>
  dplyr::filter(bcr %in% c("can11", "can13"), is.finite(obs_on_coalition_mean),
                obs_on_coalition_mean >= 100) |>
  dplyr::left_join(flags, by = "subbasin") |>
  dplyr::mutate(cf_obs = bf_on_coalition_mean / obs_on_coalition_mean,
                ratio_hi_lo = n_highhf / n_lowhf) |>
  dplyr::select(species, bcr, subbasin, HYBAS_ID, obs_on_coalition_mean, bf_on_coalition_mean,
                cf_obs, n_lowhf, n_highhf, ratio_hi_lo, frac_outside_aoa)

# verdict ------------------------------------------------------------------------------

verdict <- t1 |>
  dplyr::group_by(species, bcr) |>
  dplyr::summarise(
    knn_same_vs_all = sum(knn_same, na.rm = TRUE) / sum(knn_all[!is.na(knn_same)]),
    n_sub_same = sum(!is.na(knn_same)),
    frac_obs_outside_aoa = sum(obs_out_aoa) / sum(obs),
    bart_check_ratio = (sum(bart_check) / sum(bart_check_ref)) / (sum(obs_check) / sum(obs_check_ref)),
    bart_chk_obs = sum(bart_check) / sum(obs_check), bart_new_obs = sum(bart_new_check) / sum(obs_check),
    bart_qm = sum(bart * qm_ratio),
    obs = sum(obs), d0 = sum(d0), bart = sum(bart), knn_all = sum(knn_all),
    bart_obs_in_aoa  = sum(bart_in_aoa) / sum(obs_in_aoa),
    bart_obs_out_aoa = sum(bart_out_aoa) / sum(obs_out_aoa),
    knn_obs_in_aoa   = sum(knn_in_aoa) / sum(obs_in_aoa),
    knn_obs_out_aoa  = sum(knn_out_aoa) / sum(obs_out_aoa),
    .groups = "drop") |>
  dplyr::mutate(bart_obs = bart / obs, knn_obs = knn_all / obs, bart_vs_knn = bart / knn_all,
                qm_obs = bart_qm / obs, qm_vs_bart = bart_qm / bart)

top_gap <- t2 |> dplyr::filter(type != "total") |>
  dplyr::group_by(species, bcr, var) |> dplyr::summarise(gap = sum(delta_gap), veg = sum(delta_veg), .groups = "drop")

verdict$label <- vapply(seq_len(nrow(verdict)), function(r) {
  v  <- verdict[r, ]
  tg <- top_gap |> dplyr::filter(species == v$species, bcr == v$bcr) |> dplyr::arrange(gap) |> head(3)
  ex <- t3 |> dplyr::filter(species == v$species, bcr == v$bcr, var %in% tg$var, out_fp_bart > cut)
  lab <- c(
    if (isTRUE(v$bart_vs_knn < 1 - cut)) sprintf("BART artefact: cf %.0f%% below kNN from the same training land (drivers: %s)",
                                         100 * (1 - v$bart_vs_knn), paste(tg$var, collapse = ", ")),
    if (isTRUE(v$n_sub_same > 0 && abs(v$knn_same_vs_all - 1) > cut))
      sprintf("cross-BCR training land: same-BCR donors give %.2fx", v$knn_same_vs_all),
    if (nrow(ex)) sprintf("bird-model extrapolation on %s", paste(ex$var, collapse = ", ")),
    if (isTRUE(v$qm_vs_bart > 1 + cut))
      sprintf("backfill misses the training land's zeros: mapped onto it, cf is %.2fx BART's", v$qm_vs_bart),
    if (isTRUE(abs(v$bart_obs_out_aoa / v$bart_obs_in_aoa - 1) > cut))
      sprintf("abiotic extrapolation: cf/obs %.2f outside the AOA vs %.2f inside",
              v$bart_obs_out_aoa, v$bart_obs_in_aoa))
  if (!length(lab)) "real signal, given the premise that remnant land stands in for pre-industrial land"
  else paste(lab, collapse = "; ")
}, character(1))

write.csv(t4, file.path(out_dir, "t4_across_subbasins.csv"), row.names = FALSE)
write.csv(verdict, file.path(out_dir, "verdict.csv"), row.names = FALSE)

message("\n== T4: Spearman of footprint cf/obs across subbasins (obs >= 100 birds) ==")
print(as.data.frame(t4 |> dplyr::group_by(species, bcr) |> dplyr::summarise(
  n = dplyr::n(),
  rho_ratio = stats::cor(cf_obs, ratio_hi_lo, method = "spearman", use = "complete.obs"),
  rho_aoa   = stats::cor(cf_obs, frac_outside_aoa, method = "spearman", use = "complete.obs"),
  .groups = "drop")), digits = 3)
message("\n== Top drivers of the BART - kNN gap (birds) ==")
print(as.data.frame(top_gap |> dplyr::group_by(species, bcr) |> dplyr::slice_min(gap, n = 4)), digits = 3)
if (nzchar(bf_new_dir)) {
  message("\n== Candidate backfill vs production (direct footprint, share of observed birds kept) ==")
  print(as.data.frame(verdict |> dplyr::select(species, bcr, bart_chk_obs, bart_new_obs, knn_obs)), digits = 3)
}
message("\n== Verdict ==")
print(as.data.frame(verdict), digits = 3)
