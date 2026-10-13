# Does 08A's `_isNA` flag make the backfilled broadleaf cover treed? (2026-10-12)
#
# 08A fills a predictor's missing values with the training median and adds a `<col>_isNA` flag.
# A biotic predictor is backfilled before the covariates that use it, and the backfill is never
# NA, so its flag is 0 on every backfilled pixel. SCANFIheightcv_* is NA on ~95% of subbasin 56's
# training land (no or almost no trees), so the models after it (prcC, prcD) may read every
# backfilled pixel as one of the few treed training pixels.
#
# Replays 08A's continuous cascade up to TARGET in one subbasin with production's seeds and
# settings. Training never depends on the backfill, so each covariate is fitted ONCE:
#   A  production: backfilled predictors never NA (gbart's own x.test, checked against the
#      production backfill tif, which proves the replay)
#   B  modelled blanks: for each biotic covariate that is missing on part of the training land,
#      a probit BART (same predictors) predicts where it is missing; backfilled pixels with
#      P(missing) > 0.5 are set NA, so the flag reads 1 there, as on the training land
# Arms A and B are both predicted from the same fits with predict(), so they differ only in
# their inputs. Output: logs/diag_isna_flag_refit_S{s}.csv
#
# Run: SUB=56 Rscript Rscripts/misc/diag_isna_flag_refit.R  (~10 min)

suppressPackageStartupMessages({ library(BART); library(terra); library(dplyr) })
ia_dir <- getwd()
s      <- as.integer(Sys.getenv("SUB", "56"))
target <- Sys.getenv("TARGET", "SCANFIprcD_5x5")
year   <- 2020L
draw_tif <- Sys.getenv("BF_TIF", file.path(ia_dir, "cluster_logs", "grassland_diag", "backfill",
                                            sprintf("subbasin_%d_backfill.tif", s)))

# covariate classes, as 07 builds them
predictor_metadata <-
  dplyr::tibble(BAMexploreR::predictor_metadata) |>
  dplyr::filter(version == "v5") |>
  dplyr::select(predictor, definition, predictor_class) |>
  dplyr::mutate(predictor = sub("Year", "year", predictor), predictor = sub("Method", "method", predictor))
soil_covs <- tibble::tibble(
  predictor = paste0(rep(c("cec", "soc"), each = 6), "_",
                     c("0-5cm", "100-200cm", "15-30cm", "30-60cm", "5-15cm", "60-100cm"), "_mean_1000"),
  predictor_class = "Soil Properties")
actually_biotic_what <- c("Peatland_5x5", "Peatland_1km")
abiotic_vars <- predictor_metadata |>
  dplyr::filter(predictor_class %in% c("Annual Climate", "Climate Normals", "Topography", "Wetland", "Disturbance")) |>
  tibble::add_row(predictor = "CAfire", predictor_class = "Time Since Disturbance") |>
  dplyr::filter(!(predictor %in% actually_biotic_what)) |>
  dplyr::bind_rows(soil_covs)
biotic_vars <- predictor_metadata |>
  dplyr::filter(!(predictor_class %in% c(abiotic_vars$predictor_class, "Time", "Method"))) |>
  dplyr::bind_rows(tibble::tibble(predictor = actually_biotic_what, predictor_class = "Wetland"))
neworder <- readRDS(file.path(ia_dir, "data", "raw_data", "biotic_variable_hierarchy.rds"))
biotic_vars <- biotic_vars[match(neworder, biotic_vars$predictor), ]
categorical_responses <- c("ABoVE_1km", "NLCD_1km", "MODISLCC_1km", "MODISLCC_5x5", "SCANFI_1km", "VLCE_1km")

# subbasin data, as 07 + 08A build it
stack_y <- rast(file.path(ia_dir, "data", "raw_data", "covariates_mosaiced", sprintf("covariates_mosaiced_%d.tif", year)))
hf <- file.path(ia_dir, "data", "raw_data", "hirshpearson")
lowhf_mask  <- project(rast(file.path(hf, "CanHF_1km_lessthan1.tif")), stack_y, method = "near")
highhf_mask <- project(rast(file.path(hf, "CanHF_1km_morethan1.tif")), stack_y, method = "near")
subs <- project(vect(file.path(ia_dir, "data", "raw_data", "hydrobasins_masked_merged_subset.gpkg")), stack_y)
subbasin_s <- subs[s]
cov_s <- mask(crop(stack_y, subbasin_s), subbasin_s)
lowhf_mask_s  <- mask(resample(lowhf_mask,  cov_s, method = "near"), subbasin_s)
highhf_mask_s <- mask(resample(highhf_mask, cov_s, method = "near"), subbasin_s)
df_train <- as.data.frame(mask(cov_s, lowhf_mask_s), xy = TRUE, na.rm = FALSE)
df_full  <- as.data.frame(cov_s, xy = TRUE, na.rm = FALSE, cells = TRUE)
backfill_idx <- which(values(highhf_mask_s) == 1)
df_backfill  <- df_full[backfill_idx, , drop = FALSE]
xy_ok <- which(!is.na(df_full$x) & !is.na(df_full$y))
xy_center <- colMeans(df_full[xy_ok, c("x", "y")]); xy_scale <- apply(df_full[xy_ok, c("x", "y")], 2, sd)
for (d in c("df_train", "df_backfill")) {
  z <- get(d)
  z[, c("x", "y")] <- sweep(sweep(z[, c("x", "y")], 2, xy_center, "-"), 2, xy_scale, "/")
  assign(d, z)
}
footprint_cols <- intersect(names(df_backfill), abiotic_vars$predictor[abiotic_vars$predictor_class == "Disturbance"])
df_backfill[, footprint_cols] <- 0
# 08A keeps every pixel row and filters per covariate with `idx`; the training rows are those
# with any value (mask() leaves NA rows outside the low-HF mask)
lo_rows <- which(values(lowhf_mask_s, mat = FALSE) == 1)
df_train <- df_train[lo_rows, , drop = FALSE]

abiotic_cols <- intersect(names(df_train), abiotic_vars$predictor)
biotic_cols  <- intersect(names(df_train), biotic_vars$predictor)
biotic_cols  <- as.character(na.omit(biotic_cols[match(neworder, biotic_cols)]))
biotic_cols_cont <- setdiff(biotic_cols, categorical_responses)
run_cols <- biotic_cols_cont[seq_len(match(target, biotic_cols_cont))]
message(sprintf("subbasin %d: %d training, %d backfilled pixels; cascade: %s", s, nrow(df_train),
                nrow(df_backfill), paste(run_cols, collapse = ", ")))

# 08A's predictor preparation, split into a recipe (from the training rows and production's
# backfill) and its application, so arm B gets exactly the columns, fills and flags the fit saw
make_recipe <- function(tr, bf) {
  tr <- tr[, colSums(!is.na(tr)) > 0, drop = FALSE]
  col_sd <- sapply(tr, function(x) sd(as.numeric(x), na.rm = TRUE))
  tr <- tr[, !is.na(col_sd) & col_sd > 0, drop = FALSE]
  bf <- bf[, colnames(tr), drop = FALSE]
  keep <- colnames(bf)[colSums(!is.na(bf)) > 0]
  tr <- tr[, keep, drop = FALSE]; bf <- bf[, keep, drop = FALSE]
  na_cols <- setdiff(names(tr)[colSums(is.na(tr)) > 0L | colSums(is.na(bf)) > 0L], c("x", "y"))
  fills <- vapply(na_cols, function(cna) { f <- stats::median(tr[[cna]], na.rm = TRUE)
                                          if (is.finite(f)) f else 0 }, numeric(1))
  flagged <- na_cols[vapply(na_cols, function(cna) length(unique(as.integer(is.na(tr[[cna]])))) > 1L, logical(1))]
  list(cols = keep, fills = fills, flagged = flagged)
}
apply_recipe <- function(rc, df) {
  x <- df[, rc$cols, drop = FALSE]
  for (cna in names(rc$fills)) {
    if (cna %in% rc$flagged) x[[paste0(cna, "_isNA")]] <- as.integer(is.na(x[[cna]]))
    x[[cna]][is.na(x[[cna]])] <- rc$fills[[cna]]
  }
  # 08A appends each flag right after the loop's fills; it then reorders the backfill to the
  # training columns, which carry the flags in the same append order
  if (anyNA(x)) stop("NA left in the design: ", paste(names(x)[colSums(is.na(x)) > 0], collapse = ", "))
  as.matrix(x)
}
summ <- function(v, arm, x_draws = NULL, x_val = NULL) {
  if (!is.null(x_draws)) { d <- pmax(expm1(x_draws), 0)
    data.frame(var = v, arm = arm, share_na = NA_real_, share_zero = mean(d == 0),
               median = median(d), mean = mean(d))
  } else data.frame(var = v, arm = arm, share_na = mean(is.na(x_val)),
                    share_zero = mean(x_val == 0, na.rm = TRUE), median = median(x_val, na.rm = TRUE),
                    mean = mean(x_val, na.rm = TRUE))
}

bf_r <- rast(draw_tif)
if (!compareGeom(bf_r, cov_s, stopOnError = FALSE)) stop("backfill tif grid differs from the subbasin crop")
bf_A <- df_backfill; bf_B <- df_backfill
rows <- list(); check <- list()
for (b in run_cols) {
  t0 <- proc.time()[3]
  idx <- which(!is.na(df_train[[b]]))
  if (length(idx) == 0 || length(unique(df_train[[b]][idx])) < 2) {
    cv <- unique(df_train[[b]][idx]); if (length(cv) == 0) cv <- NA
    bf_A[[b]] <- cv; bf_B[[b]] <- cv
    message(b, ": constant (", cv, ")"); next
  }
  b_before  <- biotic_cols_cont[seq_len(match(b, biotic_cols_cont) - 1)]
  predictors <- c(abiotic_cols, b_before, "x", "y")
  y <- as.numeric(df_train[[b]][idx])
  set.seed(abs(as.integer(sprintf("%d%03d", s, which(biotic_cols == b)))))
  holdout_idx <- sample(seq_len(length(y)), size = round(0.10 * length(y)))
  tr_all <- df_train[idx, predictors, drop = FALSE]
  # 08A imputes on all training rows, THEN deploy_gbart splits off the holdout
  rc  <- make_recipe(tr_all, bf_A[, predictors, drop = FALSE])
  x_tr <- apply_recipe(rc, tr_all)[-holdout_idx, , drop = FALSE]
  x_A  <- apply_recipe(rc, bf_A)
  if (length(unique(log1p(y[-holdout_idx]))) <= 2) stop(b, ": degenerate path not replicated")
  fit <- BART::gbart(x.train = x_tr, y.train = log1p(y[-holdout_idx]), x.test = x_A, type = "wbart",
                     k = 3, ntree = 50L, ndpost = 700L, nskip = 300L, sparse = TRUE,
                     sigest = sd(log1p(y[-holdout_idx])), sigdf = 3, printevery = 10000L)
  draw_idx <- sample(nrow(fit$yhat.test), 100L)

  # replay check: production's stored draws at the backfilled pixels
  lyr  <- paste0(b, "_draw_", sprintf("%03d", 1:100))
  prod <- terra::values(bf_r[[lyr]], mat = TRUE)[df_backfill$cell, , drop = FALSE]
  mine <- t(fit$yhat.test[draw_idx, , drop = FALSE])
  check[[b]] <- data.frame(var = b, max_abs_diff = max(abs(prod - mine), na.rm = TRUE),
                           na_mismatch = sum(is.na(prod) != is.na(mine)))

  cols  <- colnames(fit$varcount)
  pr_A  <- predict(fit, newdata = x_A[, cols, drop = FALSE])
  bf_A[[b]] <- colMeans(expm1(fit$yhat.test))          # production's cascade value
  x_B   <- apply_recipe(rc, bf_B)
  pr_B  <- predict(fit, newdata = x_B[, cols, drop = FALSE])
  bf_B[[b]] <- colMeans(expm1(pr_B))
  rows[[length(rows) + 1]] <- summ(b, "A_production", x_draws = pr_A[draw_idx, ])
  rows[[length(rows) + 1]] <- summ(b, "B_modelled_blanks", x_draws = pr_B[draw_idx, ])

  # arm B: where is b missing? probit BART on ALL training rows, same predictors
  miss <- as.integer(is.na(df_train[[b]]))
  if (length(unique(miss)) > 1L) {
    rcm <- make_recipe(df_train[, predictors, drop = FALSE], bf_B[, predictors, drop = FALSE])
    set.seed(abs(as.integer(sprintf("%d%03d", s, which(biotic_cols == b)))) + 1L)
    pf <- BART::pbart(x.train = apply_recipe(rcm, df_train), y.train = miss,
                      x.test = apply_recipe(rcm, bf_B), ntree = 50L, ndpost = 700L, nskip = 300L,
                      sparse = TRUE, printevery = 10000L)
    p_miss <- pf$prob.test.mean
    bf_B[[b]][p_miss > 0.5] <- NA_real_
    rows[[length(rows) + 1]] <- data.frame(var = b, arm = "B_modelled_blanks", share_na = mean(p_miss > 0.5),
                                           share_zero = NA, median = NA, mean = NA)
  }
  rows[[length(rows) + 1]] <- summ(b, "training_land", x_val = df_train[[b]])
  message(sprintf("%s: %.0fs, replay max |diff| %.2g", b, proc.time()[3] - t0, check[[b]]$max_abs_diff))
}

res <- dplyr::bind_rows(rows)
chk <- dplyr::bind_rows(check)
print(chk, row.names = FALSE)
if (any(chk$max_abs_diff > 1e-4 | chk$na_mismatch > 0)) warning("replay does NOT reproduce production's draws")
res <- res |> dplyr::group_by(var, arm) |>
  dplyr::summarise(dplyr::across(c(share_na, share_zero, median, mean), ~ dplyr::first(stats::na.omit(.x))),
                   .groups = "drop") |>
  dplyr::arrange(match(var, run_cols), arm)
print(as.data.frame(res |> dplyr::mutate(dplyr::across(where(is.numeric), ~ signif(.x, 3)))), row.names = FALSE)
dir.create(file.path(ia_dir, "logs"), showWarnings = FALSE)
write.csv(res, file.path(ia_dir, "logs", sprintf("diag_isna_flag_refit_S%d.csv", s)), row.names = FALSE)
write.csv(chk, file.path(ia_dir, "logs", sprintf("diag_isna_flag_refit_S%d_replay.csv", s)), row.names = FALSE)
