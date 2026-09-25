# How much do the footprint covariates (V5 "Disturbance" class) move the BART backfill?
# One subbasin at a time, a few key vegetation covariates, 08A's training set (low-HF pixels,
# abiotic predictors + scaled x, y; preceding biotics omitted for simplicity, in all variants):
#   A_obs : predictors incl. footprint, backfilled with the pixels' OBSERVED footprint (08A today)
#   A_0   : same model, footprint covariates set to 0 at the backfilled pixels
#   C     : footprint covariates dropped from the predictors
# Reported on the raw scale (expm1 of each log1p draw, as 12F uses them), averaged over pixels.
suppressPackageStartupMessages({ library(terra); library(BART) })
REPO <- "C:/Users/mannf/Drive/boreal_avian_modelling_project/ImpactAssessment"; setwd(REPO)
set.seed(1)
subs <- as.integer(strsplit(Sys.getenv("SUBS", "57,98"), ",")[[1]])
resp <- c("SCANFIprcD_1km", "SCANFIheight_1km", "SCANFIbiomass_5x5")
pm <- BAMexploreR::predictor_metadata; pm <- pm[pm$version == "v5", ]
soil <- c("cec_0-5cm_mean_1000", "cec_100-200cm_mean_1000", "cec_15-30cm_mean_1000", "cec_30-60cm_mean_1000",
          "cec_5-15cm_mean_1000", "cec_60-100cm_mean_1000", "soc_0-5cm_mean_1000", "soc_100-200cm_mean_1000",
          "soc_15-30cm_mean_1000", "soc_30-60cm_mean_1000", "soc_5-15cm_mean_1000", "soc_60-100cm_mean_1000")
abio <- c(pm$predictor[pm$predictor_class %in% c("Annual Climate", "Climate Normals", "Topography", "Wetland", "Disturbance") &
                      !pm$predictor %in% c("Peatland_5x5", "Peatland_1km")], "CAfire", soil)
dist <- pm$predictor[pm$predictor_class == "Disturbance"]
cov  <- rast("data/raw_data/covariates_mosaiced/covariates_mosaiced_2020.tif")
abio <- intersect(abio, names(cov))
gp   <- vect("data/raw_data/hydrobasins_masked_merged_subset.gpkg")
lo_m <- rast("data/raw_data/hirshpearson/CanHF_1km_lessthan1.tif"); hi_m <- rast("data/raw_data/hirshpearson/CanHF_1km_morethan1.tif")
out <- list()
for (s in subs) {
  sub <- gp[s]; cs <- mask(crop(cov[[c(abio, resp)]], sub), sub)
  lo <- values(mask(resample(lo_m, cs, method = "near"), sub), mat = FALSE) == 1
  hi <- values(mask(resample(hi_m, cs, method = "near"), sub), mat = FALSE) == 1
  df <- as.data.frame(cs, xy = TRUE, na.rm = FALSE)
  ok <- !is.na(df$x); df$x <- (df$x - mean(df$x[ok])) / sd(df$x[ok]); df$y <- (df$y - mean(df$y[ok])) / sd(df$y[ok])
  li <- which(lo %in% TRUE); hiI <- which(hi %in% TRUE)
  if (length(hiI) > 20000) hiI <- sort(sample(hiI, 20000))
  for (b in resp) {
    tr <- li[!is.na(df[[b]][li])]
    if (length(tr) < 50L) { message(s, " ", b, " | only ", length(tr), " training pixels - skipped"); next }
    P  <- c(abio, "x", "y")
    Xtr <- df[tr, P, drop = FALSE]; Xhi <- df[hiI, P, drop = FALSE]
    keep <- vapply(P, function(v) sum(!is.na(Xtr[[v]])) > 1 && isTRUE(sd(Xtr[[v]], na.rm = TRUE) > 0) && any(!is.na(Xhi[[v]])), NA)
    P <- P[keep]; Xtr <- Xtr[, P, drop = FALSE]; Xhi <- Xhi[, P, drop = FALSE]
    for (v in P) { m <- median(Xtr[[v]], na.rm = TRUE); Xtr[[v]][is.na(Xtr[[v]])] <- m; Xhi[[v]][is.na(Xhi[[v]])] <- m }
    dP <- intersect(dist, P)
    y  <- log1p(df[[b]][tr])
    fa <- tryCatch(gbart(as.matrix(Xtr), y, ndpost = 200, nskip = 200, printevery = 1e6),
                   error = function(e) { message(s, " ", b, " | gbart failed: ", conditionMessage(e),
                                                 " (n_train ", length(tr), ", p ", ncol(Xtr), ")"); NULL })
    if (is.null(fa)) next
    X0 <- Xhi; for (v in dP) X0[[v]] <- 0
    pa  <- predict(fa, as.matrix(Xhi)); p0 <- predict(fa, as.matrix(X0))
    Pn  <- setdiff(P, dist)
    fc  <- gbart(as.matrix(Xtr[, Pn]), y, ndpost = 200, nskip = 200, printevery = 1e6)
    pc  <- predict(fc, as.matrix(Xhi[, Pn]))
    rng <- sapply(dP, function(v) sprintf("%s train [%.2g,%.2g] backfill median %.2g", v, min(Xtr[[v]]), max(Xtr[[v]]), median(Xhi[[v]])))
    vc  <- fa$varcount.mean; share <- sum(vc[names(vc) %in% dP]) / sum(vc)
    out[[length(out) + 1]] <- data.frame(subbasin = s, var = b, n_train = length(tr), n_backfill = length(hiI),
      train_obs = mean(df[[b]][tr]), backfill_obs = mean(df[[b]][hiI], na.rm = TRUE),
      A_obs = mean(expm1(pa)), A_0 = mean(expm1(p0)), C = mean(expm1(pc)),
      dist_split_share = round(share, 3), dist_in_model = paste(dP, collapse = "+"))
    message(s, " ", b, " | ", paste(rng, collapse = " ; "))
    saveRDS(do.call(rbind, out), "cluster_logs/sanity/diag_extreme_bcr_bart_footprint.rds")
  }
}
out <- do.call(rbind, out); out[, 5:9] <- signif(out[, 5:9], 4); print(out, row.names = FALSE)
saveRDS(out, "cluster_logs/sanity/diag_extreme_bcr_bart_footprint.rds")
