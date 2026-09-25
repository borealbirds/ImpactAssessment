# Ecological sanity check of C2's extreme BCR impacts.
# Part 1 (densities): per species x BCR, weighted observed density (birds/km2) on the
#   footprint (cid 256 superset), on low-HF land (the BART training pixels) and on the rest,
#   against the backfilled density on the footprint (tables' bf_on / footprint area).
# Part 2 (vegetation, BCRs whose backfill mosaic was pulled): what the backfill puts on the
#   footprint, against what is there now and what the low-HF land holds.
suppressPackageStartupMessages({ library(terra); library(dplyr) })
REPO <- "C:/Users/mannf/Drive/boreal_avian_modelling_project/ImpactAssessment"; setwd(REPO)
source("Rscripts/12E_shapley_utils.R")
secs  <- canonical_sectors()
gpkg  <- vect("data/raw_data/hydrobasins_masked_merged_subset.gpkg")
covs  <- rast("data/raw_data/covariates_mosaiced/covariates_mosaiced_2020.tif")
cases <- read.table(text = "
sp   bcr
CAWA can11
CAWA can13
CAWA can60
CAWA can61
OVEN can11
OVEN can13
OVEN can12
OVEN can80
OVEN can81
OVEN can60
OVEN can61", header = TRUE, stringsAsFactors = FALSE)

grid_for <- function(sp, bcr) {
  tmpl <- rast(sprintf("data/derived_data/predictions/%s/%s/2020/observed_mean.tif", sp, bcr))
  w    <- rast(sprintf("cluster_logs/sanity/weights/%s/%s/weight.tif", sp, bcr))
  stopifnot(isTRUE(compareGeom(w, tmpl, stopOnError = FALSE)))
  prj  <- function(f) values(project(rast(file.path("data/raw_data/hirshpearson", f)), tmpl, method = "near"), mat = FALSE)
  P    <- sapply(secs, function(s) { v <- prj(paste0(s, ".tif")); v[is.na(v)] <- 0; v })   # pressure score
  S    <- P > 0
  hi   <- prj("CanHF_1km_morethan1.tif"); lo <- prj("CanHF_1km_lessthan1.tif")
  tab  <- readRDS(sprintf("data/derived_data/density_tables/%s_2020_coalition_256.rds", sp))
  tab  <- tab[tab$bcr == bcr, ]
  z    <- values(rasterize(gpkg[tab$subbasin], tmpl, field = "first_HYBAS_ID"), mat = FALSE)
  fp   <- rowSums(S) > 0 & !is.na(hi) & hi >= 1
  low  <- !is.na(lo) & lo == 1
  tr   <- readRDS(sprintf("data/derived_data/density_tables/%s_2020_coalition_129.rds", sp))   # roads alone
  tr   <- tr[tr$bcr == bcr, ]
  list(tmpl = tmpl, w = values(w, mat = FALSE), obs = values(tmpl, mat = FALSE), S = S, P = P, tr = tr,
       fp = fp & !is.na(z), low = low & !is.na(z), other = !is.na(z) & !fp & !low, tab = tab, z = z)
}

dens <- list(); grids <- list()
for (r in seq_len(nrow(cases))) {
  sp <- cases$sp[r]; bcr <- cases$bcr[r]
  g  <- grid_for(sp, bcr); if (sp == "OVEN" && bcr %in% c("can11", "can12", "can13")) grids[[paste(sp, bcr)]] <- g
  birds <- g$w * g$obs * 100
  A <- function(m) sum(g$w[m], na.rm = TRUE); N <- function(m) sum(birds[m], na.rm = TRUE)
  bf <- sum(g$tab$bf_on_coalition_mean); oo <- sum(g$tab$obs_on_coalition_mean, na.rm = TRUE)
  roads_only <- g$fp & g$S[, "roads"] & rowSums(g$S) == 1
  dens[[r]] <- data.frame(
    sp, bcr,
    area_km2 = round(A(g$fp | g$low | g$other)),
    fp_share = round(A(g$fp) / A(g$fp | g$low | g$other), 3),
    low_share = round(A(g$low) / A(g$fp | g$low | g$other), 3),
    roads_only_share_of_fp = round(A(roads_only) / A(g$fp), 3),
    ag_share_of_fp = round(A(g$fp & (g$S[, "crop"] | g$S[, "pasture"])) / A(g$fp), 3),
    D_low = round(N(g$low) / A(g$low), 2),
    D_other = round(N(g$other) / A(g$other), 2),
    D_fp_obs = round(N(g$fp) / A(g$fp), 2),
    D_fp_bf = round(bf / A(g$fp), 2),
    bf_over_low = round(bf / A(g$fp) / (N(g$low) / A(g$low)), 2),
    tab_vs_raster_obs_fp = round(oo / N(g$fp), 4),
    impact_pct = round(100 * (bf - oo) / sum(g$tab$obs_total_mean), 1),
    fp_maxP_lt1 = round(A(g$fp & apply(g$P, 1, max) < 1) / A(g$fp), 3),
    fp_maxP_lt4 = round(A(g$fp & apply(g$P, 1, max) < 4) / A(g$fp), 3),
    roads_only_medP = round(median(g$P[roads_only, "roads"]), 2),
    D_roads_obs = round(sum(g$tr$obs_on_coalition_mean, na.rm = TRUE) / A(g$fp & g$S[, "roads"]), 2),
    D_roads_bf  = round(sum(g$tr$bf_on_coalition_mean) / A(g$fp & g$S[, "roads"]), 2))
  message(sp, " ", bcr, " done")
}
dens <- bind_rows(dens)
print(dens, row.names = FALSE)
saveRDS(list(dens = dens), "cluster_logs/sanity/diag_extreme_bcr_density_vegetation.rds")

# ---- Part 2: vegetation on the footprint ----
cat_vars  <- c("MODISLCC_1km", "SCANFI_1km", "VLCE_1km")
cont_vars <- c("SCANFIheight_1km", "SCANFIprcD_1km", "SCANFIprcC_1km", "SCANFIbiomass_5x5")
veg <- list()
for (bcr in c("can11", "can12", "can13")) {
  key <- grep(paste0(" ", bcr, "$"), names(grids), value = TRUE)[1]; g <- grids[[key]]
  bfm <- rast(sprintf("cluster_logs/sanity/mosaics/%s_backfilled.tif", bcr))
  stopifnot(isTRUE(compareGeom(bfm, g$tmpl, stopOnError = FALSE)))
  obs_l <- function(v) values(project(crop(covs[[v]], project(ext(g$tmpl), crs(g$tmpl), crs(covs)) ), g$tmpl, method = "near"), mat = FALSE)
  cells_fp <- which(g$fp & g$w > 0); cells_lo <- which(g$low & g$w > 0)
  wf <- g$w[cells_fp]; wl <- g$w[cells_lo]
  for (v in cat_vars) {
    o  <- obs_l(v); b <- values(bfm[[v]], mat = FALSE)
    sh <- function(x, w) { t <- tapply(w, x, sum); round(100 * t / sum(w[!is.na(x)]), 1) }
    cls <- sort(unique(c(o[cells_fp], b[cells_fp], o[cells_lo])))
    veg[[length(veg) + 1]] <- data.frame(bcr, var = v, class = cls,
      fp_obs = as.numeric(sh(o[cells_fp], wf)[as.character(cls)]),
      fp_bf  = as.numeric(sh(b[cells_fp], wf)[as.character(cls)]),
      low_obs = as.numeric(sh(o[cells_lo], wl)[as.character(cls)]))
  }
  for (v in cont_vars) {
    o   <- obs_l(v)
    lyr <- paste0(v, "_draw_", sprintf("%03d", 1:20))
    lyr <- intersect(lyr, names(bfm)); if (length(lyr) == 0) next
    d   <- terra::extract(bfm[[lyr]], cells_fp); d <- expm1(as.matrix(d)); d[!is.finite(d)] <- NA
    bm  <- rowMeans(d, na.rm = TRUE); ok <- is.finite(bm)
    wm  <- function(x, w) { k <- is.finite(x); sum(x[k] * w[k]) / sum(w[k]) }
    veg[[length(veg) + 1]] <- data.frame(bcr, var = v, class = NA,
      fp_obs = wm(o[cells_fp], wf), fp_bf = wm(bm, wf), low_obs = wm(o[cells_lo], wl),
      bf_coverage = round(mean(ok), 3))
  }
  message(bcr, " vegetation done")
}
veg <- bind_rows(veg)
print(veg, row.names = FALSE)
saveRDS(list(dens = dens, veg = veg), "cluster_logs/sanity/diag_extreme_bcr_density_vegetation.rds")
