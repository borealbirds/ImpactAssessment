# Split C2's footprint impact into (a) the bird model's DIRECT response to the footprint
# covariates 12F zeroes on the bf side (CanHF, canroad, CCNL: observed vegetation, footprint
# covariates set to 0) and (b) the VEGETATION effect (backfilled biotic covariates):
#   impact = [N(d0) - N(obs)] + [N(bf) - N(d0)]
# N(obs) from observed_bootstraps.tif, N(d0) re-predicted here with 12F's exact path (12G walk,
# per-model factor levels, NaN -> NA, both caps, weight), N(bf) from the C1 tables.
# Every 4th bootstrap (8 of 32) to keep it cheap; obs and d0 use the same 8.
suppressPackageStartupMessages({ library(terra); library(dplyr); library(gbm) })
REPO <- "C:/Users/mannf/Drive/boreal_avian_modelling_project/ImpactAssessment"; setwd(REPO)
SAN  <- file.path(REPO, "cluster_logs/sanity")
source("Rscripts/12E_shapley_utils.R"); source("Rscripts/12F_predict_species_all_coalitions.R")
source("Rscripts/12G_gbm_tree_walk.R"); Rcpp::sourceCpp("Rscripts/12G_gbm_tree_walk.cpp")
load("data/raw_data/SpeciesPredictionTruncationValues.Rdata")
categorical_responses <- c("ABoVE_1km", "NLCD_1km", "MODISLCC_1km", "MODISLCC_5x5", "SCANFI_1km", "VLCE_1km")
pm <- BAMexploreR::predictor_metadata; pm <- pm[pm$version == "v5", ]
dist_vars <- pm$predictor[pm$predictor_class == "Disturbance"]
secs <- canonical_sectors(); gpkg <- vect("data/raw_data/hydrobasins_masked_merged_subset.gpkg")
cases <- data.frame(
  sp    = c("CAWA", "CAWA", "OVEN", "OVEN", "OVEN"),
  bcr   = c("can11", "can13", "can11", "can12", "can13"),
  rdata = c(file.path(SAN, "models/CAWA_can11.Rdata"), file.path(SAN, "models/CAWA_can13.Rdata"),
            file.path(SAN, "models/OVEN_can11.Rdata"), file.path(SAN, "models/OVEN_can12.Rdata"),
            file.path(SAN, "models/OVEN_can13.Rdata")),
  stack = c(file.path(SAN, "stacks/can11_2020.tif"), file.path(SAN, "stacks/can13_2020.tif"),
            file.path(SAN, "stacks/can11_2020.tif"), file.path(SAN, "stacks/can12_2020.tif"),
            file.path(SAN, "stacks/can13_2020.tif")), stringsAsFactors = FALSE)
if (nzchar(Sys.getenv("ONLY"))) cases <- cases[paste(cases$sp, cases$bcr) %in% strsplit(Sys.getenv("ONLY"), ",")[[1]], ]
boots <- seq(1, 32, by = 4)
# coalitions reported: all 8, roads alone, crop + pasture
cids <- c(all = 256, roads = sectors_to_coalition_id("roads", secs),
          agric = sectors_to_coalition_id(c("crop", "pasture"), secs))

out <- list()
for (r in seq_len(nrow(cases))) {
  sp <- cases$sp[r]; bcr <- cases$bcr[r]; t0 <- Sys.time()
  stack_obs <- rast(cases$stack[r])
  prj <- function(f) values(project(rast(file.path("data/raw_data/hirshpearson", f)), stack_obs[[1]], method = "near"), mat = FALSE)
  S   <- sapply(secs, function(s) { v <- prj(paste0(s, ".tif")) > 0; v[is.na(v)] <- FALSE; v })
  hi  <- prj("CanHF_1km_morethan1.tif")
  tab <- readRDS(sprintf("data/derived_data/density_tables/%s_2020_coalition_256.rds", sp)); tab <- tab[tab$bcr == bcr, ]
  z   <- values(rasterize(gpkg[tab$subbasin], stack_obs[[1]], field = "first_HYBAS_ID"), mat = FALSE)
  w   <- values(rast(file.path(SAN, "weights", sp, bcr, "weight.tif")), mat = FALSE)
  sup <- which(rowSums(S) > 0 & !is.na(hi) & hi >= 1 & !is.na(z) & is.finite(w) & w > 0)
  mem <- sapply(cids, function(cid) Reduce(`|`, lapply(coalition_id_to_sectors(cid, secs), function(s) S[sup, s])))
  e <- new.env(); load(cases$rdata[r], envir = e); b.list <- e$b.list; rm(e)
  vars <- b.list[[1]]$var.names; cat_v <- intersect(vars, categorical_responses); d_v <- intersect(dist_vars, vars)
  X <- values(stack_obs[[unique(unlist(lapply(b.list, `[[`, "var.names")))]])[sup, , drop = FALSE]
  X[is.nan(X)] <- NA_real_; X <- as.data.frame(X, check.names = FALSE)
  qsp <- q.out[q.out$spp == sp, ]$densmax; q99 <- readRDS(sprintf("data/derived_data/predictions/%s/truncation_params.rds", sp))[[paste0(bcr, "_2020")]]$q99
  ob  <- rast(sprintf("data/derived_data/predictions/%s/%s/2020/observed_bootstraps.tif", sp, bcr))
  ri  <- sapply(boots, function(i) { s <- summary(b.list[[i]], plotit = FALSE); sum(s$rel.inf[s$var %in% d_v]) })
  res <- matrix(NA_real_, length(boots), 2 * length(cids), dimnames = list(NULL, c(paste0("obs_", names(cids)), paste0("d0_", names(cids)))))
  for (bi in seq_along(boots)) {
    i <- boots[bi]; m <- b.list[[i]]
    Xd <- X; for (v in d_v) Xd[[v]] <- 0
    Xd <- as_model_factors(Xd, m, cat_v)
    pd <- pmin(pmin(gbm_predict_design(m, gbm_design(m, Xd)), qsp), q99) * w[sup] * 100
    po <- values(ob[[i]], mat = FALSE)[sup] * w[sup] * 100; po[is.na(po)] <- 0
    res[bi, ] <- c(colSums(po * mem), colSums(pd * mem))
  }
  bf <- sapply(cids, function(cid) { t <- readRDS(sprintf("data/derived_data/density_tables/%s_2020_coalition_%d.rds", sp, cid)); sum(t$bf_on_coalition_mean[t$bcr == bcr]) })
  oo <- sapply(cids, function(cid) { t <- readRDS(sprintf("data/derived_data/density_tables/%s_2020_coalition_%d.rds", sp, cid)); sum(t$obs_on_coalition_mean[t$bcr == bcr], na.rm = TRUE) })
  cm <- colMeans(res)
  for (k in names(cids)) out[[length(out) + 1]] <- data.frame(sp, bcr, coalition = k,
    obs_total = round(sum(tab$obs_total_mean)), obs_on_tab = round(oo[k]), obs_on_8boot = round(cm[paste0("obs_", k)]),
    d0 = round(cm[paste0("d0_", k)]), bf = round(bf[k]),
    direct = round(cm[paste0("d0_", k)] - cm[paste0("obs_", k)]), vegetation = round((bf[k] - oo[k]) - (cm[paste0("d0_", k)] - cm[paste0("obs_", k)])),
    total = round(bf[k] - oo[k]), dist_vars_in_model = paste(d_v, collapse = "+"),
    dist_relinf_pct = round(mean(ri), 1), row.names = NULL)
  message(sp, " ", bcr, ": ", length(sup), " weight>0 footprint pixels, ", round(difftime(Sys.time(), t0, units = "mins"), 1), " min")
  rm(b.list, X, S); gc()
}
out <- bind_rows(out); print(out, row.names = FALSE)
saveRDS(out, file.path(SAN, "diag_extreme_bcr_decomposition.rds"))
