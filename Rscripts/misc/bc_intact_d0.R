# ---
# title: Intact land's density at footprint 0, BC
# author: Mannfred Boehm
# ---
#
# bc_intact_benchmark.R compares backfilled density with intact (low-HF) land at the same
# elevation. bf has every footprint covariate at exactly 0, but only 16-32% of low-HF land has
# CanHF_1km == 0, and the bird models step down between a little footprint and none (TODO.md,
# VATH). Observed density on intact land is therefore not comparable with bf. This predicts
# intact land the way bf is predicted, its own vegetation with the footprint covariates at 0, so
# that vegetation is what differs. Same path as bc_d0_pixels.R, whose d0 reproduces 12F's.
#
# Low-HF pixels of each BC (subbasin, BCR) unit, at most n_cap per 100 m elevation band, all 32
# bootstraps. n_land keeps each band's true area, since the sample is stratified.
#
# Runs locally from the RProject root, after bc_d0_tables.R. Cached per species x BCR.

library(terra)
library(gbm)
library(dplyr)

cc    <- FALSE
local <- TRUE
ia_dir  <- getwd()
nm_root <- "G:/Shared drives/BAM_NationalModels5"
year    <- 2020L
bin_m   <- 100
n_cap   <- 40L
categorical_responses <- c("ABoVE_1km", "NLCD_1km", "MODISLCC_1km", "MODISLCC_5x5",
                           "SCANFI_1km", "VLCE_1km")

pred    <- file.path(ia_dir, "data", "derived_data", "predictions")
out_dir <- file.path(ia_dir, "data", "derived_data", "bc_d0_pixels", "intact")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
cache_of <- function(sp, bcr) file.path(out_dir, sprintf("%s_%s.rds", sp, bcr))

source(file.path(ia_dir, "Rscripts", "12F_predict_species_all_coalitions.R"))  # as_model_factors()
source(file.path(ia_dir, "Rscripts", "12G_gbm_tree_walk.R"))
Rcpp::sourceCpp(file.path(ia_dir, "Rscripts", "12G_gbm_tree_walk.cpp"))

pm   <- tibble::tibble(BAMexploreR::predictor_metadata)
dist <- pm$predictor[pm$version == "v5" & pm$predictor_class == "Disturbance"]
load(file.path(ia_dir, "data", "raw_data", "SpeciesPredictionTruncationValues.Rdata"))

units <- readRDS(file.path(ia_dir, "data", "derived_data", "bc_d0_direct_pixels.rds")) |>
  distinct(species, bcr, sid)
dem  <- rast(file.path(ia_dir, "data", "raw_data", "dem_bc_1km.tif"))
low  <- rast(file.path(ia_dir, "data", "raw_data", "hirshpearson", "CanHF_1km_lessthan1.tif"))
subs <- vect(file.path(ia_dir, "data", "raw_data", "hydrobasins_masked_merged_subset.gpkg"))
subs$sid <- seq_len(nrow(subs))

for (bcr in sort(unique(units$bcr))) {
  sps <- unique(units$species[units$bcr == bcr])
  sps <- sps[!file.exists(cache_of(sps, bcr))]
  if (length(sps) == 0L) next

  stack_tmp <- file.path(tempdir(), sprintf("%s_%d.tif", bcr, year))
  file.copy(file.path(nm_root, "gis", "stacks", basename(stack_tmp)), stack_tmp, overwrite = TRUE)
  stk  <- rast(stack_tmp)
  tmpl <- stk[[1]]

  sids <- unique(units$sid[units$bcr == bcr])
  zone <- values(rasterize(subs[sids, ], tmpl, field = "sid"), mat = FALSE)
  lowv <- values(project(low, tmpl, method = "near"), mat = FALSE)
  elev <- values(project(dem, tmpl, method = "near"), mat = FALSE)
  land <- which(!is.na(zone) & !is.na(lowv) & lowv == 1 & !is.na(elev))
  bin  <- floor(elev[land] / bin_m)
  key  <- paste(zone[land], bin)

  set.seed(1)                                            # one sample per BCR, shared by species
  pick <- unlist(lapply(split(land, key), function(x) if (length(x) > n_cap) sample(x, n_cap) else x),
                 use.names = FALSE)

  for (sp in sps) {
    t0 <- Sys.time()
    w_all <- values(rast(file.path(pred, sp, bcr, year, "weight.tif")), mat = FALSE)
    o_all <- values(rast(file.path(pred, sp, bcr, year, "observed_mean.tif")), mat = FALSE)
    n_land <- tibble(key = key, in_w = w_all[land] > 0) |>
      group_by(key) |> summarise(n_land = sum(in_w), .groups = "drop")

    cells <- pick[w_all[pick] > 0]
    tp  <- readRDS(file.path(pred, sp, "truncation_params.rds"))[[paste(bcr, year, sep = "_")]]
    qsp <- q.out$densmax[q.out$spp == sp]
    stopifnot(length(qsp) == 1L, is.finite(tp$q99), isTRUE(all.equal(tp$densmax, qsp)))

    e <- new.env(parent = emptyenv())
    load(file.path(nm_root, "output", "06_bootstraps", sp, paste0(sp, "_", bcr, ".Rdata")), envir = e)
    b <- e$b.list; rm(e); invisible(gc(verbose = FALSE))
    vars   <- b[[1]]$var.names
    dist_b <- intersect(dist, vars)
    cat_b  <- intersect(vars, categorical_responses)

    X <- values(stk[[vars]])[cells, , drop = FALSE]
    X[is.nan(X)] <- NA_real_
    X <- as.data.frame(X, check.names = FALSE)
    for (v in dist_b) X[[v]] <- 0

    d0 <- matrix(0, length(cells), length(b))
    for (i in seq_along(b)) {
      Xi <- as_model_factors(X, b[[i]], cat_b)
      D  <- gbm_design(b[[i]], Xi)
      if (i == 1L) gbm_check_fast(b[[i]], Xi, D, label = paste(sp, bcr))
      d0[, i] <- pmin(pmin(gbm_predict_design(b[[i]], D), qsp), tp$q99)
    }
    w  <- w_all[cells]
    fr <- data.frame(species = sp, bcr = bcr, cell = cells, sid = zone[cells],
                     bin = floor(elev[cells] / bin_m), w = w,
                     d0 = rowMeans(d0) * w * 100, obs = o_all[cells] * w * 100)
    saveRDS(list(frame = fr, n_land = n_land), cache_of(sp, bcr))
    message(format(Sys.time(), "%H:%M"), " ", sp, " ", bcr, ": ", nrow(fr), " intact px | obs ",
            round(mean(fr$obs, na.rm = TRUE), 2), " vs footprint-0 ", round(mean(fr$d0), 2),
            " birds/km2 | ", format(round(Sys.time() - t0, 1)))
    rm(b, X, d0, fr); invisible(gc(verbose = FALSE))
  }
  unlink(stack_tmp)
}
message("INTACT_D0_DONE")
