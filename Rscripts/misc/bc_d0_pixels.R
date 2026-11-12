# ---
# title: Per-pixel d0 for BC's footprint, six western species
# author: Mannfred Boehm
# ---
#
# 12F never stores d0 per pixel: pixel_change/ holds obs and cf, and cf is bf on a coalition's
# direct footprint but d0 on the rest of it. So on direct pixels we have bf and not d0. This
# re-predicts d0 (observed vegetation, footprint covariates at 0) through 12F's own path, for
# every kept footprint pixel in a subbasin that touches BC, with all 32 bootstraps.
#
# Two checks come free. On footprint pixels outside the direct footprint, cf IS 12F's d0, so
# the re-prediction must reproduce it pixel by pixel. On direct pixels, its subbasin sums must
# equal bf - shapley_vegetation (bc_d0_tables.R).
#
# Runs locally from the RProject root. Reads each species x BCR's bootstrap models and the BCR
# stack from G: by exact file name, and from Fir (pulled by Globus): pixel_change/, weight.tif,
# truncation_params.rds. Cached per species x BCR; ONLY=VATH:can9 runs one pair.

library(terra)
library(gbm)
library(dplyr)

cc    <- FALSE
local <- TRUE
ia_dir  <- getwd()
nm_root <- "G:/Shared drives/BAM_NationalModels5"
year    <- 2020L

species <- c("VATH", "TOWA", "SWTH", "RCKI", "RBNU", "OCWA")
bcrs    <- c("can9", "can5", "can10", "can40", "can60", "can61")   # every subunit in BC
sectors <- c("built", "crop", "dam_and_associated_reservoir", "mines", "oil_gas", "pasture",
             "rail", "roads")
categorical_responses <- c("ABoVE_1km", "NLCD_1km", "MODISLCC_1km", "MODISLCC_5x5",
                           "SCANFI_1km", "VLCE_1km")

out_dir <- file.path(ia_dir, "data", "derived_data", "bc_d0_pixels")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
cache_of <- function(sp, bcr) file.path(out_dir, sprintf("%s_%s.rds", sp, bcr))

source(file.path(ia_dir, "Rscripts", "12F_predict_species_all_coalitions.R"))  # as_model_factors()
source(file.path(ia_dir, "Rscripts", "12G_gbm_tree_walk.R"))
Rcpp::sourceCpp(file.path(ia_dir, "Rscripts", "12G_gbm_tree_walk.cpp"))

pm   <- tibble::tibble(BAMexploreR::predictor_metadata)
dist <- pm$predictor[pm$version == "v5" & pm$predictor_class == "Disturbance"]
load(file.path(ia_dir, "data", "raw_data", "SpeciesPredictionTruncationValues.Rdata"))

hp   <- file.path(ia_dir, "data", "raw_data", "hirshpearson")
bc   <- project(vect(file.path(ia_dir, "data", "raw_data", "bc_boundary.gpkg")), "EPSG:5072")
subs <- vect(file.path(ia_dir, "data", "raw_data", "hydrobasins_masked_merged_subset.gpkg"))
subs$sid <- seq_len(nrow(subs))
ss   <- read.csv(file.path(ia_dir, "data", "derived_data", "sector_effects", "shapley_subbasin.csv"))
ss   <- unique(ss[ss$year == year, c("species", "bcr", "subbasin", "HYBAS_ID")])

# 14B's rows only: a species x BCR BAM withheld is not reported, so it is not worked on here
pairs <- unique(ss[ss$species %in% species & ss$bcr %in% bcrs, c("species", "bcr")])
only  <- Sys.getenv("ONLY", "")
if (nzchar(only)) {
  o     <- strsplit(only, ":")[[1]]
  pairs <- pairs[pairs$species == o[1] & pairs$bcr == o[2], ]
}
pairs <- pairs[!file.exists(cache_of(pairs$species, pairs$bcr)), ]
message(nrow(pairs), " species x BCR to do")

for (bcr in intersect(bcrs, pairs$bcr)) {
  t_bcr <- Sys.time()
  stack_tmp <- file.path(tempdir(), sprintf("%s_%d.tif", bcr, year))
  file.copy(file.path(nm_root, "gis", "stacks", basename(stack_tmp)), stack_tmp, overwrite = TRUE)
  stk  <- rast(stack_tmp)
  tmpl <- stk[[1]]

  # as 12F: each sector's rasters by nearest neighbour onto the stack grid
  direct_any <- logical(ncell(tmpl))
  for (sec in sectors) {
    s <- values(project(rast(file.path(hp, paste0(sec, ".tif"))), tmpl, method = "near"), mat = FALSE)
    d <- values(project(rast(file.path(hp, paste0(sec, "_direct.tif"))), tmpl, method = "near"), mat = FALSE)
    direct_any <- direct_any | (!is.na(d) & d == 1 & !is.na(s) & s > 0)
  }
  sid_bcr <- unique(ss$subbasin[ss$bcr == bcr])
  zone <- values(rasterize(subs[sid_bcr, ], tmpl, field = "sid"), mat = FALSE)
  inbc <- values(rasterize(bc, tmpl, background = 0), mat = FALSE) == 1

  for (sp in pairs$species[pairs$bcr == bcr]) {
    t0 <- Sys.time()
    px <- rast(file.path(ia_dir, "data", "derived_data", "density_tables", "pixel_change",
                         sprintf("%s_%d_%s.tif", sp, year, bcr)))
    stopifnot(compareGeom(px, tmpl))
    pv    <- values(px)
    cells <- which(!is.na(pv[, "obs"]))
    z     <- zone[cells]
    cells <- cells[z %in% unique(z[inbc[cells]])]        # kept pixels of subbasins touching BC
    if (length(cells) == 0L) { message(sp, " ", bcr, ": no kept pixels in BC"); next }

    w_r <- rast(file.path(ia_dir, "data", "derived_data", "predictions", sp, bcr, year, "weight.tif"))
    stopifnot(names(w_r)[1] == "weight_v3_touches", compareGeom(w_r, tmpl))
    w   <- values(w_r, mat = FALSE)[cells]
    stopifnot(all(is.finite(w)))

    tp   <- readRDS(file.path(ia_dir, "data", "derived_data", "predictions", sp,
                              "truncation_params.rds"))[[paste(bcr, year, sep = "_")]]
    qsp  <- q.out$densmax[q.out$spp == sp]
    stopifnot(length(qsp) == 1L, is.finite(tp$q99), isTRUE(all.equal(tp$densmax, qsp)))

    e <- new.env(parent = emptyenv())
    load(file.path(nm_root, "output", "06_bootstraps", sp, paste0(sp, "_", bcr, ".Rdata")), envir = e)
    b <- e$b.list; rm(e); invisible(gc(verbose = FALSE))
    vars   <- b[[1]]$var.names
    stopifnot(all(vapply(b, function(m) identical(m$var.names, vars), logical(1L))))
    dist_b <- intersect(dist, vars)
    cat_b  <- intersect(vars, categorical_responses)

    # 12F's order exactly: values() (not extract(), which returns labels for factor layers),
    # NaN -> NA_real_ so gbm takes the missing branch, then the footprint covariates to 0
    X <- values(stk[[vars]])[cells, , drop = FALSE]
    X[is.nan(X)] <- NA_real_
    X <- as.data.frame(X, check.names = FALSE)
    for (v in dist_b) X[[v]] <- 0

    rows <- which(w != 0)                                # weight 0 stays 0, as 12F
    d0   <- matrix(0, length(cells), length(b))
    for (i in seq_along(b)) {
      Xi <- as_model_factors(X[rows, , drop = FALSE], b[[i]], cat_b)
      D  <- gbm_design(b[[i]], Xi)
      if (i == 1L) gbm_check_fast(b[[i]], Xi, D, label = paste(sp, bcr))
      d0[rows, i] <- pmin(pmin(gbm_predict_design(b[[i]], D), qsp), tp$q99)
    }
    d0 <- d0 * w * 100                                   # capped, weighted, birds per pixel
    rm(b, X); invisible(gc(verbose = FALSE))

    fr <- data.frame(cell = cells, x = xFromCell(tmpl, cells), y = yFromCell(tmpl, cells),
                     sid = zone[cells], in_bc = inbc[cells], direct = direct_any[cells],
                     w = w, obs = pv[cells, "obs"], cf = pv[cells, "cf"], d0 = rowMeans(d0))

    # check: off the direct footprint, cf is 12F's d0
    nd  <- !fr$direct
    gap <- abs(fr$d0[nd] - fr$cf[nd])
    chk <- c(n_nondirect = sum(nd), max_abs = if (any(nd)) max(gap) else NA_real_,
             sum_d0 = sum(fr$d0[nd]), sum_cf = sum(fr$cf[nd]),
             frac_px_off = if (any(nd)) mean(gap > 1e-4 * pmax(fr$cf[nd], 1e-3)) else NA_real_)

    keep <- fr$direct                                    # per-bootstrap d0 where bf exists
    saveRDS(list(species = sp, bcr = bcr, frame = fr, check = chk, cells_direct = fr$cell[keep],
                 d0_boot = d0[keep, , drop = FALSE]), cache_of(sp, bcr))
    message(format(Sys.time(), "%H:%M"), " ", sp, " ", bcr, ": ", nrow(fr), " px (", sum(fr$direct),
            " direct, ", sum(fr$in_bc & fr$direct), " direct in BC) | nondirect d0 ",
            round(chk["sum_d0"]), " vs cf ", round(chk["sum_cf"]), ", max gap ",
            signif(chk["max_abs"], 3), ", px off ", signif(chk["frac_px_off"], 3), " | ",
            format(round(Sys.time() - t0, 1)))
    rm(d0, fr, pv); invisible(gc(verbose = FALSE))
  }
  unlink(stack_tmp)
  message(bcr, " done in ", format(round(Sys.time() - t_bcr, 1)))
}
message("D0_PIXELS_DONE")
