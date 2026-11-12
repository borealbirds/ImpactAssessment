# ---
# title: Western forest birds' density as the footprint covariates go from 0 to observed
# author: Mannfred Boehm
# ---
#
# The direct part of 14B's split (Open Limitation #9) is the bird model's response to the
# footprint covariates at observed vegetation: obs - d0, where d0 has every footprint
# ("Disturbance") covariate at 0. This figure draws that response as a curve: every footprint
# covariate scaled together, s x its observed value for s = 0 (d0), 0.1, ..., 1 (obs). One figure
# per BCR, one panel per species. Follows VATH's industry-positive result (TODO.md, 2026-11-03;
# figure_vath_footprint_pdp.R).
#
# Each BCR uses one random sample of n_px footprint pixels (12F's superset: any of the 8 sectors
# > 0 and CanHF >= 1, inside the BCR) for every species, so the y axis is birds per footprint
# km2 of the BCR, and (curve at observed - curve at 0) x the BCR's footprint km2 is the direct
# part. Approximations against 12F, so it matches 14B only roughly: weight = V5 range only (no
# water or data-limit mask), no q99 cap (densmax only), and the national covariate mosaic
# instead of the BCR stack.
#
# Runs locally from the RProject root. Reads each species x BCR's bootstrap models from G: by
# exact file name (150-370 MB each). Each species x BCR's curves are cached in
# output_figures/direct_part/curves/, so a rerun computes only what is missing; REDO=1
# recomputes everything.

library(terra)
library(gbm)
library(dplyr)
library(ggplot2)

cc    <- FALSE
local <- TRUE
if (!cc && local) { ia_dir <- getwd() }
nm_root <- "G:/Shared drives/BAM_NationalModels5"

species <- c(RBNU = "Red-breasted Nuthatch", SWTH = "Swainson's Thrush",
             RCKI = "Ruby-crowned Kinglet",  OCWA = "Orange-crowned Warbler",
             TOWA = "Townsend's Warbler",    VATH = "Varied Thrush")
bcr_lab <- c(can10 = "BCR 10, Northern Rockies", can5 = "BCR 5, Northern Pacific Rainforest",
             can9  = "BCR 9, Great Basin")
year    <- 2020
n_boot  <- 8L       # the first n_boot of each species x BCR's 32 bootstrap models
n_px    <- 5000L
s_grid  <- seq(0, 1, by = 0.1)
sectors <- c("crop", "pasture", "built", "roads", "rail", "oil_gas", "mines",
             "dam_and_associated_reservoir")
categorical_responses <- c("ABoVE_1km", "NLCD_1km", "MODISLCC_1km", "MODISLCC_5x5",
                           "SCANFI_1km", "VLCE_1km")

out_dir   <- file.path(ia_dir, "output_figures", "direct_part")
cache_dir <- file.path(out_dir, "curves")
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
cache_of  <- function(sp, bcr) file.path(cache_dir, sprintf("%s_%s_scaling.rds", sp, bcr))
redo      <- Sys.getenv("REDO") == "1"

# curves ------------------------------------------------------------------------------
todo <- expand.grid(sp = names(species), bcr = names(bcr_lab), stringsAsFactors = FALSE)
todo <- todo[redo | !file.exists(cache_of(todo$sp, todo$bcr)), ]

if (nrow(todo) > 0) {
  source(file.path(ia_dir, "Rscripts", "12F_predict_species_all_coalitions.R"))  # as_model_factors()
  source(file.path(ia_dir, "Rscripts", "12G_gbm_tree_walk.R"))
  Rcpp::sourceCpp(file.path(ia_dir, "Rscripts", "12G_gbm_tree_walk.cpp"))

  pm   <- tibble::tibble(BAMexploreR::predictor_metadata)
  dist <- pm$predictor[pm$version == "v5" & pm$predictor_class == "Disturbance"]
  load(file.path(ia_dir, "data", "raw_data", "SpeciesPredictionTruncationValues.Rdata"))

  mos  <- rast(file.path(ia_dir, "data", "raw_data", "covariates_mosaiced",
                         sprintf("covariates_mosaiced_%d.tif", year)))
  bcrs <- vect(file.path(ia_dir, "data", "raw_data", "Regions", "BAM_BCR_NationalModel_Unbuffered.shp"))
  hp   <- file.path(ia_dir, "data", "raw_data", "hirshpearson")
  sb   <- read.csv(file.path(ia_dir, "data", "derived_data", "sector_effects", "shapley_bcr.csv"))

  for (bcr in unique(todo$bcr)) {
    t0    <- Sys.time()
    poly  <- bcrs[bcrs$country == "can" & bcrs$subUnit == as.integer(sub("can", "", bcr))]
    tmpl  <- crop(mos[[1]], poly, snap = "out")
    inbcr <- values(rasterize(poly, tmpl, touches = TRUE, background = 0), mat = FALSE) == 1
    hf    <- values(project(rast(file.path(hp, "CanHF_1km_morethan1.tif")), tmpl, method = "near"), mat = FALSE)
    fp    <- Reduce(`|`, lapply(sectors, function(s) {
      v <- values(project(rast(file.path(hp, paste0(s, ".tif"))), tmpl, method = "near"), mat = FALSE)
      !is.na(v) & v > 0
    }))
    cells <- which(fp & !is.na(hf) & hf >= 1 & inbcr)
    n_fp  <- length(cells)
    set.seed(1)
    cells <- if (n_fp > n_px) sort(sample(cells, n_px)) else cells
    X_all <- terra::extract(mos, cellFromXY(mos, xyFromCell(tmpl, cells)))
    for (j in seq_along(X_all)) if (is.double(X_all[[j]])) X_all[[j]][is.nan(X_all[[j]])] <- NA_real_
    message(bcr, ": ", format(n_fp, big.mark = ","), " footprint pixels, ", length(cells), " sampled")

    for (sp in todo$sp[todo$bcr == bcr]) {
      range_r <- rast(file.path(ia_dir, "data", "raw_data", "v5_gis", "ranges", paste0(sp, ".tif")))
      w       <- values(classify(project(range_r, tmpl, method = "bilinear"), cbind(NA, 0)), mat = FALSE)[cells]
      densmax <- q.out$densmax[q.out$spp == sp]
      stopifnot(length(densmax) == 1L)
      keep    <- w > 0

      e <- new.env()
      load(file.path(nm_root, "output", "06_bootstraps", sp, paste0(sp, "_", bcr, ".Rdata")), envir = e)
      b <- e$b.list[seq_len(min(n_boot, length(e$b.list)))]
      rm(e); gc(verbose = FALSE)
      vars   <- unique(unlist(lapply(b, `[[`, "var.names")))
      X      <- X_all[keep, vars, drop = FALSE]
      dist_b <- intersect(dist, vars)
      cat_b  <- intersect(vars, categorical_responses)

      scl <- vector("list", length(b))
      for (i in seq_along(b)) {
        m <- b[[i]]
        if (!any(keep)) { scl[[i]] <- data.frame(boot = i, x = s_grid, dens = 0); next }
        Xf <- as_model_factors(X, m, cat_b)
        D  <- gbm_design(m, Xf)
        if (i == 1L) gbm_check_fast(m, Xf, D, label = paste(sp, bcr))
        scl[[i]] <- data.frame(boot = i, x = s_grid, dens = vapply(s_grid, function(s) {
          Ds <- D
          Ds[, dist_b] <- if (s == 0) 0 else D[, dist_b] * s   # s = 0 is d0, NA included, as 12F
          sum(pmin(gbm_predict_design(m, Ds), densmax) * w[keep] * 100) / length(cells)
        }, numeric(1L)))
      }
      scl    <- bind_rows(scl)
      direct <- with(scl, (dens[x == 1] - dens[x == 0]) * n_fp)   # birds, one per bootstrap
      saveRDS(list(species = sp, bcr = bcr, n_fp = n_fp, n_px = length(cells), n_in_range = sum(keep),
                   footprint_vars = dist_b, scl = scl, direct = direct,
                   direct_14B = -sum(sb$shapley_direct_mean[sb$species == sp & sb$bcr == bcr & sb$year == year])),
              cache_of(sp, bcr))
      message(sp, " ", bcr, " done ", format(round(Sys.time() - t0, 1)),
              " | direct ", round(mean(direct)), " vs 14B ",
              round(-sum(sb$shapley_direct_mean[sb$species == sp & sb$bcr == bcr & sb$year == year])))
      rm(b, X); gc(verbose = FALSE)
    }
  }
}

g   <- expand.grid(sp = names(species), bcr = names(bcr_lab), stringsAsFactors = FALSE)
res <- Map(function(sp, bcr) readRDS(cache_of(sp, bcr)), g$sp, g$bcr)
curves <- bind_rows(lapply(res, function(r) mutate(r$scl, species = r$species, bcr = r$bcr)))
write.csv(curves, file.path(out_dir, "footprint_scaling_curves.csv"), row.names = FALSE)
write.csv(bind_rows(lapply(res, function(r) data.frame(
  species = r$species, bcr = r$bcr, footprint_km2 = r$n_fp, sampled = r$n_px, sampled_in_range = r$n_in_range,
  direct_mean = round(mean(r$direct)), direct_q05 = round(quantile(r$direct, 0.05)),
  direct_q95 = round(quantile(r$direct, 0.95)), direct_14B = round(r$direct_14B)))),
  file.path(out_dir, "footprint_scaling_direct.csv"), row.names = FALSE)

# figures -----------------------------------------------------------------------------
blue <- "#2a78d6"
ink  <- "#3d3d3a"
theme_set(theme_minimal(base_size = 11) +
            theme(panel.grid.minor = element_blank(),
                  panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
                  strip.text = element_text(size = 11, face = "bold", colour = ink, hjust = 0),
                  axis.title = element_text(colour = ink), axis.text = element_text(colour = "grey35"),
                  panel.spacing = unit(1.2, "lines"),
                  plot.margin = margin(8, 24, 8, 8),
                  plot.title = element_text(face = "bold", size = 13, colour = ink)))

for (bcr in names(bcr_lab)) {
  d <- curves |> filter(bcr == !!bcr) |>
    mutate(name = factor(.env$species[.data$species], levels = .env$species))
  m <- d |> group_by(name, x) |> summarise(dens = mean(dens), .groups = "drop")
  fig <- ggplot() +
    geom_line(data = d, aes(x, dens, group = boot), colour = "grey60", alpha = 0.45, linewidth = 0.35) +
    geom_line(data = m, aes(x, dens), colour = blue, linewidth = 1) +
    geom_point(data = filter(m, x %in% c(0, 1)), aes(x, dens), colour = blue, fill = "white",
               shape = 21, size = 2.4, stroke = 1) +
    facet_wrap(~ name, ncol = 3, scales = "free_y") +
    scale_x_continuous(breaks = seq(0, 1, 0.25), labels = c("0", "0.25", "0.5", "0.75", "observed"),
                       expand = expansion(mult = 0.03)) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.08))) +
    labs(title = bcr_lab[[bcr]],
         x = "Footprint covariates as a fraction of their observed values",
         y = "Birds per footprint km\u00b2")
  f <- file.path(out_dir, sprintf("footprint_scaling_%s.png", bcr))
  ggsave(f, fig, width = 12, height = 7.5, dpi = 200, bg = "white")
  message("wrote ", f)
}
