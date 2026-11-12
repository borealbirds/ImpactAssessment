# ---
# title: VATH density against the footprint covariates, vegetation held at observed values
# author: Mannfred Boehm
# ---
#
# Why VATH gains birds from industry (TODO.md, 2026-11-03): its +332k is the direct part of
# 14B's split (+409k: the bird model's response to the footprint covariates, at observed
# vegetation), not the backfill (vegetation part -77k). This figure shows that response, one
# row per BCR that carries it:
#   left   partial dependence of the BCR's leading footprint covariate (can10 CanHF_1km,
#          can5 canroad_5x5), averaged over the BCR's footprint pixels rather than the training
#          surveys, since the direct part is evaluated there. The other footprint covariates
#          stay at their observed values. Rugs: survey deciles (the models' var.levels, median
#          over bootstraps) and footprint pixel deciles; the band is the lowest survey decile.
#   right  every footprint covariate scaled together from 0 (12F's d0) to its observed value
#          (obs). The gap between the two ends is the direct part, compared with 14B's.
# Footprint pixels as 12F's superset: any of the 8 sectors > 0 and CanHF >= 1, inside the BCR.
# Approximations against 12F, so the direct part matches 14B only roughly: weight = V5 range
# only (no water or data-limit mask), no q99 cap (densmax only), a random sample of n_px
# footprint pixels, and the national covariate mosaic instead of the BCR stack.
#
# Runs locally from the RProject root. Reads each BCR's bootstrap models from G: by exact file
# name (~350 MB for can10). About 18 min with n_boot = 16; the curves are cached, so REDO=1
# recomputes them and otherwise the cache is replotted.

library(terra)
library(gbm)
library(dplyr)
library(ggplot2)
library(patchwork)

cc    <- FALSE
local <- TRUE
if (!cc && local) { ia_dir <- getwd() }
nm_root <- "G:/Shared drives/BAM_NationalModels5"

species <- "VATH"
year    <- 2020
focal   <- c(can10 = "CanHF_1km", can5 = "canroad_5x5")
bcr_lab <- c(can10 = "BCR 10, Northern Rockies", can5 = "BCR 5, Northern Pacific Rainforest")
x_lab   <- c(CanHF_1km   = "Local human footprint score, CanHF_1km",
             canroad_5x5 = "Landscape road density, canroad_5x5")
n_boot  <- 16L      # the first n_boot of each BCR's 32 bootstrap models
n_px    <- 20000L
n_grid  <- 21L
s_grid  <- seq(0, 1, by = 0.1)
sectors <- c("crop", "pasture", "built", "roads", "rail", "oil_gas", "mines",
             "dam_and_associated_reservoir")
categorical_responses <- c("ABoVE_1km", "NLCD_1km", "MODISLCC_1km", "MODISLCC_5x5",
                           "SCANFI_1km", "VLCE_1km")

out_dir   <- file.path(ia_dir, "output_figures", "direct_part")
cache_rds <- file.path(out_dir, paste0(species, "_footprint_pdp_curves.rds"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# curves ------------------------------------------------------------------------------
if (!file.exists(cache_rds) || Sys.getenv("REDO") == "1") {
  source(file.path(ia_dir, "Rscripts", "12F_predict_species_all_coalitions.R"))  # as_model_factors()
  source(file.path(ia_dir, "Rscripts", "12G_gbm_tree_walk.R"))
  Rcpp::sourceCpp(file.path(ia_dir, "Rscripts", "12G_gbm_tree_walk.cpp"))

  pm   <- tibble::tibble(BAMexploreR::predictor_metadata)
  dist <- pm$predictor[pm$version == "v5" & pm$predictor_class == "Disturbance"]
  load(file.path(ia_dir, "data", "raw_data", "SpeciesPredictionTruncationValues.Rdata"))
  densmax <- q.out$densmax[q.out$spp == species]

  mos     <- rast(file.path(ia_dir, "data", "raw_data", "covariates_mosaiced",
                            sprintf("covariates_mosaiced_%d.tif", year)))
  bcrs    <- vect(file.path(ia_dir, "data", "raw_data", "Regions", "BAM_BCR_NationalModel_Unbuffered.shp"))
  hp      <- file.path(ia_dir, "data", "raw_data", "hirshpearson")
  range_r <- rast(file.path(ia_dir, "data", "raw_data", "v5_gis", "ranges", paste0(species, ".tif")))
  sb      <- read.csv(file.path(ia_dir, "data", "derived_data", "sector_effects", "shapley_bcr.csv"))

  res <- lapply(names(focal), function(bcr) {
    t0 <- Sys.time()
    fv    <- focal[[bcr]]
    poly  <- bcrs[bcrs$country == "can" & bcrs$subUnit == as.integer(sub("can", "", bcr))]
    tmpl  <- crop(mos[[1]], poly, snap = "out")
    inbcr <- values(rasterize(poly, tmpl, touches = TRUE, background = 0), mat = FALSE) == 1
    hf    <- values(project(rast(file.path(hp, "CanHF_1km_morethan1.tif")), tmpl, method = "near"), mat = FALSE)
    fp    <- Reduce(`|`, lapply(sectors, function(s) {
      v <- values(project(rast(file.path(hp, paste0(s, ".tif"))), tmpl, method = "near"), mat = FALSE)
      !is.na(v) & v > 0
    }))
    w_all <- values(classify(project(range_r, tmpl, method = "bilinear"), cbind(NA, 0)), mat = FALSE)
    cells <- which(fp & !is.na(hf) & hf >= 1 & inbcr & w_all > 0)
    n_fp  <- length(cells)
    set.seed(1)
    cells <- if (n_fp > n_px) sort(sample(cells, n_px)) else cells
    w     <- w_all[cells]

    e <- new.env()
    load(file.path(nm_root, "output", "06_bootstraps", species, paste0(species, "_", bcr, ".Rdata")), envir = e)
    b <- e$b.list[seq_len(min(n_boot, length(e$b.list)))]
    rm(e)
    vars <- unique(unlist(lapply(b, `[[`, "var.names")))
    X <- terra::extract(mos[[vars]], cellFromXY(mos, xyFromCell(tmpl, cells)))
    for (j in seq_along(X)) if (is.double(X[[j]])) X[[j]][is.nan(X[[j]])] <- NA_real_
    dist_b <- intersect(dist, vars)
    cat_b  <- intersect(vars, categorical_responses)

    fp_q   <- quantile(X[[fv]], seq(0.1, 0.9, 0.1), na.rm = TRUE)
    x_grid <- seq(0, quantile(X[[fv]], 0.95, na.rm = TRUE), length.out = n_grid)
    sv_q   <- apply(sapply(b, function(m) m$var.levels[[match(fv, m$var.names)]]), 1, median)

    pdp <- scl <- vector("list", length(b))
    for (i in seq_along(b)) {
      m  <- b[[i]]
      Xf <- as_model_factors(X, m, cat_b)
      D  <- gbm_design(m, Xf)
      if (i == 1L) gbm_check_fast(m, Xf, D, label = bcr)
      dens <- function(D) mean(pmin(gbm_predict_design(m, D), densmax) * w * 100)
      pdp[[i]] <- data.frame(boot = i, x = x_grid, dens = vapply(x_grid, function(g) {
        D[, fv] <- g
        dens(D)
      }, numeric(1L)))
      scl[[i]] <- data.frame(boot = i, x = s_grid, dens = vapply(s_grid, function(s) {
        D[, dist_b] <- if (s == 0) 0 else D[, dist_b] * s   # s = 0 is d0, NA included, as 12F
        dens(D)
      }, numeric(1L)))
      message(bcr, " bootstrap ", i, "/", length(b), " ", format(round(Sys.time() - t0, 1)))
    }
    scl <- bind_rows(scl)
    direct <- with(scl, (dens[x == 1] - dens[x == 0]) * n_fp)   # birds, one per bootstrap
    list(bcr = bcr, focal = fv, n_fp = n_fp, n_px = length(cells),
         pdp = bind_rows(pdp), scl = scl, survey_q = sv_q, fp_q = fp_q, direct = direct,
         direct_14B = -sum(sb$shapley_direct_mean[sb$species == species & sb$bcr == bcr & sb$year == year]))
  })
  names(res) <- names(focal)
  saveRDS(res, cache_rds)
  write.csv(bind_rows(lapply(res, function(r) bind_rows(
    mutate(r$pdp, bcr = r$bcr, panel = paste0("partial dependence: ", r$focal)),
    mutate(r$scl, bcr = r$bcr, panel = "all footprint covariates x s")))),
    file.path(out_dir, paste0(species, "_footprint_pdp_curves.csv")), row.names = FALSE)
}
res <- readRDS(cache_rds)

# figure ------------------------------------------------------------------------------
blue   <- "#2a78d6"
orange <- "#eb6834"
ink    <- "#3d3d3a"
fmt_k  <- function(x) paste0(ifelse(x >= 0, "+", "\u2212"), formatC(abs(x) / 1000, format = "f", digits = 0), "k")
theme_set(theme_minimal(base_size = 11) +
            theme(panel.grid.minor = element_blank(),
                  panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
                  plot.title = element_text(size = 11, face = "bold", colour = ink),
                  axis.title = element_text(colour = ink), axis.text = element_text(colour = "grey35")))

rows <- lapply(res, function(r) {
  mean_pdp <- r$pdp |> group_by(x) |> summarise(dens = mean(dens), .groups = "drop")
  mean_scl <- r$scl |> group_by(x) |> summarise(dens = mean(dens), .groups = "drop")
  ylim  <- range(c(0, r$pdp$dens, r$scl$dens))
  ylim[2] <- ylim[2] * 1.12
  x_max <- max(r$pdp$x)
  q10   <- r$survey_q[2]          # var.levels: 0, 10, ..., 100% quantiles
  sv_q  <- r$survey_q[2:10][r$survey_q[2:10] <= x_max]
  fp_q  <- r$fp_q[r$fp_q <= x_max]

  left <- ggplot() +
    annotate("rect", xmin = 0, xmax = min(q10, x_max), ymin = -Inf, ymax = Inf, fill = orange, alpha = 0.10) +
    annotate("text", x = min(q10, x_max) + 0.012 * x_max, y = 0.3 * ylim[2], label = "shaded: below the\nlowest survey decile",
             size = 3, colour = ink, hjust = 0, vjust = 0.5, lineheight = 0.9) +
    geom_line(data = r$pdp, aes(x, dens, group = boot), colour = "grey60", alpha = 0.35, linewidth = 0.35) +
    geom_line(data = mean_pdp, aes(x, dens), colour = blue, linewidth = 1) +
    geom_rug(data = data.frame(x = sv_q), aes(x = x), sides = "b", colour = orange,
             linewidth = 0.8, length = unit(0.035, "npc")) +
    geom_rug(data = data.frame(x = fp_q), aes(x = x), sides = "t", colour = blue,
             linewidth = 0.8, length = unit(0.035, "npc")) +
    annotate("text", x = x_max, y = ylim[1], label = "ticks: survey deciles", hjust = 1, vjust = -1.6,
             size = 2.8, colour = ink) +
    annotate("text", x = x_max, y = ylim[2], label = "ticks: footprint pixel deciles", hjust = 1, vjust = 2.2,
             size = 2.8, colour = ink) +
    coord_cartesian(xlim = c(0, x_max), ylim = ylim, expand = FALSE) +
    labs(title = paste0(bcr_lab[[r$bcr]], ": partial dependence"),
         x = x_lab[[r$focal]], y = "VATH on footprint pixels (birds/km\u00b2)")

  d_mean <- mean(r$direct)
  d_q    <- quantile(r$direct, c(0.05, 0.95))
  ends   <- mean_scl |> filter(x %in% c(0, 1)) |>
    mutate(lab = ifelse(x == 0, "footprint\ncovariates at 0", "observed"),
           hj = ifelse(x == 0, 0, 1), xl = ifelse(x == 0, 0.015, 0.985))
  right <- ggplot() +
    geom_line(data = r$scl, aes(x, dens, group = boot), colour = "grey60", alpha = 0.35, linewidth = 0.35) +
    geom_line(data = mean_scl, aes(x, dens), colour = blue, linewidth = 1) +
    geom_point(data = ends, aes(x, dens), colour = blue, fill = "white", shape = 21, size = 2.6, stroke = 1) +
    geom_text(data = ends, aes(xl, dens, label = lab, hjust = hj), vjust = 1.4, size = 2.9,
              colour = ink, lineheight = 0.9) +
    annotate("text", x = 0.02, y = ylim[2], hjust = 0, vjust = 1, size = 3, colour = ink, lineheight = 1,
             label = sprintf("observed \u2212 at 0: %s birds on %s footprint km\u00b2 (5\u201395%%: %s to %s)\n14B's direct part for %s: %s birds",
                             fmt_k(d_mean), format(r$n_fp, big.mark = ","), fmt_k(d_q[1]), fmt_k(d_q[2]),
                             r$bcr, fmt_k(r$direct_14B))) +
    scale_x_continuous(breaks = seq(0, 1, 0.25), labels = c("0", "0.25", "0.5", "0.75", "observed")) +
    coord_cartesian(xlim = c(-0.02, 1.02), ylim = ylim, expand = FALSE) +
    labs(title = paste0(bcr_lab[[r$bcr]], ": all footprint covariates together"),
         x = "Footprint covariates as a fraction of their observed values", y = NULL)
  left | right
})

fig <- wrap_plots(rows, ncol = 1) +
  plot_annotation(
    title = "Varied Thrush: modelled density against the footprint covariates, vegetation held at observed values",
    subtitle = paste0("Each grey line is one of ", max(res[[1]]$pdp$boot), " of the BCR's 32 bootstrap models; blue is their mean."),
    caption = paste0(
      "Averaged over ", format(res[[1]]$n_px, big.mark = ","), " random footprint pixels per BCR (any of the 8 sectors and CanHF \u2265 1), ",
      "each with its observed vegetation and climate. Left: the other footprint covariates stay at their observed values.\n",
      "Footprint covariates: CanHF_1km, CanHF_5x5, canroad_1km, canroad_5x5, CCNL_1km. CanHF_1km is Hirsh-Pearson's 300 m score at the cell centre.\n",
      "Weighted by the V5 range only (no water or data-limit mask) and capped at densmax only (no q99), so it differs a little from 14B."),
    theme = theme(plot.title = element_text(face = "bold", size = 13, colour = ink),
                  plot.subtitle = element_text(colour = "grey35"),
                  plot.caption = element_text(hjust = 0, size = 8.5, colour = "grey35")))
ggsave(file.path(out_dir, paste0(species, "_footprint_pdp.png")), fig, width = 13, height = 9.5, dpi = 200, bg = "white")
message("wrote ", file.path(out_dir, paste0(species, "_footprint_pdp.png")))
