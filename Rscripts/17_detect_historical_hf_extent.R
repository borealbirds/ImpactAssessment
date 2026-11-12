# ---
# title: How much of the 2020 footprint already existed in each earlier year
# author: Mannfred Boehm
# created: October 27, 2026
# ---
# 07 ran 1995-2015 on the 2020 footprint (HF_MASK_YEAR=2020), because Hirsh-Pearson exists only
# for 2020. The vegetation part of each year's impact corrects for that by itself: where the
# footprint came later, that year's SCANFI is still natural, so the backfill takes nothing away.
# The direct part does not: V5's CanHF_1km/_5x5 are the 2020 values in every year's stack, so
# setting them to 0 counts an effect on land with no industry yet. This script measures f_y, the
# share of the 2020 footprint already there in year y, and from it how far that moves the
# industry effect 12F-14B report for each year.
#
# f_y is a difference-in-differences on the year-specific covariates in V5's stacks:
#   deficit_y = x_y - (mean x_y over the intact (2020 low-HF) pixels in a win x win km window)
#   f_y       = sum of deficit_y / sum of deficit_2020, over one sector's 2020 direct footprint
# with x = SCANFI canopy height (the treed share in SCANFI's and VLCE's classes as checks; VLCE is
# annual Landsat, mapped independently of SCANFI). A pixel cleared before y carries its
# deficit in both years; one cleared after y has about none in year y, as it was still like its
# surroundings. Taking the surroundings in the same year cancels fire, budworm and SCANFI's own
# year-to-year offsets. Summing before dividing keeps single-pixel noise out of the ratio.
# Roads use canroad_1km (StatCan road files, year-matched from 2001, so 1995 and 2000 read 2001):
# road density barely changes the canopy of a 1 km cell, so its deficit is the land beside it.
#
# Sensitivity: on the share 1 - f_y of the footprint that came after y, the true year-y observed
# density is d0 (no footprint), so the gap there should be 0, not the direct part. With f_y per
# sector x BCR applied to 14B's direct part (cf - obs sign):
#   gap'_y = vegetation_y + f_y * direct_y      obs'_y = obs_y + (1 - f_y) * direct_y
# This assumes the direct effect per footprint pixel is uniform within a sector x BCR. For roads
# in 2005-2015 it is an upper bound, since canroad (part of the direct part) is year-matched there.
#
# Run locally after 14B, with every year's covariates_mosaiced_{year}.tif pulled from Fir:
#   Rscript Rscripts/17_detect_historical_hf_extent.R
# Writes data/derived_data/expanding_HF/footprint_presence_by_year.csv, footprint_sensitivity.csv
# and output_figures/footprint_history/footprint_presence_1995_2020.png.
# ---

library(terra)
suppressPackageStartupMessages(library(ggplot2))
terraOptions(memfrac = 0.5)

ia_dir  <- getwd()   # the local RProject
raw_dir <- file.path(ia_dir, "data", "raw_data")
hp_dir  <- file.path(raw_dir, "hirshpearson")
se_dir  <- file.path(ia_dir, "data", "derived_data", "sector_effects")
out_dir <- file.path(ia_dir, "data", "derived_data", "expanding_HF")
fig_dir <- file.path(ia_dir, "output_figures", "footprint_history")
dir.create(out_dir, showWarnings = FALSE)
dir.create(fig_dir, showWarnings = FALSE)

years      <- seq(1995, 2020, 5)
sectors    <- c("built", "crop", "pasture", "roads", "rail", "dam_and_associated_reservoir",
                "mines", "oil_gas")
indicators <- c(height = "Canopy height (SCANFI)", scanfi_treed = "Treed (SCANFI class)",
                vlce_treed = "Treed (VLCE class)", canroad = "Road density (canroad)")
win        <- 11    # km; the window the intact surroundings are taken from
min_intact <- 10    # intact pixels a window needs
min_pix    <- 50    # a sector x BCR with fewer footprint pixels, or with no canopy deficit in
min_def    <- 1     # 2020 (mean > -1 m), takes its sector's national f

stack_path <- \(y) file.path(raw_dir, "covariates_mosaiced", sprintf("covariates_mosaiced_%d.tif", y))
if (length(miss <- years[!file.exists(stack_path(years))]))
  stop("no covariate mosaic for ", paste(miss, collapse = ", "), " - pull it from Fir")
stacks <- setNames(lapply(years, \(y)
  rast(stack_path(y))[[c("SCANFIheight_1km", "SCANFI_1km", "VLCE_1km", "canroad_1km")]]), years)
grid   <- stacks[[1]][[1]]

# SCANFI height is never NA on mapped land (canopy closure is NA on 24-81% of it, so it cannot
# serve). Treed = SCANFI classes 5-7 (broadleaf, conifer, mixed); VLCE (annual Landsat) 81 treed
# wetland, 210 conifer, 220 broadleaf, 230 mixedwood.
indicator <- function(s, ind) switch(ind,
  height       = s[["SCANFIheight_1km"]],
  scanfi_treed = mask(s[["SCANFI_1km"]] %in% 5:7, s[["SCANFI_1km"]]),
  vlce_treed   = mask(s[["VLCE_1km"]] %in% c(81, 210, 220, 230), s[["VLCE_1km"]]),
  canroad      = s[["canroad_1km"]])

intact <- crop(rast(file.path(hp_dir, "CanHF_1km_lessthan1.tif")), grid)
hf     <- crop(rast(file.path(hp_dir, "CanHF_1km_morethan1.tif")), grid)

bcrs <- vect(file.path(raw_dir, "Regions", "BAM_BCR_NationalModel_Unbuffered.shp"))
bcrs <- bcrs[bcrs$country == "can", ]
bcrs$code <- paste0(bcrs$country, bcrs$subUnit)
bcr_r <- rasterize(project(bcrs, grid), grid, field = "code")


# ----------------------------------------------
# each year's deficit against the intact surroundings

w <- matrix(1, win, win)
deficit <- function(x) {
  xi <- mask(x, intact)
  n  <- focal(!is.na(xi), w, fun = "sum", na.rm = TRUE)
  m  <- focal(xi, w, fun = "sum", na.rm = TRUE) / n
  x - ifel(n >= min_intact, m, NA)
}

def <- list()
for (ind in names(indicators)) {
  f_def <- file.path(out_dir, paste0("deficit_", ind, ".tif"))   # ~30 min to rebuild all four
  if (file.exists(f_def) && all(file.mtime(f_def) > file.mtime(stack_path(years)))) {
    def[[ind]] <- rast(f_def)
    next
  }
  d <- rast(lapply(as.character(years), \(y) {
    message(ind, " ", y)
    x <- indicator(stacks[[y]], ind)
    if (ind == "canroad") x else deficit(x)   # roads: the density itself
  }))
  names(d) <- years
  def[[ind]] <- mask(d, sum(is.na(d)) == 0, maskvalues = 0,   # the same pixels every year
                     filename = f_def, overwrite = TRUE, datatype = "FLT4S")
}


# ----------------------------------------------
# f_y per indicator x sector x BCR: summed deficit in year y over the 2020 direct footprint

zone_sums <- function(x, zone) {
  s <- zonal(x, zone, fun = "sum", na.rm = TRUE)
  n <- zonal(!is.na(x[[1]]), zone, fun = "sum", na.rm = TRUE)
  s <- data.frame(bcr = s[[1]], n_pixels = n[[2]], s[, -1], check.names = FALSE)
  s[s$n_pixels > 0, ]   # a BCR with none of the sector's pixels left
}
long <- function(s, ind, sec) {
  do.call(rbind, lapply(as.character(years), \(y) data.frame(
    indicator = ind, sector = sec, bcr = s$bcr, year = as.integer(y), n_pixels = s$n_pixels,
    mean_value = s[[y]] / s$n_pixels, share_present = s[[y]] / s[[as.character(max(years))]])))
}

rows <- list()
for (sec in c(sectors, "intact")) {
  message("zonal sums: ", sec)
  fp <- if (sec == "intact") intact else
    mask(crop(rast(file.path(hp_dir, paste0(sec, "_direct.tif"))), grid), hf) == 1
  zone <- mask(bcr_r, fp, maskvalues = c(0, NA))
  for (ind in names(def)) {
    if (ind == "canroad" && sec != "roads") next
    s <- zone_sums(def[[ind]], zone)
    nat <- data.frame(bcr = "all", n_pixels = sum(s$n_pixels),
                      t(colSums(s[, -(1:2)], na.rm = TRUE)), check.names = FALSE)
    rows[[paste(sec, ind)]] <- long(rbind(s, nat), ind, sec)
  }
}
pres <- do.call(rbind, rows)
write.csv(pres, file.path(out_dir, "footprint_presence_by_year.csv"), row.names = FALSE)

wide <- function(x, v) reshape(x[, c("indicator", "sector", "year", v)], direction = "wide",
                               idvar = c("indicator", "sector"), timevar = "year")
cat("\nShare of the 2020 direct footprint present in each year (national):\n")
print(wide(pres[pres$bcr == "all" & pres$sector != "intact", ], "share_present"),
      digits = 3, row.names = FALSE)
cat("\nMean deficit of the intact pixels (should be about 0):\n")   # their share is noise / noise
print(wide(pres[pres$bcr == "all" & pres$sector == "intact", ], "mean_value"), digits = 2, row.names = FALSE)


# ----------------------------------------------
# sensitivity: the direct part scaled by f_y, per sector x BCR

f_ind <- setNames(ifelse(sectors == "roads", "canroad", "height"), sectors)
f <- do.call(rbind, lapply(sectors, \(sec) {
  p   <- pres[pres$sector == sec & pres$indicator == f_ind[[sec]], ]
  nat <- p[p$bcr == "all", c("year", "share_present")]
  p   <- p[p$bcr != "all", ]
  d20 <- p$mean_value[match(paste(p$bcr, max(years)), paste(p$bcr, p$year))]
  weak <- p$n_pixels < min_pix | (f_ind[[sec]] == "height" & d20 > -min_def)
  p$f <- ifelse(weak | !is.finite(p$share_present),
                nat$share_present[match(p$year, nat$year)], p$share_present)
  rbind(p[, c("sector", "bcr", "year", "f")],
        data.frame(sector = sec, bcr = "all", year = nat$year, f = nat$share_present))
}))
f$f <- pmin(pmax(f$f, 0), 1)   # a deficit larger than 2020's (regrowth since) counts as all present

sp <- trimws(readLines(file.path(ia_dir, "Rscripts", "species_113.txt")))
sp <- sp[nzchar(sp) & !startsWith(sp, "#")]
sh <- do.call(rbind, lapply(years, \(y) {
  x <- read.csv(file.path(se_dir, paste0("shapley_bcr", if (y == 2020) "" else paste0("_", y), ".csv")))
  x$year <- y
  x[x$species %in% sp, ]
}))
if (length(setdiff(sp, sh$species))) stop("14B's BCR tables lack ", paste(setdiff(sp, sh$species), collapse = ", "))
key  <- paste(f$sector, f$bcr, f$year)
sh$f <- f$f[match(paste(sh$sector, sh$bcr, sh$year), key)]
no_px <- is.na(sh$f)   # a BCR with none of the sector's direct footprint takes its national f
sh$f[no_px] <- f$f[match(paste(sh$sector[no_px], "all", sh$year[no_px]), key)]

sens <- do.call(rbind, lapply(years, \(y) {
  x  <- sh[sh$year == y, ]
  ob <- x[!duplicated(paste(x$species, x$bcr)), ]
  data.frame(year         = y,
             observed     = sum(ob$obs_population),
             gap          = sum(x$shapley_mean),          # 14B's sign: cf - obs
             direct       = sum(x$shapley_direct_mean),
             vegetation   = sum(x$shapley_vegetation_mean),
             gap_adj      = sum(x$shapley_vegetation_mean + x$f * x$shapley_direct_mean),
             observed_adj = sum(ob$obs_population) + sum((1 - x$f) * x$shapley_direct_mean))
}))
write.csv(sens, file.path(out_dir, "footprint_sensitivity.csv"), row.names = FALSE)
cat("\nThe 113 species (14B's sign, cf - obs), M birds:\n")
print(cbind(year = sens$year, round(sens[, -1] / 1e6, 1)), row.names = FALSE)
y0 <- sens$year == min(years); y1 <- sens$year == max(years)
cat(sprintf("\nChange in the gap %d -> %d: %+.1f M as run, %+.1f M with f_y\n", min(years), max(years),
            (sens$gap[y1] - sens$gap[y0]) / 1e6, (sens$gap_adj[y1] - sens$gap_adj[y0]) / 1e6))
cat(sprintf("Observed change %d -> %d: %+.1f M as run, %+.1f M with f_y\n", min(years), max(years),
            (sens$observed[y1] - sens$observed[y0]) / 1e6,
            (sens$observed_adj[y1] - sens$observed_adj[y0]) / 1e6))


# ----------------------------------------------
# figure: national share present, one panel per sector

nat <- pres[pres$bcr == "all" & pres$sector != "intact", ]
nat$sector    <- factor(nat$sector, levels = sectors)
nat$indicator <- factor(nat$indicator, levels = names(indicators), labels = indicators)
p <- ggplot(nat, aes(year, share_present, linetype = indicator, colour = indicator)) +
  geom_hline(yintercept = 1, colour = "grey80", linewidth = 0.3) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.3) +
  facet_wrap(~ sector, nrow = 2) +
  scale_colour_manual(values = c("black", "grey55", "#E69F00", "#0072B2"), name = NULL) +
  scale_linetype_manual(values = c("solid", "dashed", "dashed", "solid"), name = NULL) +
  scale_x_continuous(breaks = years[c(1, 3, 5)]) +
  labs(x = NULL, y = "Share of the 2020 direct footprint present") +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")
ggsave(file.path(fig_dir, "footprint_presence_1995_2020.png"), p, width = 8, height = 4.8, dpi = 300)
