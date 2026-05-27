suppressPackageStartupMessages(library(terra))
bf_dir  <- file.path("data","derived_data","predictions_coalitions","9","CAWA")
obs_dir <- file.path("data","derived_data","predictions","CAWA")
bcrs <- list.files(bf_dir)
for (bcr in bcrs) {
  bf  <- rast(file.path(bf_dir,  bcr, "2020", "backfilled_mean.tif"))
  obs <- rast(file.path(obs_dir, bcr, "2020", "observed_mean.tif"))
  imp <- bf - obs
  v   <- values(imp, mat = FALSE)
  nz  <- v[is.finite(v) & abs(v) > .Machine$double.eps]
  ov  <- values(obs, mat = FALSE)
  pos <- sum(nz > 0); neg <- sum(nz < 0)
  cat(sprintf("%s  CRS=%-8s  n_nz=%6d  pos=%5d  neg=%5d  max+=%.4f  max-=%.4f  obs_max=%.4f\n",
              bcr, sub(".*EPSG:", "EPSG:", crs(imp, describe = TRUE)$code),
              length(nz), pos, neg,
              if (pos) max(nz[nz > 0]) else NA_real_,
              if (neg) min(nz[nz < 0]) else NA_real_,
              max(ov, na.rm = TRUE)))
}
