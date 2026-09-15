# Read-only: measure per-band NA fraction in the FOUNDATION PROJECT source stacks
# (G:/.../gis/stacks/{bcr}_{year}.tif) to see whether FFP/MCMT etc. are sparse
# in the upstream source, not just in our mosaic. Run locally with G: mounted.

suppressPackageStartupMessages(library(terra))
terraOptions(progress = 0)

stacks_dir <- "G:/Shared drives/BAM_NationalModels5/gis/stacks"
out_csv    <- "C:/Users/mannf/Drive/boreal_avian_modelling_project/ImpactAssessment/logs/source_stack_climate.csv"

# pick a few BCRs spanning bad->ok coverage from rootcause_bcr_rollup.csv
bcrs <- c("can80", "can10", "can14", "can61", "can60")
year <- 2020

climate_bands <- c("FFP_1km","MCMT_1km","SHM_1km","MSP_1km","EXT_1km","CMD_1km",
                   "bFFP_1km","eFFP_1km","PPTwt_1km","TD_1km","NFFD_1km","EMT_1km",
                   "AHM_1km","RH_1km")

res <- list()
for (b in bcrs) {
  f <- file.path(stacks_dir, paste0(b, "_", year, ".tif"))
  if (!file.exists(f)) { cat("MISSING", f, "\n"); next }
  r <- rast(f)
  nm <- names(r)
  ncell_total <- ncell(r)
  nn <- terra::global(r, fun = "notNA")
  df <- data.frame(bcr = b, band = rownames(nn),
                   notNA = nn[[1]],
                   na_frac = 1 - nn[[1]]/ncell_total,
                   row.names = NULL)
  # reference: a "full" land band's notNA (use the max as the land-mask baseline)
  df$frac_of_max <- df$notNA / max(df$notNA)
  res[[b]] <- df
  cat("\n===", b, "(", nlyr(r), "bands,", ncell_total, "cells ) ===\n")
  sub <- df[df$band %in% climate_bands, ]
  sub <- sub[order(-sub$na_frac), ]
  print(sub, digits = 4, row.names = FALSE)
  cat("max notNA (land-mask baseline):", max(df$notNA),
      " min notNA:", min(df$notNA), "\n")
}

all <- do.call(rbind, res)
write.csv(all, out_csv, row.names = FALSE)
cat("\nwrote", out_csv, "\n")
