# Read-only check: is the climate-normal NA in covariates_mosaiced a SPARSE SOURCE BAND
# (NA over the whole grid) or a SPATIAL MASK coinciding with high-HF pixels?
# Run locally on the project machine. Writes logs/climate_source_sparsity.csv.

suppressPackageStartupMessages(library(terra))
terraOptions(progress = 0)

ia_dir   <- "C:/Users/mannf/Drive/boreal_avian_modelling_project/ImpactAssessment"
cov_path <- file.path(ia_dir, "data/raw_data/covariates_mosaiced/covariates_mosaiced_2020.tif")
hf_path  <- file.path(ia_dir, "data/raw_data/hirshpearson/CanHF_1km_morethan1.tif")
out_csv  <- file.path(ia_dir, "logs/climate_source_sparsity.csv")

r  <- rast(cov_path)
nm <- names(r)
ncell_total <- ncell(r)

cat("covariates_mosaiced_2020.tif:", nlyr(r), "bands,", ncell_total, "cells\n")
cat("CRS:", crs(r, describe = TRUE)$code, "\n")
cat("bands:\n"); print(nm)

# whole-grid non-NA count per band (single pass)
nn <- terra::global(r, fun = "notNA")
df <- data.frame(
  band            = rownames(nn),
  notNA_wholegrid = nn[[1]],
  na_frac_wholegrid = 1 - nn[[1]] / ncell_total,
  row.names = NULL
)

# high-HF subset: align morethan1 mask to covariate grid, measure NA there
hf_na <- rep(NA_real_, nrow(df)); hf_n <- NA_integer_
if (file.exists(hf_path)) {
  hf <- rast(hf_path)
  if (!compareGeom(hf, r[[1]], stopOnError = FALSE)) {
    hf <- project(hf, r[[1]], method = "near")
  }
  hf_idx <- which(values(hf, mat = FALSE) == 1)
  hf_n <- length(hf_idx)
  cat("high-HF pixels (morethan1==1) on covariate grid:", hf_n, "\n")
  if (hf_n > 0) {
    vv <- r[hf_idx]                       # data.frame: hf_n rows x nlyr cols
    hf_na <- 1 - colSums(!is.na(vv)) / hf_n
    hf_na <- as.numeric(hf_na[df$band])
  }
}
df$na_frac_highHF <- hf_na
df$n_highHF       <- hf_n

df <- df[order(-df$na_frac_wholegrid), ]
write.csv(df, out_csv, row.names = FALSE)
cat("\n=== per-band NA fraction (sorted by whole-grid NA) ===\n")
print(df, digits = 4, row.names = FALSE)
cat("\nwrote", out_csv, "\n")
