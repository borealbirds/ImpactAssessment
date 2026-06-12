# ---
# title: Diagnostic — why is CAWA oil_gas negative in can60?
# author: Mannfred Boehm
# ---
# Compares observed vs backfilled biotic covariates at can60 oil_gas footprint
# pixels and joins each covariate's gbm importance from CAWA's can60 bootstrap
# models, so the dominant driver of the negative impact is visible.
#
#   module load StdEnv/2023 r/4.4.0
#   cd /home/mannfred/scratch/impact_assessment/Rscripts && Rscript diag_can60_oil_gas.R
# ---

suppressPackageStartupMessages({
  library(terra)
  library(gbm)
  library(dplyr)
  library(tibble)
})

ia_dir   <- "/home/mannfred/scratch/impact_assessment"
nm_root  <- "/home/mannfred/projects/def-ecknight/NationalModels"
species  <- "CAWA"
year     <- 2020
bcr_code <- "can60"
sector   <- "oil_gas"
hirsh_dir <- file.path(ia_dir, "data", "raw_data", "hirshpearson")

cat("=== loading rasters ===\n")
obs <- terra::rast(file.path(nm_root, "gis", "stacks", paste0(bcr_code, "_", year, ".tif")))
bf  <- terra::rast(file.path(ia_dir, "data", "derived_data", "bart_models_mosaics",
                             year, paste0(bcr_code, "_backfilled.tif")))
og  <- terra::project(terra::rast(file.path(hirsh_dir, paste0(sector, ".tif"))),
                      obs[[1]], method = "near")
chf <- terra::project(terra::rast(file.path(hirsh_dir, "CanHF_1km_morethan1.tif")),
                      obs[[1]], method = "near")

cell_idx <- which(
  terra::values(og,  mat = FALSE) > 0 &
  terra::values(chf, mat = FALSE) >= 1
)
cat("oil_gas footprint pixels in", bcr_code, ":", length(cell_idx), "\n")

cat("=== continuous biotic covariates: obs vs bf (draw 1) ===\n")
draws  <- grep("_draw_001$", names(bf), value = TRUE)
cv     <- sub("_draw_001$", "", draws)
cv_obs <- intersect(cv, names(obs))

cont_diff <- lapply(cv_obs, function(v) {
  o <- terra::values(obs[[v]], mat = FALSE)[cell_idx]
  b <- expm1(terra::values(bf[[paste0(v, "_draw_001")]], mat = FALSE)[cell_idx])
  tibble(
    covariate = v,
    type      = "continuous",
    n_obs     = sum(is.finite(o)),
    n_bf      = sum(is.finite(b)),
    obs_med   = median(o, na.rm = TRUE),
    bf_med    = median(b, na.rm = TRUE),
    obs_mean  = mean(o,   na.rm = TRUE),
    bf_mean   = mean(b,   na.rm = TRUE)
  )
}) |> bind_rows() |>
  mutate(pct_diff = ifelse(obs_med != 0, 100 * (bf_med - obs_med) / obs_med, NA_real_))

cat("--- continuous covs (n =", nrow(cont_diff), ") ---\n")
print(cont_diff |> arrange(pct_diff), n = Inf)

cat("=== categorical biotic covariates: obs class vs bf class (mode/_mean) ===\n")
cat_layers <- setdiff(names(bf), draws)
cat_layers <- cat_layers[!grepl("_(sd|q025|q975)$", cat_layers)]
cat_in_obs <- intersect(sub("_mean$", "", cat_layers), names(obs))

cat_diff <- lapply(cat_in_obs, function(v) {
  o    <- terra::values(obs[[v]], mat = FALSE)[cell_idx]
  blyr <- if (v %in% names(bf)) v else paste0(v, "_mean")
  if (!blyr %in% names(bf)) return(NULL)
  b <- terra::values(bf[[blyr]], mat = FALSE)[cell_idx]
  agree <- sum(o == b, na.rm = TRUE) / sum(is.finite(o) & is.finite(b))
  tibble(covariate = v, type = "categorical",
         n_obs = sum(is.finite(o)), n_bf = sum(is.finite(b)),
         agreement = agree)
}) |> bind_rows()

cat("--- categorical covs (n =", nrow(cat_diff), ") ---\n")
print(cat_diff |> arrange(agreement), n = Inf)

cat("=== pulling CAWA can60 gbm importance ===\n")
rdata_path <- file.path(nm_root, "output", "06_bootstraps", species,
                        paste0(species, "_", bcr_code, ".Rdata"))
e <- new.env(parent = emptyenv())
load(rdata_path, envir = e)
b.list <- e$b.list

imp_long <- lapply(seq_along(b.list), function(i) {
  m  <- b.list[[i]]
  ri <- summary(m, plotit = FALSE, n.trees = m$n.trees)
  tibble(boot = i, covariate = as.character(ri$var), rel_inf = ri$rel.inf)
}) |> bind_rows()

imp_mean <- imp_long |>
  group_by(covariate) |>
  summarise(rel_inf = mean(rel_inf), .groups = "drop") |>
  arrange(desc(rel_inf))

cat("=== top 20 CAWA covariates by mean rel_inf (can60) ===\n")
print(imp_mean |> head(20), n = 20)

cat("\n=== diff joined to importance, sorted by rel_inf ===\n")
diff_all <- bind_rows(
  cont_diff |> select(covariate, type, obs_med, bf_med, pct_diff),
  cat_diff  |> mutate(obs_med = NA_real_, bf_med = NA_real_, pct_diff = NA_real_) |>
    select(covariate, type, obs_med, bf_med, pct_diff, agreement)
)

joined <- imp_mean |>
  left_join(diff_all, by = "covariate") |>
  arrange(desc(rel_inf))

print(joined |> head(30), n = 30)

cat("\n=== continuous covs sorted by absolute bf-obs gap (median scale) ===\n")
print(cont_diff |> arrange(pct_diff), n = 20)

out_csv <- file.path(ia_dir, "logs", "diag_can60_oil_gas.csv")
dir.create(dirname(out_csv), showWarnings = FALSE, recursive = TRUE)
write.csv(joined, out_csv, row.names = FALSE)
cat("\nwrote", out_csv, "\n")
