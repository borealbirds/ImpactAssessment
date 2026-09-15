# Decisive test: at pixels where V5 source StandardGreenup is NA, did V5 still
# produce an observed prediction? If yes -> gbm tolerates the NA and 12C's
# complete.cases gate is over-dropping (option: relax gate). If no -> V5 also
# masked them (option: accept the gap).
suppressMessages({library(terra)})
terra::terraOptions(progress = 0)

st  <- "G:/Shared drives/BAM_NationalModels5/gis/stacks"
pred<- "C:/Users/mannf/Drive/boreal_avian_modelling_project/ImpactAssessment/data/derived_data/predictions"

for (bcr in c("can80","can60","can11")) {
  src <- terra::rast(file.path(st, paste0(bcr, "_2020.tif")))
  if (!("StandardGreenup_1km" %in% names(src))) { cat(bcr,"no greenup\n"); next }
  gnu <- src[["StandardGreenup_1km"]]
  obs <- terra::rast(file.path(pred, "CAWA", bcr, "2020", "observed_bootstraps.tif"))[[1]]

  # align obs to source grid
  if (!terra::compareGeom(gnu, obs, stopOnError=FALSE, messages=FALSE)) {
    obs <- terra::resample(obs, gnu, method="near")
  }
  na_gnu <- is.na(gnu)
  has_obs <- !is.na(obs)

  n_na_gnu        <- terra::global(na_gnu, "sum", na.rm=TRUE)[1,1]
  n_na_gnu_w_obs  <- terra::global(na_gnu & has_obs, "sum", na.rm=TRUE)[1,1]
  n_obs_tot       <- terra::global(has_obs, "sum", na.rm=TRUE)[1,1]
  n_obs_gnu_na    <- n_na_gnu_w_obs
  cat(sprintf("\n%s: greenup-NA pixels = %d; of those WITH an observed prediction = %d (%.1f%%)\n",
              bcr, n_na_gnu, n_na_gnu_w_obs,
              if (n_na_gnu>0) 100*n_na_gnu_w_obs/n_na_gnu else NA))
  cat(sprintf("   total observed-predicted pixels = %d; share of them that have greenup NA = %.1f%%\n",
              n_obs_tot, if (n_obs_tot>0) 100*n_obs_gnu_na/n_obs_tot else NA))
}
