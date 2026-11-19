# ---
# title: Impact Assessment: set the order in which biotic covariates are backfilled
# author: Mannfred Boehm
# ---

# each continuous covariate's BART model (08A) takes every covariate above it as a predictor,
# so backfilled values cascade down this list. the land-cover classes come last: they are
# backfilled with mbart and never predict the continuous covariates.
# the two phenology layers are V5 "Annual Climate" covariates, so 07 uses them as predictors and
# never backfills them; they are listed because the production file (2026-02-26) has them.
# ("LFheigth" is V5's spelling)

neworder <-
  c("Peatland_1km", "Peatland_5x5",
    "SCANFIclosure_1km", "SCANFIclosure_5x5",
    "SCANFIbiomass_1km", "SCANFIbiomass_5x5",
    "SCANFIheight_1km", "SCANFIheight_5x5",
    "SCANFIheightcv_1km", "SCANFIheightcv_5x5",
    "SCANFIprcC_1km", "SCANFIprcC_5x5",
    "SCANFIprcD_1km", "SCANFIprcD_5x5",
    "SCANFIBalsamFir_1km", "SCANFIBalsamFir_5x5",
    "SCANFIBlackSpruce_1km", "SCANFIBlackSpruce_5x5",
    "SCANFIDouglasFir_1km", "SCANFIDouglasFir_5x5",
    "SCANFIJackPine_1km", "SCANFIJackPine_5x5",
    "SCANFILodgepolePine_1km", "SCANFILodgepolePine_5x5",
    "SCANFIPonderosaPine_1km", "SCANFIPonderosaPine_5x5",
    "SCANFITamarack_1km", "SCANFITamarack_5x5",
    "SCANFIWhiteRedPine_1km", "SCANFIWhiteRedPine_5x5",
    "LFbiomass_1km", "LFbiomass_5x5",
    "ONRLbiomass_1km", "ONRLbiomass_5x5",
    "ETHheight_1km", "ETHheight_5x5",
    "ETHheightcv_1km", "ETHheightcv_5x5",
    "LFheigth_1km", "LFheigth_5x5",
    "LFheigthcv_1km", "LFheigthcv_5x5",
    "RAPTCC_1km", "RAPTCC_5x5",
    "LFcrownclosure_1km", "LFcrownclosure_5x5",
    "StandardDormancy_1km", "StandardGreenup_1km",
    "VLCE_1km", "MODISLCC_1km", "MODISLCC_5x5", "ABoVE_1km", "SCANFI_1km", "NLCD_1km")

saveRDS(neworder, file.path("data", "raw_data", "biotic_variable_hierarchy.rds"))
