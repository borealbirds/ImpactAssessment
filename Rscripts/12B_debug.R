# ---
# title: Diagnostic wrapper for 12B/12C segfault — single job, one BCR, CHECKPOINT messages
# author: Mannfred Boehm
# ---
# Run as a single SLURM job (not array) via 12B_debug.sh.
# Targets species=BANS, coalition_id=2 (coalition=built) — the exact case from the crash log.
# Prints "CHECKPOINT: ..." before and after every terra call and large allocation.
# The last CHECKPOINT printed before the crash identifies the failing operation.
#
# Key differences from production 12B + 12C:
#   - Unique per-job terra temp dir to prevent concurrent GDAL conflicts
#   - Inline predict_species_bcr_debug() rather than sourcing 12C
#   - rdata_files limited to first BCR only
#   - n_scen capped at 2 to reduce runtime
# ---

suppressPackageStartupMessages({
  library(BAMexploreR)
  library(gbm)
  library(terra)
  library(tidyverse)
})


# set paths (mirrors 12B exactly) -----------------------------------------------

nm_root <- "/home/mannfred/projects/def-ecknight/NationalModels"

cc    <- TRUE
local <- FALSE

if (cc)            { ia_dir <- "/home/mannfred/scratch/impact_assessment" }
if (!cc && local)  { ia_dir <- getwd() }
if (!cc && !local) { ia_dir <- file.path("G:/Shared drives/BAM_NationalModels5", "data", "Extras", "sandbox_data", "impactassessment_sandbox") }


# unique terra temp dir — prevents concurrent GDAL file conflicts
.terra_tmp <- file.path(
  Sys.getenv("SLURM_TMPDIR", unset = tempdir()),
  paste0("terra_debug_",
         Sys.getenv("SLURM_JOB_ID", unset = "local"), "_1")
)
dir.create(.terra_tmp, recursive = TRUE, showWarnings = FALSE)
terra::terraOptions(tempdir = .terra_tmp)
message("CHECKPOINT: terra tempdir = ", .terra_tmp)


# import data (mirrors 12B exactly) ---------------------------------------------

bam_boundary         <- terra::vect(file.path(ia_dir, "data", "raw_data", "Regions", "BAM_BCR_NationalModel_Unbuffered.shp"))
all_subbasins_subset <- terra::vect(file.path(ia_dir, "data", "raw_data", "hydrobasins_masked_merged_subset.gpkg"))

bcr_subbasins_ref <- {
  hits <- terra::relate(all_subbasins_subset, bam_boundary, relation = "intersects")
  ij   <- which(hits, arr.ind = TRUE)
  tibble(
    sub_index = ij[, 1],
    HYBAS_ID  = all_subbasins_subset$first_HYBAS_ID[ij[, 1]],
    bcr_label = paste(bam_boundary$country[ij[, 2]], bam_boundary$subUnit[ij[, 2]], sep = "_"),
    bcr_code  = gsub("_", "", bcr_label))
}

categorical_responses <- c("ABoVE_1km", "NLCD_1km", "MODISLCC_1km", "MODISLCC_5x5", "SCANFI_1km", "VLCE_1km")

predictor_metadata <-
  dplyr::tibble(BAMexploreR::predictor_metadata) |>
  dplyr::filter(version == "v5") |>
  dplyr::select(predictor, definition, predictor_class) |>
  dplyr::mutate(dplyr::across('predictor', stringr::str_replace, 'Year', 'year')) |>
  dplyr::mutate(dplyr::across('predictor', stringr::str_replace, 'Method', 'method'))

actually_biotic_what <- c("Peatland_5x5", "Peatland_1km")
actually_biotic_df   <- tibble::tibble(predictor = actually_biotic_what,
                                       predictor_class = c("Wetland", "Wetland"))

abiotic_vars <-
  predictor_metadata |>
  dplyr::filter(predictor_class %in% c("Annual Climate", "Climate Normals", "Topography", "Wetland", "Disturbance", "Time", "Method")) |>
  dplyr::filter(!(predictor %in% actually_biotic_what))

biotic_continuous_vars <-
  predictor_metadata |>
  dplyr::filter(!(predictor_class %in% abiotic_vars$predictor_class)) |>
  dplyr::bind_rows(actually_biotic_df) |>
  dplyr::filter(!predictor %in% categorical_responses) |>
  dplyr::pull(predictor)

industry_rasters <- list.files(file.path(ia_dir, "data", "raw_data", "hirshpearson"),
                               pattern = "\\.tif$", full.names = TRUE)
industry_names   <- sub(".tif", "", basename(industry_rasters))
industry_stack   <- setNames(lapply(industry_rasters, terra::rast), industry_names)

disturbance_vars <-
  dplyr::tibble(BAMexploreR::predictor_metadata) |>
  dplyr::filter(version == "v5") |>
  dplyr::select(predictor, definition, predictor_class) |>
  dplyr::filter(predictor_class == "Disturbance")

source(file.path(ia_dir, "Rscripts", "12E_shapley_utils.R"))
# Note: we do NOT source 12C — the debug function is defined inline below


# hardcoded parameters ----------------------------------------------------------

species      <- "BANS"
coalition_id <- 2L
year         <- 2020
hirsh_dir    <- file.path(ia_dir, "data", "raw_data", "hirshpearson")

coalition <- coalition_id_to_sectors(coalition_id, canonical_sectors())
message("CHECKPOINT: species=", species,
        " coalition_id=", coalition_id,
        " coalition=", paste(coalition, collapse = "+"))


# inline debug function ---------------------------------------------------------

predict_species_bcr_debug <- function(species, year, all_subbasins_subset,
                                      coalition, coalition_id, hirsh_dir) {
  
  coalition_label <- paste(coalition, collapse = "+")
  message("CHECKPOINT: entered predict_species_bcr_debug")
  
  message("CHECKPOINT: about to list.files for rdata_files")
  rdata_files <- list.files(file.path(nm_root, "output/06_bootstraps", species),
                            pattern = "can.*\\.Rdata$", full.names = TRUE)
  message("CHECKPOINT: found ", length(rdata_files), " rdata files; limiting to 1")
  if (length(rdata_files) == 0) stop("No rdata files found for species=", species)
  rdata_files <- rdata_files[1]
  
  lapply(rdata_files, function(rdata_path) {
    on.exit({ gc(); message("CHECKPOINT: on.exit gc done") })
    
    message("CHECKPOINT: about to load ", basename(rdata_path))
    e <- new.env(parent = emptyenv())
    load(rdata_path, envir = e)
    if (!exists("b.list", envir = e)) stop("b.list not found in ", basename(rdata_path))
    b.list   <- e$b.list
    bcr_code <- attr(b.list[[1]], "bcr")
    message("CHECKPOINT: load done; bcr=", bcr_code,
            "; mem=", round(sum(gc()[, 2])), " MB")
    
    sub_ids <-
      bcr_subbasins_ref |>
      dplyr::filter(bcr_code == !!bcr_code) |>
      dplyr::pull(sub_index) |>
      unique()
    message("CHECKPOINT: sub_ids found; n=", length(sub_ids))
    
    if (length(sub_ids) == 0) {
      message("CHECKPOINT: no subbasins in this BCR — skipping")
      return(NULL)
    }
    
    message("CHECKPOINT: about to terra::rast for stack_obs")
    stack_obs <- terra::rast(file.path(nm_root, "gis/stacks",
                                       paste0(bcr_code, "_", year, ".tif")))
    message("CHECKPOINT: stack_obs loaded; nlyr=", terra::nlyr(stack_obs),
            "; mem=", round(sum(gc()[, 2])), " MB")
    
    message("CHECKPOINT: about to terra::rast for CanHF")
    canHF_raw <- terra::rast(file.path(hirsh_dir, "CanHF_1km_morethan1.tif"))
    message("CHECKPOINT: canHF loaded; about to terra::project (can take 30s)")
    canHF_r <- terra::project(canHF_raw, stack_obs, method = "near")
    rm(canHF_raw); gc()
    message("CHECKPOINT: canHF projected; mem=", round(sum(gc()[, 2])), " MB")
    
    message("CHECKPOINT: about to per-sector terra::project loop (", length(coalition), " sectors)")
    union_mask <- terra::rast(stack_obs[[1]]); terra::values(union_mask) <- 0L
    for (sec in coalition) {
      message("CHECKPOINT: projecting sector=", sec)
      sec_r <- terra::project(terra::rast(file.path(hirsh_dir, paste0(sec, ".tif"))),
                              stack_obs, method = "near")
      union_mask <- terra::ifel(sec_r > 0, 1L, union_mask)
      message("CHECKPOINT: sector=", sec, " done")
    }
    sector_mask <- terra::ifel((union_mask == 1L) & (canHF_r >= 1), 1, NA)
    rm(canHF_r, union_mask); gc()
    message("CHECKPOINT: coalition mask built; active_px=",
            terra::global(sector_mask, "notNA")[[1]],
            "; mem=", round(sum(gc()[, 2])), " MB")
    
    bf_mosaic_path <- file.path(ia_dir, "data", "derived_data", "bart_models_mosaics",
                                year, paste0(bcr_code, "_backfilled.tif"))
    if (!file.exists(bf_mosaic_path)) {
      message("CHECKPOINT: no backfilled mosaic found for ", bcr_code, " — skipping BCR")
      return(NULL)
    }
    message("CHECKPOINT: about to terra::rast for bf_mosaic")
    stack_bf <- terra::rast(bf_mosaic_path)
    message("CHECKPOINT: bf_mosaic loaded; nlyr=", terra::nlyr(stack_bf),
            "; mem=", round(sum(gc()[, 2])), " MB")
    
    draw_layer_names <- names(stack_bf)[grep("_draw_[0-9]{3}$", names(stack_bf))]
    draw_covs        <- unique(sub("_draw_[0-9]{3}$", "", draw_layer_names))
    n_draws          <- 100L
    n_scen           <- 2L  # capped at 2 for diagnostic (production uses 100)
    
    model_vars_shared  <- b.list[[1]]$var.names
    cat_vars_shared    <- intersect(model_vars_shared, categorical_responses)
    dist_shared        <- intersect(disturbance_vars$predictor, model_vars_shared)
    biotic_cont_shared <- intersect(setdiff(model_vars_shared, categorical_responses),
                                    biotic_continuous_vars)
    
    message("CHECKPOINT: about to terra::rasterize subbasin polygons")
    subbasin_zone_r <- terra::rasterize(all_subbasins_subset[sub_ids, ],
                                        stack_obs[[1]], field = "first_HYBAS_ID")
    message("CHECKPOINT: rasterize done; mem=", round(sum(gc()[, 2])), " MB")
    
    sector_cell_idx <- which(!is.na(terra::values(sector_mask, mat = FALSE)))
    sector_zones    <- terra::values(subbasin_zone_r, mat = FALSE)[sector_cell_idx]
    message("CHECKPOINT: sector_cell_idx computed; coalition_pixels=", length(sector_cell_idx))
    
    message("CHECKPOINT: about to terra::values(stack_obs) — large alloc")
    obs_all_vals <- terra::values(stack_obs)
    X_obs_sector <- as.data.frame(obs_all_vals[sector_cell_idx, , drop = FALSE],
                                  check.names = FALSE)
    rm(obs_all_vals); gc()
    message("CHECKPOINT: obs values extracted; dim=",
            nrow(X_obs_sector), "x", ncol(X_obs_sector),
            "; mem=", round(sum(gc()[, 2])), " MB")
    
    message("CHECKPOINT: about to extract BART draw values (", length(draw_covs),
            " covs x ", n_draws, " draws)")
    draw_vals_sector <- setNames(
      lapply(draw_covs, function(v) {
        sapply(seq_len(n_draws), function(d) {
          lyr_vals <- terra::values(stack_bf[[paste0(v, "_draw_", sprintf("%03d", d))]],
                                    mat = FALSE)
          vals <- expm1(lyr_vals[sector_cell_idx])
          vals[!is.finite(vals)] <- 0
          pmax(vals, 0)
        })
      }),
      draw_covs
    )
    message("CHECKPOINT: BART draws extracted; mem=", round(sum(gc()[, 2])), " MB")
    
    cat_vals_sector <- setNames(
      lapply(cat_vars_shared, function(v)
        terra::values(stack_bf[[v]], mat = FALSE)[sector_cell_idx]),
      cat_vars_shared
    )
    message("CHECKPOINT: categorical values extracted")
    
    obs_dir       <- file.path(ia_dir, "data", "derived_data", "predictions",
                               species, bcr_code, year)
    obs_boot_path <- file.path(obs_dir, "observed_bootstraps.tif")
    if (!file.exists(obs_boot_path)) stop("observed_bootstraps.tif not found: ", obs_boot_path)
    
    message("CHECKPOINT: about to terra::rast(obs_boot_path)")
    obs_boot_stack <- terra::rast(obs_boot_path)
    n_boot  <- terra::nlyr(obs_boot_stack)
    obs_preds <- lapply(seq_len(n_boot), function(i) obs_boot_stack[[i]])
    rm(obs_boot_stack)
    message("CHECKPOINT: obs_boot_stack loaded; n_boot=", n_boot,
            "; mem=", round(sum(gc()[, 2])), " MB")
    
    bf_preds <- vector("list", n_boot)
    for (i in seq_len(n_boot)) bf_preds[[i]] <- vector("list", n_scen)
    
    for (i in seq_along(b.list)) {
      model <- b.list[[i]]
      for (k in seq_len(n_scen)) {
        set.seed((sum(utf8ToInt(paste0(species, bcr_code))) + i * 1000L + k) %%
                   .Machine$integer.max)
        chosen <- sample(n_draws, 1)
        X_k <- X_obs_sector
        for (v in draw_covs)       { if (v %in% names(X_k)) X_k[[v]] <- draw_vals_sector[[v]][, chosen] }
        for (v in cat_vars_shared) { if (v %in% names(X_k)) X_k[[v]] <- cat_vals_sector[[v]] }
        for (v in dist_shared)     { if (v %in% names(X_k)) X_k[[v]] <- 0 }
        
        message("CHECKPOINT: about to gbm::predict.gbm; boot=", i, " scen=", k)
        bf_preds[[i]][[k]] <- gbm::predict.gbm(model, X_k,
                                               n.trees = model$n.trees,
                                               type = "response")
        message("CHECKPOINT: predict.gbm done; boot=", i, " scen=", k)
      }
      gc()
    }
    
    message("CHECKPOINT: about to agg() — population aggregation")
    zone_ids  <- sort(unique(terra::values(subbasin_zone_r, na.rm = TRUE)))
    hybas_ids <- all_subbasins_subset$first_HYBAS_ID[sub_ids]
    idx       <- match(hybas_ids, zone_ids)
    valid_px         <- !is.na(sector_zones)
    sector_zones_flt <- sector_zones[valid_px]
    
    for (i in seq_len(n_boot)) {
      obs_r     <- obs_preds[[i]] * 100
      obs_on_fp <- terra::mask(obs_r, sector_mask)
      terra::zonal(obs_r,     subbasin_zone_r, "sum", na.rm = TRUE)
      terra::zonal(obs_on_fp, subbasin_zone_r, "sum", na.rm = TRUE)
      for (k in seq_len(n_scen)) {
        bf_vals    <- bf_preds[[i]][[k]] * 100
        tapply(bf_vals[valid_px], sector_zones_flt, sum, na.rm = TRUE)
      }
    }
    
    message("CHECKPOINT: agg done; returning result")
    tibble(species = species, subbasin = sub_ids, bcr = bcr_code,
           coalition_id = coalition_id)
  })
}


# run ---------------------------------------------------------------------------

message("CHECKPOINT: calling predict_species_bcr_debug")
res <- predict_species_bcr_debug(
  species              = species,
  year                 = year,
  all_subbasins_subset = all_subbasins_subset,
  coalition            = coalition,
  coalition_id         = coalition_id,
  hirsh_dir            = hirsh_dir
)
message("CHECKPOINT: completed without crash")
