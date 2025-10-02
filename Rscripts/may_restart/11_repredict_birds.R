# ---
# title: Impact Assessment: re-predict bird densities after backfilling
# author: Mannfred Boehm
# created: September 28, 2025
# ---

library(furrr)
library(terra)
library(tidyterra)
library(tidyverse)


# set root path
root <- "G:/Shared drives/BAM_NationalModels5"
ia_dir <- file.path(root, "data", "Extras", "sandbox_data", "impactassessment_sandbox")



# import data ------------------------------------------------------

# import BCR boundaries
bam_boundary <- terra::vect(file.path(root, "Regions", "BAM_BCR_NationalModel_Unbuffered.shp")) 

# import merged + subsetted subbasins
all_subbasins_subset <- vect(file.path(ia_dir, "hydrobasins_masked_merged_subset.gpkg"))

# import low hf layer
lowhf_mask <- terra::rast(file.path(ia_dir, "CanHF_1km_lessthan1.tif"))



# helper functions -------------------------------------------------

# gets directories for all species density models for a given BCR
.models_for_bcr <- function(models_root, bcr_id) {
  Sys.glob(file.path(models_root, "*", paste0("*_", bcr_id, ".Rdata"))) |>
    tibble(path = _) |>
    mutate(spec_code = basename(dirname(path)), bcr_id = bcr_id)
}

# loads gbm model lists (lists of 32 bootstraps)
.load_gbm_list <- function(path) {
  e <- new.env(parent = emptyenv()); load(path, envir = e)
  nm <- ls(e)
  hit <- purrr::detect(nm, ~ {
    x <- get(.x, envir = e); is.list(x) && length(x) > 0 && inherits(x[[1]], "gbm")
  })
  if (is.null(hit)) stop("no list of gbm models found in: ", path)
  get(hit, envir = e)
}

# finds best trees per bootstrap
.best_trees <- function(m) {
  bt <- tryCatch(m$gbm.call$best.trees, error = function(...) NULL)
  if (!is.null(bt)) return(bt)
  nt <- tryCatch(m$n.trees, error = function(...) NULL)
  if (is.null(nt)) stop("cannot determine n.trees for model.")
  nt
}

# 
.align_df <- function(df, varnames) {
  miss <- setdiff(varnames, colnames(df))
  if (length(miss)) stop("missing predictors: ", paste(miss, collapse = ", "))
  df[, varnames, drop = FALSE]
}





# define density prediction function -------------------------------
predict_bird_density_per_bcr_year <- function(
    year,
    bcr,
       
) {
  
  # PART I: import and format observed and backfilled covariate data ----------------------
  # subset bam_boundary to the current BCR
  bcr_poly  <- bam_boundary[which(paste("country", "subUnit", sep="_") == bcr)]
  
  # find which subbasins are in the current BCR
  subbasin_hits <- unique(intersect(all_subbasins_subset["HYBAS_ID"], bcr_poly)$HYBAS_ID)
  
  # import backfilled subbasins (only those in current BCR) 
  backfill_path_fn <- function(hits) {file.path(ia_dir, "xgboost_models", as.character(year), sprintf("subbasin=%03d", hits), sprintf("backfilled_stack_subbasin-%03d.tif", hits))}
  backfill_paths_table <- data.frame(HYBAS_ID = subbasin_hits, path = vapply(subbasin_hits, backfill_path_fn, FUN.VALUE = character(1)))
  backfilled_stacks <- setNames(lapply(backfill_paths_table$path, terra::rast), backfill_paths_table$HYBAS_ID)

  # convert `backfilled_stacks` list to a single raster stack
  vars <- sort(unique(unlist(lapply(backfilled_stacks, names)))) # get possible layer names
  mosaic_one <- function(v) {  
    rlist1 <- lapply(backfilled_stacks, \(r) if (v %in% names(r)) r[[v]] else NULL) # if a covariate is in a stack, extract that layer
    rlist2 <- purrr::compact(rlist) # remove NULL layers
    Reduce(terra::cover, rlist2) # fill NAs in A with B
  }
  flat_bf <- terra::rast(lapply(vars, mosaic_one)); names(flat_bf) <- vars
  flat_bf <- terra::mask(flat_bf, bcr_poly) # constrain to current BCR polygon
  
  # import industry footprint cropped to current BCR
  industry_bcr <- terra::crop(terra::vect(combined_poly_path), bcr_poly)
  
  # import and crop observed subbasins to industry pixels in the current BCR
  stack_y <- 
    terra::rast(file.path(root, "gis", "stacks", paste0(bcr, "_", year, ".tif"))) |> 
    terra::crop(x = _, y = industry_bcr) |> 
    terra::mask(x = _, mask = industry_bcr)
    
  
  # PART II: find species with models in the current BCR----------------------
  
  
  
  # Precompute per-subbasin footprint cell coordinates (and areas)
  sub_list <- vector("list", length = nrow(sub_in_bcr))
  names(sub_list) <- as.character(sub_in_bcr[[subbasin_id_col]])
  
  for (i in seq_along(sub_list)) {
    sid <- as.integer(sub_in_bcr[[subbasin_id_col]][i])
    bf_path <- backfill_stack_fn(sid)
    if (!file.exists(bf_path)) { sub_list[[i]] <- NULL; next }
    
    r_bf <- rast(bf_path)
    # industry mask for this subbasin on the backfilled grid
    ind_sub <- crop(ind_bcr, sub_in_bcr[i, ])
    if (nrow(ind_sub) == 0) { sub_list[[i]] <- NULL; next }
    mask_sub <- rasterize(ind_sub, r_bf[[1]], field = 1, background = NA, touches = TRUE)
    cells <- which(!is.na(values(mask_sub)))
    if (!length(cells)) { sub_list[[i]] <- NULL; next }
    
    xy <- xyFromCell(r_bf, cells)
    area_vec <- if (use_cell_area) terra::cellSize(r_bf[[1]], unit = "ha")[cells] else rep(1, length(cells))
    sub_list[[i]] <- list(sid = sid, bf_path = bf_path, xy = xy, area = area_vec)
  }
  sub_list <- purrr::compact(sub_list)
  if (!length(sub_list)) return(list(per_species = tibble(), per_subbasin = tibble()))
  
  # Species that actually have models for this BCR
  mdl_tbl <- .models_for_bcr(models_root, bcr_id)
  if (!nrow(mdl_tbl)) return(list(per_species = tibble(), per_subbasin = tibble()))
  
  # Summaries per species (load model once; stream through subbasins)
  per_species <- pmap_dfr(mdl_tbl, function(path, spec_code, bcr_id) {
    gbms <- .load_gbm_list(path)
    varn <- gbms[[1]]$var.names
    nB   <- length(gbms)
    tot_obs <- numeric(nB); tot_bkf <- numeric(nB)
    
    for (sl in sub_list) {
      # observed predictors from mosaic at the subbasin footprint cells
      X_obs <- terra::extract(rs_mosaic[[varn]], sl$xy, ID = FALSE) |>
        as.data.frame() |>
        .align_df(varn)
      
      # backfilled predictors from the subbasin stack at the same cells
      r_bf <- rast(sl$bf_path)
      X_bkf <- terra::extract(r_bf[[varn]], sl$xy, ID = FALSE) |>
        as.data.frame() |>
        .align_df(varn)
      
      # accumulate per-bootstrap totals
      for (b in seq_len(nB)) {
        nt <- .best_trees(gbms[[b]])
        p_obs <- stats::predict(gbms[[b]], newdata = X_obs, type = "response", n.trees = nt)
        p_bkf <- stats::predict(gbms[[b]], newdata = X_bkf, type = "response", n.trees = nt)
        tot_obs[b] <- tot_obs[b] + sum(p_obs * sl$area, na.rm = TRUE)
        tot_bkf[b] <- tot_bkf[b] + sum(p_bkf * sl$area, na.rm = TRUE)
      }
    }
    
    tibble(
      spec_code = spec_code,
      bcr_id = bcr_id,
      year = year,
      observed_mean = mean(tot_obs), observed_sd = sd(tot_obs),
      backfilled_mean = mean(tot_bkf), backfilled_sd = sd(tot_bkf),
      delta_mean = mean(tot_obs - tot_bkf), delta_sd = sd(tot_obs - tot_bkf),
      ratio_mean = mean(tot_obs / tot_bkf), ratio_sd = sd(tot_obs / tot_bkf)
    )
  })
  
  # Count footprint cells per subbasin (handy for QA/QC)
  per_subbasin <- tibble(
    year = year,
    bcr_id = bcr_id,
    sub_id = vapply(sub_list, `[[`, integer(1), "sid"),
    n_cells = vapply(sub_list, function(z) nrow(z$xy), integer(1))
  )
  
  if (!is.null(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    readr::write_csv(per_species,  file.path(out_dir, sprintf("species_summary_%s_%s.csv", bcr_id, year)))
    readr::write_csv(per_subbasin, file.path(out_dir, sprintf("footprint_cells_per_subbasin_%s_%s.csv", bcr_id, year)))
  }
  
  list(per_species = per_species, per_subbasin = per_subbasin)
}


# predict bird densities  ------------------------------------------
# For each subbasin x year × species:
# observed: predict abundance on raster stack with industry footprint
# backfilled: predict abundance on raster stack with footprint replaced by intact vegetation



