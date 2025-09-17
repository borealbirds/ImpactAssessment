# ---
# title: Impact Assessment: mosaic BCR covariate stacks and crop to subbasin extent
# author: Mannfred Boehm
# created: May 7, 2025
# ---

library(BAMexploreR)
library(tidyverse)
library(terra)

# set root path
root <- "G:/Shared drives/BAM_NationalModels5"
ia_dir <- file.path(root, "data", "Extras", "sandbox_data", "impactassessment_sandbox")

# import 2020 subbasin multi-polygon (all other years are a subset of this area)
all_subbasins <- vect(file.path(ia_dir, "hydrobasins_subset_2020.gpkg"))

# identify BCRs needed for mosaicing (19 BCRs)
bcrs_needed <- BAMexploreR::bam_get_bcr(version = "v5", ext = all_subbasins)
  
# define temporal scope
years <- seq(from = 1990, to = 2020, by = 5)

# -----------------------------------------------------
# for every year, import BCR covariate stacks and mosaic


build_mosaics_by_year <- function(
    root,
    years,
    all_subbasins,
    bcrs_needed = bcrs_needed,
    categorical_predictors = c("ABoVE_1km","method","NLCD_1km","MODISLCC_1km",
                               "MODISLCC_5x5","SCANFI_1km","VLCE_1km"),
    stacks_dir = file.path(root, "gis", "stacks"),
    outdir     = file.path(ia_dir)
){
  
  # define output directory and file names
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  out_files <- setNames(vector("list", length(years)), as.character(years))
  target_crs <- terra::crs(all_subbasins)
  
  # for every year, find all BCR stacks
  for (y in years) {

    files_y <- file.path(stacks_dir, paste0(bcrs_needed, "_", y, ".tif"))
    files_y <- files_y[file.exists(files_y)]
    
    # import relevant stacks for current year and crop 
    stacks <- lapply(files_y, function(f) {
      
      # import stack f
      r <- terra::rast(f)
      
      # guard tiles with missing crs
      if (is.na(terra::crs(r)) || terra::crs(r) == "") terra::crs(r) <- target_crs
      
      # project the bounding box of the all_subbasins (fast)
      # then crop the current stack (less memory downstream)
      subb_local_ext <- terra::project(
        terra::ext(all_subbasins),
        from = target_crs,
        to   = terra::crs(r)
      )
     
      terra::crop(r, subb_local_ext)  # extent crop only
    })
    
    # index covariate names from the current stack
    all_vars <- unique(unlist(lapply(stacks, names)))
  
    # progress bar for per-variable loop
    pb   <- utils::txtProgressBar(min = 0, max = length(all_vars), style = 3)
    prog <- 0L
    
    # create containers for mosaics to go into
    mos_layers <- vector("list", length(all_vars))
    names(mos_layers) <- all_vars
    tmpl_5072 <- NULL
    
    # nested loop over layers within stacks: nm is a single covariate (layer) name
    for (nm in all_vars) {
      
      # for every stack, if nm (covariate) exists pull it as a single SpatRaster
      lyr_list <- lapply(stacks, function(r) if (nm %in% names(r)) r[[nm]] else NULL)
      
      # drop NULLs so we're left with only the stacks that actually have nm layer
      lyr_list <- Filter(Negate(is.null), lyr_list)
      
      # if no layers had that variable, skip mosaicing this variable and move on
      if (length(lyr_list) == 0L) next
      
      # check if nm layer is categorical
      is_cat <- nm %in% categorical_predictors
      if (is_cat) lyr_list <- lapply(lyr_list, terra::as.factor)
      
      # if layer exists in more than one BCR, mosaic 
      mos <- if (length(lyr_list) == 1L) {
        lyr_list[[1]]
      } else {
        do.call(
          terra::mosaic,
          c(lyr_list, list(fun = if (is_cat) "modal" else "mean"))
        )
      }
      
      # ensure categorical predictors remain factors after mosaicing
      names(mos) <- nm
      if (is_cat) mos <- terra::as.factor(mos)
      
      
      # project current covariate (nm) to 5072; first one becomes the template grid
      if (is.null(tmpl_5072)) {
        mos_5072 <- terra::project(
          mos,
          target_crs,
          method   = if (is_cat) "near" else "bilinear",
          filename = tempfile(fileext = ".tif"), overwrite = TRUE # use tmp file to avoid storing in RAM
        )
        tmpl_5072 <- mos_5072[[1]]
      } else {
        mos_5072 <- terra::project(
          mos,
          tmpl_5072,
          method   = if (is_cat) "near" else "bilinear",
          filename = tempfile(fileext = ".tif"), overwrite = TRUE
        )
      }
      
      # ensure categorical layers are factors
      if (is_cat) mos_5072 <- terra::as.factor(mos_5072)
      
      # store
      mos_layers[[nm]] <- mos_5072
      
      # update progress bar
      prog <- prog + 1L
      if (prog %% 5L == 0L || prog == length(all_vars)) utils::setTxtProgressBar(pb, prog)
      
      } # close nested loop
    
    # finish progress bar
    close(pb)
    
    # re-stack mosaiced layers
    mos_layers <- Filter(Negate(is.null), mos_layers)
    if (length(mos_layers) == 0L) { 
      message("No variables mosaicked for ", y); 
      next 
    }
    r_mos <- if (length(mos_layers) == 1L) mos_layers[[1]] else terra::rast(mos_layers)
    
    
    # CAfire: add time-since-disturbance layer for this year ---
    caf_path <- file.path(ia_dir, paste0("CAfire_", y, "_masked.tif"))
   
     if (file.exists(caf_path)) {
       
      caf <- terra::rast(caf_path)
      
      if (is.na(terra::crs(caf)) || terra::crs(caf) == "") terra::crs(caf) <- target_crs
      
      caf_aligned <- terra::project(
        caf, r_mos, method = "bilinear",
        filename = tempfile(fileext = ".tif"), overwrite = TRUE
      )
      names(caf_aligned) <- "CAfire"
      
      if ("CAfire" %in% names(r_mos)) r_mos <- r_mos[[setdiff(names(r_mos), "CAfire")]]
      r_mos <- c(r_mos, caf_aligned)
      
    } else {
      
      message("Year ", y, ": CAfire not found at ", caf_path, " (skipping)")
      
    }
    
    # crop and mask to subbasin extent
    r_out <- 
      terra::crop(x = r_mos, y = all_subbasins) |> 
      terra::mask(x = _, mask = all_subbasins)
    
    # write out
    out_path <- file.path(outdir, paste0("covariates_mosaiced_", y, ".tif"))
    terra::writeRaster(r_out, out_path, overwrite = TRUE,
                       wopt = list(gdal = c("COMPRESS=LZW", "ZLEVEL=9")))
    message("writing ", out_path)
    out_files[[as.character(y)]] <- out_path
  }
  
  # returns without printing
  invisible(out_files)
}

build_mosaics_by_year(root, 
                      years = years,
                      bcrs_needed = bcrs_needed, 
                      all_subbasins = all_subbasins)






# add CAfire to mosaiced stacks
years <- seq(1990, 2020, by=5)

add_cafire_to_mosaics <- function(ia_dir, years){
  
  for (y in years) {
    
    r_mos   <- terra::rast(file.path(ia_dir,  paste0("covariates_mosaiced_", y, ".tif")))           # EPSG:5072, cropped/masked already
    cafire  <- terra::rast(file.path(ia_dir,  paste0("CAfire_", y, "_masked.tif")))
    
    # reproject/resample CAfire to match r_mos geometry 
    cafire_aligned <- terra::project(cafire, r_mos, method = "bilinear")
    names(cafire_aligned) <- "CAfire"
    
    r_out <- c(r_mos, cafire_aligned)
    
    terra::writeRaster(
      r_out, paste0("covariates_mosaiced_", y, ".tif"), overwrite = TRUE,
      wopt = list(gdal = c("COMPRESS=LZW", "ZLEVEL=9")))
    
    message("added CAfire to ", paste0("covariates_mosaiced_", y, ".tif"))
  }
  
  invisible(TRUE)
}

add_cafire_to_mosaics(ia_dir, years)



