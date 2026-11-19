# ---
# title: Impact Assessment: mosaic BCR covariate stacks and crop to subbasin extent
# author: Mannfred Boehm
# created: May 7, 2025
# ---

# per year: Elly's 19 BCR covariate stacks, mosaicked layer by layer onto one EPSG:5072 grid,
# plus CAfire (02) and the soil layers (01), cropped to the subbasins (05)
#   -> data/raw_data/covariates_mosaiced/covariates_mosaiced_{year}.tif
# runs on Fir, one array task per year: bash 06_build_covariate_stacks.sh 1995 2000

library(BAMexploreR)
library(terra)

terraOptions(memfrac = 0.60)   # project() returns different values under memory pressure

# terra sizes its blocks from the node's free RAM, not the job's allocation (see 11)
slurm_mem_mb <- suppressWarnings(as.numeric(Sys.getenv("SLURM_MEM_PER_NODE")))
if (is.finite(slurm_mem_mb)) terraOptions(memmax = 0.3 * slurm_mem_mb / 1024)

nm_root <- "/home/mannfred/projects/def-ecknight/NationalModels"
raw_dir <- "/home/mannfred/scratch/impact_assessment/data/raw_data"
out_dir <- file.path(raw_dir, "covariates_mosaiced")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

year <- as.integer(Sys.getenv("YEARS", "2020"))   # the .sh passes one year per task

all_subbasins_subset <- terra::vect(file.path(raw_dir, "hydrobasins_masked_merged_subset.gpkg"))
target_crs <- terra::crs(all_subbasins_subset)

categorical_predictors <- c("ABoVE_1km", "method", "NLCD_1km", "MODISLCC_1km", "MODISLCC_5x5",
                            "SCANFI_1km", "VLCE_1km")



# -----------------------------------------------------
# mosaic the BCR stacks, one covariate at a time

bcrs_needed <- BAMexploreR::bam_get_bcr(version = "v5", ext = all_subbasins_subset)
files_y <- file.path(nm_root, "gis", "stacks", paste0(bcrs_needed, "_", year, ".tif"))
if (!all(file.exists(files_y))) stop("no BCR stack at ", paste(files_y[!file.exists(files_y)], collapse = ", "))
message("year ", year, ": mosaicking ", length(files_y), " BCR stacks")

# each stack cropped to the subbasins' bounding box
stacks <- lapply(files_y, function(f) {
  r <- terra::rast(f)
  if (is.na(terra::crs(r)) || terra::crs(r) == "") terra::crs(r) <- target_crs
  terra::crop(r, terra::project(terra::ext(all_subbasins_subset), from = target_crs, to = terra::crs(r)))
})

all_vars   <- unique(unlist(lapply(stacks, names)))
mos_layers <- setNames(vector("list", length(all_vars)), all_vars)
tmpl_5072  <- NULL

for (nm in all_vars) {

  # this covariate from every stack that has it
  lyr_list <- Filter(Negate(is.null), lapply(stacks, function(r) if (nm %in% names(r)) r[[nm]] else NULL))

  is_cat <- nm %in% categorical_predictors
  if (is_cat) lyr_list <- lapply(lyr_list, terra::as.factor)

  # where BCRs overlap: the modal class, or the mean
  mos <- if (length(lyr_list) == 1L) lyr_list[[1]] else
    do.call(terra::mosaic, c(lyr_list, list(fun = if (is_cat) "modal" else "mean")))
  names(mos) <- nm
  if (is_cat) mos <- terra::as.factor(mos)

  # to EPSG:5072, written to a temp file to keep it out of RAM; the first covariate sets the grid
  mos_5072 <- terra::project(mos, if (is.null(tmpl_5072)) target_crs else tmpl_5072,
                             method = if (is_cat) "near" else "bilinear",
                             filename = tempfile(fileext = ".tif"), overwrite = TRUE)
  if (is.null(tmpl_5072)) tmpl_5072 <- mos_5072[[1]]
  if (is_cat) mos_5072 <- terra::as.factor(mos_5072)

  mos_layers[[nm]] <- mos_5072
}

r_mos <- terra::rast(mos_layers)



# -----------------------------------------------------
# add CAfire and soil, crop to the subbasins

caf_path <- file.path(raw_dir, "CAfire", paste0("CAfire_", year, "_masked.tif"))
if (!file.exists(caf_path)) stop("no CAfire layer at ", caf_path, " (built by 02)")
caf_aligned <- terra::project(terra::rast(caf_path), r_mos, method = "bilinear",
                              filename = tempfile(fileext = ".tif"), overwrite = TRUE)
names(caf_aligned) <- "CAfire"
r_mos <- c(r_mos[[setdiff(names(r_mos), "CAfire")]], caf_aligned)   # replaces any CAfire layer a stack carries

soil_aligned <-
  terra::project(terra::rast(file.path(raw_dir, "isric_soil", "isric_soil_covariates_masked.tif")),
                 r_mos, method = "bilinear", filename = tempfile(fileext = ".tif"), overwrite = TRUE) |>
  terra::crop(r_mos) |>
  terra::mask(r_mos[[1]])
r_mos <- c(r_mos, soil_aligned)

r_out <-
  terra::crop(r_mos, all_subbasins_subset) |>
  terra::mask(all_subbasins_subset)

# 07 reads covariates by name, so every year needs 2020's layers, in 2020's order, on its grid
ref_path <- file.path(out_dir, "covariates_mosaiced_2020.tif")
if (file.exists(ref_path)) {
  ref <- terra::rast(ref_path)
  if (!identical(names(r_out), names(ref)))
    stop("layers differ from ", ref_path, ". only here: ", paste(setdiff(names(r_out), names(ref)), collapse = ", "),
         "; only there: ", paste(setdiff(names(ref), names(r_out)), collapse = ", "))
  terra::compareGeom(r_out, ref, stopOnError = TRUE)
}

out_path <- file.path(out_dir, paste0("covariates_mosaiced_", year, ".tif"))
terra::writeRaster(r_out, out_path, overwrite = TRUE, wopt = list(gdal = c("COMPRESS=LZW", "ZLEVEL=9")))
message("wrote ", out_path)
