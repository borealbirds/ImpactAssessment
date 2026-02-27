# ---
# title: Impact Assessment: re-predict bird densities after backfilling
# author: Mannfred Boehm
# created: September 28, 2025
# ---

suppressPackageStartupMessages({
  library(BAMexploreR)
  library(gbm)
  library(terra)
  library(tidyverse)
})


# set paths ------------------------------------------------------

# colleague's NationalModels directory — paths here must not change
nm_root <- "/home/mannfred/projects/def-ecknight/NationalModels"

cc    <- TRUE    # TRUE = Compute Canada cluster
local <- FALSE   # TRUE = local RProject machine

if (cc)            { ia_dir <- "/home/mannfred/scratch/impact_assessment" }
if (!cc && local)  { ia_dir <- getwd() }
if (!cc && !local) { ia_dir <- file.path("G:/Shared drives/BAM_NationalModels5", "data", "Extras", "sandbox_data", "impactassessment_sandbox") }


# import data ------------------------------------------------------

# import BCR boundaries
bam_boundary <- terra::vect(file.path(ia_dir, "data", "raw_data", "Regions", "BAM_BCR_NationalModel_Unbuffered.shp"))

# import merged + subsetted subbasins
all_subbasins_subset <- terra::vect(file.path(ia_dir, "data", "raw_data", "hydrobasins_masked_merged_subset.gpkg"))

# ------------------------------------------------------
# create a reference table for which subbasins are in which BCRs 
# some subbasins will be in multiple BCRs, and that's OK

bcr_subbasins_ref <- {
  
    # logical matrix: rows=subbasins, cols=BCRs
    hits <- terra::relate(centroids(all_subbasins_subset), bam_boundary, relation = "intersects")
    
    # row/col indexes of TRUE
    ij <- which(hits, arr.ind = TRUE)
    
    tibble(
       sub_index = ij[, 1],  
       HYBAS_ID  = all_subbasins_subset$first_HYBAS_ID[ij[, 1]],
       bcr_label = paste(bam_boundary$country[ij[, 2]],
                    bam_boundary$subUnit[ij[, 2]], sep = "_"),
       bcr_code  = gsub("_", "", bcr_label))
}

# ------------------------------------------------------
# define helper function:
# mosaic (backfilled) subbasin stacks into a single BCR-wide stack

mosaic_backfilled_stacks <- function(sub_ids, year) {
  
  paths <- file.path(
    file.path(ia_dir, "data", "derived_data", "bart_models", year),
    paste0("subbasin_", sub_ids),
    paste0("subbasin_", sub_ids, "_backfill.tif")
  )
  
  paths <- paths[file.exists(paths)]
  if (length(paths) == 0) return(NULL)
  
  stacks <- lapply(paths, terra::rast)
  
  ref <- stacks[[1]]
  stacks <- lapply(stacks, \(r)
                   if (!terra::compareGeom(r, ref, stopOnError = FALSE))
                     terra::resample(r, ref, method = "near") else r
  )
  
  vars <- sort(unique(unlist(lapply(stacks, names))))
  out <- terra::rast(lapply(vars, function(v) {
  # get only stacks that actually have this variable
  available <- lapply(stacks, function(r) if (v %in% names(r)) r[[v]] else NULL)
  available <- Filter(Negate(is.null), available)
  
  if (length(available) == 0) return(NULL)  # skip variable entirely
  Reduce(terra::cover, available)
  }))
  
  out <- out[[!sapply(out, is.null)]]
  names(out) <- vars[!sapply(out, is.null)]
  out
}

make_counterfactual_stack <- function(stack_obs, mu_stack, sd_stack, z) {

    out <- stack_obs

    # overwrite continuous biotic covariates
    for (v in names(mu_stack)) {
      out[[v]] <- mu_stack[[v]] + z * sd_stack[[v]]
    }

    out
}



# define covariate types -------------------------------------------------------------

# define categorical covariates
categorical_responses <- c("ABoVE_1km", "NLCD_1km","MODISLCC_1km", "MODISLCC_5x5","SCANFI_1km","VLCE_1km")

predictor_metadata <-
  dplyr::tibble(BAMexploreR::predictor_metadata) |>
  dplyr::filter(version == "v5") |>
  dplyr::select(predictor, definition, predictor_class) |>
  dplyr::mutate(dplyr::across('predictor', stringr::str_replace, 'Year', 'year')) |>
  dplyr::mutate(dplyr::across('predictor', stringr::str_replace, 'Method','method'))


# convert some abiotic variables to biotic variables
actually_biotic_what <- c("Peatland_5x5", "Peatland_1km")
actually_biotic_df <- tibble::tibble(predictor = actually_biotic_what, predictor_class = c("Wetland", "Wetland"))

# define abiotic variables
abiotic_vars <-
  predictor_metadata |>
  dplyr::filter(predictor_class %in% c("Annual Climate", "Climate Normals", "Topography", "Wetland", "Disturbance", "Time", "Method")) |>
  dplyr::filter(!(predictor %in% actually_biotic_what))

# define biotic variables
biotic_continuous_vars <-
  predictor_metadata |>
  dplyr::filter(!(predictor_class %in% abiotic_vars$predictor_class)) |>
  dplyr::bind_rows(actually_biotic_df) |>
  dplyr::filter(!predictor %in% categorical_responses) |>
  dplyr::pull(predictor)

# import human footprint rasters ------------------------------------------------------

industry_rasters <- list.files(file.path(ia_dir, "data", "raw_data", "hirshpearson"), pattern = "\\.tif$", full.names = TRUE)
industry_names <- sub(".tif", "", basename(industry_rasters))
industry_stack <- setNames(lapply(industry_rasters, terra::rast), industry_names)

# define disturbance variables, which we'll set to zero when re-predicting
disturbance_vars <-
  dplyr::tibble(BAMexploreR::predictor_metadata) |>
  dplyr::filter(version == "v5") |>
  dplyr::select(predictor, definition, predictor_class) |>
  dplyr::filter(predictor_class == "Disturbance")


# import density prediction script -------------------------------------------------

source(file.path(ia_dir, "Rscripts", "12B_predict_species_bcr.R"))



# run one species on one core ------------------------------------------------------

species_vec <- c("CAWA")
# species_vec <- sort(list.dirs(file.path(nm_root, "output/06_bootstraps"), full.names = FALSE, recursive = FALSE))[4]
year <- 2020

# get species index from SLURM
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
species <- species_vec[task_id]
message("running species: ", species)

res <- predict_species_bcr(species, year = year, all_subbasins_subset = all_subbasins_subset)
dir.create(file.path(ia_dir, "data", "derived_data", "density_tables"), showWarnings = FALSE)
saveRDS(res, file = file.path(ia_dir, "data", "derived_data", "density_tables", paste0(species, "_", year, ".rds")))
message(Sys.time(), " nice.")

