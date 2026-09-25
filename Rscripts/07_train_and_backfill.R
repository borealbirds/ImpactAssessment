# ---
# title: Impact Assessment: train models per subbasin and backfill industry footprints
# author: Mannfred Boehm
# created: August 7, 2025
# ---

#1. attach packages ----------------------------------------------
print("* attaching packages on master *")
library(BART)
library(BAMexploreR)
library(terra)
library(tidyverse)


#2. define local or cluster --------------------------------------
test  <- FALSE
cc    <- TRUE   # TRUE = Compute Canada cluster
local <- FALSE  # TRUE = local RProject machine (overrides Google Drive path)

# if working on cluster, extract the array task index from the SLURM script
args <- commandArgs(trailingOnly = TRUE)
if (cc && length(args) == 0) {
  stop("no task index supplied")
}
task_id <- if (cc) as.integer(args[1]) else 1L

# years to backfill: YEARS, comma-separated, from the submitting shell's environment
# (default 2020), e.g.
#   bash 07_submit_backfill_years.sh 2015 2020
# which exports YEARS=2015,2020 and submits 07 with --array sized to the years (and 11 after
# it). Not --export=ALL,YEARS=2015,2020: sbatch splits --export on commas. Tasks come in blocks of one year: task t runs year YEARS[(t - 1) %/% n_sub + 1]
# on subbasin (t - 1) %% n_sub + 1, with n_sub = 674 subbasins.
years <- suppressWarnings(as.integer(strsplit(trimws(Sys.getenv("YEARS", "2020")), "[, ]+")[[1]]))
if (length(years) == 0L || anyNA(years))
  stop("YEARS must be comma-separated years, e.g. YEARS=2015,2020; got '", Sys.getenv("YEARS"), "'")

#3. set root path ------------------------------------------------
print("* setting root file path *")

if (cc)            { ia_dir <- "/home/mannfred/scratch/impact_assessment" }
if (!cc && local)  { ia_dir <- getwd() }
if (!cc && !local) { ia_dir <- file.path("G:/Shared drives/BAM_NationalModels5", "data", "Extras", "sandbox_data", "impactassessment_sandbox") }

print(ia_dir)

#4. define model covariates --------------------------------------
# define predictor and response variables
# store predictor metadata as a reference
predictor_metadata <-
  dplyr::tibble(BAMexploreR::predictor_metadata) |>
  dplyr::filter(version == "v5") |>
  dplyr::select(predictor, definition, predictor_class) |>
  dplyr::mutate(dplyr::across('predictor', stringr::str_replace, 'Year', 'year')) |>
  dplyr::mutate(dplyr::across('predictor', stringr::str_replace, 'Method','method'))

# define soil covariate names
soil_covs <- tibble::tibble(predictor = c("cec_0-5cm_mean_1000", "cec_100-200cm_mean_1000",
                                  "cec_15-30cm_mean_1000", "cec_30-60cm_mean_1000",  
                                  "cec_5-15cm_mean_1000", "cec_60-100cm_mean_1000", 
                                  "soc_0-5cm_mean_1000",  "soc_100-200cm_mean_1000",
                                  "soc_15-30cm_mean_1000", "soc_30-60cm_mean_1000",  
                                  "soc_5-15cm_mean_1000", "soc_60-100cm_mean_1000"),
                    predictor_class = rep("Soil Properties", 12))

# convert some abiotic variables to biotic variables
actually_biotic_what <- c("Peatland_5x5", "Peatland_1km")
actually_biotic_df <- tibble::tibble(predictor = actually_biotic_what, predictor_class = c("Wetland", "Wetland"))

# define abiotic variables (V5 abiotic + CAfire + soil properties)
abiotic_vars <-
  predictor_metadata |> 
  dplyr::filter(predictor_class %in% c("Annual Climate", "Climate Normals", "Topography", "Wetland", "Disturbance")) |> 
  tibble::add_row(predictor = "CAfire", predictor_class ="Time Since Disturbance") |> 
  dplyr::filter(!(predictor %in% actually_biotic_what)) |> 
  dplyr::bind_rows(soil_covs)

# define biotic variables
biotic_vars <-
  predictor_metadata |> 
  dplyr::filter(!(predictor_class %in% c(abiotic_vars$predictor_class, "Time", "Method"))) |> 
  dplyr::bind_rows(actually_biotic_df)

# re-order biotic variables 
neworder <- readRDS(file = file.path(ia_dir, "data", "raw_data", "biotic_variable_hierarchy.rds"))
biotic_vars <- biotic_vars[match(neworder, biotic_vars$predictor), ]

if(exists("biotic_vars")){
  print("* biotic_vars successfully constructed *")
}else{print("* biotic_vars not constructed *")}


#5. import helper functions ----------------------------------------
 
# logfile function to track progress
make_logger <- function(logfile) { # create a new log file
  dir.create(dirname(logfile), recursive = TRUE, showWarnings = FALSE)
  function(fmt, ...) { # write a new line
    line <- sprintf("[%s pid=%d host=%s] %s\n",
                    format(Sys.time(), "%F %T"),
                    Sys.getpid(),
                    Sys.info()[["nodename"]],
                    sprintf(fmt, ...))
    cat(line, file = logfile, append = TRUE)
  } # close new line writing function
} # close file generating function

# training and backfilling function (subbasin level)
source(file.path(ia_dir, "Rscripts", "08A_train_and_backfill_subbasin_s.R"))


#6. train models and backfill biotic features for year y -----------------------------

# which (year, subbasin) this task runs
all_subbasins_subset <- terra::vect(file.path(ia_dir, "data", "raw_data", "hydrobasins_masked_merged_subset.gpkg"))
n_sub <- nrow(all_subbasins_subset)
n_task <- n_sub * length(years)
array_max <- suppressWarnings(as.integer(Sys.getenv("SLURM_ARRAY_TASK_MAX", NA)))
if (!is.na(array_max) && array_max < n_task)
  message("WARNING: the array ends at task ", array_max, " but YEARS=", paste(years, collapse = ","),
          " needs 1-", n_task, " (", n_sub, " subbasins per year) - tasks ", array_max + 1, "-", n_task,
          " were never submitted")
if (is.na(task_id) || task_id < 1L || task_id > n_task) {
  message("task ", task_id, " is past the ", n_task, " (year, subbasin) pairs of YEARS=",
          paste(years, collapse = ","), " - nothing to do")
  quit(save = "no", status = 0)
}
year           <- years[(task_id - 1L) %/% n_sub + 1L]
subbasin_index <- (task_id - 1L) %% n_sub + 1L
print(sprintf("* task %d: year %d, subbasin %d *", task_id, year, subbasin_index))

# pre-mosaiced covariate stack for this year (built locally by 06)
stack_path <- file.path(ia_dir, "data", "raw_data", "covariates_mosaiced", sprintf("covariates_mosaiced_%d.tif", year))
if (!file.exists(stack_path))
  stop("no covariate stack for ", year, ": ", stack_path, " - build it with 06_build_covariate_stacks.R ",
       "(its multi-year line) and stage it")
stack_y <- terra::rast(stack_path)

# define categorical features
categorical_responses = c("ABoVE_1km", "NLCD_1km","MODISLCC_1km", "MODISLCC_5x5","SCANFI_1km","VLCE_1km")

# low- and high-HF masks (training vs backfilled pixels) for this year's footprint:
# hirshpearson/CanHF_1km_{lessthan1,morethan1}_{year}.tif, or the undated files for 2020.
# Another year's footprint differs from 2020's, so a year without its own masks stops here;
# HF_MASK_YEAR=<y> at submission uses year y's masks instead (e.g. HF_MASK_YEAR=2020).
mask_year <- as.integer(Sys.getenv("HF_MASK_YEAR", year))
hf_mask_path <- function(kind) {
  dated <- file.path(ia_dir, "data", "raw_data", "hirshpearson", sprintf("CanHF_1km_%s_%d.tif", kind, mask_year))
  if (file.exists(dated)) return(dated)
  if (mask_year == 2020L) return(file.path(ia_dir, "data", "raw_data", "hirshpearson", sprintf("CanHF_1km_%s.tif", kind)))
  stop("no ", kind, " footprint mask for ", mask_year, " (expected ", dated, "). To backfill ", year,
       " with the 2020 masks instead, submit with HF_MASK_YEAR=2020.")
}
print(sprintf("* footprint masks: %s, %s *", hf_mask_path("lessthan1"), hf_mask_path("morethan1")))

# import low hf layer and project to current stack
lowhf_mask <- terra::rast(hf_mask_path("lessthan1"))
lowhf_mask <- terra::project(x=lowhf_mask, y=stack_y, method = "near")

# import high hf layer and project to current stack
highhf_mask <- terra::rast(hf_mask_path("morethan1"))
highhf_mask <- terra::project(x=highhf_mask, y=stack_y, method = "near")

# project subbasin boundaries to current stack
all_subbasins_subset <- terra::project(x=all_subbasins_subset, y=stack_y)
                
backfill_results <- tryCatch(
                  train_and_backfill_subbasin_s(
                  subbasin_index = subbasin_index, 
                  year           = year,
                  stack_y        = stack_y,
                  lowhf_mask     = lowhf_mask,
                  highhf_mask    = highhf_mask,
                  abiotic_vars   = abiotic_vars, 
                  biotic_vars    = biotic_vars,
                  ia_dir         = ia_dir,
                  quiet          = FALSE,
                  neworder       = neworder,
                  categorical_responses = categorical_responses,
                  all_subbasins_subset  = all_subbasins_subset,
                  cc = cc
                ), # close train_and_backfill_subbasin_s
                
                error = function(e) {
                  message("Error in subbasin ", subbasin_index, ": ", conditionMessage(e))
                  return(list(
                    subbasin = subbasin_index,
                    error = conditionMessage(e)
                  ))
                } # close error
                
            ) # close trycatch
   

#11. stop the cluster----
#print("* stopping cluster :-)*")
#stopCluster(cl)

#12. save backfilled raster for this species x year
print("* saving raster file *")
print(backfill_results)

#if(cc){ q() }

