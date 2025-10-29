# ---
# title: Impact Assessment: train models per subbasin and backfill industry footprints
# author: Mannfred Boehm
# created: August 7, 2025
# ---

#1. attach packages ----------------------------------------------
print("* attaching packages on master *")
library(BART)
library(BAMexploreR)
library(parallel)
library(terra)
library(tidyverse)


#2. define local or cluster --------------------------------------
test <- TRUE
cc <- FALSE


#3. set number of tasks for local vs cluster ---------------------
if(cc){ n_tasks <- 32}
if(!cc | test){ n_tasks <- 4}


#4. create and register clusters ---------------------------------
# creates 32 copies of R running in parallel via 32 tasks, on one of the cluster's sockets (processors). 
# Belgua has ~965 nodes
print("* creating clusters *")
cl <- makePSOCKcluster(n_tasks, type="PSOCK")

# print number of tasks and host name for confirmation
cl



#5. set root path ------------------------------------------------
print("* setting root file path *")

if(!cc){root <- "G:/Shared drives/BAM_NationalModels5"}
if(cc){root <- "/home/mannfred/scratch"}

ia_dir <- file.path(root, "data", "Extras", "sandbox_data", "impactassessment_sandbox")

tmpcl <- clusterExport(cl, c("root"))


#6. attach packages on clusters ----------------------------------
# `clusterEvalQ` evaluates a literal expression on each cluster node. 
print("* Loading packages on workers *")
tmpcl <- clusterEvalQ(cl, library(BART))
tmpcl <- clusterEvalQ(cl, library(BAMexploreR))
tmpcl <- clusterEvalQ(cl, library(terra))
tmpcl <- clusterEvalQ(cl, library(tidyverse))



#7. define model covariates --------------------------------------
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
actually_biotic_what <- c("StandardDormancy_1km", "StandardGreenup_1km", "Peatland_5x5", "Peatland_1km")
actually_biotic_df <- tibble::tibble(predictor = actually_biotic_what, predictor_class = c("Annual Climate", "Annual Climate", "Wetland", "Wetland"))

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
neworder <- readRDS(file = file.path(ia_dir, "biotic_variable_hierarchy.rds"))
biotic_vars <- biotic_vars[match(neworder, biotic_vars$predictor), ]



#8. define data paths ----------------------------------------

# subbasins polygon
all_subbasins_subset_path <- file.path(ia_dir, "hydrobasins_masked_merged_subset.gpkg")

# high human footprint raster
highhf_mask_path <- file.path(ia_dir, "hirshpearson", "CanHF_1km_morethan1.tif")

# low human footprint raster
lowhf_mask_path <- file.path(ia_dir, "hirshpearson", "CanHF_1km_lessthan1.tif")

# covariate stack for year y
stack_y_path <- file.path(ia_dir, sprintf("covariates_mosaiced_%d.tif", 2020))

# categorical covariates
categorical_responses = c("ABoVE_1km", "NLCD_1km","MODISLCC_1km",
                          "MODISLCC_5x5","SCANFI_1km","VLCE_1km")

# training and backfilling function (subbasin level)
source(file.path(getwd(), "Rscripts", "may_restart", "08_train_and_backfill_subbasin_s.R"))


#9. define model training function ----------------------------------------

# train models per subbasin per year
train_and_backfill_per_year <- function(
    year,
    all_subbasins_subset_path,
    lowhf_mask_path,
    highhf_mask_path,
    stack_y_path,
    categorical_responses,
    quiet = FALSE){
  
 ###
 # import covariate stack, low HF layer for training, high HF layer for backfilling
 ###
 
 # import pre-mosaiced covariate stack for year_y
 stack_y <- terra::rast(stack_y_path)
 
 # ensure categoricals are factors
 cats_present <- intersect(categorical_responses, names(stack_y))
 for (nm in cats_present) stack_y[[nm]] <- terra::as.factor(stack_y[[nm]])

 # import low hf layer and project to current stack
 lowhf_mask <- terra::rast(lowhf_mask_path)
 lowhf_mask <- terra::project(x=lowhf_mask, y=stack_y, method = "near")
 
 # import high hf layer and project to current stack
 highhf_mask <- terra::rast(highhf_mask_path)
 highhf_mask <- terra::project(x=highhf_mask, y=stack_y, method = "near")
 
 # import subbasin boundaries
 all_subbasins_subset <- terra::vect(all_subbasins_subset_path)
 
 
 ###
 # train and backfill models per subbasin
 ###
 
 # apply model fitting and backfilling over all subbasins
 
 for (i in 1:length(all_subbasins_subset)) {
 
 res <- train_and_backfill_subbasin_s(
     subbasin_index        = i,
     year                  = year,
     stack_y               = stack_y,
     lowhf_mask            = lowhf_mask,
     highhf_mask           = highhf_mask,
     all_subbasins_subset  = all_subbasins_subset, 
     abiotic_vars          = abiotic_vars,
     biotic_vars           = biotic_vars,
     categorical_responses = categorical_responses,
     ia_dir                = ia_dir,
     neworder              = neworder,
     quiet                 = quiet)
 
 } # close loop over subbasins
 
 invisible(res)
 
} # close train_and_backfill_per_year()


# logfile function to track progress
make_logger <- function(logfile) {
  dir.create(dirname(logfile), recursive = TRUE, showWarnings = FALSE)
  function(fmt, ...) {
    line <- sprintf("[%s pid=%d host=%s] %s\n",
                    format(Sys.time(), "%F %T"),
                    Sys.getpid(),
                    Sys.info()[["nodename"]],
                    sprintf(fmt, ...))
    cat(line, file = logfile, append = TRUE)
  }
}


#10. export the necessary variables and functions to the cluster -------------------
print("* exporting objects and functions to cluster *")
clusterExport(cl, c("neworder", "abiotic_vars", "biotic_vars", "categorical_responses",
                    "train_and_backfill_subbasin_s", "all_subbasins_subset_path", 
                    "highhf_mask_path", "lowhf_mask_path", "stack_y_path", 
                    "ia_dir", "train_and_backfill_per_year", "make_logger"))

#11. train models and backfill biotic features for year y -----------------------------
print("* running backfilling in parallel *")

# run backfilling in parallel by subbasin
subs <- 1:length(all_subbasins_subset)  
subs <- 1:2 # for testing
# parLapply runs the backfilling for every i in `subs`
backfill_results <- 
  parLapplyLB(cl, 
              X = subs, 
              fun = function(i) train_and_backfill_subbasin_s(
                subbasin_index = i, 
                year           = 2020,
                stack_y        = terra::rast(stack_y_path),
                lowhf_mask     = terra::rast(lowhf_mask_path),
                highhf_mask    = terra::rast(highhf_mask_path),
                abiotic_vars   = abiotic_vars, 
                biotic_vars    = biotic_vars,
                ia_dir         = ia_dir,
                quiet          = FALSE,
                neworder       = neworder,
                categorical_responses = categorical_responses,
                all_subbasins_subset  = terra::vect(all_subbasins_subset_path)))
  

#12. stop the cluster----
print("* stopping cluster :-)*")
stopCluster(cl)

#13. save backfilled raster for this species x year
print("* saving raster file *")
print(backfill_results)

if(cc){ q() }

