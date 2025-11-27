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



#8. import helper functions ----------------------------------------
 
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
source(file.path(getwd(), "Rscripts", "may_restart", "08_train_and_backfill_subbasin_s.R"))



#9. export the necessary variables and functions to the cluster -------------------
print("* exporting objects and functions to cluster *")
clusterExport(cl, c("neworder", "abiotic_vars", "biotic_vars", 
                    "train_and_backfill_subbasin_s", "ia_dir", "make_logger"))


#10. train models and backfill biotic features for year y -----------------------------
print("* running backfilling in parallel *")

# run backfilling in parallel by subbasin
subs <- c(301, 240, 57, 491, 100) # for testing
# parLapply runs the backfilling for every i in `subs`
backfill_results <- 
  parLapplyLB(cl, 
              X = subs, 
              fun = function(i, year = 2020) {
                
                # load spatial objects inside of each worker to avoid "external pointer is not valid"
                library(terra)
                # import pre-mosaiced covariate stack for year_y
                stack_y <- terra::rast(file.path(ia_dir, sprintf("covariates_mosaiced_%d.tif", year)))
                
                # ensure categoricals are factors
                categorical_responses = c("ABoVE_1km", "NLCD_1km","MODISLCC_1km", "MODISLCC_5x5","SCANFI_1km","VLCE_1km")
                cats_present <- intersect(categorical_responses, names(stack_y))
                for (cat in cats_present) stack_y[[cat]] <- terra::as.factor(stack_y[[cat]])
                
                # import low hf layer and project to current stack
                lowhf_mask <- terra::rast(file.path(ia_dir, "hirshpearson", "CanHF_1km_lessthan1.tif"))
                lowhf_mask <- terra::project(x=lowhf_mask, y=stack_y, method = "near")
                
                # import high hf layer and project to current stack
                highhf_mask <- terra::rast(file.path(ia_dir, "hirshpearson", "CanHF_1km_morethan1.tif"))
                highhf_mask <- terra::project(x=highhf_mask, y=stack_y, method = "near")
                
                # import subbasin boundaries and project to current stack
                all_subbasins_subset <- terra::vect(file.path(ia_dir, "hydrobasins_masked_merged_subset.gpkg"))
                all_subbasins_subset <- terra::project(x=all_subbasins_subset, y=stack_y)
                
                tryCatch(
                  train_and_backfill_subbasin_s(
                  subbasin_index = i, 
                  year           = 2020,
                  stack_y        = stack_y,
                  lowhf_mask     = lowhf_mask,
                  highhf_mask    = highhf_mask,
                  abiotic_vars   = abiotic_vars, 
                  biotic_vars    = biotic_vars,
                  ia_dir         = ia_dir,
                  quiet          = FALSE,
                  neworder       = neworder,
                  categorical_responses = categorical_responses,
                  all_subbasins_subset  = all_subbasins_subset
                ), # close train_and_backfill_subbasin_s
                
                error = function(e) {
                  message("Error in subbasin ", i, ": ", conditionMessage(e))
                  return(list(
                    subbasin = i,
                    error = conditionMessage(e)
                  ))
                } # close error
                
            ) # close trycatch
    } # close function in parapply
)  # close parapply 

#11. stop the cluster----
print("* stopping cluster :-)*")
stopCluster(cl)

#12. save backfilled raster for this species x year
print("* saving raster file *")
print(backfill_results)

if(cc){ q() }

