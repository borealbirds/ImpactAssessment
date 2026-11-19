# ---
# title: Impact Assessment: re-predict bird densities for all 255 coalitions, one BCR per task
# author: Mannfred Boehm
# ---

# one array task = one species x BCR (row SLURM_ARRAY_TASK_ID of coalition_task_table(), 12F):
#   -> data/derived_data/density_tables/by_bcr/{species}_{year}_{bcr}.rds
# written even when 12F skips the BCR (result = NULL), so 12H can tell "skipped" from "never
# ran". 12H then merges the BCRs into the tables 14B reads.
# submit: bash 12D_repredict_all_coalitions.sh OSFL GRSP   (chains 12H)

library(gbm)
library(terra)

nm_root <- "/home/mannfred/projects/def-ecknight/NationalModels"
ia_dir  <- "/home/mannfred/scratch/impact_assessment"

# terra's temp files go to the node's local disk, one directory per task
terra_tmp <- file.path(Sys.getenv("SLURM_TMPDIR", unset = tempdir()),
                       paste0("terra_", Sys.getenv("SLURM_JOB_ID"), "_", Sys.getenv("SLURM_ARRAY_TASK_ID")))
dir.create(terra_tmp, recursive = TRUE, showWarnings = FALSE)
terra::terraOptions(tempdir = terra_tmp)

source(file.path(ia_dir, "Rscripts", "12E_shapley_utils.R"))
source(file.path(ia_dir, "Rscripts", "12F_predict_species_all_coalitions.R"))
source(file.path(ia_dir, "Rscripts", "12G_gbm_tree_walk.R"))
Rcpp::sourceCpp(file.path(ia_dir, "Rscripts", "12G_gbm_tree_walk.cpp"))



# -----------------------
# 12F's globals

load(file.path(ia_dir, "data", "raw_data", "SpeciesPredictionTruncationValues.Rdata"))   # q.out

bam_boundary         <- terra::vect(file.path(ia_dir, "data", "raw_data", "Regions", "BAM_BCR_NationalModel_Unbuffered.shp"))
all_subbasins_subset <- terra::vect(file.path(ia_dir, "data", "raw_data", "hydrobasins_masked_merged_subset.gpkg"))

# subbasin x BCR pairs that intersect (a subbasin can be in several BCRs)
ij <- which(terra::relate(all_subbasins_subset, bam_boundary, relation = "intersects"), arr.ind = TRUE)
bcr_subbasins_ref <- data.frame(sub_index = ij[, 1], bcr_code = paste0(bam_boundary$country, bam_boundary$subUnit)[ij[, 2]])

pm <- dplyr::filter(BAMexploreR::predictor_metadata, version == "v5")
categorical_responses  <- c("ABoVE_1km", "NLCD_1km", "MODISLCC_1km", "MODISLCC_5x5", "SCANFI_1km", "VLCE_1km")
footprint_vars         <- pm$predictor[pm$predictor_class == "Disturbance"]
biotic_continuous_vars <- setdiff(c(pm$predictor[pm$predictor_class %in% c("Landcover", "Veg Species", "Veg Structure")],
                                    "Peatland_1km", "Peatland_5x5"),
                                  categorical_responses)



# -----------------------
# this task's species x BCR

species_vec <- strsplit(Sys.getenv("SPECIES"), "[, ]+")[[1]]
year        <- as.integer(Sys.getenv("YEARS", "2020"))   # the .sh submits one year per job
tasks       <- coalition_task_table(species_vec, nm_root)
task_id     <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
if (task_id > nrow(tasks)) {
  message("no task ", task_id, " (the table has ", nrow(tasks), ")")
  quit(save = "no", status = 0)
}
task <- tasks[task_id, ]
message(Sys.time(), " | task ", task_id, " of ", nrow(tasks), ": ", task$species, " ", task$bcr, " ", year)

result <- predict_species_all_coalitions(task$species, year = year, rdata_path = task$rdata_path,
                                         all_subbasins_subset = all_subbasins_subset,
                                         hirsh_dir = file.path(ia_dir, "data", "raw_data", "hirshpearson"),
                                         save_arrays_ids = coalition_array_ids())

# written to a temp file and renamed, so a killed task leaves no truncated file. code_md5 lets
# 12H refuse BCRs made by different versions of the prediction code
out_dir <- file.path(ia_dir, "data", "derived_data", "density_tables", "by_bcr")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_f <- file.path(out_dir, paste0(task$species, "_", year, "_", task$bcr, ".rds"))
code_files <- file.path(ia_dir, "Rscripts", c("12E_shapley_utils.R", "12F_predict_species_all_coalitions.R",
                                              "12G_gbm_tree_walk.R", "12G_gbm_tree_walk.cpp"))
saveRDS(list(species     = task$species,
             year        = year,
             bcr         = task$bcr,
             bcr_order   = task$bcr_order,
             code_md5    = paste(unname(tools::md5sum(code_files)), collapse = ":"),
             job         = paste0(Sys.getenv("SLURM_ARRAY_JOB_ID"), "_", task_id),
             created     = Sys.time(),
             test_n_boot = as.integer(Sys.getenv("TEST_N_BOOT", "0")),
             result      = result),
        paste0(out_f, ".tmp"))
if (!file.rename(paste0(out_f, ".tmp"), out_f)) stop("could not move ", out_f, ".tmp into place")
message(Sys.time(), " | wrote ", out_f, if (is.null(result)) " (BCR skipped by 12F)")
