# ---
# title: Impact Assessment: merge 12D's per-BCR results into the coalition tables 14B reads
# author: Mannfred Boehm
# ---
# Run once every 12D task has finished:
#   sbatch --dependency=afterany:<12D job id> 12H_merge_bcr_tables.sh
# (afterany, not afterok: a failed 12D task must reach the completeness check below, which
# names the array indices to resubmit, instead of leaving this job pending forever.)
#
# Reads  density_tables/by_bcr/{species}_{year}_{bcr}.rds  (one per 12D task) and writes
#   density_tables/{species}_{year}_coalition_{cid}.rds             (cid 2..256)
#   density_tables/arrays/{species}_{year}_coalition_{cid}_arrays.rds (full coalition + 8 singletons)
#   density_tables/{species}_{year}_shapley_samples.rds             (per-sample subbasin Shapley, 14B)
# BCRs are bound in each species' 06_bootstraps list.files() order - the order the
# one-job-per-species 12D summed them in - so the national arrays' floating-point sums, and
# every table, are bit-identical to that single-job run.
#
# Honours TEST_BCR like 12D, so a smoke test merges just its own BCR(s). Pass the smoke's
# TEST_N_BOOT too: 12H refuses a merge whose TEST_N_BOOT differs from the per-BCR files'.

suppressPackageStartupMessages({ library(tidyverse) })

nm_root <- "/home/mannfred/projects/def-ecknight/NationalModels"   # must not change

cc <- TRUE; local <- FALSE
if (cc)            { ia_dir <- "/home/mannfred/scratch/impact_assessment" }
if (!cc && local)  { ia_dir <- getwd() }
if (!cc && !local) { ia_dir <- file.path("G:/Shared drives/BAM_NationalModels5", "data", "Extras", "sandbox_data", "impactassessment_sandbox") }

source(file.path(ia_dir, "Rscripts", "12E_shapley_utils.R"))
source(file.path(ia_dir, "Rscripts", "12F_predict_species_all_coalitions.R"))

species_vec <- c("CAWA", "OVEN")   # must match 12D_repredict_all_coalitions.R
year <- 2020

dt_dir  <- file.path(ia_dir, "data", "derived_data", "density_tables")
arr_dir <- file.path(dt_dir, "arrays")
dir.create(arr_dir, recursive = TRUE, showWarnings = FALSE)

tasks <- coalition_task_table(species_vec, nm_root)
tasks$file <- file.path(dt_dir, "by_bcr", paste0(tasks$species, "_", year, "_", tasks$bcr, ".rds"))
message(Sys.time(), " | expecting ", nrow(tasks), " per-BCR result(s)")

# completeness: every expected species x BCR must have a result (or a "skipped" record) ----
# The resubmit hint must carry the smoke settings: task indices are positions in the
# TEST_BCR-filtered table, so without TEST_BCR the same index is a different BCR.
test_bcr    <- Sys.getenv("TEST_BCR", "")
test_n_boot <- Sys.getenv("TEST_N_BOOT", "")
test_export <- if (nchar(test_bcr) > 0L)
  paste0(" --export=ALL,TEST_BCR=", test_bcr,
         if (nchar(test_n_boot) > 0L) paste0(",TEST_N_BOOT=", test_n_boot) else "") else ""
missing <- tasks[!file.exists(tasks$file), ]
if (nrow(missing) > 0L)
  stop(nrow(missing), " per-BCR result(s) missing: ",
       paste0(missing$species, " ", missing$bcr, collapse = ", "),
       " - resubmit with: sbatch --array=", paste(missing$task, collapse = ","), test_export,
       " 12D_repredict_all_coalitions.sh (add --mem=64G if its log shows oom_kill or dead workers)")

recs <- lapply(tasks$file, readRDS)

# provenance: all from one code version, and each file is what its name says ----------
prov <- tibble(species   = tasks$species, bcr = tasks$bcr,
               job       = vapply(recs, `[[`, character(1L), "job"),
               created   = format(do.call(c, lapply(recs, `[[`, "created"))),
               n_boot    = vapply(recs, function(r) if (is.null(r$test_n_boot) || r$test_n_boot == 0L) "all" else as.character(r$test_n_boot), character(1L)),
               skipped   = vapply(recs, function(r) is.null(r$result), logical(1L)),
               code_md5  = vapply(recs, `[[`, character(1L), "code_md5"))
message(paste(capture.output(print(select(prov, -code_md5), n = Inf)), collapse = "\n"))
bad_id <- !(vapply(recs, `[[`, character(1L), "species") == tasks$species &
              vapply(recs, `[[`, character(1L), "bcr") == tasks$bcr &
              vapply(recs, `[[`, numeric(1L), "year") == year &
              vapply(recs, function(r) is.null(r$result) || identical(r$result$bcr_code, r$bcr), logical(1L)))
if (any(bad_id))
  stop("per-BCR file(s) do not hold what their name says: ", paste(basename(tasks$file[bad_id]), collapse = ", "))
if (length(unique(prov$code_md5)) > 1L)
  stop("per-BCR results come from ", length(unique(prov$code_md5)), " different versions of ",
       "12E/12F/12G - a stale file from an earlier run would be merged with new ones. Wipe ",
       "density_tables/by_bcr/*.rds and re-run 12D, or re-run the odd tasks out.")
if (length(unique(prov$n_boot)) > 1L)
  stop("per-BCR results mix TEST_N_BOOT settings: ", paste(unique(prov$n_boot), collapse = ", "))
if (!identical(prov$n_boot[1], if (nchar(test_n_boot) > 0L && test_n_boot != "0") test_n_boot else "all"))
  stop("per-BCR results were run with TEST_N_BOOT=", prov$n_boot[1], " but this merge has TEST_N_BOOT=",
       if (nchar(test_n_boot) > 0L) test_n_boot else "(unset)",
       " - pass the same TEST_N_BOOT to 12H as to 12D, so a smoke is never merged as production")
if (any(prov$n_boot != "all"))
  message(Sys.time(), " | NOTE: these are TEST_N_BOOT=", prov$n_boot[1], " smoke results")

# merge per species, in 06_bootstraps order, and write what 14B / 15A read ----------------
for (sp in unique(tasks$species)) {
  rows <- which(tasks$species == sp)
  rows <- rows[order(tasks$bcr_order[rows])]
  res  <- combine_bcr_results(lapply(recs[rows], `[[`, "result"), coalition_array_ids())

  for (cid in names(res$tables_by_cid)) {
    tbl <- res$tables_by_cid[[cid]]
    if (is.null(tbl) || nrow(tbl) == 0) {
      message(Sys.time(), " | WARNING: empty table for ", sp, " coalition ", cid, " — not written")
      next
    }
    saveRDS(tbl, file = file.path(dt_dir, paste0(sp, "_", year, "_coalition_", cid, ".rds")))
  }
  message(Sys.time(), " | ", sp, " | wrote ", length(res$tables_by_cid), " coalition tables from ",
          sum(!prov$skipped[rows]), " BCR(s)")

  # national bootstrap x scenario arrays for 15A (full coalition + 8 singletons)
  for (cid in names(res$arrays_by_cid)) {
    arr <- res$arrays_by_cid[[cid]]
    if (is.null(arr)) {
      message(Sys.time(), " | WARNING: no arrays for ", sp, " coalition ", cid, " — not written")
      next
    }
    saveRDS(arr, file = file.path(arr_dir, paste0(sp, "_", year, "_coalition_", cid, "_arrays.rds")))
  }
  message(Sys.time(), " | ", sp, " | wrote ", length(Filter(Negate(is.null), res$arrays_by_cid)),
          " array files")

  # per-sample subbasin Shapley values for 14B's uncertainty
  if (is.null(res$shapley_samples)) {
    message(Sys.time(), " | WARNING: no Shapley samples for ", sp, " — not written")
  } else {
    saveRDS(res$shapley_samples, file = file.path(dt_dir, paste0(sp, "_", year, "_shapley_samples.rds")))
    message(Sys.time(), " | ", sp, " | wrote Shapley samples: ", dim(res$shapley_samples$phi)[1L],
            " subbasin rows x ", dim(res$shapley_samples$phi)[3L], " samples")
  }
}

message(Sys.time(), " nice.")
