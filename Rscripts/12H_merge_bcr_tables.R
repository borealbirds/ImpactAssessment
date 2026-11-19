# ---
# title: Impact Assessment: merge 12D's per-BCR results into the coalition tables 14B reads
# author: Mannfred Boehm
# ---

# chained after 12D by its .sh (afterany, so a failed 12D task reaches the completeness check
# below, which names the tasks to resubmit), with the same SPECIES and TEST_* settings. writes
#   density_tables/{species}_{year}_coalition_{cid}.rds               (cid 2..256)
#   density_tables/arrays/{species}_{year}_coalition_{cid}_arrays.rds (all sectors, each sector alone; 15A)
#   density_tables/{species}_{year}_shapley_samples.rds               (per-sample subbasin Shapley values; 14B)
# BCRs are bound in each species' 06_bootstraps order, which fixes the floating-point order of
# the national sums

nm_root <- "/home/mannfred/projects/def-ecknight/NationalModels"
ia_dir  <- "/home/mannfred/scratch/impact_assessment"

source(file.path(ia_dir, "Rscripts", "12E_shapley_utils.R"))
source(file.path(ia_dir, "Rscripts", "12F_predict_species_all_coalitions.R"))

species_vec <- strsplit(Sys.getenv("SPECIES"), "[, ]+")[[1]]
year        <- as.integer(Sys.getenv("YEARS", "2020"))
test_bcr    <- Sys.getenv("TEST_BCR", "")
test_n_boot <- Sys.getenv("TEST_N_BOOT", "")
dt_dir      <- file.path(ia_dir, "data", "derived_data", "density_tables")
dir.create(file.path(dt_dir, "arrays"), recursive = TRUE, showWarnings = FALSE)

tasks <- coalition_task_table(species_vec, nm_root)
tasks$file <- file.path(dt_dir, "by_bcr", paste0(tasks$species, "_", year, "_", tasks$bcr, ".rds"))
message(Sys.time(), " | expecting ", nrow(tasks), " per-BCR results for ", year)

# every species x BCR must have a file. task indices are positions in the TEST_BCR-filtered
# table, so the resubmit line carries the smoke settings
missing <- tasks[!file.exists(tasks$file), ]
if (nrow(missing) > 0L)
  stop(nrow(missing), " per-BCR results missing: ", paste(missing$species, missing$bcr, collapse = ", "),
       " - resubmit with: YEARS=", year, " SPECIES=", paste(species_vec, collapse = ","),
       " sbatch --array=", paste(missing$task, collapse = ","),
       if (nchar(test_bcr) > 0L) paste0(" --export=ALL,TEST_BCR=", test_bcr,
                                        if (nchar(test_n_boot) > 0L) paste0(",TEST_N_BOOT=", test_n_boot)),
       " 12D_repredict_all_coalitions.sh (add --mem=64G if its log shows oom_kill or dead workers)")

recs <- lapply(tasks$file, readRDS)

# provenance: each file holds what its name says, all from one code version, one TEST_N_BOOT
n_boot <- vapply(recs, function(r) if (r$test_n_boot == 0L) "all" else as.character(r$test_n_boot), character(1L))
prov <- data.frame(species = tasks$species, bcr = tasks$bcr,
                   job     = vapply(recs, `[[`, character(1L), "job"),
                   created = format(do.call(c, lapply(recs, `[[`, "created"))),
                   n_boot  = n_boot,
                   skipped = vapply(recs, function(r) is.null(r$result), logical(1L)))
message(paste(capture.output(print(prov, row.names = FALSE)), collapse = "\n"))

bad_id <- !(vapply(recs, `[[`, character(1L), "species") == tasks$species &
              vapply(recs, `[[`, character(1L), "bcr") == tasks$bcr &
              vapply(recs, `[[`, numeric(1L), "year") == year &
              vapply(recs, function(r) is.null(r$result) || identical(r$result$bcr_code, r$bcr), logical(1L)))
if (any(bad_id)) stop("per-BCR files that do not hold what their name says: ", paste(basename(tasks$file[bad_id]), collapse = ", "))
if (length(unique(vapply(recs, `[[`, character(1L), "code_md5"))) > 1L)
  stop("per-BCR results come from different versions of 12E/12F/12G: a stale file would be merged with new ",
       "ones. rerun 12D for the species (or the odd tasks out)")
if (length(unique(n_boot)) > 1L || n_boot[1] != (if (nchar(test_n_boot) > 0L && test_n_boot != "0") test_n_boot else "all"))
  stop("TEST_N_BOOT differs between the per-BCR files (", paste(unique(n_boot), collapse = ", "), ") or from this ",
       "merge's (", test_n_boot, "): pass 12H the same TEST_N_BOOT as 12D, so a smoke is never merged as production")

# merge per species, in 06_bootstraps order
for (sp in unique(tasks$species)) {
  rows <- which(tasks$species == sp)
  rows <- rows[order(tasks$bcr_order[rows])]
  res  <- combine_bcr_results(lapply(recs[rows], `[[`, "result"), coalition_array_ids())

  for (cid in names(res$tables_by_cid)) {
    tbl <- res$tables_by_cid[[cid]]
    if (is.null(tbl) || nrow(tbl) == 0) {
      message(Sys.time(), " | WARNING: empty table for ", sp, " coalition ", cid, " - not written")
      next
    }
    saveRDS(tbl, file.path(dt_dir, paste0(sp, "_", year, "_coalition_", cid, ".rds")))
  }
  for (cid in names(res$arrays_by_cid)) {
    if (is.null(res$arrays_by_cid[[cid]])) next
    saveRDS(res$arrays_by_cid[[cid]], file.path(dt_dir, "arrays", paste0(sp, "_", year, "_coalition_", cid, "_arrays.rds")))
  }
  if (is.null(res$shapley_samples)) {
    message(Sys.time(), " | WARNING: no Shapley samples for ", sp, " - not written")
  } else {
    saveRDS(res$shapley_samples, file.path(dt_dir, paste0(sp, "_", year, "_shapley_samples.rds")))
  }
  message(Sys.time(), " | ", sp, " | merged ", sum(!prov$skipped[rows]), " BCR(s)")
}
