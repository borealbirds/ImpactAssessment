# ---
# title: testing the effects of an offset mismatch on density estimates
# author: Mannfred Boehm + claude opus5.5 
# created: November 11, 2026
# ---
#
# for a given `species_vec` and `bcr` this script refits V5's first `N_BOOT` bootstraps of in two treatments:
# with and without the offset misalignment. 
#   v5      - offsets paired by position, as V5 did (should reproduce V5's stored models)
#   matched - offsets matched by survey id 
#
# both treatments keep V5's surveys per bootstrap (bootlist), gbm seed (set.seed(i)), number of trees,
# learning rate and covariate set. note that the hyperparameters come from V5's tuning, which ran with the
# same mismatch, so the matched treatment isolates the bootstrap fit, not a full rerun of V5. 
#
# for each treatment we run several scripts from the IA pipeline: 
# density estimatation, truncation, weighting, range + water masking, conversion to 1km2 cells
# 
# this script runs on the Fir cluster. 
# everything is written under "/home/mannfred/scratch/impact_assessment/offset_test/"
# PR1_offset_refit_observed.sh submits two chained steps:
# (1) STEP=run performs one task per species x treatment, 
# 8 cores (one per bootstrap, where I decided that N_BOOT = 8)
# (2) STEP=summary generates summary tables
# 
#
# outputs (offset_test/summary/):
#   PR1_{bcr}_population.csv  one row per species x treatment: BCR population (mean, 5th/95th
#                             percentiles over bootstraps), q99, plus production's 32-bootstrap row
#   PR1_{bcr}_effect.csv      matched - v5 per species: population, paired by bootstrap, and how
#                             the mean density map moved (cells weighted by weight.tif)
#   PR1_{bcr}_offsets.csv     per species x bootstrap: how wrong V5's offsets are, and whether the
#                             v5 arm reproduced V5's stored trees
#   PR1_{bcr}_raw_check.csv   per species x arm x bootstrap: raw prediction vs V5's own
#                             (v5 arm: should be 0 difference; matched arm: what the bug moved)


suppressPackageStartupMessages({
  library(dplyr); library(terra); library(gbm)
})

if (!nzchar(Sys.getenv("SLURM_JOB_ID")))
  stop("runs on Fir only: cd Rscripts/misc/v5_upstream && bash PR1_offset_refit_observed.sh")

nm_root  <- "/home/mannfred/projects/def-ecknight/NationalModels"
prod_dir <- "/home/mannfred/scratch/impact_assessment"
test_dir <- file.path(prod_dir, "offset_test")
year     <- 2020L
arms     <- c("v5", "matched") #arms is claude's word for "treatment"
categorical_responses <- c("ABoVE_1km", "NLCD_1km", "MODISLCC_1km", "MODISLCC_5x5", "SCANFI_1km", "VLCE_1km")

step        <- Sys.getenv("STEP")

# also ran VATH, RCKI, TOWA for can5
species_vec <- strsplit(Sys.getenv("SPECIES", "OVEN,SWTH,REVI"), "[, ]+")[[1]] 
bcr         <- Sys.getenv("BCR", "can81")
n_boot      <- as.integer(Sys.getenv("N_BOOT", "8"))
n_cores     <- max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")))
if (!step %in% c("run", "summary")) stop("STEP must be run or summary; got '", step, "'")
if (length(species_vec) == 0L || is.na(n_boot) || n_boot < 1L) stop("bad SPECIES or N_BOOT")

tasks <- expand.grid(arm = arms, species = species_vec, stringsAsFactors = FALSE)
if (step == "run") {
  k <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
  if (is.na(k) || k < 1L || k > nrow(tasks))
    stop("array task ", k, " has no species x arm (", nrow(tasks), " tasks)")
  species <- tasks$species[k]
  arm     <- tasks$arm[k]
}
message(Sys.time(), " | STEP=", step, if (step == "run") paste0(" species=", species, " arm=", arm),
        " bcr=", bcr, " n_boot=", n_boot, " cores=", n_cores)

fit_path    <- function(sp, a) file.path(test_dir, "fits", a, paste0(sp, "_", bcr, ".Rdata"))
log_path    <- function(sp, a) file.path(test_dir, "fits", a, paste0(sp, "_", bcr, "_fitlog.rds"))
stack_path  <- function(sp, a) file.path(test_dir, "predictions", a, paste0(sp, "_", bcr, "_", year, ".tif"))
res_path    <- function(sp, a) file.path(test_dir, "results", paste0(sp, "_", a, "_", bcr, "_", year, ".rds"))
stored_path <- function(sp) file.path(nm_root, "output/06_bootstraps", sp, paste0(sp, "_", bcr, ".Rdata"))
weight_path <- function(sp) file.path(prod_dir, "data/derived_data/predictions", sp, bcr, year, "weight.tif")
load_blist  <- function(f) { e <- new.env(parent = emptyenv()); load(f, envir = e); e$b.list } # load a boostrap file (.Rdata) into a new empty environment and returns only `b.list``

# if a job is killed partway through a write, it leaves only a _tmp file behind
save_atomic <- function(write_fun, f) {
  dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
  tmp <- sub("(\\.[A-Za-z]+)$", "_tmp\\1", f)
  write_fun(tmp)
  if (!file.rename(tmp, f)) stop("could not move ", tmp, " into place")
}

# run on the results of a parallel mclapply over bootstraps
# looks for NULL or try-error
check_workers <- function(res, what) {
  died <- vapply(res, is.null, logical(1L))
  if (any(died)) stop(what, ": ", sum(died), " worker(s) died without a result (out of memory? check the log for oom_kill)")
  bad <- vapply(res, inherits, logical(1L), "try-error")
  if (any(bad)) stop(what, ": ", sum(bad), " worker(s) failed: ", as.character(res[[which(bad)[1L]]]))
}

# 12C's weight on this grid (range x not-water x inside data limit x inside the BCR), as 12F reads it
read_weight <- function(sp, grid) {
  w <- terra::rast(weight_path(sp))
  if (!identical(names(w)[1], "weight_v3_touches"))
    stop(weight_path(sp), " is stale (band '", names(w)[1], "'): rebuild it with 12C")
  if (!isTRUE(terra::compareGeom(w, grid, stopOnError = FALSE))) w <- terra::resample(w, grid, method = "near")
  w
}

# birds over the BCR, per layer: density (birds/ha) x 100 per 1 km2 cell, weighted by w
birds <- function(stack, w) unname(terra::global(stack * w, "sum", na.rm = TRUE)[, 1]) * 100

# ------------------------------------------------------------------
# STEP run, part 1: refit bootstraps 1-8 with V5's brt_boot() from 06.Bootstrap.R
# brt_boot() finds surveys with ids in a bootstrap's sample,
# pulls bird counts, covariates, and offsets, then fits a model with gbm()

if (step == "run" && !file.exists(fit_path(species, arm))) {
  suppressPackageStartupMessages(library(Matrix))            # V5's `bird` is a sparse Matrix

  # V5's stored models: their settings, and the reference the v5 arm should reproduce
  stored <- load_blist(stored_path(species))
  
  # check number of bootstraps for current species
  if (length(stored) < n_boot || any(vapply(stored[seq_len(n_boot)], is.null, logical(1L))))
    stop(species, " ", bcr, ": V5's b.list has fewer than ", n_boot, " fitted bootstraps")
  
  # keeps only the first 8 bootstraps and drops the other 24
  stored <- stored[seq_len(n_boot)]; invisible(gc())

  # use this list of vars to choose covariate columns
  vars  <- stored[[1]]$var.names

  # extract hyperparameters from each of the 8 models
  trees <- unique(vapply(stored, function(m) as.numeric(m$n.trees), numeric(1L)))
  lr    <- unique(vapply(stored, function(m) m$shrinkage, numeric(1L)))

  # check settings for every `stored` model
  same  <- vapply(stored, function(m) identical(m$var.names, vars) && m$interaction.depth == 3 &&
                    identical(m$distribution$name, "poisson") && m$bag.fraction == 0.5 &&
                    m$n.minobsinnode == 10, logical(1L))
  if (length(trees) != 1L || length(lr) != 1L || !all(same))
    stop(species, " ", bcr, ": V5's bootstraps do not share one setting (trees ", paste(trees, collapse = "/"),
         ", shrinkage ", paste(lr, collapse = "/"), ") or differ from brt_boot()'s gbm call")
  message(Sys.time(), " | V5 settings: ", trees, " trees, shrinkage ", lr, ", ", length(vars), " variables")

  # get survey, bird count, covariate and offset tables used to rebuild each bootstrap
  data_f <- file.path(prod_dir, "data", "raw_data", "04_NM5.0_data_stratify.Rdata")
  if (!file.exists(data_f))
    stop("V5's training data is missing: Globus G:/Shared drives/BAM_NationalModels5/data/",
         "04_NM5.0_data_stratify.Rdata to ", data_f)
  d <- new.env(); load(data_f, envir = d)
  message(Sys.time(), " | loaded ", data_f)
  bcrlist <- d$bcrlist 
  bootlist <- d$bootlist 
  bird <- d$bird
  cov <- d$cov
  visit <- d$visit
  offsets <- d$offsets
  rm(d); invisible(gc())
  
  # stage data to match V5's conventions
  bcr.i <- bcr
  spp.i <- species
  covsnew.i <- data.frame(var = vars) # one-column data frame of covariate names

  # `inputs` is a list with 8 elements (one per bootstrap)
  # each element is a list with three sub-elements: 
  # $dat is a df of survey data with bird counts, year, and var.names as columns
  # $off is the aligned or misaligned offset for a corresponding survey in $dat
  # $log stores diagnostic info
  inputs <- lapply(seq_len(n_boot), function(i) {

    # much of the following is verbatim from brt_boot() steps 6-11 (06.Bootstrap.R )
    visit.i <- bcrlist[bcrlist[,bcr.i] > 0, c("id", bcr.i)] |>
      dplyr::filter(id %in% bootlist[[i + 2]])

    bird.i <- bird[as.character(visit.i$id), spp.i]
    cov.i <- cov[cov$id %in% visit.i$id, colnames(cov) %in% covsnew.i$var]
    year.i <- visit[visit$id %in% visit.i$id, "year"]
    dat.i <- cbind(bird.i, year.i, cov.i) |>
      dplyr::rename(count = bird.i)

    # returns offsets in offsets' own row order, not in visit.i's order (wrong)
    off.i <- offsets[offsets$id %in% visit.i$id, spp.i] 
      
    # fix alignment: give each survey's its `own` offset by matching `id`s
    off_id <- offsets$id[offsets$id %in% visit.i$id]
    own    <- offsets[[spp.i]][match(visit.i$id, offsets$id)] # looks up each offset up by survey id instead
    
    if (!identical(cov$id[cov$id %in% visit.i$id], visit.i$id) ||
        !identical(visit$id[visit$id %in% visit.i$id], visit.i$id))
      stop("bootstrap ", i, ": cov or visit is not in bcrlist's order, so V5's cbind() is wrong too")
    if (!identical(names(dat.i)[-1], vars))
      stop("bootstrap ", i, ": training columns differ from V5's var.names")
    if (!setequal(off_id, visit.i$id) || anyNA(own) || length(own) != nrow(dat.i))
      stop("bootstrap ", i, ": offsets do not cover the bootstrap's surveys")
    other <- off_id != visit.i$id
    
    # both arms run in parallel because the job array has one task per species × arm
    list(dat = dat.i,
         off = if (arm == "v5") off.i else own, #
         log = data.frame(species = spp.i, bcr = bcr.i, boot = i, n_rows = nrow(dat.i),
                          n_detected = sum(dat.i$count > 0),
                          rows_with_other_surveys_offset = sum(other),
                          pct_rows_with_other_surveys_offset = round(100 * mean(other), 2),
                          cor_v5_offset_vs_own = round(cor(off.i, own), 4),
                          mean_abs_offset_error = round(mean(abs(off.i - own)), 4),
                          sd_own_offset = round(sd(own), 4)))
  })
  rm(bcrlist, bootlist, bird, cov, visit, offsets); invisible(gc())
  print(do.call(rbind, lapply(inputs, `[[`, "log")))
  
  # fit models using `gbm`
  # brt_boot() steps 13-14 (verbatim fit call and attributes from 06.Bootstrap.R)
  fit_one <- function(dat.i, off.i, i) {
    t0 <- proc.time()
    set.seed(i)
    b.i <- gbm::gbm(dat.i$count ~ . + offset(off.i),
                    data = dat.i[, -1],
                    n.trees = trees,
                    interaction.depth = 3,
                    shrinkage = lr,
                    distribution = "poisson",
                    keep.data = FALSE,
                    n.cores = 1)
    attr(b.i, "bcr")    <- bcr
    attr(b.i, "spp")    <- species
    attr(b.i, "time")   <- (proc.time() - t0)[3]
    attr(b.i, "n")      <- nrow(dat.i)
    attr(b.i, "ncount") <- sum(dat.i$count > 0)
    attr(b.i, "offset_pairing") <- arm
    b.i
  }
  message(Sys.time(), " | fitting ", n_boot, " bootstraps on ", n_cores, " core(s)")
  b.list <- parallel::mclapply(seq_len(n_boot), function(i) fit_one(inputs[[i]]$dat, inputs[[i]]$off, i),
                               mc.cores = n_cores, mc.preschedule = FALSE)
  check_workers(b.list, paste(species, arm, "fit"))

  fitlog <- do.call(rbind, lapply(seq_len(n_boot), function(i) cbind(
    arm = arm, inputs[[i]]$log,
    trees_identical_to_v5_stored = identical(b.list[[i]]$trees, stored[[i]]$trees) &&
      identical(b.list[[i]]$c.split, stored[[i]]$c.split),
    initF = b.list[[i]]$initF, initF_v5_stored = stored[[i]]$initF,
    fit_minutes = round(unname(attr(b.list[[i]], "time")) / 60, 1))))
  
  print(fitlog[, c("arm", "boot", "trees_identical_to_v5_stored", "initF", "initF_v5_stored", "fit_minutes")])
  
  if (arm == "v5" && !all(fitlog$trees_identical_to_v5_stored))
    message(Sys.time(), " | NOTE: the v5 arm did not reproduce V5's stored trees for bootstrap(s) ",
            paste(fitlog$boot[!fitlog$trees_identical_to_v5_stored], collapse = ", "),
            "; raw_check reports how far its predictions are from V5's")
  
  # save_atomic writes to a temporary filename and 
  # renames it to the final one only when the write has finished
  save_atomic(function(f) save(b.list, file = f), fit_path(species, arm))
  save_atomic(function(f) saveRDS(fitlog, f), log_path(species, arm))
  message(Sys.time(), " | wrote ", fit_path(species, arm))
  rm(b.list, stored, inputs); invisible(gc())
}


# ------------------------------------------------------------------
# STEP run, part 2: density predictions over the BCR, 
# truncated and weighted as V5 publishes them

if (step == "run" && !file.exists(res_path(species, arm))) {
  b.list <- load_blist(fit_path(species, arm))
  if (length(b.list) != n_boot || !identical(attr(b.list[[1]], "bcr"), bcr) ||
      !identical(attr(b.list[[1]], "offset_pairing"), arm))
    stop(fit_path(species, arm), " does not hold ", n_boot, " ", arm, " fits for ", bcr)

  .terra_tmp <- file.path(Sys.getenv("SLURM_TMPDIR", unset = tempdir()),
                          paste0("terra_", Sys.getenv("SLURM_JOB_ID"), "_", Sys.getenv("SLURM_ARRAY_TASK_ID")))
  dir.create(.terra_tmp, recursive = TRUE, showWarnings = FALSE)
  terra::terraOptions(tempdir = .terra_tmp)

  categorical_responses <- c("ABoVE_1km", "NLCD_1km", "MODISLCC_1km", "MODISLCC_5x5", "SCANFI_1km", "VLCE_1km")  # as 12D
  load(file.path(prod_dir, "data", "raw_data", "SpeciesPredictionTruncationValues.Rdata"))   # q.out
  source(file.path(prod_dir, "Rscripts", "12B_v5_truncate.R"))  # v5_truncate()
  source(file.path(prod_dir, "Rscripts", "12F_predict_species_all_coalitions.R"))  # as_model_factors()
  source(file.path(prod_dir, "Rscripts", "12G_gbm_tree_walk.R"))
  
  # custom C++ version of gbm::predict that is much faster
  Rcpp::sourceCpp(file.path(prod_dir, "Rscripts", "12G_gbm_tree_walk.cpp"))
  densmax <- q.out$densmax[q.out$spp == species]
  if (length(densmax) != 1L || !is.finite(densmax)) stop("no usable densmax for ", species)

  stack_obs <- terra::rast(file.path(nm_root, "gis/stacks", paste0(bcr, "_", year, ".tif")))
  
  # import V5's predictions to validate our training/prediction scheme 
  v5_raw    <- terra::rast(file.path(nm_root, "output/07_predictions", species,
                                     paste0(species, "_", bcr, "_", year, ".tif")))
 
  # ensure the two rasters have the same extent, nrows, ncols, and CRS
  if (!isTRUE(terra::compareGeom(v5_raw, stack_obs, stopOnError = FALSE)) || terra::nlyr(v5_raw) < n_boot)
    stop("V5's raw prediction is not on the ", bcr, " stack grid or has fewer than ", n_boot, " layers")
  
  # V5's original predictions, cells, and covariate rows
  E     <- terra::values(v5_raw[[seq_len(n_boot)]], mat = TRUE)
  cells <- which(rowSums(!is.na(E)) > 0L)
  E     <- E[cells, , drop = FALSE]
  vars  <- unique(unlist(lapply(b.list, `[[`, "var.names")))
  if (length(setdiff(vars, names(stack_obs))))
    stop("model variables missing from the stack: ", paste(setdiff(vars, names(stack_obs)), collapse = ", "))
  X <- terra::values(stack_obs[[vars]], mat = TRUE)[cells, , drop = FALSE]
  X[is.nan(X)] <- NA_real_                       # gbm sends only NA_real_ down its missing branch
  X <- as.data.frame(X, check.names = FALSE)
  cat_vars <- intersect(b.list[[1]]$var.names, categorical_responses)  

  #  chunking keeps each worker's peak at roughly 250,000 rows
  chunks   <- split(seq_along(cells), ceiling(seq_along(cells) / 250000L))
  message(Sys.time(), " | predicting ", length(cells), " cells x ", n_boot, " bootstraps")
  
  raw <- parallel::mclapply(seq_len(n_boot), function(i) {
    m <- b.list[[i]]
    p <- rep(NA_real_, length(cells))
    for (j in seq_along(chunks)) {
      r <- chunks[[j]][!is.na(E[chunks[[j]], i])]
      if (length(r) == 0L) next
      Xf <- as_model_factors(X[r, , drop = FALSE], m, cat_vars)
      
      # gbm_* functions are custom made and meant to mirror gbm::predict
      # gbm_design builds the design matrix that predict.gbm would build from the data frame Xf for model m
      D  <- gbm_design(m, Xf)
      if (j == 1L) gbm_check_fast(m, Xf, D, label = paste(species, bcr, arm, "bootstrap", i))
      
      # replaces gbm::predict with a C++ version that gives bit-identical results and runs 4x faster
      p[r] <- gbm_predict_design(m, D)
      rm(Xf, D); invisible(gc(full = FALSE))
    }
    p
  }, mc.cores = n_cores, mc.preschedule = FALSE)
  check_workers(raw, paste(species, arm, "observed prediction"))
  P <- do.call(cbind, raw); rm(raw, X, b.list); invisible(gc())

  # raw_check is a table that checks this arm's re-predictions against V5's own original predictions,
  # with one row per bootstrap. its job is to confirm that our prediction path reproduces V5 faithfully
  raw_check <- do.call(rbind, lapply(seq_len(n_boot), function(i) {
    ok <- !is.na(E[, i])
    rd <- abs(P[ok, i] - E[ok, i]) / pmax(abs(E[ok, i]), 1e-30)
    data.frame(species = species, arm = arm, boot = i, cells = sum(ok),
               na_where_v5_predicted = sum(is.na(P[ok, i])),
               max_rel_diff_vs_v5 = signif(max(rd, na.rm = TRUE), 3),
               pct_cells_off_v5 = round(100 * mean(rd > 1e-5, na.rm = TRUE), 3),
               raw_sum_ratio_vs_v5 = round(sum(P[ok, i], na.rm = TRUE) / sum(E[ok, i]), 5))
  }))
  print(raw_check)

  # V5's 10.Truncate.R steps 5-6, as 12A (IA pipeline) applies them: densmax, then q99 from this arm's own mean
  r_raw <- terra::rast(v5_raw, nlyrs = n_boot)
  vals  <- matrix(NA_real_, terra::ncell(r_raw), n_boot)
  vals[cells, ] <- P
  terra::values(r_raw) <- vals
  rm(vals, P, E, stack_obs); invisible(gc())
  tr <- v5_truncate(r_raw, spp = species, bcr = bcr, densmax = densmax,
                    q99 = NULL, project_to = NULL, apply_masks = FALSE)
  if (is.null(tr)) stop("all-NA bootstrap mean, no q99")
  save_atomic(function(f) terra::writeRaster(tr$stack, f, overwrite = TRUE,
                                             wopt = list(gdal = c("INTERLEAVE=BAND", "COMPRESS=DEFLATE", "BIGTIFF=YES"))),
              stack_path(species, arm))

  w   <- read_weight(species, v5_raw)
  pop <- birds(terra::rast(stack_path(species, arm)), w)
  message(Sys.time(), " | ", species, " ", bcr, " ", arm, ": ", round(mean(pop)), " birds (bootstraps ",
          paste(round(pop), collapse = ", "), "); q99 ", signif(tr$q99, 6))

  # production's observed side (V5's stored models, 32 bootstraps, q99 from all 32), same sum
  pop_prod <- NULL; q99_prod <- NA_real_
  if (arm == "v5") {
    pf <- file.path(prod_dir, "data/derived_data/predictions", species, bcr, year, "observed_bootstraps.tif")
    if (file.exists(pf)) pop_prod <- birds(terra::rast(pf), w)
    q99_prod <- tryCatch(readRDS(file.path(prod_dir, "data/derived_data/predictions", species,
                                           "truncation_params.rds"))[[paste(bcr, year, sep = "_")]]$q99,
                         error = function(e) NA_real_)
  }
  save_atomic(function(f) saveRDS(list(species = species, arm = arm, bcr = bcr, year = year, n_boot = n_boot,
                                       densmax = densmax, q99 = tr$q99, q99_production = q99_prod,
                                       population = pop, population_production = pop_prod,
                                       raw_check = raw_check), f),
              res_path(species, arm))
  message(Sys.time(), " | wrote ", res_path(species, arm))
}

# ------------------------------------------------------------------
# STEP summary: population per arm, and matched - v5 paired by bootstrap
if (step == "summary") {
  q05 <- function(x) unname(quantile(x, 0.05, na.rm = TRUE))
  q95 <- function(x) unname(quantile(x, 0.95, na.rm = TRUE))
  out_dir <- file.path(test_dir, "summary")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  ok   <- file.exists(res_path(tasks$species, tasks$arm))
  if (any(!ok)) message(Sys.time(), " | MISSING results (see offset_test/logs/): ",
                        paste(tasks$species[!ok], tasks$arm[!ok], collapse = ", "))
  have <- tasks[ok, ]

  pop_rows <- list(); eff_rows <- list(); raw_rows <- list()
  for (sp in unique(have$species)) {
    X <- setNames(lapply(have$arm[have$species == sp], function(a) readRDS(res_path(sp, a))),
                  have$arm[have$species == sp])
    for (a in names(X)) {
      p <- X[[a]]$population
      pop_rows[[paste(sp, a)]] <- data.frame(
        species = sp, bcr = bcr, arm = a, n_boot = length(p), q99 = signif(X[[a]]$q99, 6),
        population = round(mean(p)), population_q05 = round(q05(p)), population_q95 = round(q95(p)),
        population_min = round(min(p)), population_max = round(max(p)))
      raw_rows[[paste(sp, a)]] <- X[[a]]$raw_check
    }
    pp <- X$v5$population_production
    if (!is.null(pp)) {
      pop_rows[[paste(sp, "production")]] <- data.frame(
        species = sp, bcr = bcr, arm = "production (V5's stored models)", n_boot = length(pp),
        q99 = signif(X$v5$q99_production, 6),
        population = round(mean(pp)), population_q05 = round(q05(pp)), population_q95 = round(q95(pp)),
        population_min = round(min(pp)), population_max = round(max(pp)),
        population_first_n_boot = round(mean(pp[seq_len(min(n_boot, length(pp)))])))
    }
    if (all(arms %in% names(X))) {
      a <- X$v5$population; b <- X$matched$population
      if (length(a) != length(b)) stop(sp, ": the arms have different bootstrap counts")
      d <- b - a                                           # paired: same surveys and seed
      # mean density maps (over bootstraps), on the cells the weight keeps
      w  <- read_weight(sp, terra::rast(stack_path(sp, "v5")))
      wv <- terra::values(w, mat = FALSE)
      mv <- terra::values(terra::app(terra::rast(stack_path(sp, "v5")), mean, na.rm = TRUE), mat = FALSE)
      mm <- terra::values(terra::app(terra::rast(stack_path(sp, "matched")), mean, na.rm = TRUE), mat = FALSE)
      k  <- which(is.finite(wv) & wv > 0 & is.finite(mv) & is.finite(mm))
      ratio <- mm[k] / mv[k]
      eff_rows[[sp]] <- data.frame(
        species = sp, bcr = bcr,
        population_v5 = round(mean(a)), population_matched = round(mean(b)),
        difference = round(mean(d)), difference_q05 = round(q05(d)), difference_q95 = round(q95(d)),
        difference_min = round(min(d)), difference_max = round(max(d)),
        difference_pct_of_v5 = round(100 * mean(d) / mean(a), 2),
        bootstraps_higher_when_matched = sum(d > 0),
        q99_v5 = signif(X$v5$q99, 6), q99_matched = signif(X$matched$q99, 6),
        map_cells = length(k),
        map_cor = round(cor(mv[k], mm[k]), 4),
        map_median_ratio = round(median(ratio[is.finite(ratio)]), 4),
        map_pct_cells_up_10pct = round(100 * mean(ratio > 1.1, na.rm = TRUE), 1),
        map_pct_cells_down_10pct = round(100 * mean(ratio < 1 / 1.1, na.rm = TRUE), 1))
    }
  }
  pop_tab <- dplyr::bind_rows(pop_rows)
  eff_tab <- dplyr::bind_rows(eff_rows)
  raw_tab <- dplyr::bind_rows(raw_rows)
  off_tab <- dplyr::bind_rows(lapply(seq_len(nrow(tasks)), function(t) {
    f <- log_path(tasks$species[t], tasks$arm[t]); if (file.exists(f)) readRDS(f)
  }))
  write.csv(pop_tab, file.path(out_dir, paste0("PR1_", bcr, "_population.csv")), row.names = FALSE)
  write.csv(eff_tab, file.path(out_dir, paste0("PR1_", bcr, "_effect.csv")), row.names = FALSE)
  write.csv(raw_tab, file.path(out_dir, paste0("PR1_", bcr, "_raw_check.csv")), row.names = FALSE)
  write.csv(off_tab, file.path(out_dir, paste0("PR1_", bcr, "_offsets.csv")), row.names = FALSE)

  options(width = 220)
  cat("\n== offsets in V5's training rows (v5 arm's fit log) ==\n")
  if (nrow(off_tab))
    print(off_tab[off_tab$arm == "v5", c("species", "boot", "n_rows", "pct_rows_with_other_surveys_offset",
                                         "cor_v5_offset_vs_own", "mean_abs_offset_error",
                                         "trees_identical_to_v5_stored", "fit_minutes")], row.names = FALSE)
  cat("\n== raw prediction vs V5's own, per arm (over bootstraps) ==\n")
  if (nrow(raw_tab))
    print(as.data.frame(raw_tab |> group_by(species, arm) |>
      summarise(max_rel_diff = max(max_rel_diff_vs_v5), pct_cells_off = mean(pct_cells_off_v5),
                raw_sum_ratio = round(mean(raw_sum_ratio_vs_v5), 4), .groups = "drop")), row.names = FALSE)
  cat("\n== BCR population per arm (birds) ==\n")
  print(pop_tab, row.names = FALSE)
  cat("\n== matched - v5 ==\n")
  print(eff_tab, row.names = FALSE)
  message(Sys.time(), " | wrote ", out_dir)
}
