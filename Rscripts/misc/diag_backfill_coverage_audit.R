# ---
# title: Backfill-mosaic coverage audit (Open Limitation #5 diagnosis)
# author: Mannfred Boehm
# ---
# READ-ONLY diagnostic. Locates WHERE in the pipeline the high-HF coalition pixels
# disappear before 12C's complete.cases() gate, and WHICH covariate(s) drive the
# failures. Writes five coverage_stage*.csv to logs/. Does not modify any pipeline
# state. See plan: "Diagnose & Locate the Source of Backfill-Mosaic Pixel Dropout".
#
# Denominator framing: script 05 subsets to 674 subbasins that ALL contain high-HF
# pixels, and 12C's superset (super_idx) is EXACTLY the high-HF footprint pixels
# (any sector > 0 AND CanHF >= 1). So every coalition pixel belongs, by construction,
# to a backfilled subbasin — dropout against that denominator is the anomaly we hunt,
# not the expected "only some subbasins are backfilled" sparsity.
#
# Stages (see plan):
#   0  per-subbasin / per-BCR high-HF denominator
#   1  subbasin-raster production audit (07/08): which subbasins wrote no raster + why
#   2  per-subbasin layer completeness / degeneracy (08 output)
#   3  BCR mosaic completeness (11 output)
#   4  complete.cases NA decomposition (replicates 12C:228-242) — the smoking gun
#
# Run on a Fir login/compute node:
#   Rscript --vanilla Rscripts/misc/diag_backfill_coverage_audit.R
# or sbatch the companion .sh.

suppressPackageStartupMessages({
  library(BAMexploreR); library(gbm); library(terra); library(tidyverse)
})
terra::terraOptions(progress = 0)

# ---- paths / globals (mirror 12B_repredict_all_coalitions.R) ------------------
nm_root <- "/home/mannfred/projects/def-ecknight/NationalModels"   # must not change

cc <- TRUE; local <- FALSE
if (cc)            { ia_dir <- "/home/mannfred/scratch/impact_assessment" }
if (!cc && local)  { ia_dir <- getwd() }
if (!cc && !local) { ia_dir <- file.path("G:/Shared drives/BAM_NationalModels5", "data", "Extras", "sandbox_data", "impactassessment_sandbox") }

year        <- 2020
species_vec <- c("CAWA", "OVEN")
hirsh_dir   <- file.path(ia_dir, "data", "raw_data", "hirshpearson")
sb_dir      <- file.path(ia_dir, "data", "derived_data", "bart_models", year)
mo_dir      <- file.path(ia_dir, "data", "derived_data", "bart_models_mosaics", year)
log_dir     <- file.path(ia_dir, "logs")
out_dir     <- log_dir
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

.terra_tmp <- file.path(Sys.getenv("SLURM_TMPDIR", unset = tempdir()),
                        paste0("terra_diag_", Sys.getenv("SLURM_JOB_ID", unset = "local")))
dir.create(.terra_tmp, recursive = TRUE, showWarnings = FALSE)
terra::terraOptions(tempdir = .terra_tmp)

# the only covariates 12C's complete.cases() actually checks are the BRT model_vars.
load(file.path(ia_dir, "data", "raw_data", "SpeciesPredictionTruncationValues.Rdata"))
bam_boundary         <- terra::vect(file.path(ia_dir, "data", "raw_data", "Regions", "BAM_BCR_NationalModel_Unbuffered.shp"))
all_subbasins_subset <- terra::vect(file.path(ia_dir, "data", "raw_data", "hydrobasins_masked_merged_subset.gpkg"))

bcr_subbasins_ref <- {
  hits <- terra::relate(all_subbasins_subset, bam_boundary, relation = "intersects")
  ij   <- which(hits, arr.ind = TRUE)
  tibble(sub_index = ij[, 1],
         HYBAS_ID  = all_subbasins_subset$first_HYBAS_ID[ij[, 1]],
         bcr_label = paste(bam_boundary$country[ij[, 2]], bam_boundary$subUnit[ij[, 2]], sep = "_"),
         bcr_code  = gsub("_", "", bcr_label))
}

categorical_responses <- c("ABoVE_1km", "NLCD_1km", "MODISLCC_1km", "MODISLCC_5x5", "SCANFI_1km", "VLCE_1km")
predictor_metadata <-
  dplyr::tibble(BAMexploreR::predictor_metadata) |>
  dplyr::filter(version == "v5") |>
  dplyr::select(predictor, definition, predictor_class) |>
  dplyr::mutate(dplyr::across('predictor', stringr::str_replace, 'Year', 'year')) |>
  dplyr::mutate(dplyr::across('predictor', stringr::str_replace, 'Method', 'method'))
actually_biotic_what <- c("Peatland_5x5", "Peatland_1km")
actually_biotic_df   <- tibble::tibble(predictor = actually_biotic_what, predictor_class = c("Wetland", "Wetland"))
abiotic_vars <- predictor_metadata |>
  dplyr::filter(predictor_class %in% c("Annual Climate", "Climate Normals", "Topography", "Wetland", "Disturbance", "Time", "Method")) |>
  dplyr::filter(!(predictor %in% actually_biotic_what))
biotic_continuous_vars <- predictor_metadata |>
  dplyr::filter(!(predictor_class %in% abiotic_vars$predictor_class)) |>
  dplyr::bind_rows(actually_biotic_df) |>
  dplyr::filter(!predictor %in% categorical_responses) |>
  dplyr::pull(predictor)
disturbance_vars <- dplyr::tibble(BAMexploreR::predictor_metadata) |>
  dplyr::filter(version == "v5") |>
  dplyr::select(predictor, definition, predictor_class) |>
  dplyr::filter(predictor_class == "Disturbance")

source(file.path(ia_dir, "Rscripts", "12E_shapley_utils.R"))
sectors   <- canonical_sectors()
n_sectors <- length(sectors)

# =============================================================================
# STAGE 1 (global, species-independent): subbasin-raster production audit
# =============================================================================
# Which of the 674 subbasins wrote a *_backfill.tif, and for those that didn't,
# what error did 07 log? 07_train_and_backfill.R:144 emits
#   "Error in subbasin <i>: <conditionMessage>"
# The 08B_deploy_mbart.R next-in-function bug surfaces as "no loop for break/next".
message(Sys.time(), " | STAGE 1: subbasin raster production audit")
n_sub_total <- 674L
present_vec <- vapply(seq_len(n_sub_total), function(i)
  file.exists(file.path(sb_dir, paste0("subbasin_", i),
                        paste0("subbasin_", i, "_backfill.tif"))), logical(1))

# scrape every file in logs/ for the per-subbasin error lines
err_map <- setNames(rep(NA_character_, n_sub_total), seq_len(n_sub_total))
log_files <- list.files(log_dir, full.names = TRUE, recursive = TRUE)
log_files <- log_files[!grepl("coverage_stage", log_files)]
for (lf in log_files) {
  ln <- tryCatch(readLines(lf, warn = FALSE), error = function(e) character(0))
  hit <- grep("Error in subbasin [0-9]+:", ln, value = TRUE)
  for (h in hit) {
    m <- regmatches(h, regexec("Error in subbasin ([0-9]+): (.*)$", h))[[1]]
    if (length(m) == 3L) {
      si <- as.integer(m[2]); msg <- trimws(m[3])
      if (si >= 1L && si <= n_sub_total) err_map[as.character(si)] <- msg
    }
  }
}
stage1 <- tibble(
  subbasin       = seq_len(n_sub_total),
  raster_present = present_vec,
  error_string   = unname(err_map[as.character(seq_len(n_sub_total))]),
  next_bug       = grepl("no loop for break/next", err_map[as.character(seq_len(n_sub_total))], fixed = TRUE)
)
write_csv(stage1, file.path(out_dir, "coverage_stage1_subbasins.csv"))
message(Sys.time(), " | STAGE 1: ", sum(!present_vec), " of ", n_sub_total,
        " subbasins have NO raster; ", sum(stage1$next_bug, na.rm = TRUE),
        " match the mbart next-in-function bug")

# =============================================================================
# STAGE 2 (global, species-independent): per-subbasin layer completeness/degeneracy
# =============================================================================
# For each present subbasin raster, for each continuous backfill covariate, record
# how many _draw_* layers exist, whether the draws are FLAT (degenerate gbart emits
# identical draws -> cause #3) and whether draw_001 is all-NA at the footprint.
# Flatness is detected from 3 probe draws (001/050/100): identical => flat.
message(Sys.time(), " | STAGE 2: per-subbasin layer completeness / degeneracy")
probe_draws <- c(1L, 50L, 100L)
stage2_rows <- vector("list", 0L)
for (i in which(present_vec)) {
  sf <- file.path(sb_dir, paste0("subbasin_", i), paste0("subbasin_", i, "_backfill.tif"))
  r  <- tryCatch(terra::rast(sf), error = function(e) NULL)
  if (is.null(r)) next
  lyr   <- names(r)
  dcovs <- unique(sub("_draw_[0-9]{3}$", "", lyr[grep("_draw_[0-9]{3}$", lyr)]))
  for (v in dcovs) {
    dl_all <- lyr[grepl(paste0("^", v, "_draw_[0-9]{3}$"), lyr)]
    probe  <- intersect(paste0(v, "_draw_", sprintf("%03d", probe_draws)), dl_all)
    flat <- NA; allNA <- NA
    if (length(probe) >= 1L) {
      vals <- terra::values(r[[probe]], mat = TRUE)
      keep <- stats::complete.cases(vals)          # footprint pixels (non-NA draws)
      allNA <- !any(rowSums(!is.na(vals)) > 0)
      if (any(keep) && ncol(vals) >= 2L) {
        vk   <- vals[keep, , drop = FALSE]
        flat <- max(abs(vk - vk[, 1L])) < 1e-9     # all probe draws identical
      } else if (any(keep)) {
        flat <- NA                                  # only one probe draw present
      }
    }
    stage2_rows[[length(stage2_rows) + 1L]] <- tibble(
      subbasin = i, covariate = v, n_draw_layers = length(dl_all),
      flat = flat, allNA = allNA)
  }
  if (i %% 50L == 0L) message(Sys.time(), "   ...stage2 subbasin ", i)
}
stage2 <- if (length(stage2_rows)) bind_rows(stage2_rows) else
  tibble(subbasin = integer(), covariate = character(),
         n_draw_layers = integer(), flat = logical(), allNA = logical())
write_csv(stage2, file.path(out_dir, "coverage_stage2_layers.csv"))
message(Sys.time(), " | STAGE 2: ", nrow(stage2), " subbasin x covariate rows")

# =============================================================================
# STAGES 0, 3, 4 (per species x BCR): denominator, mosaic completeness, NA decomp
# =============================================================================
stage0_rows <- vector("list", 0L)
stage3_rows <- vector("list", 0L)
stage4_rows <- vector("list", 0L)

for (species in species_vec) {
  rdata_files <- list.files(file.path(nm_root, "output/06_bootstraps", species),
                            pattern = "can.*\\.Rdata$", full.names = TRUE)
  message(Sys.time(), " | ", species, " | ", length(rdata_files), " BCR models")

  for (rdata_path in rdata_files) {
    res <- tryCatch({
      e <- new.env(parent = emptyenv()); load(rdata_path, envir = e)
      b.list <- e$b.list; rm(e)
      bcr_code <- attr(b.list[[1]], "bcr")
      message(Sys.time(), " | ", species, " ", bcr_code, " | auditing")

      sub_ids <- bcr_subbasins_ref |>
        dplyr::filter(bcr_code == !!bcr_code) |>
        dplyr::pull(sub_index) |> unique()
      # NB: inside tryCatch's expr (a function body) so `next`/`break` are illegal
      # (the very bug we are diagnosing) — bail with return(NULL).
      if (length(sub_ids) == 0) { message("   no subbasins — skip"); return(NULL) }

      stack_obs <- terra::rast(file.path(nm_root, "gis/stacks", paste0(bcr_code, "_", year, ".tif")))

      # ---- SUPERSET mask: any sector > 0 AND CanHF >= 1 (identical to 12C) ----
      canHF_r <- terra::project(
        terra::rast(file.path(hirsh_dir, "CanHF_1km_morethan1.tif")), stack_obs, method = "near")
      sec_rasters <- setNames(vector("list", n_sectors), sectors)
      union_mask  <- terra::rast(stack_obs[[1]]); terra::values(union_mask) <- 0L
      for (sec in sectors) {
        sec_r <- terra::project(
          terra::rast(file.path(hirsh_dir, paste0(sec, ".tif"))), stack_obs, method = "near")
        sec_rasters[[sec]] <- sec_r
        union_mask <- terra::ifel(sec_r > 0, 1L, union_mask)
      }
      super_mask <- terra::ifel((union_mask == 1L) & (canHF_r >= 1), 1, NA)
      rm(union_mask); gc()
      n_super_px <- terra::global(super_mask, "notNA")[[1]]
      if (n_super_px == 0) { message("   superset empty — skip"); return(NULL) }

      # zones (subbasin) at superset pixels — identical to 12C:139-148
      subbasin_zone_r <- terra::rasterize(all_subbasins_subset[sub_ids, ], stack_obs[[1]],
                                          field = "first_HYBAS_ID")
      super_idx   <- which(!is.na(terra::values(super_mask, mat = FALSE)))
      super_zones <- terra::values(subbasin_zone_r, mat = FALSE)[super_idx]

      # ---- STAGE 0: per-subbasin high-HF denominator for this BCR ----
      zone_tab  <- table(super_zones)
      hybas_map <- bcr_subbasins_ref |>
        dplyr::filter(sub_index %in% sub_ids) |>
        dplyr::distinct(sub_index, HYBAS_ID)
      stage0_rows[[length(stage0_rows) + 1L]] <- tibble(
        species = species, bcr = bcr_code,
        HYBAS_ID = as.numeric(names(zone_tab)),
        n_highhf_px = as.integer(zone_tab)) |>
        dplyr::left_join(hybas_map, by = "HYBAS_ID")

      # ---- model_vars infrastructure (identical to 12C:122-137) ----
      model_vars_shared  <- b.list[[1]]$var.names
      cat_vars_shared    <- intersect(model_vars_shared, categorical_responses)
      dist_shared        <- intersect(disturbance_vars$predictor, model_vars_shared)
      cat_levels_shared  <- setNames(
        lapply(cat_vars_shared, function(v)
          as.character(b.list[[1]]$var.levels[[match(v, model_vars_shared)]])), cat_vars_shared)

      # ---- load mosaic ----
      bf_mosaic_path <- file.path(mo_dir, paste0(bcr_code, "_backfilled.tif"))
      mosaic_exists  <- file.exists(bf_mosaic_path)
      stack_bf <- if (mosaic_exists) terra::rast(bf_mosaic_path) else NULL
      bf_names <- if (mosaic_exists) names(stack_bf) else character(0)
      draw_layer_names <- bf_names[grep("_draw_[0-9]{3}$", bf_names)]
      draw_covs        <- unique(sub("_draw_[0-9]{3}$", "", draw_layer_names))
      n_draws          <- 100L

      # continuous backfill covariates the BRT model actually consumes
      req_cont <- intersect(setdiff(model_vars_shared, categorical_responses), biotic_continuous_vars)

      # ---- STAGE 3: mosaic completeness over the high-HF superset ----
      for (v in union(req_cont, draw_covs)) {
        in_mosaic <- v %in% draw_covs
        frac_nonNA <- NA_real_
        if (in_mosaic) {
          l1 <- paste0(v, "_draw_001")
          if (l1 %in% bf_names) {
            vv <- terra::values(stack_bf[[l1]], mat = FALSE)[super_idx]
            frac_nonNA <- mean(!is.na(vv))
          }
        } else {
          frac_nonNA <- 0
        }
        stage3_rows[[length(stage3_rows) + 1L]] <- tibble(
          species = species, bcr = bcr_code, covariate = v,
          required_by_model = v %in% req_cont, in_mosaic = in_mosaic,
          frac_highhf_px_nonNA = frac_nonNA)
      }

      # ---- STAGE 4: complete.cases NA decomposition (replicates 12C:177-242) ----
      # Build ONE column per model_var with its FINAL value: observed-origin from
      # the stack, categorical/draw from the mosaic, disturbance = 0. 12C reads the
      # whole stack then overwrites cat/draw/dist; restricting the stack read to the
      # observed-origin columns yields an identical X_rep[, model_vars] but avoids
      # OOM on large BCRs.
      if (mosaic_exists) {
        obs_origin <- setdiff(model_vars_shared, c(cat_vars_shared, draw_covs, dist_shared))
        obs_have   <- intersect(obs_origin, names(stack_obs))

        n_px  <- length(super_idx)
        X_rep <- as.data.frame(matrix(NA_real_, nrow = n_px, ncol = length(model_vars_shared)),
                               check.names = FALSE)
        names(X_rep) <- model_vars_shared
        if (length(obs_have) > 0L) {
          obs_vals <- terra::values(stack_obs[[obs_have]])[super_idx, , drop = FALSE]
          for (v in obs_have) X_rep[[v]] <- obs_vals[, v]
          rm(obs_vals); gc()
        }

        # categorical backfill values (12C:193-203) + factor conversion (12C:232-236)
        for (v in cat_vars_shared) {
          if (!v %in% names(X_rep)) next
          lyr_name <- if (v %in% bf_names) v else paste0(v, "_mean")
          cv <- if (lyr_name %in% bf_names)
            terra::values(stack_bf[[lyr_name]], mat = FALSE)[super_idx] else NA_integer_
          lvls <- cat_levels_shared[[v]]
          X_rep[[v]] <- if (!is.null(lvls) && length(lvls) > 0L)
            factor(as.character(cv), levels = lvls) else cv
        }
        # disturbance -> 0 (12C:231): never NA
        for (v in dist_shared) if (v %in% names(X_rep)) X_rep[[v]] <- 0
        # continuous draws, first scenario (12C:182-190, 237): expm1, non-finite -> NA
        for (v in draw_covs) {
          if (!v %in% names(X_rep)) next
          l1 <- paste0(v, "_draw_001")
          if (!l1 %in% bf_names) { X_rep[[v]] <- NA_real_; next }
          d1 <- expm1(terra::values(stack_bf[[l1]], mat = FALSE)[super_idx])
          d1[!is.finite(d1)] <- NA_real_
          X_rep[[v]] <- pmax(d1, 0)
        }

        Xm     <- X_rep[, model_vars_shared, drop = FALSE]
        na_mat <- is.na(Xm)
        row_na <- rowSums(na_mat)
        complete_frac <- mean(row_na == 0)
        message(Sys.time(), " | ", species, " ", bcr_code,
                " | complete superset pixels: ", sum(row_na == 0), " / ", length(super_idx),
                " (", round(100 * complete_frac, 1), "%)  [match against 12C log]")

        origin_of <- function(v) {
          if (v %in% cat_vars_shared)      "categorical"
          else if (v %in% draw_covs)       "draw"
          else if (v %in% dist_shared)     "disturbance"
          else                             "observed"
        }
        for (v in model_vars_shared) {
          na_v  <- na_mat[, v]
          sole  <- na_v & (row_na == 1L)
          stage4_rows[[length(stage4_rows) + 1L]] <- tibble(
            species = species, bcr = bcr_code, covariate = v,
            origin = origin_of(v),
            na_frac_highhf   = mean(na_v),
            sole_blocker_frac = mean(sole))
        }
        rm(X_rep, Xm, na_mat); gc()
      }

      rm(sec_rasters, canHF_r, super_mask, subbasin_zone_r); if (!is.null(stack_bf)) rm(stack_bf)
      gc(); invisible(NULL)
    }, error = function(e) {
      message(Sys.time(), " | ERROR auditing ", basename(rdata_path), ": ", conditionMessage(e))
      NULL
    })
  }
}

stage0 <- if (length(stage0_rows)) bind_rows(stage0_rows) else tibble()
stage3 <- if (length(stage3_rows)) bind_rows(stage3_rows) else tibble()
stage4 <- if (length(stage4_rows)) bind_rows(stage4_rows) else tibble()
write_csv(stage0, file.path(out_dir, "coverage_stage0_denominator.csv"))
write_csv(stage3, file.path(out_dir, "coverage_stage3_mosaic.csv"))
write_csv(stage4, file.path(out_dir, "coverage_stage4_na_decomp.csv"))

# ---- localizer joins (printed summary) ---------------------------------------
# Stage-1 share: fraction of high-HF pixels that sit in subbasins with NO raster.
if (nrow(stage0)) {
  miss_sub <- stage1$subbasin[!stage1$raster_present]
  loc <- stage0 |>
    dplyr::group_by(species, bcr) |>
    dplyr::summarise(
      n_highhf = sum(n_highhf_px, na.rm = TRUE),
      n_highhf_in_missing = sum(n_highhf_px[sub_index %in% miss_sub], na.rm = TRUE),
      .groups = "drop") |>
    dplyr::mutate(frac_in_missing_subbasins = n_highhf_in_missing / n_highhf)
  write_csv(loc, file.path(out_dir, "coverage_localizer_bcr.csv"))
  message(Sys.time(), " | wrote coverage_localizer_bcr.csv")
}

message(Sys.time(), " | DONE. CSVs in ", out_dir)
