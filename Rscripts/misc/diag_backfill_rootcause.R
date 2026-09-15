# ---
# title: Backfill ROOT-CAUSE audit (Open Limitation #5, stage 2)
# author: Mannfred Boehm
# ---
# READ-ONLY diagnostic. The first audit (diag_backfill_coverage_audit.R) LOCATED the
# 12C complete.cases dropout to a "block coverage gap": within subbasins that DO have
# backfill rasters, all draw + categorical layers go NA together on a BCR-specific
# fraction of high-HF pixels. This script measures the cause on the INPUT side, on the
# native 07/08 grid (covariates_mosaiced_{year}.tif), by replicating 08A's predictor
# construction (08A:29-78) WITHOUT running BART, and decomposing where the NA originates.
#
# Mechanism under test (confirmed by code reading, see plan):
#   * 08A writes a value at EVERY high-HF pixel (08A:58,230) — the gap is not unwritten px.
#   * 08A:151 keeps a df_backfill_bart row for every high-HF px; rows with any NA abiotic
#     predictor return NaN from BART (08B_gbart:77,95,101).
#   * 08A:124-132: each later continuous covariate uses earlier (already-backfilled) ones
#     as predictors -> a NaN cascades down the whole hierarchy (uniform block-missingness).
#   * 08A:101-121 / 08B_*:const branches emit a constant (NA if no training obs, 08A:105)
#     layer -> covariate-specific cap (e.g. near-absent SCANFI species).
#
# Outputs (to logs/):
#   rootcause_abiotic_na_by_cov.csv   per subbasin x abiotic covariate: na_frac at high-HF px
#   rootcause_cascade_by_subbasin.csv per subbasin: abiotic-NA seed + cascade first/last frac
#   rootcause_rare_response.csv       per biotic covariate: n subbasins hitting const / const=NA
#   rootcause_bcr_rollup.csv          per BCR: backfillable frac vs first audit's Stage-3 ref
#
# Run: Rscript --vanilla <abs path>/Rscripts/misc/diag_backfill_rootcause.R
# Smoke: DIAG_SUBBASINS=1,2,3 Rscript --vanilla ...  (restrict to those subbasin indices)

suppressPackageStartupMessages({
  library(BAMexploreR); library(terra); library(tidyverse)
})
terra::terraOptions(progress = 0)

# ---- execution context (mirror 07_train_and_backfill.R) -----------------------
cc    <- as.logical(Sys.getenv("DIAG_CC", unset = "TRUE"))
local <- as.logical(Sys.getenv("DIAG_LOCAL", unset = "FALSE"))
if (cc)            { ia_dir <- "/home/mannfred/scratch/impact_assessment" }
if (!cc && local)  { ia_dir <- getwd() }
if (!cc && !local) { ia_dir <- file.path("G:/Shared drives/BAM_NationalModels5", "data", "Extras", "sandbox_data", "impactassessment_sandbox") }
message("ia_dir = ", ia_dir)

year     <- 2020
log_dir  <- file.path(ia_dir, "logs"); dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
out_dir  <- log_dir

.terra_tmp <- file.path(Sys.getenv("SLURM_TMPDIR", unset = tempdir()),
                        paste0("terra_rootcause_", Sys.getenv("SLURM_JOB_ID", unset = "local")))
dir.create(.terra_tmp, recursive = TRUE, showWarnings = FALSE)
terra::terraOptions(tempdir = .terra_tmp)

# ---- covariate metadata (verbatim from 07:38-77) ------------------------------
predictor_metadata <-
  dplyr::tibble(BAMexploreR::predictor_metadata) |>
  dplyr::filter(version == "v5") |>
  dplyr::select(predictor, definition, predictor_class) |>
  dplyr::mutate(dplyr::across('predictor', stringr::str_replace, 'Year', 'year')) |>
  dplyr::mutate(dplyr::across('predictor', stringr::str_replace, 'Method', 'method'))

soil_covs <- tibble::tibble(predictor = c("cec_0-5cm_mean_1000", "cec_100-200cm_mean_1000",
                                  "cec_15-30cm_mean_1000", "cec_30-60cm_mean_1000",
                                  "cec_5-15cm_mean_1000", "cec_60-100cm_mean_1000",
                                  "soc_0-5cm_mean_1000",  "soc_100-200cm_mean_1000",
                                  "soc_15-30cm_mean_1000", "soc_30-60cm_mean_1000",
                                  "soc_5-15cm_mean_1000", "soc_60-100cm_mean_1000"),
                    predictor_class = rep("Soil Properties", 12))

actually_biotic_what <- c("Peatland_5x5", "Peatland_1km")
actually_biotic_df   <- tibble::tibble(predictor = actually_biotic_what, predictor_class = c("Wetland", "Wetland"))

abiotic_vars <-
  predictor_metadata |>
  dplyr::filter(predictor_class %in% c("Annual Climate", "Climate Normals", "Topography", "Wetland", "Disturbance")) |>
  tibble::add_row(predictor = "CAfire", predictor_class = "Time Since Disturbance") |>
  dplyr::filter(!(predictor %in% actually_biotic_what)) |>
  dplyr::bind_rows(soil_covs)

biotic_vars <-
  predictor_metadata |>
  dplyr::filter(!(predictor_class %in% c(abiotic_vars$predictor_class, "Time", "Method"))) |>
  dplyr::bind_rows(actually_biotic_df)

neworder <- readRDS(file = file.path(ia_dir, "data", "raw_data", "biotic_variable_hierarchy.rds"))
biotic_vars <- biotic_vars[match(neworder, biotic_vars$predictor), ]

categorical_responses <- c("ABoVE_1km", "NLCD_1km", "MODISLCC_1km", "MODISLCC_5x5", "SCANFI_1km", "VLCE_1km")

# class lookup for a given covariate name (abiotic only)
abiotic_class_of <- setNames(abiotic_vars$predictor_class, abiotic_vars$predictor)

# ---- spatial inputs (mirror 07:108-124) ---------------------------------------
stack_y <- terra::rast(file.path(ia_dir, "data", "raw_data", "covariates_mosaiced", sprintf("covariates_mosaiced_%d.tif", year)))

lowhf_mask  <- terra::rast(file.path(ia_dir, "data", "raw_data", "hirshpearson", "CanHF_1km_lessthan1.tif"))
lowhf_mask  <- terra::project(x = lowhf_mask, y = stack_y, method = "near")
highhf_mask <- terra::rast(file.path(ia_dir, "data", "raw_data", "hirshpearson", "CanHF_1km_morethan1.tif"))
highhf_mask <- terra::project(x = highhf_mask, y = stack_y, method = "near")

all_subbasins_subset <- terra::vect(file.path(ia_dir, "data", "raw_data", "hydrobasins_masked_merged_subset.gpkg"))
all_subbasins_subset <- terra::project(x = all_subbasins_subset, y = stack_y)
n_sub_total <- length(all_subbasins_subset)

# ---- subbasin -> BCR map (for Stage D rollup; mirror first audit) -------------
bcr_subbasins_ref <- tryCatch({
  bam_boundary <- terra::vect(file.path(ia_dir, "data", "raw_data", "Regions", "BAM_BCR_NationalModel_Unbuffered.shp"))
  bam_boundary <- terra::project(bam_boundary, all_subbasins_subset)
  hits <- terra::relate(all_subbasins_subset, bam_boundary, relation = "intersects")
  ij   <- which(hits, arr.ind = TRUE)
  tibble(sub_index = ij[, 1],
         bcr_code  = gsub("_", "", paste(bam_boundary$country[ij[, 2]], bam_boundary$subUnit[ij[, 2]], sep = "_"))) |>
    dplyr::distinct()
}, error = function(e) { message("bcr ref failed: ", conditionMessage(e)); tibble(sub_index = integer(), bcr_code = character()) })

# ---- subbasin selection (smoke mode) ------------------------------------------
sel_env <- Sys.getenv("DIAG_SUBBASINS", unset = "")
sub_indices <- if (nzchar(sel_env)) as.integer(strsplit(sel_env, ",")[[1]]) else seq_len(n_sub_total)
message("auditing ", length(sub_indices), " of ", n_sub_total, " subbasins")

# =============================================================================
# per-subbasin pass: replicate 08A:29-78 (no BART) and decompose NA
# =============================================================================
abioticA <- vector("list", 0L)   # Stage A rows
cascadeB <- vector("list", 0L)   # Stage B rows
rareC    <- vector("list", 0L)   # Stage C rows (per subbasin x biotic cov)

for (i in sub_indices) {
  ok <- tryCatch({
    subbasin_s <- all_subbasins_subset[i]
    cov_s <- stack_y |> terra::crop(y = subbasin_s) |> terra::mask(mask = subbasin_s)

    lowhf_mask_s  <- terra::resample(lowhf_mask,  cov_s, method = "near") |> terra::mask(mask = subbasin_s)
    highhf_mask_s <- terra::resample(highhf_mask, cov_s, method = "near") |> terra::mask(mask = subbasin_s)

    cov_train_s <- terra::mask(x = cov_s, mask = lowhf_mask_s)
    df_train    <- terra::as.data.frame(cov_train_s, xy = TRUE, na.rm = FALSE)

    df_full      <- terra::as.data.frame(cov_s, xy = TRUE, na.rm = FALSE, cells = TRUE)
    backfill_idx <- which(terra::values(highhf_mask_s) == 1)
    df_backfill  <- df_full[backfill_idx, , drop = FALSE]
    n_bf <- nrow(df_backfill)
    if (n_bf == 0L) return(TRUE)

    abiotic_cols <- intersect(names(df_train), abiotic_vars$predictor)
    biotic_cols  <- intersect(names(df_train), biotic_vars$predictor)
    biotic_cols  <- na.omit(biotic_cols[match(neworder, biotic_cols)])
    biotic_cols_cont <- setdiff(biotic_cols, categorical_responses)

    # ---- STAGE A: abiotic NA at high-HF (backfill) pixels ----
    for (a in abiotic_cols) {
      col <- df_backfill[[a]]
      abioticA[[length(abioticA) + 1L]] <- tibble(
        subbasin = i, covariate = a,
        predictor_class = unname(abiotic_class_of[a]),
        n_backfill_px = n_bf,
        na_frac_backfill = mean(is.na(col)),
        all_na = all(is.na(col)))
    }

    # ---- STAGE B: abiotic seed + hierarchy cascade (matches 08A col drops) ----
    # 08A:155 drops abiotic predictors that are ALL-NA over the backfill set; only
    # partially-present columns constrain coverage.
    abiotic_present <- abiotic_cols[vapply(abiotic_cols,
                          function(a) any(!is.na(df_backfill[[a]])), logical(1))]
    abiotic_ok <- if (length(abiotic_present))
      stats::complete.cases(df_backfill[, abiotic_present, drop = FALSE]) else rep(TRUE, n_bf)

    # dominant NA-contributing abiotic class among present columns
    dom_class <- NA_character_; dom_frac <- 0
    if (length(abiotic_present)) {
      cls <- unique(unname(abiotic_class_of[abiotic_present]))
      for (cl in cls) {
        cc_cols <- abiotic_present[abiotic_class_of[abiotic_present] == cl]
        cl_na <- rowSums(is.na(df_backfill[, cc_cols, drop = FALSE])) > 0
        f <- mean(cl_na)
        if (f >= dom_frac) { dom_frac <- f; dom_class <- cl }
      }
    }

    # cascade: bf_ok[[b]] = abiotic_ok & all preceding continuous biotic backfillable
    frac_first <- NA_real_; frac_last <- NA_real_
    if (length(biotic_cols_cont)) {
      bf_ok <- vector("list", length(biotic_cols_cont)); names(bf_ok) <- biotic_cols_cont
      for (b in biotic_cols_cont) {
        before <- biotic_cols_cont[seq_len(match(b, biotic_cols_cont) - 1L)]
        before_ok <- if (length(before)) Reduce(`&`, bf_ok[before]) else rep(TRUE, n_bf)
        bf_ok[[b]] <- abiotic_ok & before_ok
      }
      frac_first <- mean(bf_ok[[biotic_cols_cont[1]]])
      frac_last  <- mean(bf_ok[[biotic_cols_cont[length(biotic_cols_cont)]]])
    }

    cascadeB[[length(cascadeB) + 1L]] <- tibble(
      subbasin = i, n_backfill_px = n_bf,
      n_abiotic_cols = length(abiotic_cols), n_abiotic_present = length(abiotic_present),
      frac_any_abiotic_na = mean(!abiotic_ok),
      dominant_na_class = dom_class, dominant_na_frac = dom_frac,
      frac_backfillable_first = frac_first, frac_backfillable_last = frac_last)

    # ---- STAGE C: rare-covariate constant / const=NA (08A:101-105) ----
    # uses the low-HF TRAINING data, exactly as 08A's `idx`.
    for (b in biotic_cols) {
      tv <- df_train[[b]]
      idx_tr <- which(!is.na(tv))
      n_obs    <- length(idx_tr)
      n_unique <- if (n_obs) length(unique(tv[idx_tr])) else 0L
      rareC[[length(rareC) + 1L]] <- tibble(
        subbasin = i, covariate = b,
        kind = if (b %in% categorical_responses) "categorical" else "continuous",
        n_obs_train = n_obs, n_unique_train = n_unique,
        is_const = (n_obs == 0L || n_unique < 2L),
        is_const_NA = (n_obs == 0L))
    }

    rm(cov_s, cov_train_s, df_train, df_full, df_backfill); gc(); TRUE
  }, error = function(e) { message("subbasin ", i, " error: ", conditionMessage(e)); FALSE })

  if (i %% 25L == 0L) message(Sys.time(), "  ...subbasin ", i)
}

# =============================================================================
# assemble + write outputs
# =============================================================================
stageA <- if (length(abioticA)) bind_rows(abioticA) else tibble()
stageB <- if (length(cascadeB)) bind_rows(cascadeB) else tibble()
stageC_raw <- if (length(rareC)) bind_rows(rareC) else tibble()

write_csv(stageA, file.path(out_dir, "rootcause_abiotic_na_by_cov.csv"))
write_csv(stageB, file.path(out_dir, "rootcause_cascade_by_subbasin.csv"))

# Stage C aggregate: per covariate, how many subbasins hit const / const=NA
stageC <- if (nrow(stageC_raw)) {
  stageC_raw |>
    dplyr::group_by(covariate, kind) |>
    dplyr::summarise(
      n_subbasins_present  = sum(n_obs_train > 0L),
      n_subbasins_const    = sum(is_const),
      n_subbasins_const_NA = sum(is_const_NA),
      .groups = "drop") |>
    dplyr::arrange(dplyr::desc(n_subbasins_const_NA), dplyr::desc(n_subbasins_const))
} else tibble()
write_csv(stageC, file.path(out_dir, "rootcause_rare_response.csv"))

# Stage D rollup: subbasin cascade -> BCR (px-weighted), vs first audit Stage 3
stageD <- tibble()
if (nrow(stageB) && nrow(bcr_subbasins_ref)) {
  roll <- stageB |>
    dplyr::inner_join(bcr_subbasins_ref, by = c("subbasin" = "sub_index")) |>
    dplyr::group_by(bcr_code) |>
    dplyr::summarise(
      n_subbasins = dplyr::n(),
      # weighted mean MUST precede the n_backfill_px summary: dplyr evaluates summarise
      # args sequentially, so referencing n_backfill_px here keeps the group vector.
      frac_backfillable_last = {
        keep <- !is.na(frac_backfillable_last) & !is.na(n_backfill_px)
        if (any(keep)) sum(frac_backfillable_last[keep] * n_backfill_px[keep]) / sum(n_backfill_px[keep]) else NA_real_
      },
      n_backfill_px = sum(n_backfill_px, na.rm = TRUE),
      .groups = "drop")

  s3_path <- file.path(log_dir, "coverage_stage3_mosaic.csv")
  if (file.exists(s3_path)) {
    s3 <- readr::read_csv(s3_path, show_col_types = FALSE) |>
      dplyr::filter(required_by_model) |>
      dplyr::group_by(bcr) |>
      dplyr::summarise(stage3_nonNA_ref = mean(frac_highhf_px_nonNA, na.rm = TRUE), .groups = "drop")
    roll <- roll |>
      dplyr::left_join(s3, by = c("bcr_code" = "bcr")) |>
      dplyr::mutate(match_delta = frac_backfillable_last - stage3_nonNA_ref)
  }
  stageD <- roll
}
write_csv(stageD, file.path(out_dir, "rootcause_bcr_rollup.csv"))

message(Sys.time(), " | DONE. rootcause_*.csv in ", out_dir)
