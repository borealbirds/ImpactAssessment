# ---
# title: Shapley value sector attribution of bird population impacts
# author: Mannfred Boehm
# ---
# Computes exact Shapley values for each sector's contribution to the total
# industrial footprint impact on bird populations.  Shapley values sum exactly
# to the total HF impact v(N) = cf(all sectors) - observed.
#
# Architecture:
#   - Reads 12H's per-sample subbasin Shapley values ({species}_{year}_shapley_samples.rds:
#     one value per subbasin x sector x (bootstrap, scenario) sample, computed in 12F).
#   - Sums the samples bottom-up: subbasin -> BCR -> national, sample by sample, and
#     summarises each level over the samples (mean, SD, 5th/95th percentiles).
#   - Cross-checks the sample means against the Shapley values of the coalition density
#     tables' means (both are exact; they agree to floating-point error).
#   - Annotates subbasins with 10C's abiotic extrapolation flags.
#
# Why samples, not table SDs: every subbasin of a BCR is predicted by the same 32 bird
# models, and v(S) and v(S + j) come from the same samples, so both are strongly
# correlated. Propagating the tables' SDs as if independent understated v(S) SDs 1.1-2.7x
# and gave small sectors the large coalitions' noise (C2, 2026-09-24; see TODO.md).
#
# Outputs (data/derived_data/sector_effects/):
#   shapley_subbasin.csv  — per-sector Shapley values at each subbasin
#   shapley_bcr.csv       — aggregated to BCR
#   shapley_national.csv  — aggregated to national
#
# Release filter: BCRs whose models BAM withheld (currently CAWA can40) are
# dropped here, not upstream — see `withheld_models` below.
# ---

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

# ---- Execution context -------------------------------------------------------

cc    <- FALSE
local <- TRUE

if (cc)            { ia_dir <- "/home/mannfred/scratch/impact_assessment" }
if (!cc && local)  { ia_dir <- getwd() }
if (!cc && !local) { ia_dir <- file.path("G:/Shared drives/BAM_NationalModels5", "data", "Extras",
                                          "sandbox_data", "impactassessment_sandbox") }

# ---- Source utilities --------------------------------------------------------

source(file.path(ia_dir, "Rscripts", "12E_shapley_utils.R"))

# ---- Paths -------------------------------------------------------------------

dt_dir     <- file.path(ia_dir, "data/derived_data/density_tables")
basin_path <- file.path(ia_dir, "data/raw_data/hydrobasins_masked_merged_subset.gpkg")
flag_path  <- file.path(ia_dir, "data/derived_data/rds_files/extrapolation_flags.csv")
out_dir    <- file.path(ia_dir, "data/derived_data/sector_effects")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Load hydrobasins --------------------------------------------------------

hydrobasins <- terra::vect(basin_path)

# ---- BAM model release filter ------------------------------------------------
# Our BCR discovery reads 06_bootstraps/{spp}/can*.Rdata, which returns every model
# BAM FIT — a superset of the models BAM RELEASED. review/ModelReleaseDecisions.xlsx
# ("remove" tab) withholds CAWA can40 on AUC. 12A/12D deliberately still produce it
# (the products stay a complete record of what we ran); the release filter belongs
# here, at the point where numbers are reported.
#
# Set DROP_WITHHELD <- FALSE to keep them and inspect the difference.
#
# Verified 2026-09-11 against the workbook: the "remove" tab has 662 rows, of which
# exactly one touches our two species (CAWA can40, AUC = 1; OVEN has none). At the
# planned ~60-species scale, replace this literal with a read of the workbook —
# 14B runs locally, so G: is reachable:
#   readxl::read_excel(file.path(nm_root, "review", "ModelReleaseDecisions.xlsx"),
#                      sheet = "remove") |> dplyr::select(species = spp, bcr = region)

DROP_WITHHELD <- TRUE

withheld_models <- data.frame(
  species = "CAWA",
  bcr     = "can40",
  reason  = "withheld by BAM (review/ModelReleaseDecisions.xlsx, 'remove' tab; AUC)",
  stringsAsFactors = FALSE
)

# ---- Load extrapolation flags (optional) -------------------------------------
# 10C, one row per subbasin (gpkg row index = 12D's sub_index): the share of its
# backfilled pixels outside the training pixels' area of applicability, and the flag
# (> 0.5) built from it. A subbasin 10C could not assess (< 10 training pixels) is NA.

extrap_flags <- NULL
if (file.exists(flag_path)) {
  extrap_flags <- read.csv(flag_path, stringsAsFactors = FALSE)
  if (!"frac_outside_aoa" %in% names(extrap_flags))
    stop(flag_path, " predates the 2026-09-25 AOA rewrite of 10C - re-run 10C")
  message("Loaded extrapolation flags for ", nrow(extrap_flags), " subbasins (",
          sum(extrap_flags$flag, na.rm = TRUE), " flagged)")
  extrap_flags <- extrap_flags[, c("subbasin", "flag", "frac_outside_aoa")]
  names(extrap_flags)[names(extrap_flags) == "flag"] <- "extrapolation_flag"
}

# ---- Discover available species x year combinations --------------------------

dt_files <- list.files(dt_dir, pattern = "^.*_coalition_[0-9]+\\.rds$", full.names = TRUE)
if (length(dt_files) == 0) stop("No coalition density tables found in ", dt_dir)

# parse filenames: {species}_{year}_coalition_{id}.rds
parsed <- regmatches(basename(dt_files),
  regexec("^(.+)_([0-9]{4})_coalition_([0-9]+)\\.rds$", basename(dt_files)))
parsed <- do.call(rbind, lapply(parsed, function(x) x[2:4]))
dt_index <- data.frame(
  path         = dt_files,
  species      = parsed[, 1],
  year         = parsed[, 2],
  coalition_id = as.integer(parsed[, 3]),
  stringsAsFactors = FALSE
)

species_years <- unique(dt_index[, c("species", "year")])
message("Found ", nrow(species_years), " species x year combinations, ",
        nrow(dt_index), " coalition files total")

sectors   <- canonical_sectors()
n_sectors <- length(sectors)
n_coal    <- n_coalitions(sectors)
W_shap    <- shapley_weight_matrix(sectors)          # phi = W %*% v (12E)

# summaries over samples of a [unit x sector x sample] array -> [unit x sector] each
summarise_samples <- function(x) list(
  mean = apply(x, c(1, 2), mean),
  sd   = apply(x, c(1, 2), sd),
  q05  = apply(x, c(1, 2), quantile, 0.05, names = FALSE),
  q95  = apply(x, c(1, 2), quantile, 0.95, names = FALSE))

# one row per unit x sector (sector fastest), as the output CSVs are laid out
long <- function(m) as.vector(t(m))

# ---- Main processing: one species x year at a time ---------------------------

shapley_sub_rows <- vector("list", nrow(species_years))
shapley_bcr_rows <- vector("list", nrow(species_years))
shapley_nat_rows <- vector("list", nrow(species_years))

for (sy in seq_len(nrow(species_years))) {

  sp  <- species_years$species[sy]
  yr  <- species_years$year[sy]

  message("\n=== ", sp, " ", yr, " ===")

  # find all coalition files for this species x year
  idx <- dt_index$species == sp & dt_index$year == yr
  avail_ids <- dt_index$coalition_id[idx]
  avail_paths <- dt_index$path[idx]

  # all 2^N - 1 non-empty coalitions are needed (ID 1 = empty, v = 0); the cross-check
  # against the samples below fails on any that are missing
  expected_ids <- 2:n_coal
  missing_ids  <- setdiff(expected_ids, avail_ids)
  if (length(missing_ids) > 0)
    stop(sp, " ", yr, ": missing ", length(missing_ids), " coalition table(s): ",
         paste(head(missing_ids, 10), collapse = ", "), if (length(missing_ids) > 10) "...")

  # read all coalition density tables into a list keyed by coalition_id,
  # dropping any BCR whose model BAM withheld (see withheld_models above)
  drop_bcrs <- withheld_models$bcr[withheld_models$species == sp]
  n_dropped <- 0L
  dt_list <- setNames(
    lapply(seq_along(avail_ids), function(i) {
      d <- readRDS(avail_paths[i])
      if (DROP_WITHHELD && length(drop_bcrs) > 0L && "bcr" %in% names(d)) {
        hit <- d$bcr %in% drop_bcrs
        n_dropped <<- n_dropped + sum(hit)
        d <- d[!hit, , drop = FALSE]
      }
      d
    }),
    as.character(avail_ids)
  )
  if (length(drop_bcrs) > 0L) {
    if (DROP_WITHHELD) {
      message("  release filter: dropped ", n_dropped, " subbasin-rows in BCR(s) ",
              paste(drop_bcrs, collapse = ", "), " — ",
              paste(unique(withheld_models$reason[withheld_models$species == sp]),
                    collapse = "; "))
      if (n_dropped == 0L)
        message("  NOTE: no rows matched ", paste(drop_bcrs, collapse = ", "),
                " — either 12D did not produce it, or the BCR code has changed.")
    } else {
      message("  release filter DISABLED — withheld BCR(s) ",
              paste(drop_bcrs, collapse = ", "), " are INCLUDED in these numbers.")
    }
  }

  # get the set of all BCRs x subbasins across coalitions
  all_rows <- bind_rows(dt_list, .id = "coal_id")
  bcr_sub_ref <- distinct(all_rows, bcr, subbasin)
  sub_keys   <- paste(bcr_sub_ref$bcr, bcr_sub_ref$subbasin, sep = "::")
  bcr_lookup <- setNames(bcr_sub_ref$bcr, sub_keys)
  sub_lookup <- setNames(bcr_sub_ref$subbasin, sub_keys)

  # ---- v(S) = cf(S) - obs = bf_on_coalition - obs_on_coalition, from the tables' means ----
  # (impact of removing coalition S: positive means more birds without S)
  v_mean_mat     <- matrix(0, nrow = length(sub_keys), ncol = n_coal,
                           dimnames = list(sub_keys, as.character(1:n_coal)))
  obs_total_mean <- setNames(numeric(length(sub_keys)), sub_keys)

  for (cid_str in names(dt_list)) {
    dt   <- dt_list[[cid_str]]
    keys <- paste(dt$bcr, dt$subbasin, sep = "::")

    # 12F marks a subbasin with no kept pixels for this coalition by a NaN obs_on mean
    # and an NA (not NaN) obs_on sd; nothing is backfilled there, so v(S) = 0.
    # Keying the zero on that marker, and stopping on any other NA, keeps a real NA
    # from being zeroed silently (an is.nan() test on the sd let NA through instead).
    empty <- is.nan(dt$obs_on_coalition_mean)
    if (any(dt$bf_on_coalition_mean[empty] != 0))
      stop(sp, " coalition ", cid_str, ": no observed pixels but bf_on != 0 at ",
           paste(keys[empty & dt$bf_on_coalition_mean != 0], collapse = ", "))
    v <- dt$bf_on_coalition_mean - dt$obs_on_coalition_mean
    v[empty] <- 0
    if (anyNA(v))
      stop(sp, " coalition ", cid_str, ": NA impact with observed pixels present at ",
           paste(keys[is.na(v)], collapse = ", "))
    v_mean_mat[keys, cid_str] <- v

    # obs_total is coalition-free, except that 12F zeroes it in a coalition with no
    # footprint in the BCR (kept for bit-identity with the retired per-coalition path)
    obs_total_mean[keys] <- pmax(obs_total_mean[keys], dt$obs_total_mean, na.rm = TRUE)
  }

  # Shapley values of the mean v(S): the cross-check for the samples below
  phi_tab <- v_mean_mat %*% t(W_shap)
  for (k in head(sub_keys, 3L)) {
    ref <- compute_shapley(setNames(v_mean_mat[k, ], colnames(v_mean_mat)), sectors)
    if (max(abs(ref - phi_tab[k, ])) > 1e-9 * max(1, abs(v_mean_mat[k, ])))
      stop("shapley_weight_matrix() disagrees with compute_shapley() at ", k)
  }

  # ---- per-sample subbasin Shapley values (12F -> 12H) ----
  ss_path <- file.path(dt_dir, paste0(sp, "_", yr, "_shapley_samples.rds"))
  if (!file.exists(ss_path))
    stop("missing ", ss_path, " - it is written by 12H from 12D's per-BCR files. 12F saves ",
         "per-sample Shapley values only since 2026-09-25, so older tables need a 12D + 12H rerun.")
  ss      <- readRDS(ss_path)
  ss_keys <- paste(ss$bcr, ss$subbasin, sep = "::")
  keep    <- !(DROP_WITHHELD & ss$bcr %in% drop_bcrs)
  if (!setequal(ss_keys[keep], sub_keys) || anyDuplicated(ss_keys[keep]))
    stop(sp, " ", yr, ": the Shapley samples and the coalition tables cover different ",
         "subbasins - they come from different 12H merges")
  rows <- which(keep)[match(sub_keys, ss_keys[keep])]
  phi  <- ss$phi[rows, , , drop = FALSE]                   # [subbasin x sector x sample]
  if (!identical(dimnames(phi)[[2L]], sectors))
    stop(sp, " ", yr, ": Shapley samples carry sectors ", paste(dimnames(phi)[[2L]], collapse = ", "))
  n_samp <- dim(phi)[3L]
  message("  ", length(sub_keys), " subbasin rows x ", n_samp, " samples (",
          ss$n_boot, " bootstraps x ", ss$n_scen, " scenarios)")

  # the samples must reproduce the tables: mean of per-sample Shapley = Shapley of means
  st_sub <- summarise_samples(phi)
  gap    <- max(abs(st_sub$mean - phi_tab))
  if (gap > 1e-8 * max(1, abs(v_mean_mat)))
    stop(sp, " ", yr, ": mean of the Shapley samples is ", signif(gap, 3), " off the Shapley ",
         "values of the tables' means - samples and tables are from different runs")
  message("  sample means match the tables' Shapley values (max abs diff ", signif(gap, 3), ")")

  # total HF impact v(N) per sample = sum over sectors (efficiency holds per sample)
  tot_sub <- apply(phi, c(1, 3), sum)                      # [subbasin x sample]

  # ---- subbasin table ----
  obs_pop <- unname(obs_total_mean)
  shapley_sub_df <- data.frame(
    species            = sp,
    bcr                = rep(unname(bcr_lookup), each = n_sectors),
    year               = yr,
    subbasin           = rep(unname(sub_lookup), each = n_sectors),
    HYBAS_ID           = rep(hydrobasins$first_HYBAS_ID[as.integer(sub_lookup)], each = n_sectors),
    obs_population     = rep(round(obs_pop), each = n_sectors),
    total_HF_impact    = rep(round(rowMeans(tot_sub)), each = n_sectors),
    total_HF_impact_sd = rep(round(apply(tot_sub, 1, sd)), each = n_sectors),
    sector             = rep(sectors, times = length(sub_keys)),
    shapley_mean       = round(long(st_sub$mean), 2),
    shapley_sd         = round(long(st_sub$sd), 2),
    shapley_q05        = round(long(st_sub$q05), 2),
    shapley_q95        = round(long(st_sub$q95), 2),
    shapley_pct        = round(long(st_sub$mean / obs_pop * 100), 4),
    shapley_check      = rep(round(rowSums(st_sub$mean), 2), each = n_sectors),  # = v(N)
    stringsAsFactors   = FALSE
  )
  if (!is.null(extrap_flags))
    shapley_sub_df <- left_join(shapley_sub_df, extrap_flags, by = "subbasin")
  shapley_sub_rows[[sy]] <- shapley_sub_df

  # ---- BCR: sum the subbasin samples within each BCR, sample by sample ----
  phi_bcr <- rowsum(matrix(phi, nrow = length(sub_keys)), unname(bcr_lookup), reorder = FALSE)
  bcrs    <- rownames(phi_bcr)
  phi_bcr <- array(phi_bcr, c(length(bcrs), n_sectors, n_samp))
  st_bcr  <- summarise_samples(phi_bcr)
  tot_bcr <- apply(phi_bcr, c(1, 3), sum)
  bcr_of  <- unname(bcr_lookup)
  obs_bcr <- vapply(bcrs, function(b) sum(obs_pop[bcr_of == b]), numeric(1L))
  n_sub_b <- vapply(bcrs, function(b) sum(bcr_of == b), integer(1L))
  # extrapolation: how many subbasins are flagged, and how much of the impact they carry
  flagged <- if (is.null(extrap_flags)) rep(NA, length(sub_keys)) else
    extrap_flags$extrapolation_flag[match(as.integer(sub_lookup), extrap_flags$subbasin)]
  imp_sub <- rowMeans(tot_sub)
  n_flg_b <- vapply(bcrs, function(b) if (is.null(extrap_flags)) NA_integer_ else
    sum(flagged[bcr_of == b], na.rm = TRUE), integer(1L))
  imp_flg_b <- vapply(bcrs, function(b) if (is.null(extrap_flags)) NA_real_ else
    sum(imp_sub[bcr_of == b & flagged %in% TRUE]), numeric(1L))
  shapley_bcr_rows[[sy]] <- data.frame(
    species            = sp,
    bcr                = rep(bcrs, each = n_sectors),
    year               = yr,
    sector             = rep(sectors, times = length(bcrs)),
    obs_population     = rep(round(obs_bcr), each = n_sectors),
    total_HF_impact    = rep(round(rowMeans(tot_bcr)), each = n_sectors),
    total_HF_impact_sd = rep(round(apply(tot_bcr, 1, sd)), each = n_sectors),
    shapley_mean       = round(long(st_bcr$mean), 2),
    shapley_sd         = round(long(st_bcr$sd), 2),
    shapley_q05        = round(long(st_bcr$q05), 2),
    shapley_q95        = round(long(st_bcr$q95), 2),
    shapley_pct        = round(long(st_bcr$mean / obs_bcr * 100), 4),
    n_subbasins        = rep(n_sub_b, each = n_sectors),
    n_flagged          = rep(n_flg_b, each = n_sectors),
    total_HF_impact_flagged = rep(round(imp_flg_b), each = n_sectors),
    stringsAsFactors   = FALSE
  )

  # ---- national: sum the BCR samples, sample by sample ----
  # Bootstrap i is paired across BCRs, as V5's national estimate pairs them (mosaic of
  # bootstrap i); that is right whether or not BCRs' bootstrap i share a resample.
  phi_nat <- array(apply(phi_bcr, c(2, 3), sum), c(1L, n_sectors, n_samp))
  st_nat  <- summarise_samples(phi_nat)
  tot_nat <- apply(phi_nat, c(1, 3), sum)
  shapley_nat_rows[[sy]] <- data.frame(
    species            = sp,
    year               = yr,
    sector             = sectors,
    obs_population     = round(sum(obs_pop)),
    total_HF_impact    = round(mean(tot_nat)),
    total_HF_impact_sd = round(sd(tot_nat)),
    total_HF_impact_q05 = round(quantile(tot_nat, 0.05, names = FALSE)),
    total_HF_impact_q95 = round(quantile(tot_nat, 0.95, names = FALSE)),
    shapley_mean       = round(long(st_nat$mean), 2),
    shapley_sd         = round(long(st_nat$sd), 2),
    shapley_q05        = round(long(st_nat$q05), 2),
    shapley_q95        = round(long(st_nat$q95), 2),
    shapley_pct        = round(long(st_nat$mean / sum(obs_pop) * 100), 4),
    n_bcrs             = length(bcrs),
    n_flagged          = if (is.null(extrap_flags)) NA_integer_ else sum(flagged, na.rm = TRUE),
    total_HF_impact_flagged = round(sum(imp_flg_b)),
    stringsAsFactors   = FALSE
  )

  message(sprintf("  additivity: sum(phi) = %.3f, v(N) = %.3f (per-sample max |residual| %.2e)",
                  sum(st_nat$mean), mean(tot_nat),
                  max(abs(colSums(phi_nat[1, , ]) - tot_nat[1, ]))))
  message("  ", sp, " ", yr, ": ", nrow(shapley_sub_df), " subbasin rows, ",
          nrow(shapley_bcr_rows[[sy]]), " BCR rows")
}

# ---- Assemble and write ------------------------------------------------------

shapley_sub_all  <- bind_rows(shapley_sub_rows)
shapley_bcr_all  <- bind_rows(shapley_bcr_rows)
shapley_national <- bind_rows(shapley_nat_rows)

write.csv(shapley_sub_all, file.path(out_dir, "shapley_subbasin.csv"), row.names = FALSE)
write.csv(shapley_bcr_all, file.path(out_dir, "shapley_bcr.csv"),      row.names = FALSE)
write.csv(shapley_national, file.path(out_dir, "shapley_national.csv"), row.names = FALSE)

message("\nWrote 3 tables to ", out_dir)
message("  shapley_subbasin.csv : ", nrow(shapley_sub_all), " rows")
message("  shapley_bcr.csv      : ", nrow(shapley_bcr_all), " rows")
message("  shapley_national.csv : ", nrow(shapley_national), " rows")
