# ---
# title: Shapley value sector attribution of bird population impacts
# author: Mannfred Boehm
# ---
# Computes exact Shapley values for each sector's contribution to the total
# industrial footprint impact on bird populations.  Shapley values sum exactly
# to the total HF impact v(N) = cf(all sectors) - observed.
#
# Architecture:
#   - Reads coalition density tables (one per coalition x species x year)
#     produced by 12A/12B.
#   - For each coalition, computes v(S) = cf(S) - obs at subbasin level.
#   - Applies the Shapley formula per sector per subbasin.
#   - Aggregates bottom-up: subbasin -> BCR -> national.
#   - Optionally filters by abiotic extrapolation flags.
#
# Outputs (data/derived_data/sector_effects/):
#   shapley_subbasin.csv  — per-sector Shapley values at each subbasin
#   shapley_bcr.csv       — aggregated to BCR
#   shapley_national.csv  — aggregated to national
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

# ---- Load extrapolation flags (optional) -------------------------------------

extrap_flags <- NULL
if (file.exists(flag_path)) {
  extrap_flags <- read.csv(flag_path, stringsAsFactors = FALSE)
  message("Loaded extrapolation flags for ", nrow(extrap_flags), " subbasins (",
          sum(extrap_flags$flag, na.rm = TRUE), " flagged)")
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

sectors <- canonical_sectors()
n_sectors <- length(sectors)
n_coal    <- n_coalitions(sectors)

# ---- Main processing: one species x year at a time ---------------------------

shapley_sub_rows <- vector("list", nrow(species_years))
shapley_bcr_rows <- vector("list", nrow(species_years))

for (sy in seq_len(nrow(species_years))) {

  sp  <- species_years$species[sy]
  yr  <- species_years$year[sy]

  message("\n=== ", sp, " ", yr, " ===")

  # find all coalition files for this species x year
  idx <- dt_index$species == sp & dt_index$year == yr
  avail_ids <- dt_index$coalition_id[idx]
  avail_paths <- dt_index$path[idx]

  # check completeness (need all 2^N - 1 non-empty coalitions; ID 1 = empty, v=0)
  expected_ids <- 2:n_coal
  missing_ids  <- setdiff(expected_ids, avail_ids)
  if (length(missing_ids) > 0) {
    message("  WARNING: missing ", length(missing_ids), " coalitions: ",
            paste(head(missing_ids, 10), collapse = ", "),
            if (length(missing_ids) > 10) "..." else "")
    message("  Shapley values will be approximate (missing coalitions treated as v=0)")
  }

  # read all coalition density tables into a list keyed by coalition_id
  dt_list <- setNames(
    lapply(seq_along(avail_ids), function(i) readRDS(avail_paths[i])),
    as.character(avail_ids)
  )

  # get the set of all BCRs x subbasins across coalitions
  all_rows <- bind_rows(dt_list, .id = "coal_id")
  bcr_sub_ref <- distinct(all_rows, bcr, subbasin)

  # ---- Compute v(S) = cf(S) - obs per subbasin per coalition ----
  # v(S) = (obs - obs_on_coalition + bf_on_coalition) - obs
  #       = bf_on_coalition - obs_on_coalition
  # (impact of removing coalition S: positive means more birds without S)

  # build a matrix: rows = subbasin key, columns = coalition IDs, values = v(S) mean
  # and a parallel matrix for v(S) sd
  sub_keys <- paste(bcr_sub_ref$bcr, bcr_sub_ref$subbasin, sep = "::")
  v_mean_mat <- matrix(0, nrow = length(sub_keys), ncol = n_coal,
                        dimnames = list(sub_keys, as.character(1:n_coal)))
  v_sd_mat   <- matrix(0, nrow = length(sub_keys), ncol = n_coal,
                        dimnames = list(sub_keys, as.character(1:n_coal)))

  # also store obs_total per subbasin (from any coalition; should be identical)
  obs_total_mean <- setNames(numeric(length(sub_keys)), sub_keys)
  obs_total_sd   <- setNames(numeric(length(sub_keys)), sub_keys)
  bcr_lookup     <- setNames(bcr_sub_ref$bcr, sub_keys)
  sub_lookup     <- setNames(bcr_sub_ref$subbasin, sub_keys)

  for (cid_str in names(dt_list)) {
    dt <- dt_list[[cid_str]]
    cid <- as.integer(cid_str)

    for (r in seq_len(nrow(dt))) {
      key <- paste(dt$bcr[r], dt$subbasin[r], sep = "::")
      if (!key %in% sub_keys) next

      impact_mean <- dt$bf_on_coalition_mean[r] - dt$obs_on_coalition_mean[r]
      impact_sd   <- sqrt(dt$bf_on_coalition_sd[r]^2 + dt$obs_on_coalition_sd[r]^2)
      if (is.nan(impact_mean)) impact_mean <- 0
      if (is.nan(impact_sd))   impact_sd   <- 0

      v_mean_mat[key, as.character(cid)] <- impact_mean
      v_sd_mat[key, as.character(cid)]   <- impact_sd

      # store obs_total (same across coalitions; take max to handle rounding)
      if (obs_total_mean[key] == 0) {
        obs_total_mean[key] <- dt$obs_total_mean[r]
        obs_total_sd[key]   <- dt$obs_total_sd[r]
      }
    }
  }

  # coalition 1 (empty set) has v = 0 by definition — already initialized

  # ---- Compute Shapley values per subbasin ----

  shapley_sub <- vector("list", length(sub_keys))

  for (si in seq_along(sub_keys)) {
    key <- sub_keys[si]

    # named vector of v(S) for all coalitions
    v_vec <- setNames(v_mean_mat[key, ], as.character(1:n_coal))

    # compute Shapley values
    phi <- compute_shapley(v_vec, sectors)

    # Shapley SD: propagate from coalition SDs using the Shapley weights
    # phi_j = sum_S w(S) * [v(S+j) - v(S)]
    # Var(phi_j) ≈ sum_S w(S)^2 * [Var(v(S+j)) + Var(v(S))]
    phi_sd <- setNames(numeric(n_sectors), sectors)
    for (j in seq_len(n_sectors)) {
      others <- sectors[-j]
      var_j <- 0
      for (bits in 0:(2^(n_sectors - 1) - 1)) {
        S      <- others[which(as.logical(intToBits(bits)[1:(n_sectors - 1)]))]
        s_size <- length(S)
        w <- factorial(s_size) * factorial(n_sectors - s_size - 1L) / factorial(n_sectors)

        id_with    <- sectors_to_coalition_id(c(S, sectors[j]), sectors)
        id_without <- sectors_to_coalition_id(S, sectors)

        var_j <- var_j + w^2 * (v_sd_mat[key, as.character(id_with)]^2 +
                                 v_sd_mat[key, as.character(id_without)]^2)
      }
      phi_sd[j] <- sqrt(var_j)
    }

    # full coalition impact (total HF impact at this subbasin)
    full_id <- n_coal
    v_full  <- v_mean_mat[key, as.character(full_id)]

    obs_pop <- obs_total_mean[key]

    shapley_sub[[si]] <- data.frame(
      species         = sp,
      bcr             = bcr_lookup[key],
      year            = yr,
      subbasin        = sub_lookup[key],
      HYBAS_ID        = hydrobasins$first_HYBAS_ID[as.integer(sub_lookup[key])],
      obs_population  = round(obs_pop),
      total_HF_impact = round(v_full),
      sector          = sectors,
      shapley_mean    = round(phi, 2),
      shapley_sd      = round(phi_sd, 2),
      shapley_pct     = round(phi / obs_pop * 100, 4),
      shapley_check   = round(sum(phi), 2),  # should equal v_full
      stringsAsFactors = FALSE
    )
  }

  shapley_sub_df <- bind_rows(shapley_sub)

  # attach extrapolation flag if available
  if (!is.null(extrap_flags)) {
    shapley_sub_df <- left_join(
      shapley_sub_df,
      extrap_flags[, c("subbasin", "flag", "ks_max", "mahal_exceedance")],
      by = "subbasin"
    )
    names(shapley_sub_df)[names(shapley_sub_df) == "flag"] <- "extrapolation_flag"
  }

  shapley_sub_rows[[sy]] <- shapley_sub_df

  # ---- Aggregate to BCR (sum of subbasin Shapley values) ----

  shapley_bcr <- shapley_sub_df |>
    group_by(species, bcr, year, sector) |>
    summarise(
      obs_population  = round(sum(obs_population)),
      total_HF_impact = round(sum(total_HF_impact)),
      shapley_mean    = round(sum(shapley_mean), 2),
      shapley_sd      = round(sqrt(sum(shapley_sd^2)), 2),
      shapley_pct     = round(sum(shapley_mean) / sum(obs_population) * 100, 4),
      n_subbasins     = n(),
      n_flagged       = if ("extrapolation_flag" %in% names(shapley_sub_df))
                          sum(extrapolation_flag, na.rm = TRUE) else NA_integer_,
      .groups = "drop"
    )

  shapley_bcr_rows[[sy]] <- shapley_bcr

  message("  ", sp, " ", yr, ": ", nrow(shapley_sub_df), " subbasin rows, ",
          nrow(shapley_bcr), " BCR rows")
}

# ---- Assemble final tables ---------------------------------------------------

shapley_sub_all <- bind_rows(shapley_sub_rows)
shapley_bcr_all <- bind_rows(shapley_bcr_rows)

# ---- National scale (sum of BCR Shapley values) -----------------------------

shapley_national <- shapley_bcr_all |>
  group_by(species, year, sector) |>
  summarise(
    obs_population  = round(sum(obs_population)),
    total_HF_impact = round(sum(total_HF_impact)),
    shapley_mean    = round(sum(shapley_mean), 2),
    shapley_sd      = round(sqrt(sum(shapley_sd^2)), 2),
    shapley_pct     = round(sum(shapley_mean) / sum(obs_population) * 100, 4),
    n_bcrs          = n(),
    .groups = "drop"
  )

# ---- Verify additivity -------------------------------------------------------

check <- shapley_national |>
  group_by(species, year) |>
  summarise(
    sum_shapley     = sum(shapley_mean),
    total_HF_impact = first(total_HF_impact),
    residual        = sum_shapley - total_HF_impact,
    .groups = "drop"
  )

message("\n=== Shapley additivity check ===")
for (i in seq_len(nrow(check))) {
  message(sprintf("  %s %s: sum(phi) = %.1f, v(N) = %.1f, residual = %.1f",
                  check$species[i], check$year[i],
                  check$sum_shapley[i], check$total_HF_impact[i], check$residual[i]))
}

# ---- Write outputs -----------------------------------------------------------

write.csv(shapley_sub_all, file.path(out_dir, "shapley_subbasin.csv"), row.names = FALSE)
write.csv(shapley_bcr_all, file.path(out_dir, "shapley_bcr.csv"),      row.names = FALSE)
write.csv(shapley_national, file.path(out_dir, "shapley_national.csv"), row.names = FALSE)

message("\nWrote 3 tables to ", out_dir)
message("  shapley_subbasin.csv : ", nrow(shapley_sub_all), " rows")
message("  shapley_bcr.csv      : ", nrow(shapley_bcr_all), " rows")
message("  shapley_national.csv : ", nrow(shapley_national), " rows")
