# ---
# title: Shapley value sector attribution of bird population impacts
# author: Mannfred Boehm
# ---

# each sector's exact Shapley value at subbasin, BCR and national scale, from 12H's per-sample
# subbasin Shapley values ({species}_{year}_shapley_samples.rds), summed bottom-up sample by
# sample (a BCR's subbasins share its 32 bird models, so table SDs propagated as independent
# understated v(S) SDs 1.1-2.7x). each value splits into a direct part (the bird model's
# response to footprint covariates at 0 on the observed vegetation) and a vegetation part.
# sign: v(S) = counterfactual - observed (positive = more birds without industry; 14C flips it)
#   -> data/derived_data/sector_effects/shapley_{subbasin,bcr,national}.csv   (_{year} for years other than 2020)
# BAM's release filter is applied here: every Canadian row of the "remove" tab of
# data/raw_data/ModelReleaseDecisions.xlsx (BAM's, gitignored) is dropped from the numbers.
# on Fir, after 12H: bash 14B_sector_attribution.sh <12H job id>   (runs 14C, chains 18)

library(dplyr)

ia_dir <- if (nzchar(Sys.getenv("SLURM_JOB_ID"))) "/home/mannfred/scratch/impact_assessment" else getwd()
dt_dir  <- file.path(ia_dir, "data/derived_data/density_tables")
out_dir <- file.path(ia_dir, "data/derived_data/sector_effects")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

source(file.path(ia_dir, "Rscripts", "12E_shapley_utils.R"))

year <- as.integer(Sys.getenv("YEARS", "2020"))   # the .sh submits one year per job
ysfx <- if (year == 2020L) "" else paste0("_", year)

hydrobasins <- terra::vect(file.path(ia_dir, "data/raw_data/hydrobasins_masked_merged_subset.gpkg"))

# BAM's release filter: the BCRs whose models BAM withheld, with the reasons (its 0/1 columns)
remove_tab <- readxl::read_excel(file.path(ia_dir, "data/raw_data/ModelReleaseDecisions.xlsx"), sheet = "remove")
remove_tab <- remove_tab[grepl("^can", remove_tab$region), ]
flag_cols  <- setdiff(names(remove_tab), c("spp", "region", "mosaic"))
withheld_models <- data.frame(
  species = remove_tab$spp,
  bcr     = remove_tab$region,
  reason  = paste0(remove_tab$region, " withheld by BAM (",
                   apply(as.matrix(remove_tab[flag_cols]) == 1, 1, function(x) paste(flag_cols[which(x)], collapse = ", ")), ")"))

# 10C's extrapolation flags, one row per subbasin (its gpkg row = 12D's sub_index; NA where
# 10C had < 10 training pixels)
extrap_flags <- read.csv(file.path(ia_dir, "data/derived_data/rds_files/extrapolation_flags.csv"))
extrap_flags <- extrap_flags[, c("subbasin", "flag", "frac_outside_aoa")]
names(extrap_flags)[2] <- "extrapolation_flag"

# this year's coalition tables: {species}_{year}_coalition_{id}.rds
dt_files <- list.files(dt_dir, pattern = paste0("^.*_", year, "_coalition_[0-9]+\\.rds$"), full.names = TRUE)
if (length(dt_files) == 0) stop("no coalition tables for ", year, " in ", dt_dir)
parsed   <- do.call(rbind, regmatches(basename(dt_files), regexec("^(.+)_([0-9]{4})_coalition_([0-9]+)\\.rds$", basename(dt_files))))
dt_index <- data.frame(path = dt_files, species = parsed[, 2], year = parsed[, 3], coalition_id = as.integer(parsed[, 4]))
species_years <- unique(dt_index[, c("species", "year")])
message(nrow(species_years), " species, ", nrow(dt_index), " coalition tables")

sectors   <- canonical_sectors()
n_sectors <- length(sectors)
n_coal    <- 2L^n_sectors
W_shap    <- shapley_weight_matrix(sectors)   # phi = W %*% v

# mean, SD and 5th/95th percentiles over samples: [unit x sector x sample] -> [unit x sector] each
summarise_samples <- function(x) list(
  mean = apply(x, c(1, 2), mean),
  sd   = apply(x, c(1, 2), sd),
  q05  = apply(x, c(1, 2), quantile, 0.05, names = FALSE),
  q95  = apply(x, c(1, 2), quantile, 0.95, names = FALSE))

# one row per unit x sector (sector fastest), the CSVs' layout
long <- function(m) as.vector(t(m))

# direct and vegetation parts of a [unit x sector x sample] array, from the direct part per
# bootstrap xd [unit x sector x bootstrap]; boot_of maps each sample to its bootstrap
split_samples <- function(x, xd, boot_of) {
  xd_s <- xd[, , boot_of, drop = FALSE]
  xv   <- x - xd_s
  list(direct_mean = apply(xd_s, c(1, 2), mean), direct_sd = apply(xd_s, c(1, 2), sd),
       veg_mean    = apply(xv,   c(1, 2), mean), veg_sd    = apply(xv,   c(1, 2), sd))
}

# sum a [unit x sector x k] array over the units of each group -> [group x sector x k]
sum_units <- function(x, g) {
  r <- rowsum(matrix(x, nrow = dim(x)[1L]), g, reorder = FALSE)
  array(r, c(nrow(r), dim(x)[2L], dim(x)[3L]), dimnames = list(rownames(r), dimnames(x)[[2L]], NULL))
}



# ---------------------------------------------------
# one species at a time

shapley_sub_rows <- shapley_bcr_rows <- shapley_nat_rows <- vector("list", nrow(species_years))

for (sy in seq_len(nrow(species_years))) {

  sp <- species_years$species[sy]
  yr <- species_years$year[sy]
  message("\n=== ", sp, " ", yr, " ===")

  idx <- dt_index$species == sp & dt_index$year == yr
  missing_ids <- setdiff(2:n_coal, dt_index$coalition_id[idx])
  if (length(missing_ids) > 0)
    stop(sp, " ", yr, ": ", length(missing_ids), " coalition tables missing: ", paste(head(missing_ids, 10), collapse = ", "))

  # the coalition tables, minus the BCRs BAM withheld
  drop_bcrs <- withheld_models$bcr[withheld_models$species == sp]
  dt_list <- setNames(lapply(dt_index$path[idx], function(f) {
    d <- readRDS(f)
    d[!d$bcr %in% drop_bcrs, , drop = FALSE]
  }), as.character(dt_index$coalition_id[idx]))
  if (length(drop_bcrs) > 0L)
    message("  release filter: ", paste(unique(withheld_models$reason[withheld_models$species == sp]), collapse = "; "))

  # BAM withheld every BCR of some species (BGGN, BWWA, DUNL, GWWA, NESP, SPGR, YBCU, YTVI)
  if (sum(vapply(dt_list, nrow, integer(1))) == 0L) {
    message("  no subbasin rows left after the release filter: ", sp, " is left out")
    next
  }

  bcr_sub_ref <- distinct(bind_rows(dt_list), bcr, subbasin)
  sub_keys    <- paste(bcr_sub_ref$bcr, bcr_sub_ref$subbasin, sep = "::")
  bcr_lookup  <- setNames(bcr_sub_ref$bcr, sub_keys)
  sub_lookup  <- setNames(bcr_sub_ref$subbasin, sub_keys)

  # ---- v(S) = bf_on_coalition - obs_on_coalition, from the tables' means ----
  v_mean_mat     <- matrix(0, nrow = length(sub_keys), ncol = n_coal, dimnames = list(sub_keys, as.character(1:n_coal)))
  obs_total_mean <- setNames(numeric(length(sub_keys)), sub_keys)

  for (cid_str in names(dt_list)) {
    dt   <- dt_list[[cid_str]]
    keys <- paste(dt$bcr, dt$subbasin, sep = "::")

    # 12F marks a subbasin with no kept pixels by a NaN obs_on mean (its sd is NA, not NaN);
    # nothing is backfilled there, so v(S) = 0. any other NA is an error
    empty <- is.nan(dt$obs_on_coalition_mean)
    if (any(dt$bf_on_coalition_mean[empty] != 0))
      stop(sp, " coalition ", cid_str, ": no observed pixels but bf_on != 0 at ", paste(keys[empty & dt$bf_on_coalition_mean != 0], collapse = ", "))
    v <- dt$bf_on_coalition_mean - dt$obs_on_coalition_mean
    v[empty] <- 0
    if (anyNA(v)) stop(sp, " coalition ", cid_str, ": NA impact with observed pixels present at ", paste(keys[is.na(v)], collapse = ", "))
    v_mean_mat[keys, cid_str] <- v

    # obs_total is the same in every coalition, except 0 in one with no footprint in the BCR (12F)
    obs_total_mean[keys] <- pmax(obs_total_mean[keys], dt$obs_total_mean, na.rm = TRUE)
  }

  # Shapley values of the mean v(S), checked against the textbook formula, then against the samples
  phi_tab <- v_mean_mat %*% t(W_shap)
  for (k in head(sub_keys, 3L)) {
    ref <- compute_shapley(setNames(v_mean_mat[k, ], colnames(v_mean_mat)), sectors)
    if (max(abs(ref - phi_tab[k, ])) > 1e-9 * max(1, abs(v_mean_mat[k, ])))
      stop("shapley_weight_matrix() disagrees with compute_shapley() at ", k)
  }

  # ---- the per-sample subbasin Shapley values (12F -> 12H) ----
  ss      <- readRDS(file.path(dt_dir, paste0(sp, "_", yr, "_shapley_samples.rds")))
  ss_keys <- paste(ss$bcr, ss$subbasin, sep = "::")
  keep    <- !ss$bcr %in% drop_bcrs
  if (!setequal(ss_keys[keep], sub_keys) || anyDuplicated(ss_keys[keep]))
    stop(sp, " ", yr, ": the Shapley samples and the coalition tables cover different subbasins (different 12H merges?)")
  rows <- which(keep)[match(sub_keys, ss_keys[keep])]
  phi  <- ss$phi[rows, , , drop = FALSE]               # [subbasin x sector x sample]
  if (!identical(dimnames(phi)[[2L]], sectors)) stop(sp, " ", yr, ": the Shapley samples carry other sectors")
  phi_d   <- ss$phi_direct[rows, , , drop = FALSE]     # [subbasin x sector x bootstrap]
  boot_of <- rep(seq_len(ss$n_boot), each = ss$n_scen) # samples run scenario-within-bootstrap
  if (length(boot_of) != dim(phi)[3L] || dim(phi_d)[3L] != ss$n_boot)
    stop(sp, " ", yr, ": the Shapley samples and their direct part disagree on the sample grid")

  # the mean of the per-sample Shapley values must be the Shapley value of the mean v(S)
  st_sub <- summarise_samples(phi)
  gap    <- max(abs(st_sub$mean - phi_tab))
  if (gap > 1e-8 * max(1, abs(v_mean_mat)))
    stop(sp, " ", yr, ": the Shapley samples are ", signif(gap, 3), " off the tables' Shapley values (different runs?)")

  tot_sub <- apply(phi, c(1, 3), sum)                  # v(all sectors) per subbasin and sample
  sp_sub  <- split_samples(phi, phi_d, boot_of)
  obs_pop <- unname(obs_total_mean)
  flagged <- extrap_flags$extrapolation_flag[match(as.integer(sub_lookup), extrap_flags$subbasin)]

  # ---- subbasins ----
  shapley_sub_rows[[sy]] <- left_join(data.frame(
    species                 = sp,
    bcr                     = rep(unname(bcr_lookup), each = n_sectors),
    year                    = yr,
    subbasin                = rep(unname(sub_lookup), each = n_sectors),
    HYBAS_ID                = rep(hydrobasins$first_HYBAS_ID[as.integer(sub_lookup)], each = n_sectors),
    obs_population          = rep(round(obs_pop), each = n_sectors),
    total_HF_impact         = rep(round(rowMeans(tot_sub)), each = n_sectors),
    total_HF_impact_sd      = rep(round(apply(tot_sub, 1, sd)), each = n_sectors),
    sector                  = rep(sectors, times = length(sub_keys)),
    shapley_mean            = round(long(st_sub$mean), 2),
    shapley_sd              = round(long(st_sub$sd), 2),
    shapley_q05             = round(long(st_sub$q05), 2),
    shapley_q95             = round(long(st_sub$q95), 2),
    shapley_pct             = round(long(st_sub$mean / obs_pop * 100), 4),
    shapley_check           = rep(round(rowSums(st_sub$mean), 2), each = n_sectors),   # = v(all sectors)
    shapley_direct_mean     = round(long(sp_sub$direct_mean), 2),
    shapley_direct_sd       = round(long(sp_sub$direct_sd), 2),
    shapley_vegetation_mean = round(long(sp_sub$veg_mean), 2),
    shapley_vegetation_sd   = round(long(sp_sub$veg_sd), 2)),
    extrap_flags, by = "subbasin")

  # ---- BCRs: subbasin samples summed within each BCR ----
  bcr_of    <- unname(bcr_lookup)
  phi_bcr   <- sum_units(phi, bcr_of)
  phi_d_bcr <- sum_units(phi_d, bcr_of)
  bcrs      <- rownames(phi_bcr)
  st_bcr    <- summarise_samples(phi_bcr)
  sp_bcr    <- split_samples(phi_bcr, phi_d_bcr, boot_of)
  tot_bcr   <- apply(phi_bcr, c(1, 3), sum)
  obs_bcr   <- vapply(bcrs, function(b) sum(obs_pop[bcr_of == b]), numeric(1L))
  # flagged subbasins, and the impact they carry
  imp_sub   <- rowMeans(tot_sub)
  imp_flg_b <- vapply(bcrs, function(b) sum(imp_sub[bcr_of == b & flagged %in% TRUE]), numeric(1L))
  shapley_bcr_rows[[sy]] <- data.frame(
    species                 = sp,
    bcr                     = rep(bcrs, each = n_sectors),
    year                    = yr,
    sector                  = rep(sectors, times = length(bcrs)),
    obs_population          = rep(round(obs_bcr), each = n_sectors),
    total_HF_impact         = rep(round(rowMeans(tot_bcr)), each = n_sectors),
    total_HF_impact_sd      = rep(round(apply(tot_bcr, 1, sd)), each = n_sectors),
    shapley_mean            = round(long(st_bcr$mean), 2),
    shapley_sd              = round(long(st_bcr$sd), 2),
    shapley_q05             = round(long(st_bcr$q05), 2),
    shapley_q95             = round(long(st_bcr$q95), 2),
    shapley_pct             = round(long(st_bcr$mean / obs_bcr * 100), 4),
    n_subbasins             = rep(vapply(bcrs, function(b) sum(bcr_of == b), integer(1L)), each = n_sectors),
    n_flagged               = rep(vapply(bcrs, function(b) sum(flagged[bcr_of == b], na.rm = TRUE), integer(1L)), each = n_sectors),
    total_HF_impact_flagged = rep(round(imp_flg_b), each = n_sectors),
    shapley_direct_mean     = round(long(sp_bcr$direct_mean), 2),
    shapley_direct_sd       = round(long(sp_bcr$direct_sd), 2),
    shapley_vegetation_mean = round(long(sp_bcr$veg_mean), 2),
    shapley_vegetation_sd   = round(long(sp_bcr$veg_sd), 2))

  # ---- national: BCR samples summed, bootstrap i paired across BCRs (as V5's national mosaic pairs them) ----
  phi_nat <- sum_units(phi_bcr, rep("national", length(bcrs)))
  st_nat  <- summarise_samples(phi_nat)
  sp_nat  <- split_samples(phi_nat, sum_units(phi_d_bcr, rep("national", length(bcrs))), boot_of)
  tot_nat <- apply(phi_nat, c(1, 3), sum)
  shapley_nat_rows[[sy]] <- data.frame(
    species                 = sp,
    year                    = yr,
    sector                  = sectors,
    obs_population          = round(sum(obs_pop)),
    total_HF_impact         = round(mean(tot_nat)),
    total_HF_impact_sd      = round(sd(tot_nat)),
    total_HF_impact_q05     = round(quantile(tot_nat, 0.05, names = FALSE)),
    total_HF_impact_q95     = round(quantile(tot_nat, 0.95, names = FALSE)),
    shapley_mean            = round(long(st_nat$mean), 2),
    shapley_sd              = round(long(st_nat$sd), 2),
    shapley_q05             = round(long(st_nat$q05), 2),
    shapley_q95             = round(long(st_nat$q95), 2),
    shapley_pct             = round(long(st_nat$mean / sum(obs_pop) * 100), 4),
    n_bcrs                  = length(bcrs),
    n_flagged               = sum(flagged, na.rm = TRUE),
    total_HF_impact_flagged = round(sum(imp_flg_b)),
    shapley_direct_mean     = round(long(sp_nat$direct_mean), 2),
    shapley_direct_sd       = round(long(sp_nat$direct_sd), 2),
    shapley_vegetation_mean = round(long(sp_nat$veg_mean), 2),
    shapley_vegetation_sd   = round(long(sp_nat$veg_sd), 2))

  message(sprintf("  sum of Shapley values %.3f = v(all sectors) %.3f", sum(st_nat$mean), mean(tot_nat)))
}

out_csv <- function(name) file.path(out_dir, paste0(name, ysfx, ".csv"))
write.csv(bind_rows(shapley_sub_rows), out_csv("shapley_subbasin"), row.names = FALSE)
write.csv(bind_rows(shapley_bcr_rows), out_csv("shapley_bcr"),      row.names = FALSE)
write.csv(bind_rows(shapley_nat_rows), out_csv("shapley_national"), row.names = FALSE)
message("\nwrote shapley_{subbasin,bcr,national}", ysfx, ".csv to ", out_dir)
