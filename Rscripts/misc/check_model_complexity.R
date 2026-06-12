library(terra)

nm_root      <- "/home/mannfred/projects/def-ecknight/NationalModels"
ia_dir       <- "/home/mannfred/scratch/impact_assessment"
hirsh_dir    <- file.path(ia_dir, "data", "raw_data", "hirshpearson")
species      <- "CAWA"
coalition_id <- 2L
year         <- 2020L

# map coalition ID to sector names (mirrors 12E_shapley_utils.R::coalition_id_to_sectors)
sectors   <- sort(c("built", "crop", "dam_and_associated_reservoir",
                    "mines", "oil_gas", "pasture", "rail", "roads"))
bits      <- as.logical(intToBits(coalition_id - 1L)[seq_along(sectors)])
coalition <- sectors[bits]
cat(sprintf("Coalition %d: {%s}\n\n", coalition_id, paste(coalition, collapse = ", ")))

# load source rasters once; project to each BCR grid inside the loop
canHF_src <- terra::rast(file.path(hirsh_dir, "CanHF_1km_morethan1.tif"))
sec_srcs  <- lapply(coalition, function(s) terra::rast(file.path(hirsh_dir, paste0(s, ".tif"))))

rdata_files <- list.files(
  file.path(nm_root, "output/06_bootstraps", species),
  pattern = "can.*\\.Rdata$", full.names = TRUE)

cat(sprintf("Found %d BCR models for %s\n", length(rdata_files), species))
cat(sprintf("%-8s %8s %6s %6s %8s %12s %10s %12s\n",
            "BCR", "n.trees", "depth", "npred", "nTrain",
            "coalition_px", "pred_1(s)", "est_boot(min)"))
cat(strrep("-", 82), "\n")

results <- lapply(rdata_files, function(path) {
  bcr <- sub(paste0(".*", species, "_(.+)\\.Rdata$"), "\\1", path)

  # replicate 12C coalition mask to get exact pixel count for this BCR
  stack_ref  <- terra::rast(file.path(nm_root, "gis/stacks", paste0(bcr, "_", year, ".tif")))
  canHF_r    <- terra::project(canHF_src, stack_ref, method = "near")
  union_mask <- terra::rast(stack_ref[[1]]); terra::values(union_mask) <- 0L
  for (sec_r in sec_srcs) {
    sec_bcr    <- terra::project(sec_r, stack_ref, method = "near")
    union_mask <- terra::ifel(sec_bcr > 0, 1L, union_mask)
    rm(sec_bcr)
  }
  sector_mask <- terra::ifel((union_mask == 1L) & (canHF_r >= 1), 1, NA)
  n_px        <- as.integer(terra::global(sector_mask, "notNA")[[1]])
  rm(stack_ref, canHF_r, union_mask, sector_mask); gc()

  e <- new.env(parent = emptyenv())
  load(path, envir = e)
  b <- e$b.list[[1]]

  X_fake <- as.data.frame(
    matrix(0.5, nrow = max(n_px, 1L), ncol = length(b$var.names),
           dimnames = list(NULL, b$var.names))
  )
  t0           <- proc.time()
  suppressWarnings(gbm::predict.gbm(b, X_fake, n.trees = b$n.trees, type = "response"))
  elapsed      <- (proc.time() - t0)[["elapsed"]]
  est_boot_min <- elapsed * 100 / 60

  cat(sprintf("%-8s %8d %6d %6d %8d %12s %10.1f %12.1f\n",
              bcr, b$n.trees, b$interaction.depth, length(b$var.names), b$nTrain,
              format(n_px, big.mark = ","), elapsed, est_boot_min))

  out <- list(bcr = bcr, n.trees = b$n.trees, n_px = n_px,
              pred_sec = elapsed, est_boot_min = est_boot_min)
  rm(e, b, X_fake); gc()
  out
})

cat(strrep("-", 82), "\n")
total_est <- sum(sapply(results, `[[`, "est_boot_min"), na.rm = TRUE) * 32
cat(sprintf("\nEst. total wall time for CAWA (all BCRs, 32 boots): %.1f hours\n",
            total_est / 60))
