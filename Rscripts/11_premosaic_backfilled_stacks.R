# ---
# title: Impact Assessment: pre-mosaic backfilled subbasin stacks per BCR
# author: Mannfred Boehm
# created: March 23, 2026
# ---

# one BCR per task: the backfills (07) of every subbasin touching the BCR, mosaicked onto the
# BCR's covariate grid and masked to its polygon, so 12F reads one file per BCR
#   -> data/derived_data/bart_models_mosaics/{year}/{bcr}_backfilled.tif
# chained after 07 by 07_train_and_backfill.sh. MOSAIC_DIR=<dir> redirects the output (tests)

library(terra)

# terra sizes its blocks from the node's free RAM and can't see the Slurm allocation: on a
# 750 GB node it held whole BCR stacks in memory (can3 peaked at 316 GB), which a smaller
# request would kill. capped at 30% of the allocation (resample, merge and mask give the same
# cells in blocks)
slurm_mem_mb <- suppressWarnings(as.numeric(Sys.getenv("SLURM_MEM_PER_NODE")))
if (!is.na(slurm_mem_mb)) terra::terraOptions(memmax = 0.3 * slurm_mem_mb / 1024)

nm_root <- "/home/mannfred/projects/def-ecknight/NationalModels"
ia_dir  <- if (nzchar(Sys.getenv("SLURM_JOB_ID"))) "/home/mannfred/scratch/impact_assessment" else getwd()

bam_boundary         <- terra::vect(file.path(ia_dir, "data", "raw_data", "Regions", "BAM_BCR_NationalModel_Unbuffered.shp"))
all_subbasins_subset <- terra::vect(file.path(ia_dir, "data", "raw_data", "hydrobasins_masked_merged_subset.gpkg"))
bam_bcr_codes        <- paste0(bam_boundary$country, bam_boundary$subUnit)

# subbasin x BCR pairs that intersect (a subbasin can be in several BCRs)
ij <- which(terra::relate(all_subbasins_subset, bam_boundary, relation = "intersects"), arr.ind = TRUE)
bcr_subbasins_ref <- data.frame(sub_index = ij[, 1], bcr_code = bam_bcr_codes[ij[, 2]])

# this task's BCR: its position among the sorted Canadian BCRs (19, as 07's .sh counts them)
bcr_vec  <- sort(unique(bcr_subbasins_ref$bcr_code))
bcr_vec  <- bcr_vec[startsWith(bcr_vec, "can")]
bcr_code <- bcr_vec[as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))]
year     <- as.integer(Sys.getenv("YEARS", "2020"))   # 07's .sh passes one year per job
sub_ids  <- unique(bcr_subbasins_ref$sub_index[bcr_subbasins_ref$bcr_code == bcr_code])
message(Sys.time(), " | year ", year, " | ", bcr_code, " | subbasins=", length(sub_ids))

out_dir  <- file.path(ia_dir, "data", "derived_data", Sys.getenv("MOSAIC_DIR", "bart_models_mosaics"), year)
out_path <- file.path(out_dir, paste0(bcr_code, "_backfilled.tif"))
in_files <- file.path(ia_dir, "data", "derived_data", "bart_models", year,
                      paste0("subbasin_", sub_ids), paste0("subbasin_", sub_ids, "_backfill.tif"))
in_files <- in_files[file.exists(in_files)]

# skip only a mosaic newer than every backfill it is built from (skipping on existence alone
# once turned a whole rerun into a no-op, 2026-09-19)
if (file.exists(out_path) && length(in_files) > 0L && file.mtime(out_path) > max(file.mtime(in_files))) {
  message(Sys.time(), " | ", bcr_code, " | mosaic is newer than its backfills - skipping")
  quit(save = "no", status = 0)
}
if (length(in_files) == 0L) {
  message(Sys.time(), " | ", bcr_code, " | no subbasin backfills - skipping")
  quit(save = "no", status = 0)
}



# ----------------------------------------------
# mosaic onto the BCR's covariate grid, one covariate at a time, on disk

# the BCR's observed stack, whose grid the mosaic takes (12F reads both)
ref1    <- terra::rast(file.path(nm_root, "gis/stacks", paste0(bcr_code, "_", year, ".tif")))[[1]]
ref_ext <- terra::ext(ref1)

# the backfills that overlap the BCR grid
overlaps <- vapply(in_files, function(p) {
  e <- terra::ext(terra::rast(p))
  e$xmin < ref_ext$xmax && e$xmax > ref_ext$xmin && e$ymin < ref_ext$ymax && e$ymax > ref_ext$ymin
}, logical(1))
in_files <- in_files[overlaps]

# each subbasin's whole stack, resampled once onto its own window of the BCR grid. a window
# snapped to the grid gives the same cells as the whole grid ("near" takes the source cell
# under each target cell's centre), at a fraction of the cost: whole-grid resamples took 5 h
# for can81. INTERLEAVE=BAND, because the loop below reads one band at a time: GDAL's default
# pixel interleave made each band read scan the whole file (~23 min per covariate)
gtiff_opts <- c("INTERLEAVE=BAND", "COMPRESS=DEFLATE", "BIGTIFF=YES", "TILED=NO")
rs_dir <- tempfile("premosaic_rs_")
dir.create(rs_dir, recursive = TRUE)
sub_handles <- lapply(seq_along(in_files), function(k) {
  r     <- terra::rast(in_files[k])
  e_int <- terra::intersect(terra::ext(r), ref_ext)
  pth   <- file.path(rs_dir, sprintf("sub_%04d.tif", k))
  terra::resample(terra::crop(r, e_int), terra::crop(ref1, e_int, snap = "out"), method = "near",
                  filename = pth, overwrite = TRUE, wopt = list(gdal = gtiff_opts))
  message(Sys.time(), " | pre-resampled subbasin ", k, "/", length(in_files))
  terra::rast(pth)
})
names_per <- lapply(sub_handles, names)
vars      <- sort(unique(unlist(names_per)))

# per covariate: merge the subbasin windows (merge() keeps the first non-NA value, in subbasin
# order), extend to the BCR grid and write to its own temp file, so nothing piles up in RAM
layer_dir   <- tempfile("premosaic_layers_")
dir.create(layer_dir, recursive = TRUE)
layer_paths <- file.path(layer_dir, sprintf("layer_%05d.tif", seq_along(vars)))
for (k in seq_along(vars)) {
  layers <- lapply(which(vapply(names_per, function(n) vars[k] %in% n, logical(1))),
                   function(i) sub_handles[[i]][[vars[k]]])
  acc <- if (length(layers) == 1L) layers[[1]] else terra::merge(terra::sprc(layers))
  acc <- terra::extend(acc, ref1)
  if (!terra::compareGeom(acc, ref1, stopOnError = FALSE)) stop("merged layer ", vars[k], " is not on the BCR grid")
  terra::writeRaster(acc, layer_paths[k], overwrite = TRUE)
  rm(acc)
  gc(verbose = FALSE)
  if (k %% 50 == 0 || k == length(vars)) message(Sys.time(), " | covariate ", k, "/", length(vars))
}
stack_bf <- terra::rast(layer_paths)
names(stack_bf) <- vars



# ----------------------------------------------
# mask to the BCR polygon and write, in one streaming pass (masking, then writing, held all
# ~2500 layers in RAM and was killed at 128G). written under a temporary name and renamed when
# complete, so a killed task leaves no truncated mosaic. INTERLEAVE=BAND for 12F's band reads

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
part_path <- file.path(out_dir, paste0(bcr_code, "_backfilled_partial.tif"))
terra::mask(stack_bf, bam_boundary[bam_bcr_codes == bcr_code, ],
            filename = part_path, overwrite = TRUE, wopt = list(gdal = gtiff_opts))
if (!file.rename(part_path, out_path)) stop("could not move ", part_path, " into place")
message(Sys.time(), " | ", bcr_code, " | wrote ", out_path)
