# Fast local de-risk: surgically replace ONLY the CAfire band in the existing
# covariates_mosaiced_2020.tif with the recoded (NA->0) CAfire, replicating
# 06_build_covariate_stacks.R:159-175. All other bands are unchanged, so this is
# equivalent to a full 06 rebuild for the purpose of the cascade audit.
suppressMessages({library(terra)})
terra::terraOptions(progress = 0)

root   <- "G:/Shared drives/BAM_NationalModels5"
sbx    <- file.path(root, "data", "Extras", "sandbox_data", "impactassessment_sandbox")

# local RProject mosaic (this is what the audit reads with DIAG_LOCAL=TRUE)
mos_path <- "C:/Users/mannf/Drive/boreal_avian_modelling_project/ImpactAssessment/data/raw_data/covariates_mosaiced/covariates_mosaiced_2020.tif"
caf_path <- file.path(sbx, "CAfire", "CAfire_2020_masked.tif")   # already recoded NA->0
stopifnot(file.exists(mos_path), file.exists(caf_path))

r_mos <- terra::rast(mos_path)
cat("mosaic bands:", terra::nlyr(r_mos), " has CAfire:", ("CAfire" %in% names(r_mos)), "\n")

caf <- terra::rast(caf_path)
if (is.na(terra::crs(caf)) || terra::crs(caf) == "") terra::crs(caf) <- terra::crs(r_mos)

caf_aligned <- terra::project(caf, r_mos, method = "bilinear",
                              filename = tempfile(fileext = ".tif"), overwrite = TRUE)
names(caf_aligned) <- "CAfire"
caf_aligned <- terra::mask(caf_aligned, r_mos[[1]])   # match mosaic footprint

# NA frac of CAfire band over the mosaic footprint, before vs after
foot <- !is.na(r_mos[["CAfire"]] * 0 + 1)   # in-mosaic cells (CAfire may be all-NA-ish)
foot_n <- terra::global(!is.na(r_mos[[1]]), "sum", na.rm = TRUE)[1,1]
na_before <- terra::global(is.na(r_mos[["CAfire"]]) & !is.na(r_mos[[1]]), "sum", na.rm=TRUE)[1,1]
na_after  <- terra::global(is.na(caf_aligned)       & !is.na(r_mos[[1]]), "sum", na.rm=TRUE)[1,1]
cat(sprintf("CAfire NA within mosaic footprint: BEFORE %d/%d (%.3f)  AFTER %d/%d (%.3f)\n",
            na_before, foot_n, na_before/foot_n, na_after, foot_n, na_after/foot_n))

# rebuild stack with CAfire swapped, write in place (original preserved in sandbox root)
r_new <- r_mos[[setdiff(names(r_mos), "CAfire")]]
r_new <- c(r_new, caf_aligned)
# keep CAfire in its original position is not required (audit indexes by name)

tmp_out <- sub("\\.tif$", "_cafirefix.tif", mos_path)
terra::writeRaster(r_new, tmp_out, overwrite = TRUE,
                   wopt = list(gdal = c("COMPRESS=LZW", "ZLEVEL=9")))
cat("wrote fixed mosaic ->", tmp_out, "\n")
