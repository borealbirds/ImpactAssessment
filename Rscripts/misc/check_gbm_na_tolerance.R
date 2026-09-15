# Crux test for the "relax 12C gate" fix: does gbm::predict.gbm return FINITE
# predictions when StandardGreenup_1km / StandardDormancy_1km are NA?
suppressMessages({library(terra); library(gbm)})
terra::terraOptions(progress = 0)

st  <- "G:/Shared drives/BAM_NationalModels5/gis/stacks"
bs  <- "G:/Shared drives/BAM_NationalModels5/output/06_bootstraps/CAWA/CAWA_can80.Rdata"
e <- new.env(); load(bs, envir = e)
b.list <- get("b.list", envir = e)
model  <- b.list[[1]]
vn <- model$var.names
cat("model has", length(vn), "vars; greenup in:", "StandardGreenup_1km" %in% vn,
    "dormancy in:", "StandardDormancy_1km" %in% vn, "\n")

# build a design matrix from the source stack: take complete-case rows, then
# blank out the phenology cols to NA and compare predictions.
src <- terra::rast(file.path(st, "can80_2020.tif"))
keep_vars <- intersect(vn, names(src))
X <- as.data.frame(terra::values(src[[keep_vars]]), check.names = FALSE)
# categorical vars -> factor with the model's levels (so gbm doesn't choke on those)
cat_vars <- intersect(c("ABoVE_1km","NLCD_1km","MODISLCC_1km","MODISLCC_5x5","SCANFI_1km","VLCE_1km"), keep_vars)
for (v in cat_vars) {
  lv <- model$var.levels[[match(v, vn)]]
  X[[v]] <- factor(as.character(X[[v]]), levels = as.character(lv))
}
# rows complete on ALL model vars (so the only thing we toggle is phenology)
cc <- stats::complete.cases(X[, keep_vars, drop = FALSE])
Xc <- X[cc, , drop = FALSE]
set.seed(1); Xc <- Xc[sample(nrow(Xc), min(5000, nrow(Xc))), , drop = FALSE]
cat("test rows (fully complete):", nrow(Xc), "\n")

p_full <- gbm::predict.gbm(model, Xc, n.trees = model$n.trees, type = "response")

Xna <- Xc
for (v in intersect(c("StandardGreenup_1km","StandardDormancy_1km"), names(Xna))) Xna[[v]] <- NA_real_
p_na <- gbm::predict.gbm(model, Xna, n.trees = model$n.trees, type = "response")

cat(sprintf("p_full finite: %d/%d   p_na finite: %d/%d\n",
            sum(is.finite(p_full)), length(p_full), sum(is.finite(p_na)), length(p_na)))
cat("summary p_full:\n"); print(summary(p_full))
cat("summary p_na (phenology forced NA):\n"); print(summary(p_na))
cat(sprintf("mean abs rel diff where both finite: %.4f\n",
            mean(abs((p_na - p_full)/pmax(p_full,1e-9))[is.finite(p_full)&is.finite(p_na)])))
