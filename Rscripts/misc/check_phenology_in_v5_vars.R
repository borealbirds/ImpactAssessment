# Are StandardGreenup_1km / StandardDormancy_1km in the V5 bird-model var.names
# (= 12C model_vars_shared, which drives the complete.cases drop)? Inspect a few
# CAWA + OVEN bootstrap models across BCRs.
bs_dir <- "G:/Shared drives/BAM_NationalModels5/output/06_bootstraps"
targets <- c("StandardGreenup_1km", "StandardDormancy_1km")

inspect <- function(rdata_path) {
  e <- new.env()
  load(rdata_path, envir = e)
  # find the object holding a gbm-like model with var.names
  vn <- NULL
  for (nm in ls(e)) {
    obj <- get(nm, envir = e)
    if (is.list(obj) && !is.null(obj$var.names)) { vn <- obj$var.names; break }
    # bootstrap list: list of models
    if (is.list(obj) && length(obj) && is.list(obj[[1]]) && !is.null(obj[[1]]$var.names)) {
      vn <- obj[[1]]$var.names; break
    }
  }
  list(objs = ls(e), var.names = vn)
}

files <- c(file.path(bs_dir, "CAWA", c("CAWA_can60.Rdata", "CAWA_can11.Rdata")),
           file.path(bs_dir, "OVEN", c("OVEN_can60.Rdata")))
files <- files[file.exists(files)]

for (f in files) {
  cat("\n==== ", basename(f), " ====\n")
  r <- tryCatch(inspect(f), error = function(e) { cat("ERR:", conditionMessage(e), "\n"); NULL })
  if (is.null(r)) next
  cat("objects in file:", paste(r$objs, collapse = ", "), "\n")
  if (is.null(r$var.names)) { cat("NO var.names found\n"); next }
  cat("n var.names:", length(r$var.names), "\n")
  for (t in targets) cat(sprintf("  %-22s in var.names: %s\n", t, t %in% r$var.names))
  # show any phenology-ish names
  ph <- grep("Green|Dorman|phen|Phen", r$var.names, value = TRUE)
  cat("phenology-like var.names:", if (length(ph)) paste(ph, collapse=", ") else "(none)", "\n")
}
