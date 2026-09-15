# Is StandardGreenup/Dormancy NA already in V5 SOURCE stacks, or introduced by 06?
# Also characterize WHERE the NA falls (vs a land-cover / water band in the same stack).
suppressMessages({library(terra)})
terra::terraOptions(progress = 0)
st <- "G:/Shared drives/BAM_NationalModels5/gis/stacks"
targets <- c("StandardGreenup_1km", "StandardDormancy_1km")

for (bcr in c("can11", "can80", "can60")) {
  f <- file.path(st, paste0(bcr, "_2020.tif"))
  if (!file.exists(f)) { cat(bcr, ": no 2020 stack\n"); next }
  r <- terra::rast(f)
  cat("\n==== ", bcr, " (", terra::nlyr(r), " bands) ====\n")
  present <- intersect(targets, names(r))
  cat("phenology bands present:", paste(present, collapse=", "), "\n")
  ntot <- terra::ncell(r)
  for (t in present) {
    na <- terra::global(is.na(r[[t]]), "sum", na.rm=TRUE)[1,1]
    cat(sprintf("  %-22s source NA frac (full grid): %.3f\n", t, na/ntot))
  }
  # in-data footprint defined by a near-always-present climate band
  anchor <- intersect(c("ERAMAT_1km","AHM_1km","mTPI_1km"), names(r))[1]
  if (!is.na(anchor)) {
    foot <- !is.na(r[[anchor]])
    foot_n <- terra::global(foot, "sum", na.rm=TRUE)[1,1]
    for (t in present) {
      na_in <- terra::global(is.na(r[[t]]) & foot, "sum", na.rm=TRUE)[1,1]
      cat(sprintf("  %-22s NA frac WITHIN %s footprint: %.3f\n", t, anchor, na_in/foot_n))
    }
  }
  # cross-tab phenology-NA against a land-cover band if available
  lc <- intersect(c("VLCE_1km","SCANFI_1km","ABoVE_1km","NLCD_1km"), names(r))[1]
  if (!is.na(lc) && length(present)) {
    t <- present[1]
    na_mask <- is.na(r[[t]])
    lcv <- terra::values(r[[lc]])[,1]
    nav <- terra::values(na_mask)[,1]
    tab <- table(landcover = lcv[nav %in% TRUE], useNA="ifany")
    cat("  land-cover (", lc, ") distribution AT", t, "NA pixels (top 8):\n")
    print(head(sort(tab, decreasing=TRUE), 8))
  }
}
