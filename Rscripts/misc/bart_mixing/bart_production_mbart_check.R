# One-off: do 08B's mbart2 settings mix? Two chains at 08B's settings on SCANFI_1km in two small
# subbasins, full predictors, x.test = 200 backfill pixels with observed earlier covariates.
suppressPackageStartupMessages({ library(BART); library(terra); library(dplyr); library(posterior) })
ia_dir <- getwd(); year <- 2020L
out <- commandArgs(TRUE)[1]
say <- function(...) message(format(Sys.time(), "[%H:%M:%S] "), sprintf(...))
logp <- function(...) invisible(NULL)
source(file.path(ia_dir, "Rscripts", "08A_train_and_backfill_subbasin_s.R"))

meta <- dplyr::filter(BAMexploreR::predictor_metadata, version == "v5") |>
  dplyr::select(predictor, predictor_class) |>
  dplyr::mutate(predictor = sub("Method", "method", sub("Year", "year", predictor)))
soil <- paste0(rep(c("cec", "soc"), each = 6), "_", c("0-5", "100-200", "15-30", "30-60", "5-15", "60-100"), "cm_mean_1000")
wet <- c("Peatland_5x5", "Peatland_1km")
abiotic <- meta |> dplyr::filter(predictor_class %in% c("Annual Climate", "Climate Normals", "Topography", "Wetland", "Disturbance")) |>
  tibble::add_row(predictor = "CAfire", predictor_class = "Time Since Disturbance") |>
  dplyr::filter(!(predictor %in% wet)) |> dplyr::bind_rows(tibble::tibble(predictor = soil, predictor_class = "Soil Properties"))
biotic <- c(meta$predictor[!(meta$predictor_class %in% c(abiotic$predictor_class, "Time", "Method"))], wet)
neworder <- readRDS(file.path(ia_dir, "data", "raw_data", "biotic_variable_hierarchy.rds"))
categorical <- c("ABoVE_1km", "NLCD_1km", "MODISLCC_1km", "MODISLCC_5x5", "SCANFI_1km", "VLCE_1km")

stack_y <- rast(file.path(ia_dir, "data", "raw_data", "covariates_mosaiced", sprintf("covariates_mosaiced_%d.tif", year)))
subbasins <- project(vect(file.path(ia_dir, "data", "raw_data", "hydrobasins_masked_merged_subset.gpkg")), crs(stack_y))

rows <- list()
for (s in c(9L, 645L)) {
  poly <- subbasins[s]
  cov_s <- mask(crop(stack_y, poly), poly)
  hf <- function(kind) { r <- rast(file.path(ia_dir, "data", "raw_data", "hirshpearson", sprintf("CanHF_1km_%s.tif", kind)))
    mask(project(crop(r, project(poly, crs(r))), cov_s, method = "near"), poly) }
  low <- which(values(hf("lessthan1"), mat = FALSE) == 1); high <- which(values(hf("morethan1"), mat = FALSE) == 1)
  df <- as.data.frame(cov_s, xy = TRUE, na.rm = FALSE)
  ok <- which(!is.na(df$x)); ctr <- colMeans(df[ok, c("x", "y")]); scl <- apply(df[ok, c("x", "y")], 2, sd)
  df[, c("x", "y")] <- sweep(sweep(df[, c("x", "y")], 2, ctr, "-"), 2, scl, "/")
  train <- df[low, ]; fill <- df[high, ]
  fp <- intersect(names(fill), abiotic$predictor[abiotic$predictor_class == "Disturbance"]); fill[, fp] <- 0
  ab <- intersect(names(train), abiotic$predictor); bio <- intersect(names(train), biotic)
  bio <- as.character(na.omit(bio[match(neworder, bio)])); cont <- setdiff(bio, categorical)
  b <- "SCANFI_1km"
  idx <- which(!is.na(train[[b]]) & !(train[[b]] %in% 8))
  set.seed(s); fs <- fill[sample(nrow(fill), min(200L, nrow(fill))), , drop = FALSE]
  d <- prep_bart_design(train[idx, c(ab, cont, "x", "y"), drop = FALSE], fs, logp)
  codes <- train[[b]][idx]; present <- sort(unique(codes)); y <- match(codes, present)
  set.seed(s + 1L); ho <- sample(seq_along(y), round(0.1 * length(y)))
  keep <- setdiff(seq_along(y), ho); present2 <- sort(unique(y[keep])); ytr <- match(y[keep], present2)
  say("subbasin %d: n %d, K %d, p %d", s, length(ytr), length(present2), ncol(d$train))
  shares <- list(); pix <- list()
  for (chain in 1:3) {
    set.seed(1000L * s + chain); t0 <- proc.time()[[3]]
    capture.output(fit <- suppressWarnings(BART::mbart2(x.train = as.matrix(d$train)[keep, , drop = FALSE], y.train = ytr,
      x.test = as.matrix(d$backfill), type = "pbart", k = 3, ntree = 40L, ndpost = 500L, nskip = 150L,
      keepevery = 10L, printevery = 100000L, sparse = TRUE)))
    K <- fit$K; n_px <- nrow(d$backfill)
    P <- pnorm(fit$yhat.test)                       # [500 x (pixel-major n_px * K)]
    arr <- array(t(P), dim = c(K, n_px, nrow(P)))   # class, pixel, draw
    tot <- apply(arr, c(2, 3), sum)
    arr <- sweep(arr, c(2, 3), tot, "/")
    shares[[chain]] <- t(apply(arr, c(1, 3), mean)) # draw x class
    pix[[chain]] <- arr
    say("subbasin %d chain %d: %.0f s", s, chain, proc.time()[[3]] - t0)
  }
  nd <- 500L; prod <- 1:nd
  r_share <- sapply(seq_len(K), function(k) rhat(sapply(shares, function(m) m[prod, k])))
  e_share <- sapply(seq_len(K), function(k) ess_bulk(sapply(shares, function(m) m[prod, k])))
  big <- which(colMeans(shares[[1]]) > 0.05)
  js <- sample(n_px, min(60L, n_px))
  r_pix <- unlist(lapply(big, function(k) sapply(js, function(j) rhat(sapply(pix, function(a) a[k, j, prod])))))
  # first third against last third: is burn-in 150 enough for the class shares?
  early <- 1:167; late <- 334:500
  pool <- function(r) do.call(rbind, lapply(shares, function(m) m[r, , drop = FALSE]))
  shift <- (colMeans(pool(early)) - colMeans(pool(late))) / apply(pool(late), 2, sd)
  # between-chain effective draws of the mean class share and of pixel probability
  mns <- lapply(shares, colMeans); pv <- apply(pool(prod), 2, var)
  sq <- Reduce(`+`, lapply(combn(3, 2, simplify = FALSE), function(p) (mns[[p[1]]] - mns[[p[2]]])^2)) / 3
  rows[[length(rows) + 1L]] <- tibble(subbasin = s, n = length(ytr), K = K, classes_over_5pct = length(big),
    rhat_share_median = median(r_share[big]), rhat_share_max = max(r_share[big]),
    rhat_pixel_median = median(r_pix, na.rm = TRUE), rhat_pixel_q95 = unname(quantile(r_pix, .95, na.rm = TRUE)),
    ess_share_median = median(e_share[big]), between_chain_ess_share = median((pv / (sq / 2))[big]),
    abs_shift_early_late_sd = median(abs(shift[big])))
  print(as.data.frame(rows[[length(rows)]]), digits = 3)
  write.csv(bind_rows(rows), out, row.names = FALSE)
}
say("finished")
