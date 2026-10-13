deploy_pbart_blank <- function(
    b,
    df_train,
    train_rows,
    df_backfill,
    predictors,
    biotic_cols,
    subbasin_index,
    logp
) {

  # Where is `b` blank? Some V5 covariates are missing by design on part of the land:
  # SCANFIheightcv_* has no value where canopy height is 0 (91-98% of the training land in the
  # southern can11 subbasins), and V5's bird models send that NA down gbm's missing branch,
  # which they learned at open survey sites. The amount models (gbart/mbart) learn only from
  # the rows where `b` has a value and predict a value everywhere, so without this model the
  # backfill is never blank: every backfilled pixel looks treed to the bird model, and to the
  # later BART models in the hierarchy, whose `<b>_isNA` flag (prep_bart_design) is 0 on every
  # backfilled pixel but 1 on most of the training land.
  # A probit BART on all the training land's pixels, with `b`'s own predictors, gives
  # P(blank) per backfilled pixel; 08A blanks the backfill where it exceeds 0.5. Returns NULL
  # when `b` is never blank on the training land, and 1 everywhere when it is always blank
  # (e.g. VLCE_1km in tundra subbasins, which VLCE does not map), so 12F passes it to the
  # bird model as V5 saw that land rather than dropping the pixels.
  miss <- as.integer(is.na(df_train[[b]][train_rows]))
  if (!length(miss) || all(miss == 0L)) return(NULL)
  if (all(miss == 1L)) {
    logp("[%s] blank on all training pixels: backfill blank everywhere", b)
    return(rep(1, nrow(df_backfill)))
  }

  # a single backfill pixel (subbasins 423, 247): BART drops a one-row x.test to a vector and
  # predicts one value per predictor, so take the training share, as gbart/mbart take the mean
  if (nrow(df_backfill) == 1L) {
    logp("[%s] single backfill pixel: P(blank) = training share %.2f", b, mean(miss))
    return(mean(miss))
  }

  design <- prep_bart_design(df_train[train_rows, predictors, drop = FALSE], df_backfill, logp)

  # reproducible per (year, subbasin, covariate), and apart from the amount model's seed
  set.seed(abs(as.integer(sprintf("%d%03d", subbasin_index, which(biotic_cols == b)))) + 1000000L)

  fit <- BART::pbart(x.train = as.matrix(design$train),
                     y.train = miss,
                     x.test  = as.matrix(design$backfill),
                     ntree = 50L,
                     ndpost = 500L,
                     nskip = 300L,
                     sparse = TRUE,
                     printevery = 1000L)

  p_blank <- fit$prob.test.mean
  logp("[%s] blank on %.1f%% of training pixels; modelled blank on %.1f%% of backfilled pixels",
       b, 100 * mean(miss), 100 * mean(p_blank > 0.5))
  p_blank

} # close function
