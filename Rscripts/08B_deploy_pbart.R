# ---
# title: Impact Assessment: probit BART for where a covariate is blank (sourced by 07)
# author: Mannfred Boehm
# ---

# some V5 covariates are blank by design on part of the land: SCANFIheightcv_* has no value
# where canopy height is 0 (91-98% of the training land in southern can11 subbasins), and V5's
# bird models send that NA down gbm's missing branch, learned at open survey sites. gbart and
# mbart learn only where b has a value, so without this model every backfilled pixel would look
# treed, to the bird model and to the later BART models (their `<b>_isNA` flag would be 0 on
# every backfilled pixel but 1 on most of the training land).
# P(blank) per backfilled pixel, from a probit BART on all training pixels with b's own
# predictors; 08A blanks the backfill where it exceeds 0.5. NULL when b is never blank on the
# training land; 1 everywhere when it always is (VLCE_1km in tundra subbasins, which VLCE does
# not map), so 12F gives the bird model what V5 saw there
deploy_pbart_blank <- function(b, train, predictors, df_backfill, seed, logp) {

  miss <- as.integer(is.na(train[[b]]))
  if (!length(miss) || all(miss == 0L)) return(NULL)
  if (all(miss == 1L)) {
    logp("[%s] blank on all training pixels: backfill blank everywhere", b)
    return(rep(1, nrow(df_backfill)))
  }

  # one backfilled pixel (subbasins 423, 247): BART reduces a one-row x.test to a vector, so
  # take the training share, as gbart and mbart take the training mean
  if (nrow(df_backfill) == 1L) {
    logp("[%s] single backfill pixel: P(blank) = training share %.2f", b, mean(miss))
    return(mean(miss))
  }

  design <- prep_bart_design(train[, predictors, drop = FALSE], df_backfill, logp)

  set.seed(seed + 1000000L)   # apart from the amount model's seed
  fit <- BART::pbart(x.train    = as.matrix(design$train),
                     y.train    = miss,
                     x.test     = as.matrix(design$backfill),
                     ntree      = 50L,
                     ndpost     = 500L,
                     nskip      = 300L,
                     sparse     = TRUE,
                     printevery = 1000L)

  p_blank <- fit$prob.test.mean
  logp("[%s] blank on %.1f%% of training pixels; modelled blank on %.1f%% of backfilled pixels",
       b, 100 * mean(miss), 100 * mean(p_blank > 0.5))
  p_blank
}
