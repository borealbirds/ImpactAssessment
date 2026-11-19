# ---
# title: Impact Assessment: multinomial BART for one land-cover covariate (sourced by 07)
# author: Mannfred Boehm
# ---

# 08A calls this only when the training pixels hold at least two classes. returns the modal
# class (kept for the metrics), 100 drawn classes per pixel as layers <b>_class_001..100
# (what 12F uses), the metrics, and the holdout's actual vs predicted classes
deploy_mbart <- function(b, y, design, seed, subbasin_index, year, logp) {

  n_bf       <- nrow(design$backfill)
  class_lyrs <- paste0(b, "_class_", sprintf("%03d", seq_len(100L)))

  # mbart2 needs the classes numbered 1..K, so map the codes present to 1..K and back
  y_codes <- as.character(y)
  present <- sort(unique(y_codes), na.last = NA)
  present <- present[order(as.numeric(present))]
  inv_map <- setNames(present, seq_along(present))
  y_int   <- unname(setNames(seq_along(present), present)[y_codes])

  # no model: the modal class, and each pixel's classes drawn from the training shares
  from_training_shares <- function(n, logmsg) {
    tab   <- table(y_int)
    probs <- tab / sum(tab)
    logp("[%s] %s", b, logmsg)
    list(mode    = inv_map[names(which.max(tab))],
         maxprob = max(probs),
         entropy = -sum(probs * log(pmax(probs, 1e-12))),
         draw    = function() unname(inv_map[sample(names(tab), n, replace = TRUE, prob = probs)]))
  }

  # one backfilled pixel
  if (n_bf == 1) {
    s <- from_training_shares(1L, "single backfill pixel: skipped mbart()")
    set.seed(seed)
    layers <- list(s$mode, s$maxprob, s$entropy)
    names(layers) <- c(b, paste0(b, c("_maxprob", "_entropy")))
    return(list(backfill = s$mode, layers = c(layers, setNames(replicate(100L, s$draw(), simplify = FALSE), class_lyrs))))
  }

  set.seed(seed)
  holdout_idx <- sample(seq_len(length(y_int)), size = round(0.10 * length(y_int)))
  df_holdout  <- design$train[holdout_idx, ]
  df_train    <- design$train[-holdout_idx, ]
  y_train_int <- y_int[-holdout_idx]

  # one class left in the training split
  if (length(unique(y_train_int)) < 2) {
    s <- from_training_shares(n_bf, "degenerate after holdout split: single class")
    layers <- list(s$mode, rep(s$maxprob, n_bf), rep(s$entropy, n_bf))
    names(layers) <- c(b, paste0(b, c("_maxprob", "_entropy")))
    return(list(backfill = s$mode, layers = c(layers, setNames(replicate(100L, s$draw(), simplify = FALSE), class_lyrs))))
  }

  # mbart2 fits one probit model per class in its training rows, numbered 1..K: a class that
  # fell only in the holdout would leave a gap, so renumber over the training split. the
  # holdout keeps its true codes for the metrics
  holdout_codes <- present[y_int[holdout_idx]]
  train_codes   <- present[y_train_int]
  present       <- present[sort(unique(y_train_int))]
  y_train_int   <- match(train_codes, present)

  # latent probit draws are autocorrelated, hence keepevery = 10 (gbart keeps every draw)
  fit <- BART::mbart2(x.train    = as.matrix(df_train),
                      y.train    = y_train_int,
                      x.test     = as.matrix(design$backfill),
                      type       = "pbart",
                      k          = 3,
                      ntree      = 40L,
                      ndpost     = 500L,
                      nskip      = 150L,
                      keepevery  = 10L,
                      printevery = 350,
                      sparse     = TRUE)
  logp("mbart fit to [%s]", b)

  # class probabilities, [pixel x class]
  classes   <- present[fit$cats]
  prob_ecol <- matrix(fit$prob.test.mean, ncol = fit$K, byrow = TRUE)
  layers <- list(classes[max.col(prob_ecol, ties.method = "first")],
                 apply(prob_ecol, 1, max),
                 -rowSums(prob_ecol * log(pmax(prob_ecol, 1e-12))))
  names(layers) <- c(b, paste0(b, c("_maxprob", "_entropy")))

  # drawn classes, what 12F uses: in each of 100 posterior draws, each pixel's class is drawn
  # from that draw's class probabilities. the modal class would give every pixel its likeliest
  # class and erase minority classes (BOBO can13: mixed forest 65% vs 29% on nearby training
  # land, 16C, 2026-10-09). the probabilities are pnorm(z_k) normalised to sum to 1, not
  # mbart2's prob.test, a softmax meant for logits that roughly doubles rare classes on the
  # probit scale (can13 subbasin 98 holdout: herbs 3.1% observed, 6.1% softmax, 3.2% probit)
  draw_idx <- sample(nrow(fit$yhat.test), 100L)
  for (d in seq_along(draw_idx)) {
    P  <- stats::pnorm(matrix(fit$yhat.test[draw_idx[d], ], ncol = fit$K, byrow = TRUE))
    cp <- P
    for (j in seq_len(fit$K)[-1L]) cp[, j] <- cp[, j - 1L] + P[, j]
    u  <- stats::runif(nrow(P)) * cp[, fit$K]
    layers[[class_lyrs[d]]] <- classes[pmin(1L + rowSums(cp < u), fit$K)]
  }

  top_var <- names(sort(colMeans(do.call(rbind, fit$varprob)), decreasing = TRUE))[1]
  metrics_train <- collect_metrics_mbart(fit, df_train, train_codes, b, subbasin_index, year, top_var, present)

  # holdout, on the predictors the fit kept
  pred <- predict(fit, newdata = as.matrix(df_holdout[, colnames(fit$varprob[[1]]), drop = FALSE]))
  metrics_holdout <- collect_holdout_metrics_mbart(pred, holdout_codes, b, subbasin_index, year, top_var, present)
  prob_holdout <- matrix(pred$prob.test.mean, ncol = fit$K, byrow = TRUE)

  list(backfill  = layers[[b]],
       layers    = layers,
       metrics   = list(metrics_train, metrics_holdout),
       confusion = list(data.frame(covariate = b,
                                   subbasin  = subbasin_index,
                                   actual    = holdout_codes,
                                   predicted = classes[max.col(prob_holdout, ties.method = "first")],
                                   stringsAsFactors = FALSE)))
}
