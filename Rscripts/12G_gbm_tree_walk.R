# ---
# title: Impact Assessment: fast, bit-identical gbm prediction for 12F
# author: Mannfred Boehm
# ---
# Sourced by 12D. The compiled walk itself, gbm_tree_walk(), comes from
# 12G_gbm_tree_walk.cpp, which 12D compiles with Rcpp::sourceCpp() before 12F runs.
#
# Why: 98.7% of predict.gbm's time is gbm_pred, gbm's compiled tree walk, and that walk
# makes three R API calls per node visited. gbm_tree_walk() does the same arithmetic in
# the same order, 3.6-4.4x faster (CAWA can11, 9150 trees). R-side prep is the other 1.3%,
# which is why a rewrite in another language would not help: the cost is the tree walk.
#
# The split into gbm_design() + gbm_predict_design() also lets 12F build the design matrix
# ONCE per bootstrap and overwrite only the BART-draw columns per scenario.
#
# Bit-identity with predict.gbm is an invariant, not an aspiration: 12F calls
# gbm_check_fast() in every worker and stop()s on any difference.

# Numeric design matrix exactly as predict.gbm builds it (gbm 2.1.x/2.2.x): the model's
# columns in var.names order via model.frame, factors recoded to 0-based codes in the
# model's own level order, everything else passed through untouched.
gbm_design <- function(model, newdata) {
  if (is.null(model$Terms))
    stop("gbm_design: model has no Terms, so predict.gbm would use newdata's own columns",
         " - use gbm::predict.gbm for this model")
  x <- stats::model.frame(stats::terms(stats::reformulate(model$var.names)), newdata,
                          na.action = stats::na.pass)
  for (i in seq_len(ncol(x))) if (is.factor(x[, i]))
    x[, i] <- as.numeric(factor(x[, i], levels = model$var.levels[[i]])) - 1
  matrix(as.double(unlist(x, use.names = FALSE)), nrow = nrow(x),
         dimnames = list(NULL, model$var.names))
}

# type = "response" prediction from a gbm_design() matrix, using all model$n.trees trees.
gbm_predict_design <- function(model, D, block = 256L) {
  if (!is.null(model$num.classes) && model$num.classes != 1L)
    stop("gbm_predict_design: multi-class models are not supported")
  if (!identical(model$distribution$name, "poisson"))
    stop("gbm_predict_design: only poisson models are supported (got ",
         model$distribution$name, ")")
  if (nrow(D) == 0L) return(numeric(0))
  exp(gbm_tree_walk(D, model$trees, model$c.split, as.integer(model$var.type),
                    model$initF, as.integer(model$n.trees), as.integer(block)))
}

# The fast path must equal predict.gbm bit for bit. Checked on the first k rows.
gbm_check_fast <- function(model, newdata, D, k = 1000L, label = "") {
  k <- min(k, nrow(newdata))
  if (k == 0L) return(invisible(TRUE))
  ref  <- suppressWarnings(gbm::predict.gbm(model, newdata[seq_len(k), , drop = FALSE],
                                            n.trees = model$n.trees, type = "response"))
  fast <- gbm_predict_design(model, D[seq_len(k), , drop = FALSE])
  if (!identical(fast, ref))
    stop(label, " | gbm_tree_walk differs from gbm::predict.gbm on ", sum(fast != ref, na.rm = TRUE),
         " of ", k, " rows (max rel diff ",
         signif(max(abs(fast - ref) / abs(ref), na.rm = TRUE), 3), ") - the fast path is not ",
         "bit-identical on this platform; do not use it.")
  invisible(TRUE)
}
