# ---
# title: Impact Assessment: fast, bit-identical gbm prediction for 12F
# author: Mannfred Boehm
# ---

# 98.7% of predict.gbm's time is its compiled tree walk, which makes three R API calls per node.
# gbm_tree_walk() (12G_gbm_tree_walk.cpp, compiled by 12D) does the same arithmetic in the same
# order, 3.6-4.4x faster. splitting design (gbm_design) from prediction lets 12F build the design
# matrix once per bootstrap and swap only the BART-draw columns per scenario.
# bit-identity with predict.gbm is checked in every 12F worker (gbm_check_fast)

# the numeric design matrix predict.gbm builds (gbm 2.1.x/2.2.x): the model's columns in
# var.names order, factors as 0-based codes in the model's own levels, the rest untouched
gbm_design <- function(model, newdata) {
  if (is.null(model$Terms))
    stop("gbm_design: model has no Terms, so predict.gbm would use newdata's own columns")
  x <- stats::model.frame(stats::terms(stats::reformulate(model$var.names)), newdata, na.action = stats::na.pass)
  for (i in seq_len(ncol(x))) if (is.factor(x[, i]))
    x[, i] <- as.numeric(factor(x[, i], levels = model$var.levels[[i]])) - 1
  matrix(as.double(unlist(x, use.names = FALSE)), nrow = nrow(x), dimnames = list(NULL, model$var.names))
}

# type = "response" prediction from a gbm_design() matrix, with all of the model's trees
gbm_predict_design <- function(model, D, block = 256L) {
  if (!is.null(model$num.classes) && model$num.classes != 1L) stop("gbm_predict_design: multi-class models are not supported")
  if (!identical(model$distribution$name, "poisson")) stop("gbm_predict_design: only poisson models are supported")
  if (nrow(D) == 0L) return(numeric(0))
  exp(gbm_tree_walk(D, model$trees, model$c.split, as.integer(model$var.type),
                    model$initF, as.integer(model$n.trees), as.integer(block)))
}

# stops unless the fast path equals predict.gbm bit for bit on the first k rows
gbm_check_fast <- function(model, newdata, D, k = 1000L, label = "") {
  k <- min(k, nrow(newdata))
  if (k == 0L) return(invisible(TRUE))
  ref  <- suppressWarnings(gbm::predict.gbm(model, newdata[seq_len(k), , drop = FALSE], n.trees = model$n.trees, type = "response"))
  fast <- gbm_predict_design(model, D[seq_len(k), , drop = FALSE])
  if (!identical(fast, ref))
    stop(label, " | gbm_tree_walk differs from predict.gbm on ", sum(fast != ref, na.rm = TRUE), " of ", k, " rows")
  invisible(TRUE)
}
