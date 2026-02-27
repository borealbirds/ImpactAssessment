deploy_mbart <- function(
    b,
    df_train,
    df_train_bart,
    df_backfill,
    df_backfill_bart,
    idx,
    biotic_cols,
    subbasin_index,
    year,
    metrics,
    confusion,
    out_layers,
    logp,
    cc,
    ia_dir
) {

  # source BART metrics summary functions
  source(file.path(ia_dir, "Rscripts", "09_collect_metrics_mbart.R"))
  source(file.path(ia_dir, "Rscripts", "09_collect_holdout_metrics_mbart.R"))
  
  # BART::mbart needs the landcover classes as a consecutive sequence of integers 1..K' 
  # but actual classes could be any arbitrary non-consecutive sequence of integers 1..K
  # so, we convert 1..K -> 1..K' -> train -> backfill as 1..K' -> convert back to 1..K
  
  # the vector of original landcover codes for the selected training pixels. `idx` ensures no NAs present
  # don't use df_train_bart because it's missing `b`. 
  y_codes <- as.character(df_train[[b]][idx]) 
  
  # make sure y isn’t empty 
  if (!length(y_codes)) { logp("[%s] skip: empty y_codes", b); next } 
  
  # land cover classes that are *actually* present in the training subset
  present <- sort(unique(y_codes), na.last = NA)
  present <- present[order(as.numeric(present))]
  K <- length(present) 
  if (K == 0) { logp("[%s] skip: no classes present", b); next }
  
  # pre-allocate outputs
  b_mode    <- rep(NA_real_, nrow(df_backfill_bart))
  b_maxprob <- rep(NA_real_, nrow(df_backfill_bart))
  b_entropy <- rep(NA_real_, nrow(df_backfill_bart))
  
  # build a forward and inverse map only over present classes (size K')
  fwd_map  <- setNames(seq_along(present), present)   # original code -> 1..K'
  inv_map  <- setNames(present, seq_along(present))   # 1..K' -> original code
  
  # determine complexity of landcover classes
  if (K < 2) {
    
    # trivial: only one class in training data
    b_mode    <- rep(present, nrow(df_backfill_bart)) # the one and only landcover class present
    b_maxprob <- rep(1, nrow(df_backfill_bart)) 
    b_entropy <- rep(0, nrow(df_backfill_bart))
    
    single_code_num <- as.numeric(present)
    df_backfill[[b]] <- rep(present, nrow(df_backfill_bart))
    out_layers[[b]] <- rep(present, nrow(df_backfill_bart))
    out_layers[[paste0(b, "_maxprob")]] <- b_maxprob
    out_layers[[paste0(b, "_entropy")]] <- b_entropy
    
    logp("[%s] invariant: %s", b, present)
    
    next
    
  } else {
    
    # for K > 2, training labels mapped to 1..K'
    y_int <- unname(fwd_map[y_codes])  
    
    # check that there is >1 high HF pixel to backfill (subbasins 423,304,247,246,221,160)
    # if nrow(df_backfill_bart) == 1 then take the mean from the training pixels
    if (nrow(df_backfill_bart) == 1) {
      
      # empirical class probabilities from training data
      tab <- table(y_int)
      probs <- tab / sum(tab)
      
      # mode and uncertainty
      mode_int <- as.integer(names(which.max(tab)))
      b_mode   <- inv_map[as.character(mode_int)]
      b_maxprob <- max(probs)
      b_entropy <- -sum(probs * log(pmax(probs, 1e-12)))
      
      # populate outputs
      df_backfill[[b]] <- b_mode
      out_layers[[b]] <- b_mode
      out_layers[[paste0(b, "_maxprob")]] <- b_maxprob
      out_layers[[paste0(b, "_entropy")]] <- b_entropy
      
      logp("[%s] single backfill pixel: skipped mbart()", b)
      
      return(list(
        df_backfill = df_backfill,
        metrics     = metrics,
        confusion   = confusion,
        out_layers  = out_layers
      ))
    }
    
    # reproducible per (year, subbasin, covariate)
    set.seed(abs(as.integer(sprintf("%d%03d", subbasin_index, which(biotic_cols==b)))))
    
    # split y into 90% training and 10% holdout
    holdout_idx <- sample(seq_len(length(y_int)), size = round(0.10 * length(y_int)))
    df_holdout    <- df_train_bart[holdout_idx, ] # get holdout rows before modifying df_train_bart
    df_train_bart <- df_train_bart[-holdout_idx, ]

    # if only one class remains after holdout split, skip mbart
    y_train_int <- y_int[-holdout_idx]

    if (length(unique(y_train_int)) < 2) {

       # empirical class probabilities from full training data (before split)
       tab <- table(y_int)
       probs <- tab / sum(tab)

       mode_int <- as.integer(names(which.max(tab)))
       b_mode   <- inv_map[as.character(mode_int)]
       b_maxprob <- max(probs)
       b_entropy <- -sum(probs * log(pmax(probs, 1e-12)))

       df_backfill[[b]] <- b_mode
       out_layers[[b]] <- b_mode
       out_layers[[paste0(b, "_maxprob")]] <- rep(b_maxprob, nrow(df_backfill_bart))
       out_layers[[paste0(b, "_entropy")]] <- rep(b_entropy, nrow(df_backfill_bart))

       logp("[%s] degenerate after holdout split: single class", b)

       return(list(
       df_backfill = df_backfill,
       metrics     = metrics,
       confusion   = confusion,
       out_layers  = out_layers
       ))
    }
    
    # multinomial bart keeps every 10th posterior sample (gbart keeps every sample)
    # so mbart will take at least 10 times longer
    # gbart assumes continuous residuals -> conjugate priors -> less autocorrelation -> no thinning
    # mbart uses latent variables -> non-conjugate updates -> more autocorrelation -> needs thinning
    # mbart runs K times (e.g. 11 classes will take 11 times longer than gbart)
    fit <- BART::mbart2(x.train = as.matrix(df_train_bart),
                        y.train = y_int[-holdout_idx],
                        x.test  = as.matrix(df_backfill_bart),
                        type = "pbart",
                        k = 3,
                        ntree = 40L,
                        ndpost = 500L, #500
                        nskip = 150L, #150
                        keepevery = 10L,
                        printevery = 350,
                        sparse = TRUE)
    
    logp("mbart fit to [%s]", b)
    
    # prob_ecol is a matrix with class probabilities per backfill pixel 
    # The mbart2() output is pixel-major:
    # row = pixel, columns = class1..classK_model
    K_model    <- fit$K
    cats_model <- fit$cats
    prob_ecol <- matrix(fit$prob.test.mean, ncol = K_model, byrow = TRUE)
    colnames(prob_ecol) <- present[fit$cats]
    
    
    # predicted class index in ecological codes
    # i.e. for every pixel get the most likely landcover class
    pred_idx <- max.col(prob_ecol, ties.method = "first")
    b_mode   <- present[fit$cats][pred_idx]
    
    # uncertainty measures
    eps <- 1e-12
    b_entropy <- -rowSums(prob_ecol * log(pmax(prob_ecol, eps)))
    b_maxprob <- apply(prob_ecol, 1, max)
    
    # variable importance
    varprob_mean <- colMeans(do.call(rbind, fit$varprob))
    top_var <- names(sort(varprob_mean, decreasing = TRUE))[1]
    
    # store most likely class for rasterization
    df_backfill[[b]]                         <- b_mode
    out_layers[[b]]                          <- b_mode
    out_layers[[paste0(b, "_maxprob")]]      <- b_maxprob
    out_layers[[paste0(b, "_entropy")]]      <- b_entropy
    
    # collect metrics for categorical covariate b
    if (K >= 2){
      metrics[[length(metrics) + 1L]] <- collect_metrics_mbart(
        fit,
        y         = present[y_int[-holdout_idx]],
        X_train   = df_train_bart,
        covariate = b,
        subbasin  = subbasin_index,
        year      = year,
        top_var   = top_var,
        present   = present)
    } 
    
    # compute metrics from the holdout data
    # first, predict on holdout data (will be re-used for confusion matrix)
    df_holdout <- df_holdout[, colnames(fit$varprob[[1]]), drop = FALSE] # drop covariates not used by BART
    pred <- predict(fit, newdata = as.matrix(df_holdout))
    metrics[[length(metrics) + 1L]] <- collect_holdout_metrics_mbart(
      fit       = fit,
      pred      = pred,
      actual_holdout_class = present[y_int[holdout_idx]],
      covariate = b,
      subbasin  = subbasin_index,
      year      = year,
      top_var   = top_var,
      present   = present
    )
    
    # use predictions for confusion matrix
    prob_ecol <- matrix(pred$prob.test.mean, ncol = K_model, byrow = TRUE)
    colnames(prob_ecol) <- present[fit$cats]
    
    yhat_holdout_class <- present[max.col(prob_ecol, ties.method = "first")]
    actual_holdout_class <- present[y_int[holdout_idx]]
    
    confusion[[length(confusion) + 1L]] <- data.frame(
      covariate = b,
      subbasin  = subbasin_index,
      actual    = actual_holdout_class,
      predicted = yhat_holdout_class,
      stringsAsFactors = FALSE
    )
    
  } # close if single or multi-class
  
  list(
    df_backfill = df_backfill,
    metrics     = metrics,
    confusion   = confusion,
    out_layers  = out_layers
    )
  
} # close function
