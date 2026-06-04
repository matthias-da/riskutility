# =============================================================================
# Membership Inference Attack (MIA) via Classification
# =============================================================================
#
# PURPOSE
# -------
# Estimates membership inference risk for synthetic data by training a
# classifier to distinguish synthetic records from holdout records, then
# testing whether real training records are identified as "synthetic-like."
#
# High recall on training records indicates the synthetic data has memorised
# characteristics of its training set, posing a privacy risk.
#
# LOGIC (step by step)
# --------------------
# 1. Three datasets are required:
#    - real_data  : the records used to TRAIN the synthetic data generator
#    - synt_data  : the output of the generator
#    - hout_data  : real records from the same population, NOT used for training
#
# 2. The holdout is split 75/25 into hout_train and hout_test.
#
# 3. A classifier is trained on:
#      class 1 = synthetic data  (sampled to match hout_train size)
#      class 0 = hout_train
#    This teaches the model: "what looks like synthetic data vs. unseen real data?"
#
# 4. A balanced test set is constructed:
#      class 1 = real training records  (the ones we want to detect)
#      class 0 = hout_test
#
# 5. The classifier predicts on this test set. If a real training record is
#    classified as "synthetic-like" (class 1), the generator likely memorised it.
#
# 6. Recall on class 1 = proportion of training records successfully identified
#    = membership inference risk.
#    Baseline (no memorisation) ~ 0.5. Values >> 0.5 indicate privacy risk.
#
# 7. Steps 2-6 are repeated `num_eval_iter` times with different holdout splits
#    for stability. Results are averaged with standard errors.
#
# CLASSIFIER BACK-END
# -------------------
# The classifier is selected via the `method` argument. Currently supported:
#   - "rf" : Random forest via the {ranger} package (default)
#
# To add a new classifier in the future, implement a function with signature:
#   my_classifier(X_train, y_train, X_test, num_threads, ...)
# that returns a numeric vector of predicted P(class = 1) for each test row.
# Then register it in the `classifiers` list inside mia_classifier().
#
# =============================================================================


# -----------------------------------------------------------------------------
# Random forest back-end (ranger)
# -----------------------------------------------------------------------------
#' Train a probability random forest and return P(class=1) on test
#'
#' Internal helper for \code{\link{mia_classifier}}.
#'
#' @param X_train data.frame of training features (factors and numerics).
#' @param y_train factor with levels \code{c("0","1")}.
#' @param X_test data.frame of test features (same columns as X_train).
#' @param num_threads integer, number of threads for ranger.
#' @param num_trees integer, number of trees in the forest.
#' @return Numeric vector of length \code{nrow(X_test)} with P(class = 1).
#' @keywords internal
.mia_backend_rf <- function(X_train, y_train, X_test,
                            num_threads = 1L, num_trees = 300L) {
  if (!requireNamespace("ranger", quietly = TRUE)) {
    stop("Package 'ranger' is required for method = 'rf'. Please install it.",
         call. = FALSE)
  }

  df_train     <- X_train
  df_train$y   <- y_train

  fit <- ranger::ranger(
    y ~ .,
    data        = df_train,
    num.trees   = num_trees,
    probability = TRUE,
    num.threads = num_threads
  )

  preds  <- predict(fit, data = X_test)$predictions
  prob_1 <- preds[, "1"]
  return(prob_1)
}


# -----------------------------------------------------------------------------
# Main function
# -----------------------------------------------------------------------------
#' Membership Inference Attack metric via classification
#'
#' @param X synth_pair object, or (for the default method) a data.frame of
#'   records used to train the synthetic generator.
#' @param synt_data data.frame. Synthetic records produced by the generator.
#' @param hout_data data.frame. Holdout records from the same population,
#'   NOT used to train the generator.
#' @param Y Alias for \code{synt_data}, for naming consistency with the rest of
#'   the package. If \code{synt_data} is \code{NULL} and \code{Y} is supplied,
#'   \code{Y} is used.
#' @param holdout Alias for \code{hout_data}. If \code{hout_data} is \code{NULL}
#'   and \code{holdout} is supplied, \code{holdout} is used.
#' @param cols character vector of column names to use (optional). If NULL, all common
#' @param cat_cols  character vector of categorical column names (optional).
#'   If NULL, inferred as non-numeric columns.
#' @param num_cols  character vector of numeric column names (optional).
#'   If NULL, inferred as numeric columns.
#' @param method    character. Classification back-end to use.
#'   Currently supported: "rf" (random forest via ranger).
#' @param num_eval_iter integer. Number of evaluation iterations with different
#'   holdout splits. Default 5.
#' @param seed integer or NULL. Random seed for reproducibility.
#' @param num_trees integer. Number of trees (only used for method = "rf").
#' @param num_threads integer. Threads for parallelism. Default 1.
#' @param ... additional arguments passed to methods
#'
#' @return An object of class \code{"mia"} with components:
#' \describe{
#'   \item{precision}{Mean precision across iterations.}
#'   \item{precision_se}{Standard error of precision.}
#'   \item{recall}{Mean recall (primary privacy metric).}
#'   \item{recall_se}{Standard error of recall.}
#'   \item{macro_f1}{Mean macro F1 score.}
#'   \item{macro_f1_se}{Standard error of macro F1.}
#'   \item{privacy_pass}{Logical, TRUE if recall <= 0.55.}
#'   \item{num_eval_iter}{Number of evaluation iterations.}
#'   \item{method}{Classification method used.}
#' }
#' Values near 0.5 indicate performance close to random guessing (low risk).
#' Values > 0.5 indicate the synthetic data leaks training membership.
#' Values < 0.5 indicate performance worse than random guessing.
#' @family privacy-models
#' @author Matthias Templ
#' @export
mia_classifier <- function(X, ...) {
  UseMethod("mia_classifier")
}

#' @rdname mia_classifier
#' @export
mia_classifier.synth_pair <- function(X, ...) {
  if (!is.null(X$source) && X$source == "sdcMicro") {
    stop("mia_classifier is designed for synthetic data evaluation and is not applicable to ",
         "traditionally anonymized data (sdcMicro objects). ",
         "Use dcr(), nndr(), or ims() for distance-based privacy evaluation instead.",
         call. = FALSE)
  }
  if (is.null(X$holdout)) {
    stop("mia_classifier requires holdout data. Provide a synth_pair with holdout, ",
         "or use mia_classifier.default() directly with real_data, synt_data, hout_data.",
         call. = FALSE)
  }
  mia_classifier.default(
    X = X$original,
    synt_data = X$synthetic,
    hout_data = X$holdout,
    cols      = X$vars,
    ...
  )
}

#' @rdname mia_classifier
#' @export
mia_classifier.default <- function(X,
                           synt_data     = NULL,
                           hout_data     = NULL,
                           cols          = NULL,
                           cat_cols      = NULL,
                           num_cols      = NULL,
                           method        = "rf",
                           num_eval_iter = 5L,
                           seed          = NULL,
                           num_trees     = 300L,
                           num_threads   = 1L,
                           Y             = NULL,
                           holdout       = NULL,
                           ...) {

  # Accept the package-wide names 'Y' (synthetic) and 'holdout' as aliases for
  # the local 'synt_data' / 'hout_data'. The original names remain valid.
  if (is.null(synt_data) && !is.null(Y))       synt_data <- Y
  if (is.null(hout_data) && !is.null(holdout)) hout_data <- holdout
  if (is.null(synt_data) || is.null(hout_data)) {
    stop("Provide synthetic data (via 'Y' or 'synt_data') and holdout data ",
         "(via 'holdout' or 'hout_data').", call. = FALSE)
  }

  real_data <- X

  # --- Registry of available classifier back-ends ----------------------------
  # To add a new method, implement a function with the same signature as
  # .mia_backend_rf and add an entry here.
  classifiers <- list(
    rf = .mia_backend_rf
  )

  method <- tolower(method)
  if (!method %in% names(classifiers)) {
    stop(
      sprintf(
        "Unknown method '%s'. Available methods: %s",
        method, paste(names(classifiers), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  clf_fn <- classifiers[[method]]

  # --- Input validation ------------------------------------------------------
  if (!is.data.frame(real_data) || !is.data.frame(synt_data) || !is.data.frame(hout_data)) {
    stop("real_data, synt_data, hout_data must all be data.frames.", call. = FALSE)
  }

  num_eval_iter <- as.integer(num_eval_iter)
  num_trees     <- as.integer(num_trees)
  num_threads   <- as.integer(num_threads)

  if (is.na(num_eval_iter) || num_eval_iter < 1L) stop("num_eval_iter must be >= 1.", call. = FALSE)
  if (is.na(num_trees)     || num_trees     < 1L) stop("num_trees must be >= 1.",     call. = FALSE)
  if (is.na(num_threads)   || num_threads   < 1L) stop("num_threads must be >= 1.",   call. = FALSE)

  # --- Column alignment ------------------------------------------------------
  # Keep only columns present in all three datasets.
  common_cols <- Reduce(intersect, list(names(real_data), names(synt_data), names(hout_data)))
  if (length(common_cols) == 0L) {
    stop("No common columns across real / synthetic / holdout data.", call. = FALSE)
  }

  # If cols is specified, restrict to those columns only
  if (!is.null(cols)) {
    miss <- setdiff(cols, common_cols)
    if (length(miss)) {
      stop("cols not found in common columns: ", paste(miss, collapse = ", "), call. = FALSE)
    }
    common_cols <- cols
  }

  real_data <- as.data.frame(real_data[, common_cols, drop = FALSE])
  synt_data <- as.data.frame(synt_data[, common_cols, drop = FALSE])
  hout_data <- as.data.frame(hout_data[, common_cols, drop = FALSE])

  # --- Infer column roles if not provided ------------------------------------
  is_num <- vapply(real_data, is.numeric, logical(1))

  if (is.null(cat_cols)) {
    cat_cols <- names(real_data)[!is_num]
  } else {
    miss <- setdiff(cat_cols, common_cols)
    if (length(miss)) stop("cat_cols not found: ", paste(miss, collapse = ", "), call. = FALSE)
  }

  if (is.null(num_cols)) {
    num_cols <- names(real_data)[is_num]
  } else {
    miss <- setdiff(num_cols, common_cols)
    if (length(miss)) stop("num_cols not found: ", paste(miss, collapse = ", "), call. = FALSE)
  }

  # --- Factor encoding with shared levels ------------------------------------
  # All three datasets must share the same factor levels so the classifier
  # sees a consistent encoding across training and testing.
  for (cc in cat_cols) {
    all_levels <- unique(c(
      as.character(real_data[[cc]]),
      as.character(synt_data[[cc]]),
      as.character(hout_data[[cc]])
    ))
    real_data[[cc]] <- factor(as.character(real_data[[cc]]), levels = all_levels)
    synt_data[[cc]] <- factor(as.character(synt_data[[cc]]), levels = all_levels)
    hout_data[[cc]] <- factor(as.character(hout_data[[cc]]), levels = all_levels)
  }

  # --- Coerce numeric columns ------------------------------------------------
  for (nc in num_cols) {
    real_data[[nc]] <- suppressWarnings(as.numeric(real_data[[nc]]))
    synt_data[[nc]] <- suppressWarnings(as.numeric(synt_data[[nc]]))
    hout_data[[nc]] <- suppressWarnings(as.numeric(hout_data[[nc]]))
  }

  # --- NA handling -------------------------------------------------------------
  # No imputation is performed. The default RF back-end (ranger) handles NAs
  # natively by splitting on missingness, which preserves any signal that
  # differences in NA patterns carry (e.g. the generator may handle missingness
  # differently from the real data — that is a legitimate MIA signal).
  #
  # For factors, an explicit NA level is added so ranger can split on it.
  # Future classifier back-ends that cannot handle NAs will need their own
  # imputation strategy.

  # Flag presence of NAs so the user is aware
  na_real <- sum(is.na(real_data))
  na_synt <- sum(is.na(synt_data))
  na_hout <- sum(is.na(hout_data))
  if ((na_real + na_synt + na_hout) > 0L) {
    warning(
      sprintf(
        "NAs detected (real: %d, synthetic: %d, holdout: %d). No imputation applied; the RF back-end handles NAs natively. Other back-ends may require imputation.",
        na_real, na_synt, na_hout
      ),
      call. = FALSE
    )
  }

  # Add explicit NA factor level so ranger can split on missingness
  for (cc in cat_cols) {
    real_data[[cc]] <- addNA(real_data[[cc]])
    synt_data[[cc]] <- addNA(synt_data[[cc]])
    hout_data[[cc]] <- addNA(hout_data[[cc]])
  }

  # --- Binary classification helpers -----------------------------------------
  safe_div <- function(a, b) if (b == 0) as.numeric(NA) else a / b

  precision_bin <- function(y_true, y_pred, pos = 1L) {
    tp <- sum(y_true == pos & y_pred == pos)
    fp <- sum(y_true != pos & y_pred == pos)
    safe_div(tp, tp + fp)
  }

  recall_bin <- function(y_true, y_pred, pos = 1L) {
    tp <- sum(y_true == pos & y_pred == pos)
    fn <- sum(y_true == pos & y_pred != pos)
    safe_div(tp, tp + fn)
  }

  f1_bin <- function(prec, rec) {
    if (is.na(prec) || is.na(rec) || (prec + rec) == 0) {
      return(as.numeric(NA))
    }
    2 * prec * rec / (prec + rec)
  }

  macro_f1 <- function(y_true, y_pred) {
    p1 <- precision_bin(y_true, y_pred, 1L)
    r1 <- recall_bin(y_true, y_pred, 1L)
    p0 <- precision_bin(y_true, y_pred, 0L)
    r0 <- recall_bin(y_true, y_pred, 0L)
    mean(c(f1_bin(p1, r1), f1_bin(p0, r0)), na.rm = TRUE)
  }

  # --- Size checks -----------------------------------------------------------
  n_hout <- nrow(hout_data)
  if (n_hout < 4L)        stop("Holdout data too small for train/test split.", call. = FALSE)
  if (nrow(real_data) < 1L) stop("real_data is empty.", call. = FALSE)
  if (nrow(synt_data) < 1L) stop("synt_data is empty.", call. = FALSE)

  # --- Evaluation loop -------------------------------------------------------
  if (!is.null(seed)) set.seed(seed)

  precs <- numeric(num_eval_iter)
  recs  <- numeric(num_eval_iter)
  f1s   <- numeric(num_eval_iter)

  for (i in seq_len(num_eval_iter)) {

    # -- Split holdout 75 / 25 ------------------------------------------------
    idx     <- sample.int(n_hout)
    n_test  <- max(1L, floor(0.25 * n_hout))
    idx_test  <- idx[seq_len(n_test)]
    idx_train <- idx[-seq_len(n_test)]

    hout_train <- hout_data[idx_train, , drop = FALSE]
    hout_test  <- hout_data[idx_test,  , drop = FALSE]

    # -- Build TRAINING set: synthetic (1) vs holdout-train (0) ----------------
    # Sample synthetic data to match the holdout-train size for class balance.
    syn_n    <- nrow(synt_data)
    need_n   <- nrow(hout_train)
    syn_samp <- synt_data[
      sample.int(syn_n, size = need_n, replace = syn_n < need_n), , drop = FALSE
    ]

    X_train <- rbind(syn_samp, hout_train)
    y_train <- factor(
      c(rep.int(1L, nrow(syn_samp)), rep.int(0L, nrow(hout_train))),
      levels = c("0", "1")
    )

    # -- Build TEST set: real-training (1) vs holdout-test (0) -----------------
    # Balance to equal size by downsampling the larger group.
    n_real <- nrow(real_data)
    n_htest <- nrow(hout_test)

    if (n_real < n_htest) {
      real_samp <- real_data
      hout_samp <- hout_test[sample.int(n_htest, size = n_real, replace = FALSE), , drop = FALSE]
    } else {
      real_samp <- real_data[sample.int(n_real, size = n_htest, replace = FALSE), , drop = FALSE]
      hout_samp <- hout_test
    }

    X_test <- rbind(real_samp, hout_samp)
    y_test <- c(rep.int(1L, nrow(real_samp)), rep.int(0L, nrow(hout_samp)))

    # -- Classify using the selected back-end ----------------------------------
    prob_1 <- clf_fn(
      X_train     = X_train,
      y_train     = y_train,
      X_test      = X_test,
      num_threads = num_threads,
      num_trees   = num_trees
    )

    y_pred <- ifelse(prob_1 >= 0.5, 1L, 0L)

    # -- Record iteration metrics ----------------------------------------------
    precs[i] <- precision_bin(y_test, y_pred, 1L)
    recs[i]  <- recall_bin(y_test, y_pred, 1L)
    f1s[i]   <- macro_f1(y_test, y_pred)
  }

  # --- Aggregate across iterations -------------------------------------------
  mean_se <- function(x) {
    mu <- mean(x, na.rm = TRUE)
    se <- if (length(x) <= 1L) {
      as.numeric(NA)
    } else {
      stats::sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
    }
    list(mean = mu, se = se)
  }

  p <- mean_se(precs)
  r <- mean_se(recs)
  f <- mean_se(f1s)

  result <- list(
    precision    = p$mean,
    precision_se = p$se,
    recall       = r$mean,
    recall_se    = r$se,
    macro_f1     = f$mean,
    macro_f1_se  = f$se,
    num_eval_iter = num_eval_iter,
    method        = method
  )
  result$privacy_pass <- !is.na(r$mean) && r$mean <= 0.55
  class(result) <- "mia"
  return(result)
}


#' Print method for mia objects
#'
#' @param x an object of class "mia"
#' @param ... additional arguments (ignored)
#' @return The input object, invisibly.
#' @export
print.mia <- function(x, ...) {
  cat("Membership Inference Attack (MIA) Results\n")
  cat("==========================================\n\n")
  cat("  Method:     ", x$method, "\n")
  cat("  Iterations: ", x[["num_eval_iter"]], "\n\n")
  cat("  Precision:  ", sprintf("%.4f", x$precision),
      sprintf("(SE: %.4f)", x$precision_se), "\n")
  cat("  Recall:     ", sprintf("%.4f", x$recall),
      sprintf("(SE: %.4f)", x$recall_se), "\n")
  cat("  Macro F1:   ", sprintf("%.4f", x$macro_f1),
      sprintf("(SE: %.4f)", x$macro_f1_se), "\n\n")
  recall <- x$recall
  if (!is.na(recall)) {
    if (recall > 0.6) {
      cat("  Risk level: HIGH (recall substantially above 0.5 baseline)\n")
    } else if (recall > 0.55) {
      cat("  Risk level: MODERATE (recall moderately above 0.5 baseline)\n")
    } else {
      cat("  Risk level: LOW (recall near or below 0.5 baseline)\n")
    }
  }
  cat("  Privacy pass:", x$privacy_pass, "\n")
  invisible(x)
}


#' Summary method for mia objects
#'
#' @param object an object of class "mia"
#' @param ... additional arguments (ignored)
#' @return An object of class "summary.mia"
#' @export
summary.mia <- function(object, ...) {
  summ <- list(
    precision = object$precision,
    precision_se = object$precision_se,
    recall = object$recall,
    recall_se = object$recall_se,
    macro_f1 = object$macro_f1,
    macro_f1_se = object$macro_f1_se,
    num_eval_iter = object[["num_eval_iter"]],
    method = object[["method"]],
    recall_excess = object$recall - 0.5
  )
  class(summ) <- "summary.mia"
  summ
}


#' Print method for summary.mia objects
#'
#' @param x an object of class "summary.mia"
#' @param ... additional arguments (ignored)
#' @return The input object, invisibly.
#' @export
print.summary.mia <- function(x, ...) {
  cat("Summary: Membership Inference Attack (MIA)\n")
  cat("============================================\n\n")

  cat("Method:", x$method, "\n")
  cat("Evaluation iterations:", x$num_eval_iter, "\n\n")

  cat("Metrics (mean +/- SE):\n")
  cat("  Precision: ", sprintf("%.4f +/- %.4f", x$precision, x$precision_se), "\n")
  cat("  Recall:    ", sprintf("%.4f +/- %.4f", x$recall, x$recall_se), "\n")
  cat("  Macro F1:  ", sprintf("%.4f +/- %.4f", x$macro_f1, x$macro_f1_se), "\n\n")

  cat("Recall excess over baseline (0.5):", sprintf("%.4f", x$recall_excess), "\n")
  if (!is.na(x$recall_excess)) {
    if (x$recall_excess > 0.1) {
      cat("Interpretation: Substantial memorisation detected.\n")
    } else if (x$recall_excess > 0.05) {
      cat("Interpretation: Moderate memorisation signal.\n")
    } else if (x$recall_excess > 0) {
      cat("Interpretation: Weak memorisation signal.\n")
    } else {
      cat("Interpretation: No memorisation detected (near or below baseline).\n")
    }
  }

  invisible(x)
}


#' Plot method for mia objects
#'
#' @param x an object of class "mia"
#' @param y not used
#' @param ... additional arguments (ignored)
#' @param which integer, which plot: 1 = metric barplot (default)
#' @return No return value; called for the side effect of producing a plot.
#' @export
#' @importFrom graphics barplot arrows
plot.mia <- function(x, y = NULL, ..., which = 1) {
  if (which == 1) {
    means <- c(x$precision, x$recall, x$macro_f1)
    ses <- c(x$precision_se, x$recall_se, x$macro_f1_se)
    names(means) <- c("Precision", "Recall", "Macro F1")

    bp <- barplot(means, ylim = c(0, min(1, max(means + ses, na.rm = TRUE) * 1.2)),
                  col = c("steelblue", "coral", "seagreen"),
                  main = "MIA Classification Metrics",
                  ylab = "Score", border = NA)
    abline(h = 0.5, lty = 2, col = "grey40")

    # Add error bars
    if (!any(is.na(ses))) {
      arrows(bp, means - ses, bp, means + ses,
             angle = 90, code = 3, length = 0.05, col = "grey30")
    }
  }
  invisible(x)
}
