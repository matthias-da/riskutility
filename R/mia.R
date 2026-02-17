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
#' Train a probability random forest and return P(class=1) on test.
#' X_train data.frame of training features (factors + numerics ok).
#' y_train factor with levels c("0","1").
#' X_test  data.frame of test features (same columns as X_train).
#' num_threads integer, number of threads for ranger.
#' num_trees integer, number of trees in the forest.
#' numeric vector of length nrow(X_test) with P(class = 1).
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
#' @param real_data data.frame. Records used to train the synthetic generator.
#' @param synt_data data.frame. Synthetic records produced by the generator.
#' @param hout_data data.frame. Holdout records from the same population,
#'   NOT used to train the generator.
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
#'
#' @return Named list with MIA precision, recall, macro F1 (means and SEs).
#'   - MIA recall is the primary privacy metric.
#'   - The evaluation test set is balanced (50/50).
#'   - Values near 0.5 indicate performance close to random guessing (low risk).
#'   - Values > 0.5 indicate the synthetic data leaks training membership.
#'   - Values < 0.5 indicate performance worse than random guessing.
mia_classifier <- function(real_data,
                           synt_data,
                           hout_data,
                           cat_cols      = NULL,
                           num_cols      = NULL,
                           method        = "rf",
                           num_eval_iter = 5L,
                           seed          = NULL,
                           num_trees     = 300L,
                           num_threads   = 1L) {

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

  list(
    "MIA precision"    = p$mean,
    "MIA precision se" = p$se,
    "MIA recall"       = r$mean,
    "MIA recall se"    = r$se,
    "MIA macro F1"     = f$mean,
    "MIA macro F1 se"  = f$se
  )
}
