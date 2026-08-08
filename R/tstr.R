#' Train on Synthetic, Test on Real (TSTR) Utility Measure
#'
#' Evaluates synthetic data utility by comparing the predictive performance of
#' models trained on synthetic data versus models trained on original data.
#' A TSTR ratio near 1 indicates that the synthetic data preserves the
#' predictive structure of the original data.
#'
#' @param X A data.frame or data.table containing the original dataset.
#' @param Y A data.frame or data.table containing the synthetic/anonymized dataset.
#' @param target_var Character vector of target variable name(s). Required.
#'   If multiple targets are given, TSTR is computed for each and results are averaged.
#' @param model Character, the model type to use. One of \code{"glm"} (default),
#'   \code{"rf"} (random forest via ranger), \code{"rpart"} (CART), or
#'   \code{"gbm"} (gradient boosting). Requires the \pkg{caret} package and
#'   model-specific packages (\pkg{ranger}, \pkg{rpart}, \pkg{xgboost}).
#' @param test_fraction Numeric, fraction of original data held out for testing.
#'   Default 0.3. Ignored if \code{test_data} is provided.
#' @param test_data Optional data.frame of pre-split test data. If provided,
#'   \code{test_fraction} is ignored and the full original data is used for
#'   TRTR training.
#' @param na.rm Logical, whether to remove rows with NA values. Default TRUE.
#' @param seed Integer, random seed for reproducibility. Default NULL.
#' @param ... additional arguments passed to methods (currently unused)
#'
#' @return An object of class \code{"tstr"} containing:
#' \itemize{
#'   \item \code{tstr_ratio}: ratio of TSTR to TRTR performance (1.0 = perfect)
#'   \item \code{utility_score}: clamped ratio in \eqn{[0,1]} (higher = better)
#'   \item \code{tstr_performance}: model performance trained on synthetic data
#'   \item \code{trtr_performance}: baseline performance trained on original data
#'   \item \code{metric}: performance metric used (\code{"R2"} or \code{"AUC"})
#'   \item \code{model}: model type used
#'   \item \code{target_var}: target variable name(s)
#'   \item \code{n_train_orig}: number of original training observations
#'   \item \code{n_test}: number of test observations
#'   \item \code{n_synth}: number of synthetic observations
#'   \item \code{test_fraction}: fraction used for test split
#'   \item \code{per_target}: list of per-target results (when multiple targets)
#' }
#'
#' @details
#' The TSTR protocol evaluates whether a synthetic dataset preserves the
#' predictive relationships in the original data:
#'
#' \enumerate{
#'   \item Split original data into training (1 - \code{test_fraction}) and
#'     test (\code{test_fraction}) sets.
#'   \item \strong{TRTR} (baseline): Train on original training set, evaluate
#'     on test set.
#'   \item \strong{TSTR}: Train on full synthetic data, evaluate on the same
#'     test set.
#'   \item Compute \code{tstr_ratio = TSTR_performance / TRTR_performance}.
#' }
#'
#' For classification targets (factor/character), the AUC (Area Under the ROC
#' Curve) is computed via the Mann-Whitney U statistic. For regression targets
#' (numeric), R-squared is used.
#'
#' Models are trained without cross-validation (\code{caret::trainControl(method = "none")})
#' to ensure training uses the full training set.
#'
#' @seealso \code{\link{propscore}} for propensity score utility,
#'   \code{\link{compare_model_performance}} for cross-validated model comparison
#'
#' @references
#' Esteban, C., Hyland, S. L., and Ratsch, G. (2017). Real-valued (Medical) Time
#' Series Generation with Recurrent Conditional GANs. arXiv:1706.02633.
#'
#' @author Matthias Templ
#' @family utility
#' @export
#' @importFrom stats complete.cases as.formula predict glm
#'
#' @examples
#' set.seed(42)
#' n <- 200
#' X <- data.frame(
#'   x1 = rnorm(n),
#'   x2 = rnorm(n)
#' )
#' X$y <- 2 * X$x1 + 0.5 * X$x2 + rnorm(n, sd = 0.5)
#'
#' # Good synthetic data (same relationship)
#' Y_good <- data.frame(
#'   x1 = rnorm(n),
#'   x2 = rnorm(n)
#' )
#' Y_good$y <- 2 * Y_good$x1 + 0.5 * Y_good$x2 + rnorm(n, sd = 0.5)
#'
#' \donttest{
#' if (requireNamespace("caret", quietly = TRUE)) {
#'   result <- tstr(X, Y_good, target_var = "y", seed = 1)
#'   print(result)
#'   summary(result)
#' }
#' }
tstr <- function(X, ...) {
  UseMethod("tstr")
}

#' @rdname tstr
#' @export
tstr.synth_pair <- function(X, ...) {
  tstr.default(
    X = X$original,
    Y = X$synthetic,
    ...
  )
}

#' @rdname tstr
#' @export
tstr.default <- function(X, Y,
                         target_var = NULL,
                         model = c("glm", "rf", "rpart", "gbm"),
                         test_fraction = 0.3,
                         test_data = NULL,
                         na.rm = TRUE,
                         seed = NULL,
                         ...) {

  # --- Input validation -------------------------------------------------------
  if (is.null(target_var) || length(target_var) == 0) {
    stop("'target_var' is required.")
  }

  if (!requireNamespace("caret", quietly = TRUE)) {
    stop("Package 'caret' is required for tstr(). Please install it.")
  }

  model <- match.arg(model)

  # Map model names to caret method strings
  caret_method <- switch(model,
    glm   = "glm",
    rf    = "ranger",
    rpart = "rpart",
    gbm   = "gbm"
  )

  # Check model-specific package availability
  pkg_needed <- switch(model,
    rf    = "ranger",
    rpart = "rpart",
    gbm   = "xgboost",
    NULL
  )
  if (!is.null(pkg_needed) && !requireNamespace(pkg_needed, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required for model='%s'. Please install it.",
                 pkg_needed, model))
  }

  # Convert to data.frame
  X <- as.data.frame(X)
  Y <- as.data.frame(Y)

  # Check target vars exist
  missing_X <- setdiff(target_var, names(X))
  missing_Y <- setdiff(target_var, names(Y))
  if (length(missing_X) > 0) {
    stop(paste("Target variables missing in X:", paste(missing_X, collapse = ", ")))
  }
  if (length(missing_Y) > 0) {
    stop(paste("Target variables missing in Y:", paste(missing_Y, collapse = ", ")))
  }

  # NA handling
  if (na.rm) {
    X <- X[complete.cases(X), , drop = FALSE]
    Y <- Y[complete.cases(Y), , drop = FALSE]
    if (!is.null(test_data)) {
      test_data <- as.data.frame(test_data)
      test_data <- test_data[complete.cases(test_data), , drop = FALSE]
    }
  }

  if (nrow(X) < 10) stop("Too few complete cases in X (need at least 10).")
  if (nrow(Y) < 10) stop("Too few complete cases in Y (need at least 10).")

  # Set seed
  if (!is.null(seed)) set.seed(seed)

  # --- Split original data or use provided test data --------------------------
  if (!is.null(test_data)) {
    test_data <- as.data.frame(test_data)
    orig_train <- X
    orig_test <- test_data
    actual_test_fraction <- nrow(orig_test) / (nrow(X) + nrow(orig_test))
  } else {
    n <- nrow(X)
    n_test <- max(2L, round(n * test_fraction))
    idx_test <- sample(n, n_test)
    orig_test <- X[idx_test, , drop = FALSE]
    orig_train <- X[-idx_test, , drop = FALSE]
    actual_test_fraction <- test_fraction
  }

  # --- Per-target computation -------------------------------------------------
  per_target <- list()

  for (tv in target_var) {
    is_classification <- is.factor(X[[tv]]) || is.character(X[[tv]])

    # Ensure factor for classification
    if (is_classification) {
      orig_train[[tv]] <- as.factor(orig_train[[tv]])
      orig_test[[tv]] <- as.factor(orig_test[[tv]])
      Y[[tv]] <- as.factor(Y[[tv]])

      # Ensure consistent factor levels
      all_levels <- union(union(levels(orig_train[[tv]]), levels(orig_test[[tv]])),
                          levels(Y[[tv]]))
      orig_train[[tv]] <- factor(orig_train[[tv]], levels = all_levels)
      orig_test[[tv]] <- factor(orig_test[[tv]], levels = all_levels)
      Y[[tv]] <- factor(Y[[tv]], levels = all_levels)
    }

    # Build formula: target ~ all other common variables
    common_vars <- intersect(names(X), names(Y))
    predictor_vars <- setdiff(common_vars, tv)
    if (length(predictor_vars) == 0) {
      stop(sprintf("No predictor variables available for target '%s'.", tv))
    }
    fml <- as.formula(paste(tv, "~", paste(predictor_vars, collapse = " + ")))

    # Determine metric
    metric_name <- if (is_classification) "AUC" else "R2"

    # caret trainControl: no CV, just use full training set
    if (is_classification) {
      tc <- caret::trainControl(method = "none", classProbs = TRUE)
    } else {
      tc <- caret::trainControl(method = "none")
    }

    # Build model-specific tuneGrid (caret requires tuneGrid with method="none")
    tune_grid <- .tstr_tune_grid(model, caret_method, is_classification,
                                 ncol_predictors = length(predictor_vars))

    # GLM for classification needs family = binomial (only for 2-class)
    extra_args <- list()
    if (model == "glm" && is_classification) {
      n_levels <- length(all_levels)
      if (n_levels != 2) {
        stop(sprintf(
          "GLM classification requires exactly 2 classes, but target '%s' has %d. Use model='rf' or 'rpart'.",
          tv, n_levels))
      }
      extra_args$family <- "binomial"
    }

    if (model == "gbm") {
      extra_args$verbose <- FALSE
    }

    # --- TRTR: Train on original, test on original ---
    trtr_perf <- .tstr_train_eval(
      train_data = orig_train, test_data = orig_test,
      formula = fml, caret_method = caret_method,
      trainControl = tc, tune_grid = tune_grid,
      is_classification = is_classification,
      target_var = tv, extra_args = extra_args
    )

    # --- TSTR: Train on synthetic, test on original ---
    tstr_perf <- .tstr_train_eval(
      train_data = Y, test_data = orig_test,
      formula = fml, caret_method = caret_method,
      trainControl = tc, tune_grid = tune_grid,
      is_classification = is_classification,
      target_var = tv, extra_args = extra_args
    )

    # Compute ratio. If the train-on-real model itself has no predictive power
    # there is no transferable structure, so the ratio is undefined rather than
    # zero or one -- reporting a number here would credit (or blame) the
    # synthesiser for a relationship that does not exist in the original data.
    if (trtr_perf < 1e-10) {
      warning("tstr(): the train-on-real model has no predictive power for '",
              tv, "' (", metric_name, " = ", format(trtr_perf, digits = 3),
              "); the TSTR ratio is undefined and reported as NA.",
              call. = FALSE)
      ratio <- NA_real_
    } else {
      ratio <- tstr_perf / trtr_perf
    }

    per_target[[tv]] <- list(
      tstr_performance = tstr_perf,
      trtr_performance = trtr_perf,
      tstr_ratio = ratio,
      metric = metric_name,
      is_classification = is_classification
    )
  }

  # --- Aggregate across targets -----------------------------------------------
  ratios <- vapply(per_target, function(pt) pt$tstr_ratio, numeric(1))
  tstr_perfs <- vapply(per_target, function(pt) pt$tstr_performance, numeric(1))
  trtr_perfs <- vapply(per_target, function(pt) pt$trtr_performance, numeric(1))
  metrics <- unname(vapply(per_target, function(pt) pt$metric, character(1)))

  avg_ratio <- if (all(is.na(ratios))) NA_real_ else mean(ratios, na.rm = TRUE)
  avg_tstr <- mean(tstr_perfs)
  avg_trtr <- mean(trtr_perfs)
  metric_used <- if (length(unique(metrics)) == 1) metrics[1] else paste(unique(metrics), collapse = "/")

  utility_score <- if (is.na(avg_ratio)) NA_real_ else max(0, min(1, avg_ratio))

  result <- list(
    tstr_ratio = avg_ratio,
    utility_score = utility_score,
    tstr_performance = avg_tstr,
    trtr_performance = avg_trtr,
    metric = metric_used,
    model = model,
    target_var = target_var,
    n_train_orig = nrow(orig_train),
    n_test = nrow(orig_test),
    n_synth = nrow(Y),
    test_fraction = actual_test_fraction,
    per_target = per_target
  )

  class(result) <- "tstr"
  return(result)
}


# ---- Internal helpers --------------------------------------------------------

# Build tuneGrid for caret with method = "none"
# @param model user-facing model name
# @param caret_method caret method string
# @param is_classification logical
# @param ncol_predictors number of predictor columns
# @return data.frame tuneGrid
# @keywords internal
.tstr_tune_grid <- function(model, caret_method, is_classification, ncol_predictors) {
  switch(model,
    glm = NULL,  # glm has no tuning parameters
    rf = data.frame(
      mtry = max(1L, floor(sqrt(ncol_predictors))),
      splitrule = if (is_classification) "gini" else "variance",
      min.node.size = if (is_classification) 1L else 5L
    ),
    rpart = data.frame(cp = 0.01),
    gbm = data.frame(
      nrounds = 100L,
      max_depth = 3L,
      eta = 0.1,
      gamma = 0,
      colsample_bytree = 1,
      min_child_weight = 1,
      subsample = 1
    ),
    NULL
  )
}


# Train a model and evaluate on test data
# @param train_data data.frame for training
# @param test_data data.frame for testing
# @param formula model formula
# @param caret_method caret method string
# @param trainControl caret trainControl object
# @param tune_grid tuneGrid data.frame
# @param is_classification logical
# @param target_var character, target variable name
# @param extra_args list of extra arguments for caret::train
# @return numeric, performance metric (R2 or AUC)
# @keywords internal
.tstr_train_eval <- function(train_data, test_data, formula, caret_method,
                             trainControl, tune_grid, is_classification,
                             target_var, extra_args = list()) {

  # Build caret::train call
  train_args <- list(
    form = formula,
    data = train_data,
    method = caret_method,
    trControl = trainControl
  )
  if (!is.null(tune_grid)) {
    train_args$tuneGrid <- tune_grid
  }
  train_args <- c(train_args, extra_args)

  fit <- do.call(caret::train, train_args)

  if (is_classification) {
    # Predict probabilities for AUC
    probs <- predict(fit, newdata = test_data, type = "prob")
    truth <- test_data[[target_var]]
    auc <- .tstr_auc_mannwhitney(probs, truth)
    return(auc)
  } else {
    # Predict numeric values for R2
    preds <- predict(fit, newdata = test_data)
    truth <- test_data[[target_var]]
    r2 <- .tstr_r_squared(preds, truth)
    return(r2)
  }
}


# Compute AUC via Mann-Whitney U statistic
# For binary classification: AUC = U / (n0 * n1)
# For multiclass: macro-average of one-vs-rest AUCs
# @param probs data.frame of predicted probabilities (one column per class)
# @param truth factor of true labels
# @return numeric AUC in [0, 1]
# @keywords internal
.tstr_auc_mannwhitney <- function(probs, truth) {
  levels_truth <- levels(truth)
  n_classes <- length(levels_truth)

  if (n_classes == 2) {
    # Binary AUC
    pos_class <- levels_truth[2]
    pos_probs <- probs[[pos_class]]
    binary_truth <- as.integer(truth == pos_class)

    idx1 <- which(binary_truth == 1)
    idx0 <- which(binary_truth == 0)

    if (length(idx1) == 0 || length(idx0) == 0) return(0.5)

    n1 <- length(idx1)
    n0 <- length(idx0)
    # Mann-Whitney U
    U <- 0
    for (i in idx1) {
      U <- U + sum(pos_probs[i] > pos_probs[idx0]) +
        0.5 * sum(pos_probs[i] == pos_probs[idx0])
    }
    return(U / (n1 * n0))
  } else {
    # Macro-average one-vs-rest AUC
    aucs <- numeric(n_classes)
    for (k in seq_len(n_classes)) {
      cls <- levels_truth[k]
      if (cls %in% names(probs)) {
        pos_probs <- probs[[cls]]
      } else {
        aucs[k] <- 0.5
        next
      }
      binary_truth <- as.integer(truth == cls)
      idx1 <- which(binary_truth == 1)
      idx0 <- which(binary_truth == 0)
      if (length(idx1) == 0 || length(idx0) == 0) {
        aucs[k] <- 0.5
        next
      }
      n1 <- length(idx1)
      n0 <- length(idx0)
      U <- 0
      for (i in idx1) {
        U <- U + sum(pos_probs[i] > pos_probs[idx0]) +
          0.5 * sum(pos_probs[i] == pos_probs[idx0])
      }
      aucs[k] <- U / (n1 * n0)
    }
    return(mean(aucs))
  }
}


# Compute R-squared
# R2 = 1 - SS_res / SS_tot
# Clamps to [0, 1] to avoid negative R2 for very poor models
# @param preds numeric vector of predictions
# @param truth numeric vector of true values
# @return numeric R2 in [0, 1]
# @keywords internal
.tstr_r_squared <- function(preds, truth) {
  ss_res <- sum((truth - preds)^2)
  ss_tot <- sum((truth - mean(truth))^2)
  if (ss_tot < 1e-10) return(1)  # constant target
  r2 <- 1 - ss_res / ss_tot
  # Clamp at 0. A negative R2 means the model predicts worse than the target
  # mean; it carries no information about how much predictive structure was
  # transferred, and letting it through makes the TSTR/TRTR ratio a quotient of
  # two negative numbers (which can exceed 1 and read as high utility).
  return(max(0, r2))
}


# ---- S3 methods --------------------------------------------------------------

#' @rdname tstr
#' @param x an object of class \code{"tstr"}
#' @export
print.tstr <- function(x, ...) {
  cat("Train on Synthetic, Test on Real (TSTR)\n")
  cat("========================================\n\n")

  cat("Dataset Sizes:\n")
  cat("  Original training:", x$n_train_orig, "\n")
  cat("  Original test:    ", x$n_test, "\n")
  cat("  Synthetic:        ", x$n_synth, "\n")
  cat("  Test fraction:    ", sprintf("%.1f%%", x$test_fraction * 100), "\n\n")

  cat("Settings:\n")
  cat("  Model:    ", x$model, "\n")
  cat("  Target(s):", paste(x$target_var, collapse = ", "), "\n")
  cat("  Metric:   ", x$metric, "\n\n")

  cat("Results:\n")
  cat("  TRTR performance:", sprintf("%.4f", x$trtr_performance), "(baseline)\n")
  cat("  TSTR performance:", sprintf("%.4f", x$tstr_performance), "\n")
  cat("  TSTR ratio:      ", sprintf("%.4f", x$tstr_ratio), "\n")
  cat("  Utility score:   ", sprintf("%.4f", x$utility_score), "(higher = better)\n\n")

  if (length(x$per_target) > 1) {
    cat("Per-Target Results:\n")
    for (tv in names(x$per_target)) {
      pt <- x$per_target[[tv]]
      cat(sprintf("  %s: TRTR=%.4f, TSTR=%.4f, ratio=%.4f (%s)\n",
                  tv, pt$trtr_performance, pt$tstr_performance,
                  pt$tstr_ratio, pt$metric))
    }
    cat("\n")
  }

  cat("Interpretation:\n")
  if (is.na(x$utility_score)) {
    cat("  UNDEFINED: the train-on-real model has no predictive power, so there\n")
    cat("  is no predictive structure whose transfer could be assessed.\n")
  } else if (x$utility_score > 0.95) {
    cat("  EXCELLENT: Synthetic data preserves predictive structure very well.\n")
  } else if (x$utility_score > 0.80) {
    cat("  GOOD: Synthetic data preserves most predictive relationships.\n")
  } else if (x$utility_score > 0.50) {
    cat("  MODERATE: Some predictive information is lost in synthetic data.\n")
  } else {
    cat("  POOR: Significant predictive structure lost in synthetic data.\n")
  }

  invisible(x)
}


#' @rdname tstr
#' @param object an object of class \code{"tstr"}
#' @export
summary.tstr <- function(object, ...) {
  summ <- list(
    tstr_ratio = object$tstr_ratio,
    utility_score = object$utility_score,
    tstr_performance = object$tstr_performance,
    trtr_performance = object$trtr_performance,
    metric = object$metric,
    model = object$model,
    target_var = object$target_var,
    n_train_orig = object$n_train_orig,
    n_test = object$n_test,
    n_synth = object$n_synth,
    test_fraction = object$test_fraction,
    per_target = object$per_target,
    performance_gap = object$trtr_performance - object$tstr_performance
  )

  class(summ) <- "summary.tstr"
  return(summ)
}


#' @rdname tstr
#' @export
print.summary.tstr <- function(x, ...) {
  cat("Summary: Train on Synthetic, Test on Real (TSTR)\n")
  cat("=================================================\n\n")

  cat("Model:", x$model, " | Metric:", x$metric, "\n")
  cat("Target(s):", paste(x$target_var, collapse = ", "), "\n\n")

  cat("Performance Comparison:\n")
  cat("  TRTR (baseline):  ", sprintf("%.4f", x$trtr_performance), "\n")
  cat("  TSTR (synthetic): ", sprintf("%.4f", x$tstr_performance), "\n")
  cat("  Gap:              ", sprintf("%.4f", x$performance_gap), "\n")
  cat("  TSTR ratio:       ", sprintf("%.4f", x$tstr_ratio), "\n")
  cat("  Utility score:    ", sprintf("%.4f", x$utility_score), "\n\n")

  if (length(x$per_target) > 1) {
    cat("Per-Target Breakdown:\n")
    for (tv in names(x$per_target)) {
      pt <- x$per_target[[tv]]
      cat(sprintf("  %-15s  TRTR=%.4f  TSTR=%.4f  ratio=%.4f  (%s)\n",
                  tv, pt$trtr_performance, pt$tstr_performance,
                  pt$tstr_ratio, pt$metric))
    }
    cat("\n")
  }

  cat("Sample Sizes: train_orig =", x$n_train_orig,
      ", test =", x$n_test,
      ", synth =", x$n_synth, "\n")

  invisible(x)
}


#' @rdname tstr
#' @param y not used
#' @param which integer, which plot: 1 = grouped bar chart comparing TRTR vs TSTR
#'   performance (per target if multiple)
#' @importFrom graphics barplot legend par mtext
#' @export
plot.tstr <- function(x, y = NULL, ..., which = 1) {
  show <- rep(FALSE, 1)
  show[which] <- TRUE

  if (show[1]) {
    n_targets <- length(x$per_target)
    target_names <- names(x$per_target)

    # Build matrix: rows = TRTR/TSTR, cols = targets
    trtr_vals <- vapply(x$per_target, function(pt) pt$trtr_performance, numeric(1))
    tstr_vals <- vapply(x$per_target, function(pt) pt$tstr_performance, numeric(1))

    bar_mat <- rbind(TRTR = trtr_vals, TSTR = tstr_vals)
    colnames(bar_mat) <- target_names

    barplot(bar_mat,
            beside = TRUE,
            main = "TRTR vs TSTR Performance",
            ylab = x$metric,
            col = c("steelblue", "coral"),
            las = 1,
            ...)

    legend("topright",
           legend = c("TRTR (baseline)", "TSTR (synthetic)"),
           fill = c("steelblue", "coral"),
           bty = "n")

    mtext(sprintf("TSTR ratio = %.3f | Utility = %.3f",
                  x$tstr_ratio, x$utility_score),
          side = 3, line = 0, cex = 0.9)
  }

  invisible(x)
}
