#' Cross-Validation of Synthesizer Disclosure Risk
#'
#' Evaluates the RAPID disclosure risk of a synthesis method using
#' k-fold cross-validation.
#' In each fold the synthesizer is trained on the training split and
#' RAPID is evaluated on the held-out test split, giving a
#' distribution of risk scores that accounts for variability in both
#' the synthesis and the evaluation.
#'
#' @param original_data A data frame of original (confidential) records.
#' @param synthesizer A function that takes a data frame (and optionally
#'   a \code{seed} argument) and returns a synthetic data frame.  If
#'   \code{NULL}, \code{synthpop::syn()} is used as a default.
#' @param quasi_identifiers Character vector of quasi-identifier column
#'   names.
#' @param sensitive_attribute Character string naming the sensitive
#'   (target) variable.
#' @param k Number of folds (default 5).
#' @param stratified Logical; use stratified folds for categorical
#'   targets?  Requires the \pkg{caret} package (default \code{TRUE}).
#' @param return_details Logical; include per-fold metrics in the
#'   result?
#' @param return_all_records Logical; include per-record risk data
#'   from every fold?
#' @param seed Optional random seed.
#' @param verbose Logical; print progress messages?
#' @param ... Additional arguments passed to \code{\link{rapid}}
#'   (e.g. \code{model_type}, \code{cat_tau}, \code{num_epsilon}).
#'
#' @return An object of class \code{"rapid_cv"} with components:
#' \describe{
#'   \item{cv_summary}{List with \code{mean}, \code{sd}, \code{se},
#'     \code{median}, \code{min}, \code{max}, \code{ci_lower},
#'     \code{ci_upper} of the fold-level RAPID scores, plus
#'     aggregated model metrics.}
#'   \item{cv_details}{Data frame of per-fold results (if
#'     \code{return_details = TRUE}).}
#'   \item{fold_data}{Concatenated per-record risk data across folds
#'     (if \code{return_all_records = TRUE}).}
#'   \item{settings}{List of CV settings (\code{k}, \code{stratified},
#'     etc.).}
#'   \item{eval_method}{Evaluation method used.}
#'   \item{model_type}{Model type used.}
#'   \item{threshold}{Threshold value used.}
#' }
#'
#' @references
#' Thees, O., Mueller, N., & Templ, M. (2026). Beyond the Trade-off Curve:
#' Multivariate and Advanced Risk-Utility Maps for Evaluating Anonymized and
#' Synthetic Data. \emph{Journal of Official Statistics}.
#'
#' @seealso \code{\link{rapid}}, \code{\link{rapid_test}}
#'
#' @family rapid
#' @author Oscar Thees, Matthias Templ
#' @examples
#' \dontrun{
#' cv <- rapid_synthesizer_cv(
#'   original_data = df,
#'   synthesizer = function(data, seed = NULL) {
#'     synthpop::syn(data, m = 1, seed = seed, print.flag = FALSE)$syn
#'   },
#'   quasi_identifiers = c("age", "sex", "region"),
#'   sensitive_attribute = "income",
#'   model_type = "rf", k = 5
#' )
#' print(cv)
#' }
#'
#' @export
rapid_synthesizer_cv <- function(original_data,
                                 synthesizer = NULL,
                                 quasi_identifiers,
                                 sensitive_attribute,
                                 k = 5,
                                 stratified = TRUE,
                                 return_details = FALSE,
                                 return_all_records = FALSE,
                                 seed = 2025,
                                 verbose = TRUE,
                                 ...) {
  if (!is.null(seed)) set.seed(seed)

  # Default synthesizer: synthpop

  if (is.null(synthesizer)) {
    if (!requireNamespace("synthpop", quietly = TRUE)) {
      stop("Package 'synthpop' is required for the default synthesizer. ",
           "Install it or provide a custom synthesizer function.")
    }
    synthesizer <- function(data, seed = NULL) {
      synthpop::syn(data, m = 1, seed = seed, print.flag = FALSE)$syn
    }
  }

  stopifnot(k >= 2, nrow(original_data) >= k, is.function(synthesizer))

  target_vec <- original_data[[sensitive_attribute]]
  is_categorical <- is.factor(target_vec) || is.character(target_vec)

  # Extract parameters from ... for the result object
  dots <- list(...)
  model_type <- if ("model_type" %in% names(dots)) dots$model_type else "rf"

  eval_method <- if (is_categorical) {
    if (!is.null(dots$cat_eval_method)) dots$cat_eval_method else "RCS_marginal"
  } else {
    if (!is.null(dots$num_error_metric)) dots$num_error_metric else "symmetric"
  }

  threshold <- if (is_categorical) {
    if (!is.null(dots$cat_tau)) dots$cat_tau else NA
  } else {
    if (!is.null(dots$num_epsilon)) dots$num_epsilon else NA
  }

  # ---- Create folds ----
  if (stratified && is_categorical) {
    folds <- .create_stratified_folds(target_vec, k)
  } else {
    folds <- split(seq_len(nrow(original_data)),
                   sample(rep(seq_len(k), length.out = nrow(original_data))))
  }

  if (verbose) {
    msg <- sprintf("  RAPID Synthesizer CV (%d folds", k)
    if (stratified && is_categorical) msg <- paste0(msg, ", stratified")
    message(msg, ")")
  }

  # ---- Run CV ----
  fold_results <- vector("list", k)
  fold_data <- if (return_all_records) vector("list", k) else NULL

  for (fold_id in seq_len(k)) {
    if (verbose) message(sprintf("  Fold %d/%d: ", fold_id, k), appendLF = FALSE)

    test_idx  <- folds[[fold_id]]
    train_idx <- setdiff(seq_len(nrow(original_data)), test_idx)

    # Synthesize from training split
    if (verbose) message(sprintf("synthesizing %d  ", length(train_idx)),
                         appendLF = FALSE)
    if (!is.null(seed)) {
      synth <- synthesizer(original_data[train_idx, ], seed = seed + fold_id)
    } else {
      synth <- synthesizer(original_data[train_idx, ])
    }

    # Evaluate RAPID on test split
    if (verbose) message(sprintf("evaluating %d... ", length(test_idx)),
                         appendLF = FALSE)

    result <- rapid(
      X = original_data[test_idx, ],
      Y = synth,
      key_vars = quasi_identifiers,
      target_var = sensitive_attribute,
      return_all_records = return_all_records,
      seed = NULL, verbose = FALSE,
      ...
    )

    conf_rate <- result$rapid

    if (verbose) message(sprintf("done (risk: %.3f)", conf_rate))

    # Store per-fold metrics
    if (is_categorical) {
      fold_results[[fold_id]] <- data.frame(
        fold = fold_id,
        rapid = conf_rate,
        accuracy = result$model_metrics$accuracy,
        stringsAsFactors = FALSE
      )
    } else {
      fold_results[[fold_id]] <- data.frame(
        fold = fold_id,
        rapid = conf_rate,
        mae  = result$model_metrics$mae,
        rmse = result$model_metrics$rmse,
        rmae = result$model_metrics$rmae,
        rrmse = result$model_metrics$rrmse,
        stringsAsFactors = FALSE
      )
    }

    if (return_all_records) {
      fold_data[[fold_id]] <- result$records
      fold_data[[fold_id]]$fold <- fold_id
    }
  }

  # ---- Aggregate ----
  fold_df <- do.call(rbind, fold_results)
  rapids  <- fold_df$rapid

  cv_summary <- list(
    mean   = mean(rapids),
    sd     = stats::sd(rapids),
    se     = stats::sd(rapids) / sqrt(k),
    median = stats::median(rapids),
    min    = min(rapids),
    max    = max(rapids),
    ci_lower = mean(rapids) - stats::qnorm(0.975) * stats::sd(rapids) / sqrt(k),
    ci_upper = mean(rapids) + stats::qnorm(0.975) * stats::sd(rapids) / sqrt(k)
  )

  if (is_categorical) {
    cv_summary$mean_accuracy <- mean(fold_df$accuracy, na.rm = TRUE)
    cv_summary$sd_accuracy   <- stats::sd(fold_df$accuracy, na.rm = TRUE)
  } else {
    cv_summary$mean_mae   <- mean(fold_df$mae, na.rm = TRUE)
    cv_summary$mean_rmse  <- mean(fold_df$rmse, na.rm = TRUE)
    cv_summary$mean_rmae  <- mean(fold_df$rmae, na.rm = TRUE)
    cv_summary$mean_rrmse <- mean(fold_df$rrmse, na.rm = TRUE)
  }

  out <- list(
    cv_summary  = cv_summary,
    cv_details  = if (return_details) fold_df else NULL,
    fold_data   = if (return_all_records) do.call(rbind, fold_data) else NULL,
    settings    = list(
      k = k,
      cv_type = "synthesizer",
      stratified = stratified,
      is_categorical = is_categorical,
      n_original = nrow(original_data),
      sensitive_attribute = sensitive_attribute
    ),
    eval_method = eval_method,
    model_type  = model_type,
    threshold   = threshold
  )
  class(out) <- "rapid_cv"

  if (verbose) {
    message(sprintf("  Mean Risk: %.4f [%.4f, %.4f]",
                    cv_summary$mean, cv_summary$ci_lower, cv_summary$ci_upper))
  }

  out
}


# ---- Internal helper: stratified folds ----

#' Create stratified k-fold indices
#' @noRd
.create_stratified_folds <- function(y, k) {
  folds <- vector("list", k)
  for (cls in unique(y)) {
    class_indices <- which(y == cls)
    shuffled <- sample(class_indices)
    fold_assignment <- rep(seq_len(k), length.out = length(class_indices))
    for (fold_id in seq_len(k)) {
      folds[[fold_id]] <- c(folds[[fold_id]],
                            shuffled[fold_assignment == fold_id])
    }
  }
  lapply(folds, sample)
}


#' @rdname rapid_synthesizer_cv
#' @param x an object of class \code{"rapid_cv"}
#' @param ... additional arguments (ignored)
#' @export
print.rapid_cv <- function(x, ...) {
  cat("\n  RAPID Cross-Validation Results\n\n")
  cat(sprintf("  Evaluation method: %s\n", x$eval_method))
  cat(sprintf("  Model type:        %s\n", x$model_type))
  cat(sprintf("  K-folds:           %d\n", x$settings$k))

  if (x$settings$is_categorical) {
    cat(sprintf("  Threshold (tau):   %s\n",
                if (!is.na(x$threshold)) round(x$threshold, 3) else "default"))
  } else {
    cat(sprintf("  Threshold (epsilon): %s\n",
                if (!is.na(x$threshold)) round(x$threshold, 3) else "default"))
  }

  cat("\n  Risk Estimate:\n")
  cat(sprintf("    Mean:   %.4f\n", x$cv_summary$mean))
  cat(sprintf("    SD:     %.4f\n", x$cv_summary$sd))
  cat(sprintf("    95%% CI: [%.4f, %.4f]\n",
              x$cv_summary$ci_lower, x$cv_summary$ci_upper))
  cat("\n")
  invisible(x)
}


#' Summary method for rapid_cv objects
#'
#' @param object an object of class "rapid_cv"
#' @param ... additional arguments (ignored)
#' @return An object of class "summary.rapid_cv"
#' @export
summary.rapid_cv <- function(object, ...) {
  summ <- list(
    eval_method = object$eval_method,
    model_type = object$model_type,
    threshold = object$threshold,
    k = object$settings$k,
    stratified = object$settings$stratified,
    is_categorical = object$settings$is_categorical,
    n_original = object$settings$n_original,
    sensitive_attribute = object$settings$sensitive_attribute,
    mean = object$cv_summary$mean,
    sd = object$cv_summary$sd,
    se = object$cv_summary$se,
    median = object$cv_summary$median,
    min = object$cv_summary$min,
    max = object$cv_summary$max,
    ci_lower = object$cv_summary$ci_lower,
    ci_upper = object$cv_summary$ci_upper,
    cv_summary = object$cv_summary,
    has_details = !is.null(object$cv_details),
    has_fold_data = !is.null(object$fold_data)
  )
  class(summ) <- "summary.rapid_cv"
  summ
}


#' Print method for summary.rapid_cv objects
#'
#' @param x an object of class "summary.rapid_cv"
#' @param ... additional arguments (ignored)
#' @export
print.summary.rapid_cv <- function(x, ...) {
  cat("Summary: RAPID Cross-Validation\n")
  cat("================================\n\n")

  cat("Configuration:\n")
  cat("  Evaluation method:    ", x$eval_method, "\n")
  cat("  Model type:           ", x$model_type, "\n")
  cat("  Sensitive attribute:  ", x$sensitive_attribute, "\n")
  cat("  K-folds:              ", x$k, "\n")
  cat("  Stratified:           ", x$stratified, "\n")
  cat("  Original data size:   ", x$n_original, "\n\n")

  cat("Risk Estimate:\n")
  cat("  Mean:   ", sprintf("%.4f", x$mean), "\n")
  cat("  Median: ", sprintf("%.4f", x$median), "\n")
  cat("  SD:     ", sprintf("%.4f", x$sd), "\n")
  cat("  SE:     ", sprintf("%.4f", x$se), "\n")
  cat("  Range:  [", sprintf("%.4f", x$min), ", ",
      sprintf("%.4f", x$max), "]\n", sep = "")
  cat("  95% CI: [", sprintf("%.4f", x$ci_lower), ", ",
      sprintf("%.4f", x$ci_upper), "]\n\n", sep = "")

  if (x$is_categorical) {
    if (!is.null(x$cv_summary$mean_accuracy)) {
      cat("Model Performance:\n")
      cat("  Mean accuracy: ", sprintf("%.4f", x$cv_summary$mean_accuracy), "\n")
    }
  } else {
    if (!is.null(x$cv_summary$mean_mae)) {
      cat("Model Performance:\n")
      cat("  Mean MAE:  ", sprintf("%.4f", x$cv_summary$mean_mae), "\n")
      cat("  Mean RMSE: ", sprintf("%.4f", x$cv_summary$mean_rmse), "\n")
    }
  }

  invisible(x)
}
