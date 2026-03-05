#' Permutation Test for RAPID Disclosure Risk
#'
#' Tests whether the observed RAPID score is significantly greater than
#' expected under the null hypothesis that the sensitive attribute is
#' independent of the quasi-identifiers.  The model is fitted once on
#' the synthetic data and predictions are computed once on the original
#' quasi-identifiers; each permutation only re-evaluates the risk score
#' against shuffled labels, making the test computationally efficient.
#'
#' @param original_data A data frame of original (confidential) records.
#' @param synthetic_data A data frame of synthetic records used to train
#'   the predictive model.
#' @param quasi_identifiers Character vector of column names used as
#'   predictors.
#' @param sensitive_attribute Character string naming the sensitive
#'   (target) variable.
#' @param model_type Model to train on synthetic data.  One of
#'   \code{"rf"} (default), \code{"lm"}, \code{"cart"}, \code{"gbm"},
#'   or \code{"logit"}.
#' @param cat_tau Categorical risk threshold; see \code{\link{rapid}}.
#' @param num_epsilon Numeric error threshold; see \code{\link{rapid}}.
#' @param cat_eval_method Categorical evaluation method; see
#'   \code{\link{rapid}}.
#' @param num_error_metric Numeric error metric; see
#'   \code{\link{rapid}}.
#' @param num_epsilon_type Numeric threshold type (\code{"percentage"}
#'   or \code{"absolute"}); see \code{\link{rapid}}.
#' @param num_delta Smoothing constant for percentage metrics
#'   (default 0.01).
#' @param n_permutations Number of permutations (default 199).
#' @param alpha Significance level (default 0.05).
#' @param seed Optional random seed.
#' @param verbose Logical; print progress?
#' @param ... Additional arguments passed to the model fitting function.
#'
#' @return An object of class \code{"rapid_test"} with components:
#' \describe{
#'   \item{statistic}{Observed RAPID score.}
#'   \item{p_value}{Permutation p-value.}
#'   \item{null_distribution}{Numeric vector of permutation RAPID values.}
#'   \item{n_permutations}{Number of permutations used.}
#'   \item{alpha}{Significance level.}
#'   \item{significant}{Logical; \code{TRUE} if \code{p_value <= alpha}.}
#'   \item{null_mean}{Mean of the null distribution.}
#'   \item{null_sd}{Standard deviation of the null distribution.}
#'   \item{null_quantiles}{Named vector of null quantiles (5th, 50th, 95th percentile).}
#'   \item{observed_rapid}{The full \code{rapid} object from the
#'     unpermuted evaluation.}
#'   \item{method}{\code{"permutation"}.}
#' }
#'
#' @references
#' Thees, O., Mueller, N., & Templ, M. (2026). Beyond the Trade-off Curve:
#' Multivariate and Advanced Risk-Utility Maps for Evaluating Anonymized and
#' Synthetic Data. \emph{Journal of Official Statistics}.
#'
#' @seealso \code{\link{rapid}}, \code{\link[=confint.rapid]{confint}},
#'   \code{\link{rapid_threshold_select}}
#'
#' @family rapid
#' @examples
#' # Small runnable example with few permutations
#' set.seed(42)
#' X <- data.frame(
#'   age = sample(20:60, 80, replace = TRUE),
#'   sex = sample(c("M", "F"), 80, replace = TRUE),
#'   income = rnorm(80, 50000, 10000)
#' )
#' Y <- X
#' Y$income <- Y$income + rnorm(80, 0, 5000)
#' res <- rapid_test(X, Y,
#'                   quasi_identifiers = c("age", "sex"),
#'                   sensitive_attribute = "income",
#'                   model_type = "lm", n_permutations = 9)
#' print(res)
#'
#' \donttest{
#' # With more permutations and random forest
#' res2 <- rapid_test(X, Y,
#'                    quasi_identifiers = c("age", "sex"),
#'                    sensitive_attribute = "income",
#'                    model_type = "rf", n_permutations = 199)
#' print(res2)
#' }
#'
#' @export
rapid_test <- function(original_data, synthetic_data,
                       quasi_identifiers, sensitive_attribute,
                       model_type = "rf",
                       cat_tau = 0.3, num_epsilon = 10,
                       cat_eval_method = "RCS_marginal",
                       num_error_metric = "symmetric",
                       num_epsilon_type = "percentage",
                       num_delta = 0.01,
                       n_permutations = 199, alpha = 0.05,
                       seed = NULL, verbose = FALSE, ...) {

  if (!is.null(seed)) set.seed(seed)

  # ---- Prepare: fit model once, predict once ----
  prep <- .rapid_prepare_permutation(
    original_data = original_data,
    synthetic_data = synthetic_data,
    quasi_identifiers = quasi_identifiers,
    sensitive_attribute = sensitive_attribute,
    model_type = model_type,
    num_delta = num_delta,
    ...
  )

  # ---- Observed RAPID ----
  rapid_obs <- .compute_rapid_permuted(
    predictions = prep$predictions,
    y = prep$y_true,
    is_categorical = prep$is_categorical,
    cat_tau = cat_tau,
    num_epsilon = num_epsilon,
    cat_eval_method = cat_eval_method,
    num_error_metric = num_error_metric,
    num_epsilon_type = num_epsilon_type,
    num_delta = num_delta,
    original_data = prep$original_prepared,
    sensitive_attribute = sensitive_attribute
  )

  # ---- Also compute the full rapid object for the result ----
  observed_rapid <- rapid(
    X = original_data, Y = synthetic_data,
    key_vars = quasi_identifiers,
    target_var = sensitive_attribute,
    model_type = model_type,
    cat_tau = cat_tau, num_epsilon = num_epsilon,
    cat_eval_method = cat_eval_method,
    num_error_metric = num_error_metric,
    num_epsilon_type = num_epsilon_type,
    num_delta = num_delta,
    seed = seed, verbose = FALSE, ...
  )

  # ---- Permutation loop ----
  perm_rapids <- numeric(n_permutations)
  for (b in seq_len(n_permutations)) {
    y_perm <- sample(prep$y_true)
    perm_rapids[b] <- .compute_rapid_permuted(
      predictions = prep$predictions,
      y = y_perm,
      is_categorical = prep$is_categorical,
      cat_tau = cat_tau,
      num_epsilon = num_epsilon,
      cat_eval_method = cat_eval_method,
      num_error_metric = num_error_metric,
      num_epsilon_type = num_epsilon_type,
      num_delta = num_delta,
      original_data = prep$original_prepared,
      sensitive_attribute = sensitive_attribute
    )
    if (verbose && b %% 50 == 0) {
      message(sprintf("  Permutation %d/%d", b, n_permutations))
    }
  }

  # ---- p-value (conservative: +1 in numerator and denominator) ----
  p_value <- (1 + sum(perm_rapids >= rapid_obs)) / (1 + n_permutations)

  result <- list(
    statistic = rapid_obs,
    p_value = p_value,
    null_distribution = perm_rapids,
    n_permutations = n_permutations,
    alpha = alpha,
    significant = p_value <= alpha,
    null_mean = mean(perm_rapids),
    null_sd = stats::sd(perm_rapids),
    null_quantiles = stats::quantile(perm_rapids, probs = c(0.05, 0.50, 0.95)),
    observed_rapid = observed_rapid,
    method = "permutation"
  )
  class(result) <- "rapid_test"
  result
}


#' @rdname rapid_test
#' @param x an object of class \code{"rapid_test"}
#' @param ... additional arguments (ignored)
#' @export
print.rapid_test <- function(x, ...) {
  cat("\n  RAPID Permutation Test\n\n")
  cat(sprintf("  Observed RAPID:  %.4f\n", x$statistic))
  cat(sprintf("  p-value:         %.4f", x$p_value))
  if (x$significant) cat(" *") else cat("  ")
  cat(sprintf("  (alpha = %.2f)\n", x$alpha))
  cat(sprintf("  Permutations:    %d\n", x$n_permutations))
  cat(sprintf("  Null mean (sd):  %.4f (%.4f)\n", x$null_mean, x$null_sd))
  cat(sprintf("  Null 95%% quantile: %.4f\n", x$null_quantiles[["95%"]]))
  if (x$significant) {
    cat("  Conclusion: RAPID is significantly above chance level.\n")
  } else {
    cat("  Conclusion: No significant disclosure risk detected.\n")
  }
  cat("\n")
  invisible(x)
}


#' Summary method for rapid_test objects
#'
#' @param object an object of class "rapid_test"
#' @param ... additional arguments (ignored)
#' @return An object of class "summary.rapid_test"
#' @export
summary.rapid_test <- function(object, ...) {
  summ <- list(
    statistic = object$statistic,
    p_value = object$p_value,
    alpha = object$alpha,
    significant = object$significant,
    n_permutations = object$n_permutations,
    null_mean = object$null_mean,
    null_sd = object$null_sd,
    null_quantiles = object$null_quantiles,
    effect_size = (object$statistic - object$null_mean) / max(object$null_sd, 1e-10),
    method = object$method
  )
  class(summ) <- "summary.rapid_test"
  summ
}


#' Print method for summary.rapid_test objects
#'
#' @param x an object of class "summary.rapid_test"
#' @param ... additional arguments (ignored)
#' @export
print.summary.rapid_test <- function(x, ...) {
  cat("Summary: RAPID Permutation Test\n")
  cat("===============================\n\n")

  cat("Test Result:\n")
  cat("  Observed RAPID: ", sprintf("%.4f", x$statistic), "\n")
  cat("  p-value:        ", sprintf("%.4f", x$p_value),
      if (x$significant) " *" else "", "\n")
  cat("  Alpha:          ", sprintf("%.2f", x$alpha), "\n")
  cat("  Significant:    ", x$significant, "\n\n")

  cat("Null Distribution (", x$n_permutations, " permutations):\n", sep = "")
  cat("  Mean:           ", sprintf("%.4f", x$null_mean), "\n")
  cat("  SD:             ", sprintf("%.4f", x$null_sd), "\n")
  cat("  Quantiles:      ", paste(sprintf("%.4f", x$null_quantiles), collapse = "  "), "\n")
  cat("  Effect size:    ", sprintf("%.2f", x$effect_size), " (standardized)\n\n")

  if (x$significant) {
    cat("Conclusion: RAPID is significantly above chance level.\n")
  } else {
    cat("Conclusion: No significant disclosure risk detected.\n")
  }
  cat("\n")
  invisible(x)
}


# =====================================================================
# Internal helpers (used by rapid_test and rapid_threshold_select)
# =====================================================================

#' Prepare data for permutation-based RAPID inference
#'
#' Fits the model once on synthetic data and computes predictions on
#' original quasi-identifiers.  Returns everything needed to quickly
#' recompute RAPID under permuted labels.
#'
#' @noRd
.rapid_prepare_permutation <- function(original_data, synthetic_data,
                                       quasi_identifiers, sensitive_attribute,
                                       model_type, num_delta = 0.01, ...) {
  # Validate inputs
  stopifnot(
    is.data.frame(original_data),
    is.data.frame(synthetic_data),
    all(quasi_identifiers %in% names(original_data)),
    sensitive_attribute %in% names(original_data)
  )

  .check_rapid_deps(model_type)

  # Handle NAs
  cleaned <- .handle_sensitive_na(original_data, synthetic_data,
                                  sensitive_attribute, "constant", 0)
  original_data <- cleaned$original
  synthetic_data <- cleaned$synthetic

  original_data <- .handle_qi_na(original_data, quasi_identifiers)
  synthetic_data <- .handle_qi_na(synthetic_data, quasi_identifiers)

  # Ensure factor levels are consistent
  y_orig <- original_data[[sensitive_attribute]]
  is_categorical <- is.factor(y_orig) || is.character(y_orig)

  if (is_categorical) {
    synthetic_data <- .ensure_levels(synthetic_data, original_data,
                                     sensitive_attribute, quasi_identifiers)
    if (is.character(original_data[[sensitive_attribute]])) {
      original_data[[sensitive_attribute]] <- factor(original_data[[sensitive_attribute]])
    }
    if (is.character(synthetic_data[[sensitive_attribute]])) {
      synthetic_data[[sensitive_attribute]] <- factor(synthetic_data[[sensitive_attribute]])
    }
  }

  # Build formula and fit model
  formula <- stats::as.formula(
    paste(sensitive_attribute, "~", paste(quasi_identifiers, collapse = " + "))
  )
  fit <- .fit_rapid_model(model_type, formula, synthetic_data,
                          original_data, sensitive_attribute, ...)

  # Predict on original X
  predictions <- .predict_rapid(model_type, fit, original_data, sensitive_attribute)

  y_true <- original_data[[sensitive_attribute]]
  if (is_categorical && !is.factor(y_true)) y_true <- factor(y_true)

  list(
    predictions = predictions,
    y_true = y_true,
    is_categorical = is_categorical,
    model = fit,
    formula = formula,
    original_prepared = original_data,
    synthetic_prepared = synthetic_data
  )
}


#' Compute RAPID score from fixed predictions and (possibly permuted) labels
#'
#' This is the inner loop of the permutation test.  It takes the
#' pre-computed prediction matrix / vector and computes the proportion
#' at risk against given labels.  Each call is O(n).
#'
#' @noRd
.compute_rapid_permuted <- function(predictions, y, is_categorical,
                                    cat_tau = 0.3, num_epsilon = 10,
                                    cat_eval_method = "RCS_marginal",
                                    num_error_metric = "symmetric",
                                    num_epsilon_type = "percentage",
                                    num_delta = 0.01,
                                    original_data = NULL,
                                    sensitive_attribute = NULL) {
  if (is_categorical) {
    .compute_rapid_cat_permuted(predictions, y, cat_tau, cat_eval_method,
                                original_data, sensitive_attribute)
  } else {
    .compute_rapid_num_permuted(predictions, y, num_epsilon,
                                num_epsilon_type, num_error_metric, num_delta)
  }
}


#' Categorical RAPID from fixed probability matrix + labels
#' @noRd
.compute_rapid_cat_permuted <- function(prob_matrix, y, cat_tau,
                                        cat_eval_method, original_data,
                                        sensitive_attribute) {
  if (is.data.frame(prob_matrix)) prob_matrix <- as.matrix(prob_matrix)

  # g_i: predicted probability for the (permuted) true class
  g_i <- prob_matrix[cbind(
    seq_len(nrow(prob_matrix)),
    match(as.character(y), colnames(prob_matrix))
  )]

  if (cat_eval_method == "RCS_conditional") {
    # Class-conditional baseline (recomputed for permuted labels)
    classes <- unique(as.character(y))
    baseline_by_class <- vapply(classes, function(k) {
      idx <- which(as.character(y) == k)
      if (length(idx) == 0) return(NA_real_)
      mean(prob_matrix[idx, k], na.rm = TRUE)
    }, numeric(1))
    names(baseline_by_class) <- classes
    b_i <- baseline_by_class[as.character(y)]
    r_i <- g_i / b_i
    at_risk <- r_i > cat_tau

  } else if (cat_eval_method == "RCS_marginal") {
    # Marginal baseline from the label distribution
    marginal_freq <- prop.table(table(y))
    b_i <- as.numeric(marginal_freq[as.character(y)])
    max_improvement <- 1 - b_i
    normalized_gain <- (g_i - b_i) / max_improvement
    at_risk <- normalized_gain > cat_tau

  } else { # NCE
    ce <- -log(pmax(g_i, 1e-10))
    n_classes <- ncol(prob_matrix)
    max_entropy <- log(n_classes)
    risk_score <- 1 - ce / max_entropy
    at_risk <- risk_score > cat_tau
  }

  sum(at_risk, na.rm = TRUE) / length(at_risk)
}


#' Numeric RAPID from fixed predicted values + labels
#' @noRd
.compute_rapid_num_permuted <- function(predicted, y, num_epsilon,
                                        num_epsilon_type, num_error_metric,
                                        num_delta) {
  A <- y
  B <- predicted

  if (num_epsilon_type == "percentage") {
    error_values <- switch(num_error_metric,
      symmetric = {
        2 * abs(A - B) / (abs(A) + abs(B) + 2 * num_delta) * 100
      },
      stabilised_relative = {
        abs(A - B) / (abs(A) + num_delta) * 100
      },
      absolute = {
        stop("Percentage threshold not meaningful for absolute error metric.")
      }
    )
  } else {
    error_values <- abs(A - B)
  }

  # Low error = high risk (model predicts accurately)
  at_risk <- error_values < num_epsilon
  sum(at_risk) / length(at_risk)
}
