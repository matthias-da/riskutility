#' Permutation Test for RAPID Disclosure Risk
#'
#' Tests whether the observed RAPID is greater than expected when the sensitive
#' variable carries no information about the records' quasi-identifiers. The
#' attacker is fitted once; each permutation scores its fixed predictions
#' against randomly reassigned true values under the rule of the object (the
#' marginal baseline is that of the reassigned value), so the test needs no
#' refit.
#'
#' \code{rapid_test()} works on a fitted \code{\link{rapid}} object or on the
#' data. On the data it first calls \code{rapid()} with the arguments in
#' \code{...} and \code{n_perm = 0}, then tests that object. The scoring rule is
#' part of the object: change it with \code{\link[=update.rapid]{update}()}.
#' \code{rapid_test(rapid(X, Y, ..., seed = s, n_perm = 0), seed = s)} equals
#' \code{rapid_test(X, Y, ..., seed = s)}.
#'
#' @param X a \code{rapid} object (with stored predictions), the original data,
#'   or a \code{\link{synth_pair}}
#' @param Y data frame of synthetic data (not needed if X is a synth_pair)
#' @param key_vars character vector of quasi-identifier variable names
#' @param target_var character, name of the sensitive target variable
#' @param ... for the data methods, arguments of \code{\link{rapid}} (model
#'   and rule; the old argument names \code{original_data},
#'   \code{synthetic_data}, \code{quasi_identifiers} and
#'   \code{sensitive_attribute} are accepted). Not allowed for a \code{rapid}
#'   object.
#' @param n_permutations number of permutations (default 199)
#' @param alpha significance level (default 0.05)
#' @param seed random seed of the permutations, and of the fit for the data
#'   methods (default NULL: the random stream is left alone)
#'
#' @return An object of class \code{"rapid_test"} with components:
#' \describe{
#'   \item{statistic, rapid}{Observed RAPID of the object.}
#'   \item{p_value}{Permutation p-value, \eqn{(1 + \#\{null \ge obs\}) / (1 + B)}.}
#'   \item{null_distribution}{RAPID of each permutation.}
#'   \item{n_permutations, alpha}{Settings.}
#'   \item{significant}{Logical; \code{TRUE} if \code{p_value <= alpha}.}
#'   \item{null_mean, chance}{Mean of the null distribution (the chance level).}
#'   \item{excess}{\code{rapid - chance}.}
#'   \item{null_sd}{Standard deviation of the null distribution.}
#'   \item{null_quantiles}{Null quantiles (5, 50, 95 percent).}
#'   \item{observed_rapid}{The \code{rapid} object that was tested.}
#'   \item{method}{\code{"permutation"}.}
#' }
#'
#' @references
#' Thees, O., Mueller, N., & Templ, M. (2026). Beyond the Trade-off Curve:
#' Multivariate and Advanced Risk-Utility Maps for Evaluating Anonymized and
#' Synthetic Data. \emph{Journal of Official Statistics}.
#'
#' @seealso \code{\link{rapid}}, \code{\link[=confint.rapid]{confint}},
#'   \code{\link{rapid_threshold_select}}, \code{\link[=update.rapid]{update}}
#'
#' @family rapid
#' @author Oscar Thees, Matthias Templ
#' @examples
#' set.seed(42)
#' X <- data.frame(
#'   age = sample(20:60, 80, replace = TRUE),
#'   sex = sample(c("M", "F"), 80, replace = TRUE),
#'   income = rnorm(80, 50000, 10000)
#' )
#' X$income <- X$income + 500 * X$age
#' Y <- X
#' Y$income <- Y$income + rnorm(80, 0, 5000)
#'
#' # On a fitted object: no refit
#' r <- rapid(X, Y, key_vars = c("age", "sex"), target_var = "income",
#'            model_type = "lm")
#' res <- rapid_test(r, n_permutations = 99, seed = 1)
#' res
#' summary(res)
#'
#' # On the data: arguments in ... go to rapid()
#' rapid_test(X, Y, key_vars = c("age", "sex"), target_var = "income",
#'            model_type = "lm", n_permutations = 99, seed = 1)
#'
#' \donttest{
#' plot(res)
#' }
#'
#' @export
rapid_test <- function(X, ...) {
  if (missing(X)) {
    # Without X the data arrive under their alias names
    cl <- sys.call()
    cl[[1L]] <- rapid_test.default
    return(eval(cl, parent.frame()))
  }
  UseMethod("rapid_test")
}

#' @rdname rapid_test
#' @export
rapid_test.rapid <- function(X, n_permutations = 199, alpha = 0.05, seed = NULL, ...) {
  .rapid_check_empty_dots(list(...), .rapid_object_hint)
  .rapid_require(X, c("version", "predictions"))
  .rapid_check_test_args(n_permutations, alpha)

  if (!is.null(seed)) set.seed(seed)
  null <- .rapid_null_rates(X$rule, .rapid_null(X$attack, X$rule, n_permutations))
  statistic <- X$rapid
  p_value <- (1 + sum(null >= statistic)) / (1 + n_permutations)
  chance <- mean(null)

  structure(list(
    statistic = statistic,
    rapid = statistic,
    p_value = p_value,
    null_distribution = null,
    n_permutations = n_permutations,
    alpha = alpha,
    significant = p_value <= alpha,
    null_mean = chance,
    chance = chance,
    excess = statistic - chance,
    null_sd = if (n_permutations > 1) stats::sd(null) else NA_real_,
    null_quantiles = stats::quantile(null, probs = c(0.05, 0.50, 0.95)),
    observed_rapid = X,
    method = "permutation"
  ), class = "rapid_test")
}

#' @rdname rapid_test
#' @export
rapid_test.synth_pair <- function(X, ...) {
  # The default method unpacks the pair through rapid()
  NextMethod()
}

#' @rdname rapid_test
#' @export
rapid_test.default <- function(X, Y, key_vars, target_var, ..., n_permutations = 199,
                               alpha = 0.05, seed = NULL) {
  fit <- .rapid_fit_for(match.call(), parent.frame())
  rapid_test(fit, n_permutations = n_permutations, alpha = alpha, seed = seed)
}

#' @rdname rapid_test
#' @param x an object of class \code{"rapid_test"}
#' @param digits number of decimal places
#' @export
format.rapid_test <- function(x, digits = 3, ...) {
  f <- function(v) formatC(v, digits = digits, format = "f")
  c(
    paste0("RAPID permutation test (", x$n_permutations, " permutations)"),
    paste0("RAPID = ", f(x$rapid), ", chance level = ", f(x$chance), ", excess = ",
           f(x$excess)),
    paste0("p-value = ", formatC(x$p_value, digits = 4, format = "f"), " (alpha = ",
           format(x$alpha), "): RAPID is ",
           if (x$significant) "significantly above" else "not significantly above",
           " the chance level")
  )
}

#' @rdname rapid_test
#' @export
print.rapid_test <- function(x, ...) {
  cat(format(x, ...), sep = "\n")
  invisible(x)
}

#' @rdname rapid_test
#' @param object an object of class \code{"rapid_test"}
#' @export
summary.rapid_test <- function(object, ...) {
  structure(list(
    statistic = object$statistic,
    rapid = object$rapid,
    p_value = object$p_value,
    alpha = object$alpha,
    significant = object$significant,
    n_permutations = object$n_permutations,
    null_mean = object$null_mean,
    chance = object$chance,
    excess = object$excess,
    null_sd = object$null_sd,
    null_quantiles = object$null_quantiles,
    effect_size = (object$statistic - object$null_mean) / max(object$null_sd, 1e-10),
    rule = format(object$observed_rapid$rule),
    method = object$method
  ), class = "summary.rapid_test")
}

#' @rdname rapid_test
#' @export
format.summary.rapid_test <- function(x, ...) {
  f <- function(v) formatC(v, digits = 4, format = "f")
  c(
    "Summary: RAPID Permutation Test",
    "===============================",
    paste0("Rule:           ", x$rule),
    paste0("Observed RAPID: ", f(x$rapid)),
    paste0("p-value:        ", f(x$p_value), if (x$significant) " *"),
    paste0("Alpha:          ", format(x$alpha)),
    paste0("Null distribution (", x$n_permutations, " permutations):"),
    paste0("  Mean (chance): ", f(x$chance)),
    paste0("  SD:            ", f(x$null_sd)),
    paste0("  Quantiles:     ", paste(f(x$null_quantiles), collapse = "  "),
           " (5%, 50%, 95%)"),
    paste0("  Excess:        ", f(x$excess)),
    paste0("  Effect size:   ", formatC(x$effect_size, digits = 2, format = "f"),
           " (standardized)"),
    if (x$significant) "Conclusion: RAPID is significantly above the chance level."
    else "Conclusion: RAPID is not significantly above the chance level."
  )
}

#' @rdname rapid_test
#' @export
print.summary.rapid_test <- function(x, ...) {
  cat(format(x, ...), sep = "\n")
  invisible(x)
}

#' @rdname rapid_test
#' @param bins number of histogram bins
#' @export
plot.rapid_test <- function(x, ..., bins = 30) {
  q <- stats::quantile(x$null_distribution, 1 - x$alpha, names = FALSE)
  p <- ggplot2::ggplot(data.frame(rapid = x$null_distribution),
                       ggplot2::aes(x = .data$rapid)) +
    ggplot2::geom_histogram(bins = bins, fill = "grey70", colour = "white") +
    ggplot2::geom_vline(xintercept = q, colour = "grey30", linetype = "dashed") +
    ggplot2::geom_vline(xintercept = x$statistic, colour = "firebrick", linewidth = 1) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      x = "RAPID", y = "Permutations",
      title = "RAPID permutation test",
      subtitle = sprintf("Solid = observed (p = %.3f), dashed = null %s%% quantile",
                         x$p_value, format(100 * (1 - x$alpha)))
    )
  print(p)
  invisible(p)
}

#' @rdname rapid_test
#' @param row.names,optional not used
#' @export
as.data.frame.rapid_test <- function(x, row.names = NULL, optional = FALSE, ...) {
  data.frame(permutation = seq_along(x$null_distribution), rapid = x$null_distribution)
}


# =====================================================================
# Helpers shared with rapid_threshold_select()
# =====================================================================

.rapid_object_hint <- "The scoring rule is part of the object; change it with update()."

#' Fit the rapid object for a data method of rapid_test() or rapid_threshold_select()
#'
#' Turns the matched call of the data method into a call of `rapid()` with the
#' caller's expressions and `n_perm = 0`, so that the object stores a call that
#' can be evaluated again. `seed` stays in the call and seeds the fit.
#' @param call `match.call()` of the data method.
#' @param env The caller's frame.
#' @noRd
.rapid_fit_for <- function(call, env) {
  call <- call[!names(call) %in% c("n_permutations", "alpha", "tau_range", "epsilon_range")]
  call[[1L]] <- rapid
  call$n_perm <- 0
  eval(call, env)
}

#' RAPID of each permutation: share of flagged records per column of the null scores
#' @noRd
.rapid_null_rates <- function(rule, null) {
  flags <- .rule_flag(rule, null)
  dim(flags) <- dim(null)
  colMeans(flags)
}

#' @noRd
.rapid_check_test_args <- function(n_permutations, alpha) {
  if (!is.numeric(n_permutations) || length(n_permutations) != 1 ||
      is.na(n_permutations) || n_permutations < 1 ||
      n_permutations != round(n_permutations)) {
    .rapid_abort("'n_permutations' must be a single whole number of at least 1.", "input")
  }
  if (!is.numeric(alpha) || length(alpha) != 1 || is.na(alpha) || alpha <= 0 ||
      alpha >= 1) {
    .rapid_abort("'alpha' must be a single number between 0 and 1.", "input")
  }
  invisible(TRUE)
}


# =====================================================================
# 0.2.0 internals, used only by rapid_qi_contribution() until its rebuild
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
