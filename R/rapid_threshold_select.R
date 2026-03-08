#' Data-Driven Threshold Selection for RAPID
#'
#' Selects the largest threshold at which the observed RAPID score is
#' still significantly greater than the permutation null, following the
#' approach described in the RAPID inferential framework.  The
#' recommended threshold \eqn{\tau^*} (or \eqn{\varepsilon^*}) is the
#' maximum threshold where \code{RAPID_obs > Q_{1-alpha}(null)}.
#'
#' @inheritParams rapid_test
#' @param cat_eval_method Categorical evaluation method; see
#'   \code{\link{rapid}}.
#' @param num_error_metric Numeric error metric; see
#'   \code{\link{rapid}}.
#' @param num_epsilon_type Numeric threshold type (\code{"percentage"}
#'   or \code{"absolute"}); see \code{\link{rapid}}.
#' @param num_delta Smoothing constant for percentage metrics
#'   (default 0.01).
#' @param tau_range Numeric vector of categorical thresholds to
#'   evaluate (used when the sensitive attribute is categorical).
#' @param epsilon_range Numeric vector of numeric thresholds to
#'   evaluate (used when the sensitive attribute is continuous).
#' @param n_permutations Number of permutations (default 199).
#' @param alpha Significance level (default 0.05).
#'
#' @return An object of class \code{"rapid_threshold"} with components:
#' \describe{
#'   \item{threshold_star}{Recommended threshold, or \code{NA} if no
#'     threshold is significant.}
#'   \item{results}{Data frame with columns \code{threshold},
#'     \code{rapid_obs}, \code{null_quantile}, and
#'     \code{significant}.}
#'   \item{alpha}{Significance level used.}
#'   \item{is_categorical}{Logical; type of sensitive attribute.}
#'   \item{observed_rapid}{Full \code{rapid} object at the recommended
#'     threshold.}
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
#' # Small runnable example
#' set.seed(42)
#' X <- data.frame(
#'   age = sample(20:60, 80, replace = TRUE),
#'   sex = sample(c("M", "F"), 80, replace = TRUE),
#'   income = rnorm(80, 50000, 10000)
#' )
#' Y <- X
#' Y$income <- Y$income + rnorm(80, 0, 5000)
#' sel <- rapid_threshold_select(X, Y,
#'   quasi_identifiers = c("age", "sex"),
#'   sensitive_attribute = "income",
#'   model_type = "lm",
#'   epsilon_range = seq(5, 30, by = 5),
#'   n_permutations = 9)
#' print(sel)
#'
#' \donttest{
#' # With random forest and finer grid
#' sel2 <- rapid_threshold_select(X, Y,
#'   quasi_identifiers = c("age", "sex"),
#'   sensitive_attribute = "income",
#'   model_type = "rf", n_permutations = 199)
#' print(sel2)
#' plot(sel2)
#' }
#'
#' @export
rapid_threshold_select <- function(original_data, synthetic_data,
                                   quasi_identifiers, sensitive_attribute,
                                   model_type = "rf",
                                   cat_eval_method = "RCS_marginal",
                                   num_error_metric = "symmetric",
                                   num_epsilon_type = "percentage",
                                   num_delta = 0.01,
                                   tau_range = seq(0.05, 0.95, by = 0.05),
                                   epsilon_range = seq(1, 50, by = 2),
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

  # Choose threshold grid based on variable type
  if (prep$is_categorical) {
    thresholds <- tau_range
  } else {
    thresholds <- epsilon_range
  }
  n_thresh <- length(thresholds)

  # ---- Observed RAPID across thresholds ----
  obs_rapids <- vapply(thresholds, function(th) {
    if (prep$is_categorical) {
      .compute_rapid_permuted(
        prep$predictions, prep$y_true, TRUE,
        cat_tau = th, cat_eval_method = cat_eval_method,
        original_data = prep$original_prepared,
        sensitive_attribute = sensitive_attribute
      )
    } else {
      .compute_rapid_permuted(
        prep$predictions, prep$y_true, FALSE,
        num_epsilon = th, num_error_metric = num_error_metric,
        num_epsilon_type = num_epsilon_type, num_delta = num_delta
      )
    }
  }, numeric(1))

  # ---- Permutation null distributions across thresholds ----
  perm_matrix <- matrix(NA_real_, nrow = n_permutations, ncol = n_thresh)

  for (b in seq_len(n_permutations)) {
    y_perm <- sample(prep$y_true)
    for (j in seq_len(n_thresh)) {
      th <- thresholds[j]
      if (prep$is_categorical) {
        perm_matrix[b, j] <- .compute_rapid_permuted(
          prep$predictions, y_perm, TRUE,
          cat_tau = th, cat_eval_method = cat_eval_method,
          original_data = prep$original_prepared,
          sensitive_attribute = sensitive_attribute
        )
      } else {
        perm_matrix[b, j] <- .compute_rapid_permuted(
          prep$predictions, y_perm, FALSE,
          num_epsilon = th, num_error_metric = num_error_metric,
          num_epsilon_type = num_epsilon_type, num_delta = num_delta
        )
      }
    }
    if (verbose && b %% 50 == 0) {
      message(sprintf("  Permutation %d/%d", b, n_permutations))
    }
  }

  # ---- Compute (1-alpha)-quantile of null at each threshold ----
  null_quantiles <- apply(perm_matrix, 2, stats::quantile,
                          probs = 1 - alpha, na.rm = TRUE)
  significant <- obs_rapids > null_quantiles

  results <- data.frame(
    threshold = thresholds,
    rapid_obs = obs_rapids,
    null_quantile = null_quantiles,
    significant = significant
  )

  # ---- Recommended threshold: max threshold that is still significant ----
  if (any(significant)) {
    threshold_star <- max(thresholds[significant])
  } else {
    threshold_star <- NA_real_
  }

  # ---- Compute full rapid object at recommended threshold ----
  if (!is.na(threshold_star)) {
    rapid_args <- list(
      X = original_data, Y = synthetic_data,
      key_vars = quasi_identifiers,
      target_var = sensitive_attribute,
      model_type = model_type,
      cat_eval_method = cat_eval_method,
      num_error_metric = num_error_metric,
      num_epsilon_type = num_epsilon_type,
      num_delta = num_delta,
      seed = seed, verbose = FALSE
    )
    if (prep$is_categorical) {
      rapid_args$cat_tau <- threshold_star
    } else {
      rapid_args$num_epsilon <- threshold_star
    }
    observed_rapid <- do.call(rapid, c(rapid_args, list(...)))
  } else {
    observed_rapid <- NULL
  }

  out <- list(
    threshold_star = threshold_star,
    results = results,
    alpha = alpha,
    is_categorical = prep$is_categorical,
    observed_rapid = observed_rapid
  )
  class(out) <- "rapid_threshold"
  out
}


#' @rdname rapid_threshold_select
#' @param x an object of class \code{"rapid_threshold"}
#' @param ... additional arguments (ignored)
#' @export
print.rapid_threshold <- function(x, ...) {
  cat("\n  RAPID Data-Driven Threshold Selection\n\n")
  type_label <- if (x$is_categorical) "tau" else "epsilon"
  if (!is.na(x$threshold_star)) {
    cat(sprintf("  Recommended %s*: %.3f  (alpha = %.2f)\n",
                type_label, x$threshold_star, x$alpha))
    cat(sprintf("  RAPID at %s*:    %.4f\n",
                type_label, x$observed_rapid$rapid))
  } else {
    cat(sprintf("  No threshold significant at alpha = %.2f\n", x$alpha))
  }
  n_sig <- sum(x$results$significant)
  cat(sprintf("  Significant thresholds: %d / %d evaluated\n",
              n_sig, nrow(x$results)))
  cat("\n")
  invisible(x)
}


#' Summary method for rapid_threshold objects
#'
#' @param object an object of class "rapid_threshold"
#' @param ... additional arguments (ignored)
#' @return An object of class "summary.rapid_threshold"
#' @export
summary.rapid_threshold <- function(object, ...) {
  n_sig <- sum(object$results$significant)
  summ <- list(
    threshold_star = object$threshold_star,
    alpha = object$alpha,
    is_categorical = object$is_categorical,
    n_thresholds = nrow(object$results),
    n_significant = n_sig,
    pct_significant = 100 * n_sig / nrow(object$results),
    rapid_at_star = if (!is.na(object$threshold_star)) object$observed_rapid$rapid else NA_real_,
    threshold_range = range(object$results$threshold),
    results = object$results
  )
  class(summ) <- "summary.rapid_threshold"
  summ
}


#' Print method for summary.rapid_threshold objects
#'
#' @param x an object of class "summary.rapid_threshold"
#' @param ... additional arguments (ignored)
#' @export
print.summary.rapid_threshold <- function(x, ...) {
  cat("Summary: RAPID Threshold Selection\n")
  cat("===================================\n\n")

  type_label <- if (x$is_categorical) "tau" else "epsilon"

  cat("Threshold Grid:\n")
  cat("  Range:      [", sprintf("%.3f", x$threshold_range[1]), ", ",
      sprintf("%.3f", x$threshold_range[2]), "]\n", sep = "")
  cat("  Evaluated:  ", x$n_thresholds, "\n")
  cat("  Significant:", x$n_significant, sprintf("(%.1f%%)", x$pct_significant), "\n")
  cat("  Alpha:      ", sprintf("%.2f", x$alpha), "\n\n")

  if (!is.na(x$threshold_star)) {
    cat("Recommended ", type_label, "*: ", sprintf("%.3f", x$threshold_star), "\n", sep = "")
    cat("RAPID at ", type_label, "*:    ", sprintf("%.4f", x$rapid_at_star), "\n\n", sep = "")
  } else {
    cat("No threshold significant at alpha = ", sprintf("%.2f", x$alpha), "\n\n", sep = "")
  }

  cat("Results:\n")
  print(x$results, row.names = FALSE)
  cat("\n")
  invisible(x)
}


#' @rdname rapid_threshold_select
#' @export
plot.rapid_threshold <- function(x, ...) {
  df <- x$results

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$threshold)) +
    ggplot2::geom_line(ggplot2::aes(y = .data$rapid_obs),
                       color = "steelblue", linewidth = 1) +
    ggplot2::geom_line(ggplot2::aes(y = .data$null_quantile),
                       color = "grey50", linetype = "dashed", linewidth = 0.8) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      x = if (x$is_categorical) expression(tau) else expression(epsilon),
      y = "RAPID",
      title = "Threshold Selection",
      subtitle = sprintf("Solid = observed, dashed = null %s%% quantile",
                          round((1 - x$alpha) * 100))
    )

  if (!is.na(x$threshold_star)) {
    p <- p + ggplot2::geom_vline(xintercept = x$threshold_star,
                                  color = "firebrick", linetype = "dotted",
                                  linewidth = 0.8)
  }

  p
}
