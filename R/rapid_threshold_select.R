#' Data-Driven Threshold Selection for RAPID
#'
#' Scores the attack of a \code{rapid} object over a grid of thresholds and
#' compares the observed RAPID at each threshold with the \code{1 - alpha}
#' quantile of its permutation null distribution. The scores under permuted
#' true values are computed once, for all thresholds, and the model is not
#' fitted again.
#'
#' Where the rule flags a record when its score is above the threshold (band,
#' \code{RCS_marginal}, \code{RCS_conditional}, \code{NCE}), the grid is
#' \code{tau_range} and the recommended threshold \code{threshold_star} is the
#' largest significant one, the strictest threshold at which RAPID is still
#' above chance. Where it flags a record when its score is below the threshold
#' (point rule: the score is the error, the threshold the tolerance), the grid
#' is \code{epsilon_range} and \code{threshold_star} is the smallest
#' significant one, again the strictest. For the band rule the tolerance stays
#' that of the object and the grid runs over the gain threshold \code{num_tau}.
#'
#' On the data, \code{rapid_threshold_select()} first calls \code{rapid()} with
#' the arguments in \code{...} and \code{n_perm = 0}, as \code{\link{rapid_test}}
#' does. The grid that the rule does not use (\code{epsilon_range} for a rule
#' that flags scores above the threshold, \code{tau_range} for the point rule)
#' is ignored with a warning of class \code{"rapid_warning_ignored_argument"}
#' when it is given.
#'
#' A threshold counts as significant when the observed RAPID is above the
#' \code{1 - alpha} quantile (type 7 of \code{\link[stats]{quantile}()}) of the
#' RAPID values of the permutations. \code{\link{rapid_test}()} uses the
#' permutation p-value instead and calls RAPID significant when \eqn{p \le
#' \alpha}{p <= alpha}. The two criteria can disagree at the same threshold,
#' because the quantile is interpolated between permuted values (see
#' \code{\link{rapid_test}} for an example). Each threshold of the grid is tested
#' at level \code{alpha} on its own. There is no correction for testing many
#' thresholds, so under the null hypothesis the chance that at least one
#' threshold of the grid comes out significant can exceed \code{alpha}. Read
#' \code{threshold_star} as a description of how strict a threshold the release
#' still fails, and test the threshold you settle on with \code{rapid_test()}.
#'
#' @inheritParams rapid_test
#' @param ... for the data methods, arguments of \code{\link{rapid}} (model
#'   and rule, the old argument names are accepted). Not allowed for a
#'   \code{rapid} object. Ignored by the methods for \code{rapid_threshold}
#'   objects.
#' @param tau_range numeric vector of thresholds for rules that flag scores
#'   above the threshold (default \code{seq(0.05, 0.95, by = 0.05)})
#' @param epsilon_range numeric vector of tolerances for the point rule (a
#'   proportion for the relative error metrics, units of the target for the
#'   absolute one, default \code{seq(0.01, 0.5, by = 0.02)})
#'
#' @return An object of class \code{"rapid_threshold"} with components:
#' \describe{
#'   \item{threshold_star}{Recommended threshold, or \code{NA} if no
#'     threshold is significant.}
#'   \item{results}{Data frame with columns \code{threshold}, \code{rapid_obs},
#'     \code{chance} (mean of the null distribution), \code{null_quantile} and
#'     \code{significant} (\code{rapid_obs > null_quantile}).}
#'   \item{rapid, chance, excess}{Observed RAPID, chance level and their
#'     difference at \code{threshold_star} (\code{NA} if there is none).}
#'   \item{alpha, n_permutations}{Settings.}
#'   \item{rule}{The rule of the object whose threshold was varied.}
#'   \item{is_categorical}{Logical, the type of the sensitive variable.}
#'   \item{observed_rapid}{The object rescored at \code{threshold_star}, with
#'     \code{n_perm = 0} (\code{NULL} if there is none).}
#' }
#' \code{format()} returns the printed lines and \code{print()} the object,
#' invisibly. \code{summary()} returns an object of class
#' \code{"summary.rapid_threshold"} with the number and share of significant
#' thresholds and the table \code{results}. \code{plot()} draws the observed
#' RAPID and the null quantile against the threshold, marks
#' \code{threshold_star} and returns the ggplot object, invisibly.
#' \code{as.data.frame()} returns \code{results}.
#'
#' @references
#' Templ, M., Thees, O. and Müller, R. (2026). RAPID: Risk of Attribute
#' Prediction-Induced Disclosure. Manuscript under review at EPJ Data Science.
#'
#' @seealso \code{\link{rapid}}, \code{\link{rapid_test}}
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
#' r <- rapid(X, Y, key_vars = c("age", "sex"), target_var = "income",
#'            model_type = "lm")
#'
#' # Band rule: the grid runs over the gain threshold
#' sel <- rapid_threshold_select(r, tau_range = seq(0.1, 0.9, by = 0.1),
#'                               n_permutations = 99, seed = 1)
#' sel
#' summary(sel)
#' as.data.frame(sel)
#' plot(sel)
#'
#' # Point rule: the grid runs over the tolerance
#' sel_point <- rapid_threshold_select(update(r, num_eval_method = "point"),
#'                                     n_permutations = 99, seed = 1)
#' sel_point
#'
#' # On the data, the arguments in ... go to rapid()
#' rapid_threshold_select(X, Y, key_vars = c("age", "sex"), target_var = "income",
#'                        model_type = "lm", n_permutations = 99, seed = 1)
#'
#' @export
rapid_threshold_select <- function(X, ...) {
  if (missing(X)) {
    # Without X the data arrive under their alias names
    cl <- sys.call()
    cl[[1L]] <- rapid_threshold_select.default
    return(eval(cl, parent.frame()))
  }
  UseMethod("rapid_threshold_select")
}

#' @rdname rapid_threshold_select
#' @export
rapid_threshold_select.rapid <- function(X, tau_range = seq(0.05, 0.95, by = 0.05),
                                         epsilon_range = seq(0.01, 0.5, by = 0.02),
                                         n_permutations = 199, alpha = 0.05,
                                         seed = NULL, ...) {
  .rapid_check_empty_dots(list(...), .rapid_object_hint)
  .rapid_require(X, c("version", "predictions"))
  .rapid_check_test_args(n_permutations, alpha)
  rule <- X$rule
  above <- rule$direction == "above"
  used <- if (above) "tau_range" else "epsilon_range"
  unused <- if (above) "epsilon_range" else "tau_range"
  if (if (above) !missing(epsilon_range) else !missing(tau_range)) {
    .rapid_warn(c("Argument '", unused, "' is ignored: the rule \"", rule$method,
                  "\" flags scores ", if (above) "above" else "below", " the threshold, ",
                  "so the grid is '", used, "'."), "ignored_argument")
  }
  thresholds <- if (above) tau_range else epsilon_range
  if (!is.numeric(thresholds) || length(thresholds) == 0 || anyNA(thresholds)) {
    .rapid_abort(c("'", if (above) "tau_range" else "epsilon_range",
                   "' must be a numeric vector without missing values."), "input")
  }

  # One score matrix under permuted true values serves every threshold: the
  # score of a rule does not depend on its threshold.
  if (!is.null(seed)) set.seed(seed)
  null <- .rapid_null(X$attack, rule, n_permutations)
  per_threshold <- vapply(thresholds, function(th) {
    r <- .rule_set_threshold(rule, th)
    rates <- .rapid_null_rates(r, null)
    c(rapid_obs = mean(.rule_flag(r, X$score)), chance = mean(rates),
      null_quantile = stats::quantile(rates, 1 - alpha, names = FALSE))
  }, numeric(3))
  results <- data.frame(threshold = thresholds,
                        rapid_obs = per_threshold["rapid_obs", ],
                        chance = per_threshold["chance", ],
                        null_quantile = per_threshold["null_quantile", ])
  results$significant <- results$rapid_obs > results$null_quantile

  star <- if (!any(results$significant)) {
    NA_real_
  } else if (above) {
    max(thresholds[results$significant])
  } else {
    min(thresholds[results$significant])
  }
  at_star <- if (is.na(star)) NULL else results[match(star, thresholds), ]

  structure(list(
    threshold_star = star,
    results = results,
    rapid = if (is.na(star)) NA_real_ else at_star$rapid_obs,
    chance = if (is.na(star)) NA_real_ else at_star$chance,
    excess = if (is.na(star)) NA_real_ else at_star$rapid_obs - at_star$chance,
    alpha = alpha,
    n_permutations = n_permutations,
    rule = rule,
    is_categorical = X$is_categorical,
    observed_rapid = if (!is.na(star)) .rapid_at_threshold(X, star)
  ), class = "rapid_threshold")
}

#' @rdname rapid_threshold_select
#' @export
rapid_threshold_select.synth_pair <- function(X, ...) {
  # The default method unpacks the pair through rapid()
  NextMethod()
}

#' @rdname rapid_threshold_select
#' @export
rapid_threshold_select.default <- function(X, Y, key_vars, target_var, ...,
                                           tau_range = seq(0.05, 0.95, by = 0.05),
                                           epsilon_range = seq(0.01, 0.5, by = 0.02),
                                           n_permutations = 199, alpha = 0.05,
                                           seed = NULL) {
  fit <- .rapid_fit_for(match.call(), parent.frame())
  # Only the grids that were given, so that the object method can warn about one
  # the rule does not use
  args <- list(fit, n_permutations = n_permutations, alpha = alpha, seed = seed)
  if (!missing(tau_range)) args$tau_range <- tau_range
  if (!missing(epsilon_range)) args$epsilon_range <- epsilon_range
  do.call(rapid_threshold_select, args)
}

#' The object rescored at a threshold, with the threshold in its call
#'
#' No permutations (`n_perm = 0`): the chance level at the threshold is part of
#' the threshold selection.
#' @noRd
.rapid_at_threshold <- function(x, threshold) {
  call <- x$call
  if (x$is_categorical) {
    call$cat_tau <- threshold
  } else if (x$rule$direction == "above") {
    call$num_tau <- threshold
  } else {
    # The threshold is the tolerance, on the scale of the rule's metric
    call$num_epsilon <- threshold
    call$num_error_metric <- x$rule$metric
    call$num_epsilon_type <- NULL
  }
  call$threshold_rule <- "fixed"
  call$n_perm <- 0
  .rapid_rescore(x, .rule_set_threshold(x$rule, threshold), n_perm = 0,
                 threshold_rule = "fixed", perm_level = x$perm_level, call = call)
}

#' @rdname rapid_threshold_select
#' @param x an object of class \code{"rapid_threshold"}
#' @param digits number of decimal places
#' @export
format.rapid_threshold <- function(x, digits = 3, ...) {
  f <- function(v) formatC(v, digits = digits, format = "f")
  label <- .rapid_threshold_label(x$rule)
  c(
    paste0("RAPID threshold selection over ", nrow(x$results), " values of ", label,
           " (", x$n_permutations, " permutations, alpha = ", format(x$alpha), ")"),
    if (is.na(x$threshold_star)) {
      "No threshold is significant."
    } else {
      c(paste0("Recommended ", label, " = ", format(x$threshold_star)),
        paste0("RAPID = ", f(x$rapid), ", chance level = ", f(x$chance), ", excess = ",
               f(x$excess)))
    },
    paste0("Significant: ", sum(x$results$significant), " of ", nrow(x$results),
           " thresholds")
  )
}

#' @rdname rapid_threshold_select
#' @export
print.rapid_threshold <- function(x, ...) {
  cat(format(x, ...), sep = "\n")
  invisible(x)
}

#' @rdname rapid_threshold_select
#' @param object an object of class \code{"rapid_threshold"}
#' @export
summary.rapid_threshold <- function(object, ...) {
  n_sig <- sum(object$results$significant)
  structure(list(
    threshold_star = object$threshold_star,
    label = .rapid_threshold_label(object$rule),
    rule = format(object$rule),
    alpha = object$alpha,
    n_permutations = object$n_permutations,
    is_categorical = object$is_categorical,
    n_thresholds = nrow(object$results),
    n_significant = n_sig,
    pct_significant = 100 * n_sig / nrow(object$results),
    rapid = object$rapid,
    chance = object$chance,
    excess = object$excess,
    threshold_range = range(object$results$threshold),
    results = object$results
  ), class = "summary.rapid_threshold")
}

#' @rdname rapid_threshold_select
#' @export
format.summary.rapid_threshold <- function(x, ...) {
  f <- function(v) formatC(v, digits = 4, format = "f")
  c(
    "Summary: RAPID Threshold Selection",
    "==================================",
    paste0("Rule of the object: ", x$rule),
    paste0("Grid of ", x$label, ": [", format(x$threshold_range[1]), ", ",
           format(x$threshold_range[2]), "], ", x$n_thresholds, " values"),
    paste0("Significant: ", x$n_significant, sprintf(" (%.1f%%)", x$pct_significant),
           " at alpha = ", format(x$alpha), ", ", x$n_permutations, " permutations"),
    if (is.na(x$threshold_star)) {
      "No threshold is significant."
    } else {
      c(paste0("Recommended ", x$label, " = ", format(x$threshold_star)),
        paste0("RAPID = ", f(x$rapid), ", chance level = ", f(x$chance), ", excess = ",
               f(x$excess)))
    },
    "Results:",
    utils::capture.output(print(x$results, row.names = FALSE))
  )
}

#' @rdname rapid_threshold_select
#' @export
print.summary.rapid_threshold <- function(x, ...) {
  cat(format(x, ...), sep = "\n")
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
      x = .rapid_threshold_label(x$rule),
      y = "RAPID",
      title = "Threshold Selection",
      subtitle = sprintf("Solid = observed, dashed = null %s%% quantile",
                         format(100 * (1 - x$alpha)))
    )
  if (!is.na(x$threshold_star)) {
    p <- p + ggplot2::geom_vline(xintercept = x$threshold_star,
                                 color = "firebrick", linetype = "dotted",
                                 linewidth = 0.8)
  }
  print(p)
  invisible(p)
}

#' @rdname rapid_threshold_select
#' @param row.names,optional not used
#' @export
as.data.frame.rapid_threshold <- function(x, row.names = NULL, optional = FALSE, ...) {
  x$results
}

# Name of the varied threshold: tolerance where the score is the error
.rapid_threshold_label <- function(rule) {
  if (rule$direction == "above") "tau" else "epsilon"
}
