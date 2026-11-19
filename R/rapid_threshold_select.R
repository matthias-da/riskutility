#' Data-Driven Threshold Selection for RAPID
#'
#' Answers the question how strict a threshold a release still fails. The
#' attack of a \code{rapid} object is scored over a grid of thresholds, and at
#' each one the observed RAPID is compared with the \code{1 - alpha} quantile
#' of its permutation null distribution. \code{threshold_star} is the
#' strictest threshold at which RAPID is still above that quantile. The
#' selection is exploratory. The thresholds are tested one by one without a
#' correction for their number, and a test at \code{threshold_star} on the
#' same records is no confirmation (see Details). The scores under permuted
#' true values are computed once, for all thresholds, and the model is not
#' fitted again.
#'
#' Where a record is at risk (flagged) when its score is above the threshold
#' (band, \code{RCS_marginal}, \code{RCS_conditional}, \code{NCE}), the grid
#' is \code{tau_range} and \code{threshold_star} is the largest significant
#' threshold, the strictest one at which RAPID is still above chance. Where a
#' record is at risk when its score is below the threshold (point rule: the
#' score is the error, the threshold the tolerance), the grid is
#' \code{epsilon_range} and \code{threshold_star} is the smallest significant
#' one, again the strictest. For the band rule the tolerance stays that of the
#' object and the grid runs over the gain threshold \code{num_tau}.
#'
#' On the data, \code{rapid_threshold_select()} first calls \code{rapid()} with
#' the arguments in \code{...} and \code{n_perm = 0}, as \code{\link{rapid_test}}
#' does. The grid that the rule does not use (\code{tau_range} for the point
#' rule, \code{epsilon_range} for every other rule)
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
#' threshold of the grid comes out significant can exceed \code{alpha}.
#'
#' The selected threshold is exploratory. \code{threshold_star} says up to
#' which strictness the release still shows a RAPID above chance. The print
#' calls it "Recommended". It is no advice to report RAPID at that threshold.
#' At \code{threshold_star} RAPID can rest on few records (7 of 200 in the
#' last example). For a report, fix the threshold before looking at the
#' data, and use \code{plot()} of the selection to show how RAPID depends on
#' the threshold.
#'
#' A test with \code{rapid_test()} at \code{threshold_star} on the same object
#' is no confirmation, because the threshold was chosen by looking at the same
#' data. To test a threshold, fix it before looking at the results, or select
#' it on one data set and test it on data not used for the selection, such as
#' a holdout or another release (see the last example).
#'
#' @inheritParams rapid_test
#' @param ... for the data methods, arguments of \code{\link{rapid}} (model
#'   and rule, the old argument names are accepted). Not allowed for a
#'   \code{rapid} object. Ignored by the methods for \code{rapid_threshold}
#'   objects.
#' @param tau_range numeric vector of thresholds for rules that put records
#'   with a score above the threshold at risk (default
#'   \code{seq(0.05, 0.95, by = 0.05)}). The default suits the two gain rules
#'   (band, \code{RCS_marginal}), where a threshold is at least 0 and below 1,
#'   and \code{NCE} with its default threshold 0.5. Under
#'   \code{RCS_conditional} the score is a ratio with the default threshold 1,
#'   and the default grid ends below it. Give a grid around 1 there, for
#'   example \code{seq(0.8, 2, by = 0.1)}. The values are not checked against
#'   the range of the rule.
#' @param epsilon_range numeric vector of tolerances for the point rule
#'   (default \code{seq(0.01, 0.5, by = 0.02)}, that is 1 to 49 percent). For
#'   the relative error metrics a tolerance is a proportion. A grid with a
#'   value of 1 or more is read as percent there, as \code{num_epsilon} in
#'   \code{rapid()}, with a warning of class
#'   \code{"rapid_warning_legacy_tolerance"}. For the absolute metric a
#'   tolerance is in the units of the target, and the default grid has no
#'   meaning. Always give \code{epsilon_range} there, for example
#'   \code{seq(1000, 10000, by = 1000)} for an income.
#' @param n_permutations number of permutations for the null distribution, a
#'   whole number of at least 1 (default 199). The permutations of the
#'   selection are \code{n_permutations}, those of the fit \code{n_perm}.
#' @param alpha level of the quantile criterion, above 0 and below 1 (default
#'   0.05). Each threshold is tested at this level on its own.
#'
#' @return An object of class \code{"rapid_threshold"} with components:
#' \describe{
#'   \item{threshold_star}{Selected threshold, the strictest significant one
#'     of the grid, or \code{NA} if no threshold is significant. The print
#'     calls it "Recommended". It is exploratory, see Details.}
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
#' \code{"summary.rapid_threshold"} with the number and the percentage of
#' significant thresholds (\code{n_significant}, \code{pct_significant}), the
#' range of the grid, the values at \code{threshold_star} and the table
#' \code{results}. \code{plot()} draws the observed
#' RAPID and the null quantile against the threshold, marks
#' \code{threshold_star} and returns the ggplot object, invisibly.
#' \code{as.data.frame()} returns \code{results}.
#'
#' @references
#' Templ, M., Thees, O. and Müller, R. (2026). RAPID: Risk of Attribute
#' Prediction-Induced Disclosure. Manuscript under review at EPJ Data Science.
#'
#' @seealso \code{\link{rapid}}, \code{\link{rapid_test}}.
#'   The section "Choosing the threshold" of
#'   \code{vignette("rapid", package = "riskutility")}.
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
#' # A selected threshold is confirmed on other records: select on one half
#' # of the original records, test on the other half
#' set.seed(42)
#' X2 <- data.frame(age = sample(20:60, 400, replace = TRUE),
#'                  hours = sample(10:50, 400, replace = TRUE))
#' X2$income <- 20000 + 600 * X2$age + 300 * X2$hours + rnorm(400, 0, 5000)
#' Y2 <- X2
#' Y2$income <- Y2$income + rnorm(400, 0, 3000)
#' half <- seq_len(nrow(X2)) <= 200
#' r_select <- rapid(X2[half, ], Y2, key_vars = c("age", "hours"),
#'                   target_var = "income", model_type = "lm", seed = 1)
#' sel2 <- rapid_threshold_select(r_select, n_permutations = 199, seed = 1)
#' sel2
#' sel2$observed_rapid$n_at_risk
#' r_confirm <- rapid(X2[!half, ], Y2, key_vars = c("age", "hours"),
#'                    target_var = "income", model_type = "lm",
#'                    num_tau = sel2$threshold_star, seed = 1)
#' rapid_test(r_confirm, n_permutations = 199, seed = 2)
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
  # A relative tolerance of 1 or more is read as percent, as num_epsilon in rapid()
  if (!above && !identical(rule$metric, "absolute") && any(thresholds >= 1)) {
    .rapid_warn(c("epsilon_range has values of 1 or more and is read as percent. A relative ",
                  "tolerance is now a proportion: divide epsilon_range by 100."),
                "legacy_tolerance")
    thresholds <- thresholds / 100
  }

  # One score matrix under permuted true values serves every threshold: the
  # score of a rule does not depend on its threshold.
  if (!is.null(seed)) set.seed(seed)
  null <- .rapid_null(X$attack, rule, n_permutations)
  per_threshold <- vapply(thresholds, function(th) {
    r <- .rule_set_threshold(rule, th)
    rates <- .rapid_chance_rates(r, null)
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
  # perm_level has no meaning under a fixed threshold
  call$perm_level <- NULL
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
