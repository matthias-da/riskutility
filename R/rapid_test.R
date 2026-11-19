#' Permutation Test for RAPID Disclosure Risk
#'
#' Tests whether RAPID is above its chance level, that is, whether the attacker's
#' predictions are aligned with the true values of the original records beyond
#' what random pairing gives. The attacker is fitted once. Each permutation
#' scores its fixed predictions against randomly reassigned true values under
#' the rule of the object (the baseline is that of the reassigned value), so
#' the test needs no refit.
#'
#' \code{rapid_test()} works on a fitted \code{\link{rapid}} object or on the
#' data. On the data it first calls \code{rapid()} with the arguments in
#' \code{...} and \code{n_perm = 0}, then tests that object. A given
#' \code{n_perm} or \code{threshold_rule} is passed on instead, so that
#' \code{threshold_rule = "permutation"} works on the data as well. The scoring rule is
#' part of the object: change it with \code{\link[=update.rapid]{update}()}.
#' \code{rapid_test(rapid(X, Y, ..., seed = s, n_perm = 0), seed = s)} equals
#' \code{rapid_test(X, Y, ..., seed = s)}.
#'
#' With \eqn{B} = \code{n_permutations} permutations, \eqn{R_b} the RAPID of
#' permutation \eqn{b} and \eqn{R} the observed RAPID, the p-value is
#' \deqn{p = \frac{1 + \#\{b : R_b \ge R\}}{1 + B},}{p = (1 + #{b : R_b >= R}) /
#'   (1 + B),}
#' and RAPID counts as significant when \eqn{p \le \alpha}{p <= alpha}. The
#' null hypothesis is that the predictions carry no information about which true
#' value belongs to which record. The test holds the fitted attacker fixed, so
#' it is exact given that attacker, and a significant result says that this
#' attacker infers the sensitive value of more records than chance allows.
#'
#' The test says nothing beyond that. It holds the attacker and the release
#' fixed, so the p-value does not cover the randomness of the model fit or of
#' the synthesis. A result that is not significant does not show that the
#' release is safe. It shows that under this attacker and this rule no more
#' records are at risk than chance gives. Another attacker can do better (see
#' \code{\link{rapid_suite}}). The smallest possible p-value is
#' \eqn{1/(1 + B)}, so a significant result needs at least
#' \code{1 / alpha - 1} permutations (19 for \code{alpha = 0.05}). The test
#' covers the records of \code{X}. A holdout stored in the object takes no
#' part in it.
#'
#' The mean of the \eqn{R_b} is the chance level \eqn{\pi_0}{pi0}, and the
#' Wald statistic \eqn{z = (R - \pi_0) / \sqrt{\pi_0 (1 - \pi_0) / n}}{z = (R -
#' pi0) / sqrt(pi0 (1 - pi0) / n)} with \eqn{n} scored records gives a
#' one-sided p-value from the normal distribution. The Wald test needs a stable
#' estimate of \eqn{\pi_0}{pi0}: use \code{n_permutations >= 999} (the print
#' says so otherwise). Its variance treats every record as at risk with the same
#' probability \eqn{\pi_0}{pi0}. When records differ in that probability, as
#' they do under a forest attacker, the permutation distribution is narrower,
#' and the Wald test is conservative.
#'
#' \code{plot()} draws the null distribution with the observed RAPID and the
#' \code{1 - alpha} quantile of the null. This quantile is the criterion of
#' \code{\link{rapid_threshold_select}()}, which calls a threshold significant
#' when RAPID is above it. The two criteria can disagree at the same threshold:
#' the quantile is interpolated between permuted values, while the p-value
#' counts the permuted values at or above RAPID. With 19 permutations and
#' \code{alpha = 0.05}, for example, a RAPID equal to the largest permuted value,
#' when no other permutation reaches it, is above the interpolated 95 percent
#' quantile but has \eqn{p = 2/20 = 0.1}.
#' The verdict of \code{rapid_test()} is the one of the p-value.
#'
#' @param X a \code{rapid} object (with stored predictions), the original data,
#'   or a \code{\link{synth_pair}}
#' @param Y data frame of the released (synthetic or anonymized) data the
#'   attacker is fitted on. Not needed if \code{X} is a \code{synth_pair}.
#' @param key_vars character vector with the names of the quasi-identifiers
#' @param target_var character, name of the sensitive variable (the target)
#' @param ... for the data methods, arguments of \code{\link{rapid}} (model
#'   and rule, the old argument names \code{original_data},
#'   \code{synthetic_data}, \code{quasi_identifiers} and
#'   \code{sensitive_attribute} are accepted). Not allowed for a \code{rapid}
#'   object. Ignored by the methods for \code{rapid_test} objects.
#' @param n_permutations number of permutations of the test, a whole number
#'   of at least 1 (default 199). The p-value cannot fall below
#'   \code{1 / (1 + n_permutations)}. The permutations of the test are
#'   \code{n_permutations}, those of the fit \code{n_perm}.
#' @param alpha significance level, above 0 and below 1 (default 0.05)
#' @param seed random seed of the permutations, and of the fit for the data
#'   methods (default NULL: the random stream is left alone). With a seed,
#'   \code{set.seed()} is called and the random stream of the session is not
#'   restored afterwards. In a simulation, set the seed once before the loop
#'   and leave \code{seed = NULL}, or pass a different seed in every
#'   replication.
#'
#' @return An object of class \code{"rapid_test"} with components:
#' \describe{
#'   \item{statistic, rapid}{Observed RAPID of the object.}
#'   \item{p_value}{Permutation p-value, \eqn{(1 + \#\{null \ge obs\}) / (1 + B)}{(1 +
#'     #{null >= obs}) / (1 + B)}.}
#'   \item{null_distribution}{RAPID of each permutation.}
#'   \item{n_permutations, alpha}{Settings.}
#'   \item{significant}{Logical, \code{TRUE} if \code{p_value <= alpha}.}
#'   \item{null_mean, chance}{Mean of the null distribution (the chance level).
#'     It rests on \code{n_permutations} permutations, the chance level of the
#'     object on its \code{n_perm} (10 by default). In a report, give the
#'     chance level of the test together with its p-value.}
#'   \item{excess}{\code{rapid - chance}.}
#'   \item{wald}{List \code{(z, p_value, pi0)}: Wald test of \code{rapid} against the chance
#'     level \code{pi0 = null_mean}, one-sided \code{p_value = 1 - pnorm(z)}. \code{z} and
#'     \code{p_value} are \code{NA} when \code{pi0} is 0 or 1.}
#'   \item{null_sd}{Standard deviation of the null distribution (\code{NA}
#'     with one permutation).}
#'   \item{null_quantiles}{Null quantiles (5, 50, 95 percent).}
#'   \item{observed_rapid}{The \code{rapid} object that was tested.}
#'   \item{method}{\code{"permutation"}.}
#' }
#' \code{format()} returns the printed lines and \code{print()} the object,
#' invisibly. \code{summary()} returns an object of class
#' \code{"summary.rapid_test"} with the components above except
#' \code{null_distribution} and \code{observed_rapid}, plus \code{rule} (the
#' printed rule) and \code{effect_size},
#' \code{(rapid - null_mean) / null_sd}. \code{plot()} returns the ggplot
#' object, invisibly. \code{as.data.frame()} returns the null distribution as a
#' data frame with the columns \code{permutation} and \code{rapid}.
#'
#' @references
#' Templ, M., Thees, O. and Müller, R. (2026). RAPID: Risk of Attribute
#' Prediction-Induced Disclosure. Manuscript under review at EPJ Data Science.
#'
#' @seealso \code{\link{rapid}}, \code{\link[=confint.rapid]{confint}},
#'   \code{\link{rapid_threshold_select}}, \code{\link[=update.rapid]{update}}.
#'   The section "Testing and intervals" of
#'   \code{vignette("rapid", package = "riskutility")}.
#'
#' @family rapid
#' @author Oscar Thees, Matthias Templ
#' @examples
#' set.seed(42)
#' X <- data.frame(
#'   age = sample(20:60, 100, replace = TRUE),
#'   hours = sample(10:50, 100, replace = TRUE)
#' )
#' X$income <- 20000 + 600 * X$age + 300 * X$hours + rnorm(100, 0, 5000)
#' Y <- X
#' Y$income <- Y$income + rnorm(100, 0, 3000)
#'
#' # On a fitted object: no refit
#' r <- rapid(X, Y, key_vars = c("age", "hours"), target_var = "income",
#'            model_type = "lm")
#' res <- rapid_test(r, n_permutations = 199, seed = 1)
#' res
#' res$p_value
#' res$wald
#' summary(res)
#' head(as.data.frame(res))
#' plot(res)
#'
#' # On the data: arguments in ... go to rapid()
#' rapid_test(X, Y, key_vars = c("age", "hours"), target_var = "income",
#'            model_type = "lm", n_permutations = 199, seed = 1)
#'
#' # A release without information on income. Under the forest some records
#' # are still at risk, and the test shows that they are no more than chance
#' # gives.
#' Y0 <- Y
#' Y0$income <- sample(Y0$income)
#' rapid_test(X, Y0, key_vars = c("age", "hours"), target_var = "income",
#'            model_type = "rf", n_permutations = 199, seed = 1)
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
  null <- .rapid_chance_rates(X$rule, .rapid_null(X$attack, X$rule, n_permutations))
  statistic <- X$rapid
  p_value <- (1 + sum(null >= statistic)) / (1 + n_permutations)
  chance <- mean(null)
  wald <- .rapid_wald(statistic, chance, X$n_total)

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
    wald = wald,
    null_sd = if (n_permutations > 1) stats::sd(null) else NA_real_,
    null_quantiles = stats::quantile(null, probs = c(0.05, 0.50, 0.95)),
    observed_rapid = X,
    method = "permutation"
  ), class = "rapid_test")
}

# Wald test of RAPID against the chance level pi0 (the null mean); one-sided. Not defined
# (NA) when pi0 is 0 or 1.
.rapid_wald <- function(statistic, pi0, n) {
  z <- if (is.finite(pi0) && pi0 > 0 && pi0 < 1) {
    (statistic - pi0) / sqrt(pi0 * (1 - pi0) / n)
  } else NA_real_
  list(z = z, p_value = if (is.na(z)) NA_real_ else 1 - stats::pnorm(z), pi0 = pi0)
}

# Printed lines of the Wald test (and the hint for few permutations)
.rapid_wald_lines <- function(x, digits = 3) {
  w <- x$wald
  if (is.null(w)) return(NULL)
  out <- if (is.na(w$z)) {
    "Wald test: not defined (the chance level is 0 or 1)"
  } else {
    paste0("Wald test against the chance level: z = ", formatC(w$z, digits = digits,
           format = "f"), ", one-sided p-value = ", formatC(w$p_value, digits = 4,
           format = "f"))
  }
  if (x$n_permutations < 999) {
    out <- c(out, "Hint: use n_permutations >= 999 for the Wald test")
  }
  out
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
           " the chance level"),
    .rapid_wald_lines(x, digits)
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
    wald = object$wald,
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
    else "Conclusion: RAPID is not significantly above the chance level.",
    .rapid_wald_lines(x, 4)
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
#' caller's expressions, so that the object stores a call that can be evaluated
#' again. `seed` stays in the call and seeds the fit. The fit skips its own
#' permutations (`n_perm = 0`) unless the call gives `n_perm` or
#' `threshold_rule`.
#' @param call `match.call()` of the data method.
#' @param env The caller's frame.
#' @noRd
.rapid_fit_for <- function(call, env) {
  call <- call[!names(call) %in% c("n_permutations", "alpha", "tau_range", "epsilon_range")]
  call[[1L]] <- rapid
  if (!any(c("n_perm", "threshold_rule") %in% names(call))) call$n_perm <- 0
  eval(call, env)
}

#' @noRd
.rapid_check_test_args <- function(n_permutations, alpha) {
  if (!is.numeric(n_permutations) || length(n_permutations) != 1 ||
      !is.finite(n_permutations) || n_permutations < 1 ||
      n_permutations != round(n_permutations)) {
    .rapid_abort("'n_permutations' must be a single whole number of at least 1.", "input")
  }
  if (!is.numeric(alpha) || length(alpha) != 1 || is.na(alpha) || alpha <= 0 ||
      alpha >= 1) {
    .rapid_abort("'alpha' must be a single number between 0 and 1.", "input")
  }
  invisible(TRUE)
}
