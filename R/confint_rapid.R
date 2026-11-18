#' Confidence Intervals for RAPID Risk Estimates
#'
#' Computes confidence intervals for the RAPID disclosure risk measure
#' (proportion of records at risk). Four methods are available: Wald
#' (normal approximation), Wilson (score interval, recommended),
#' Clopper-Pearson (exact) and bootstrap (percentile).
#'
#' @param object An object of class \code{"rapid"}.
#' @param parm Unused; included for compatibility with the generic.
#' @param level Confidence level (default 0.95).
#' @param method One of \code{"wilson"} (default), \code{"wald"},
#'   \code{"clopper-pearson"}, or \code{"bootstrap"}.
#' @param n_bootstrap Number of bootstrap replicates (default 1000).
#'   Only used when \code{method = "bootstrap"}.
#' @param seed Optional random seed for the bootstrap.
#' @param ... Additional arguments (ignored).
#'
#' @details
#' The RAPID score is a sample proportion (number at risk / total).
#' Because the model is trained on separate synthetic data and applied to
#' the original records, the per-record risk indicators are conditionally
#' i.i.d. Bernoulli, making binomial confidence intervals valid.
#'
#' \describe{
#'   \item{wald}{Standard normal approximation
#'     \eqn{\hat p \pm z_{\alpha/2}\sqrt{\hat p(1-\hat p)/n}}{p +/- z * sqrt(p(1-p)/n)}.
#'     Can produce intervals outside \eqn{[0,1]}.}
#'   \item{wilson}{Score interval (Wilson, 1927). Has better coverage
#'     near 0 and 1; recommended as the default.}
#'   \item{clopper-pearson}{Exact interval from beta quantiles
#'     (Clopper and Pearson, 1934), equal to \code{binom.test()}; bounds 0 at
#'     \code{n_at_risk = 0} and 1 at \code{n_at_risk = n_total}. Conservative.}
#'   \item{bootstrap}{Resamples the \code{at_risk} indicators from
#'     \code{object$records} and returns percentile confidence limits.
#'     No model refitting is needed.}
#' }
#'
#' @return A \eqn{1 \times 2}{1x2} matrix with columns for the lower
#'   and upper confidence limits, following the standard \code{confint}
#'   convention.
#'
#' @references
#' Wilson, E.B. (1927). Probable Inference, the Law of Succession, and
#' Statistical Inference. \emph{Journal of the American Statistical Association},
#' 22(158), 209--212.
#'
#' Clopper, C.J. and Pearson, E.S. (1934). The use of confidence or fiducial limits
#' illustrated in the case of the binomial. \emph{Biometrika}, 26(4), 404--413.
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
#' r <- rapid(X, Y, key_vars = c("age", "sex"),
#'            target_var = "income", model_type = "lm")
#' confint(r)
#' confint(r, method = "wald")
#' confint(r, method = "clopper-pearson")
#'
#' \donttest{
#' # With random forest model and bootstrap CI
#' r2 <- rapid(X, Y, key_vars = c("age", "sex"),
#'             target_var = "income", model_type = "rf")
#' confint(r2)
#' confint(r2, method = "bootstrap", n_bootstrap = 500)
#' }
#'
#' @export
confint.rapid <- function(object, parm = NULL, level = 0.95,
                          method = c("wilson", "wald", "clopper-pearson", "bootstrap"),
                          n_bootstrap = 1000, seed = NULL, ...) {
  .rapid_require(object, "version")
  method <- .rapid_choice(method, c("wilson", "wald", "clopper-pearson", "bootstrap"),
                           "method")

  p <- object$rapid
  n <- object$n_total
  alpha <- 1 - level
  z <- stats::qnorm(1 - alpha / 2)

  ci <- switch(method,
    wald = {
      se <- sqrt(p * (1 - p) / n)
      c(p - z * se, p + z * se)
    },

    wilson = {
      denom <- 1 + z^2 / n
      centre <- (p + z^2 / (2 * n)) / denom
      margin <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / denom
      c(centre - margin, centre + margin)
    },

    "clopper-pearson" = {
      k <- object$n_at_risk
      c(if (k == 0) 0 else stats::qbeta(alpha / 2, k, n - k + 1),
        if (k == n) 1 else stats::qbeta(1 - alpha / 2, k + 1, n - k))
    },

    bootstrap = {
      if (!is.null(seed)) set.seed(seed)
      # Resampling the n flags with replacement gives a binomial count
      boot_props <- stats::rbinom(n_bootstrap, n, p) / n
      stats::quantile(boot_props, probs = c(alpha / 2, 1 - alpha / 2),
                      names = FALSE)
    }
  )

  ci_mat <- matrix(ci, nrow = 1,
                   dimnames = list("RAPID",
                                   paste0(100 * c(alpha / 2, 1 - alpha / 2), "%")))
  ci_mat
}
