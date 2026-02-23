#' Confidence Intervals for RAPID Risk Estimates
#'
#' Computes confidence intervals for the RAPID disclosure risk measure
#' (proportion of records at risk). Three methods are available: Wald
#' (normal approximation), Wilson (score interval, recommended), and
#' bootstrap (percentile).
#'
#' @param object An object of class \code{"rapid"}.
#' @param parm Unused; included for compatibility with the generic.
#' @param level Confidence level (default 0.95).
#' @param method One of \code{"wilson"} (default), \code{"wald"}, or
#'   \code{"bootstrap"}.
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
#'   \item{bootstrap}{Resamples the \code{at_risk} indicators from
#'     \code{object$records} and returns percentile confidence limits.
#'     No model refitting is needed.}
#' }
#'
#' @return A \eqn{1 \times 2}{1x2} matrix with columns for the lower
#'   and upper confidence limits, following the standard \code{confint}
#'   convention.
#'
#' @seealso \code{\link{rapid}}, \code{\link{rapid_test}}
#'
#' @examples
#' \dontrun{
#' r <- rapid(original, synthetic, key_vars = c("age", "sex"),
#'            target_var = "income", model_type = "rf")
#' confint(r)                          # Wilson (default)
#' confint(r, method = "wald")         # Wald
#' confint(r, method = "bootstrap")    # Bootstrap percentile
#' }
#'
#' @export
confint.rapid <- function(object, parm = NULL, level = 0.95,
                          method = c("wilson", "wald", "bootstrap"),
                          n_bootstrap = 1000, seed = NULL, ...) {
  method <- match.arg(method)

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

    bootstrap = {
      if (!is.null(seed)) set.seed(seed)
      at_risk <- object$records$at_risk
      boot_props <- vapply(seq_len(n_bootstrap), function(i) {
        mean(sample(at_risk, replace = TRUE))
      }, numeric(1))
      stats::quantile(boot_props, probs = c(alpha / 2, 1 - alpha / 2),
                      names = FALSE)
    }
  )

  ci_mat <- matrix(ci, nrow = 1,
                   dimnames = list("RAPID",
                                   paste0(100 * c(alpha / 2, 1 - alpha / 2), "%")))
  ci_mat
}
