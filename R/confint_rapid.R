#' @rdname rapid-methods
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
