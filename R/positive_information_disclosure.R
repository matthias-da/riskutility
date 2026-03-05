#' Positive Information Disclosure (PID)
#'
#' Compute the Positive Information Disclosure (PID) metric based on prior and posterior probabilities.
#' This metric measures how much an adversary’s posterior probability improves compared to the prior,
#' indicating how much additional information has been gained from an observation.
#'
#' @param prior A named numeric vector of prior probabilities. Names represent possible secrets \eqn{x^*}.
#' @param posterior A named numeric vector of posterior probabilities (same names as \code{prior}).
#'
#' @return A single numeric value representing the maximum relative improvement in belief:
#' \deqn{\max_{x^*} \left( \frac{p(x^*|y) - p(x^*)}{p(x^*)} \right)}
#'
#' @details
#' Positive Information Disclosure (PID) builds on Shannon’s perfect secrecy criterion.
#' It quantifies how much an adversary's belief about any secret \eqn{x^*} improves after observing data \eqn{y}.
#' A PID of 0 means no information gain, satisfying perfect secrecy.
#' Larger values indicate greater improvement in belief, which may pose privacy risks.
#'
#' High PID values suggest privacy leakage—especially if sensitive values receive much higher posterior probabilities.
#'
#' @references
#' Wagner, I., & Eckhoff, D. (2018). Technical Privacy Metrics: A Systematic Survey.
#' ACM Computing Surveys (CSUR), 51(3), Article 57. Section 5.2.11.
#'
#' @examples
#' # Example: Prior vs. posterior beliefs over three secrets
#' prior <- c("A" = 0.3, "B" = 0.4, "C" = 0.3)
#' posterior <- c("A" = 0.6, "B" = 0.3, "C" = 0.1)
#' positive_information_disclosure(prior, posterior)
#'
#' # Privacy evaluation: Detecting excessive belief improvement for a sensitive attribute
#' # Suppose after anonymization, the adversary can better guess sensitive conditions
#' prior <- c("depression" = 0.05, "none" = 0.95)
#' posterior <- c("depression" = 0.15, "none" = 0.85)
#' pid <- positive_information_disclosure(prior, posterior)
#' if (pid > 1) message("Privacy risk detected: PID > 1")
#'
#' @family information-theory
#' @export
positive_information_disclosure <- function(prior, posterior) {
  if (!is.numeric(prior) || !is.numeric(posterior)) {
    stop("Both 'prior' and 'posterior' must be numeric vectors.")
  }
  if (any(prior <= 0)) {
    stop("All prior probabilities must be positive.")
  }
  if (!all(names(prior) %in% names(posterior))) {
    stop("Names of 'prior' and 'posterior' must match.")
  }

  secrets <- intersect(names(prior), names(posterior))
  ratios <- (posterior[secrets] - prior[secrets]) / prior[secrets]
  max(ratios)
}
