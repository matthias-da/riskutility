#' @title Maximum Information Leakage
#'
#' @description
#' Computes the maximum information leakage from a private variable \code{X} to an observed output \code{Y}.
#' It quantifies the worst-case information gain an adversary can achieve by observing any single output value.
#'
#' @details
#' Maximum Information Leakage is defined as:
#' \deqn{privMIL = \max_{y \in Y} I(X; Y = y)}
#' where \eqn{I(X; Y = y)} is the mutual information conditioned on a single observation.
#'
#' This is a worst-case metric, assuming that the adversary selects the output value \code{y} that maximizes information gain.
#'
#' This metric is relevant in scenarios where one specific output could lead to a disproportionate privacy breach.
#' See Section 5.2.7 in Wagner and Eckhoff (2018) for a detailed discussion.
#'
#' @param joint_dist A matrix representing the joint distribution of X and Y, where rows correspond to values of X and columns to values of Y. The matrix should sum to 1.
#'
#' @return A numeric value representing the maximum information leakage (in bits).
#'
#' @references
#' Wagner, I., & Eckhoff, D. (2018). Technical Privacy Metrics: A Systematic Survey. ACM Computing Surveys (CSUR), 51(3), Article 57.
#'
#' @examples
#' # Basic example: Define joint distribution of X and Y
#' joint <- matrix(c(0.1, 0.05, 0.05,
#'                   0.2, 0.1, 0.1,
#'                   0.1, 0.15, 0.15), nrow = 3, byrow = TRUE)
#' # Normalize to sum to 1
#' joint <- joint / sum(joint)
#' max_info_leakage(joint)
#'
#' # Example use in privacy evaluation:
#' # Suppose X represents income brackets (low, medium, high)
#' # and Y is an anonymized salary category after applying noise.
#' # We want to know the worst-case information leakage from salary output Y to true income X.
#' income <- factor(c("low", "low", "medium", "medium", "high", "high"))
#' salary <- factor(c("cat1", "cat2", "cat2", "cat2", "cat3", "cat1"))
#' joint_table <- table(income, salary)
#' joint_prob <- prop.table(joint_table)
#' max_info_leakage(joint_prob)
#'
#' @family information-theory
#' @author Matthias Templ
#' @export
max_info_leakage <- function(joint_dist) {
  if (!is.matrix(joint_dist)) stop("joint_dist must be a matrix")
  if (any(joint_dist < 0)) stop("joint_dist must not contain negative values")
  if (abs(sum(joint_dist) - 1) > 1e-6) stop("joint_dist must sum to 1")

  px <- rowSums(joint_dist)
  py <- colSums(joint_dist)

  info_leak_per_y <- sapply(seq_along(py), function(j) {
    if (py[j] == 0) return(0)
    pxy <- joint_dist[, j]
    px_given_y <- pxy / py[j]
    entropy_px_given_y <- -sum(px_given_y[pxy > 0] * log2(px_given_y[pxy > 0]))
    entropy_px <- -sum(px[px > 0] * log2(px[px > 0]))
    return(entropy_px - entropy_px_given_y)
  })

  return(max(info_leak_per_y))
}
