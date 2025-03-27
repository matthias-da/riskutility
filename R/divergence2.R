#' Additional Entropy-Based Privacy Measures
#'
#' This set of functions implements various entropy-based measures used in privacy and information theory.
#' All functions assume that the input vectors \code{p} or matrices are valid probability distributions.
#'
#' \describe{
#'   \item{\code{RenyiEntropy(p, alpha = 2)}}{Computes Rényi entropy of order \code{alpha} for a probability vector \code{p}.}
#'   \item{\code{MaxEntropy(p)}}{Returns the Hartley (max) entropy, defined as the logarithm of the number of non-zero categories.}
#'   \item{\code{MinEntropy(p)}}{Computes the min-entropy, defined as the negative logarithm of the maximum probability.}
#'   \item{\code{NormalizedEntropy(p)}}{Returns the Shannon entropy normalized by the maximum entropy, resulting in a value between 0 and 1.}
#'   \item{\code{ConditionalEntropy(joint)}}{Computes the conditional entropy H(Y|X), given a joint probability matrix (rows = X, columns = Y).}
#'   \item{\code{CumulativeEntropy(x)}}{Estimates cumulative entropy based on the empirical cumulative distribution function of a numeric vector \code{x}.}
#' }
#'
#' @param p A vector of probability values (must sum to 1 for entropy-based methods).
#' @param alpha The order of Rényi entropy (must be > 0 and ≠ 1).
#' @param joint A joint probability matrix, where rows represent values of X and columns of Y.
#' @param x A numeric vector (for cumulative entropy).
#'
#' @return A numeric value representing the corresponding entropy or risk measure.
#'
#' @references
#' Rényi, A. (1961). On measures of entropy and information. In Proceedings of the 4th Berkeley Symposium on Mathematical Statistics and Probability (Vol. 1, pp. 547–561).
#'
#' Shannon, C.E. (1948). A Mathematical Theory of Communication. Bell System Technical Journal, 27(3), 379–423.
#'
#' Bender, S., Haas, A., & Klose, C. (2001). The IAB Employment Sample. Journal of Applied Social Science Studies, 121(2), 183–190.
#'
#' Ichim, D. (2009). Local neighbourhood-based record linkage and its application to the Canadian census. In Privacy in Statistical Databases (pp. 207–221). Springer.
#'
#' @keywords entropy privacy risk
#' @name EntropyMeasures
#' @aliases RenyiEntropy MaxEntropy MinEntropy NormalizedEntropy ConditionalEntropy CumulativeEntropy
#' @export
RenyiEntropy <- function(p, alpha = 2) {
  stopifnot(alpha > 0, alpha != 1)
  p <- p[p > 0]
  (1 / (1 - alpha)) * log(sum(p^alpha))
}

#' @export
MaxEntropy <- function(p) {
  log(length(p[p > 0]))
}

#' @export
MinEntropy <- function(p) {
  -log(max(p))
}

#' @export
NormalizedEntropy <- function(p) {
  p <- p[p > 0]
  H <- -sum(p * log(p))
  H_max <- log(length(p))
  H / H_max
}

#' @export
ConditionalEntropy <- function(joint) {
  p_x <- rowSums(joint)
  H_cond <- 0
  for (i in seq_len(nrow(joint))) {
    p_y_given_x <- joint[i, ] / p_x[i]
    p_y_given_x <- p_y_given_x[p_y_given_x > 0]
    H_cond <- H_cond + p_x[i] * sum(-p_y_given_x * log(p_y_given_x))
  }
  H_cond
}

#' @export
CumulativeEntropy <- function(x) {
  x <- sort(x)
  F <- ecdf(x)(x)
  -mean(log(F[F > 0]))
}
