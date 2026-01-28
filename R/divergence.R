#' Entropy measures
#'
#' Kullback-Leibler and Jensen-Shannon divergence
#'
#' @details
#' The Kullback-Leibler divergence is defined as:
#' \deqn{KLD(X,Y) = \sum X \times \log\left(\frac{X}{Y}\right)}
#'
#' If the data are multivariate, a Kullback-Leibler inspired divergence
#' is calculated through
#'
#' \deqn{KLD(X,Y) = \sum_{i,j}{X_{i,j} \left( \log\left(\frac{X_{i,j}}{Y_{i,j}}\right) +
#' \log(Y) - \log(X) \right)}}
#' with \deqn{X = \sum{X_{i,j}}} and \deqn{Y = \sum{Y_{i,j}}}
#'
#' The compositional versions of the Kullback-Leibler divergence is given
#' in the reference below.
#'
#' @param A a vector of probability densities
#' @param B a vector of probability densities
#' @return The divergence measure
#' @references J.A. Martin-Fernandez, M. Bren, C. Barcelo-Vidal,
#' V. Pawlowski-Glahn (1999). A measure of difference for compositional
#' data based on measures of divergence. In S. Lippard, A. Nass, and
#' R. Sinding-Larsen (Eds.), Proceedings of IAMG'99, Volume 1,
#' Trondheim (Norway), pp. 211-215.
#' @aliases KLDiv KLDiv_bayes JSDiv JSDiv_bayes
#' @rdname Entropy
#' @name Entropy
#' @export KLDiv
#' @export KLDiv_bayes
#' @export JSDiv
#' @export JSDiv_bayes
#'
NULL

#' @rdname Entropy
KLDiv <- function(A, B) {
  sum(A * log(A / B))
}

#' @rdname Entropy
KLDiv_bayes <- function(A, B) {
  length(A) / 2 * log(mean(A / B, na.rm = TRUE) * mean(B / A, na.rm = TRUE))
}

#' @rdname Entropy
JSDiv <- function(A, B) {
  M <- (A + B) / 2
  jsd <- 0.5 * KLDiv(A, M) + 0.5 * KLDiv(B, M)
  return(jsd)
}

#' @rdname Entropy
JSDiv_bayes <- function(A, B) {
  M <- (A + B) / 2
  jsd <- 0.5 * KLDiv_bayes(A, M) + 0.5 * KLDiv_bayes(B, M)
  return(jsd)
}

#' @rdname Entropy
#' @export
CrossEntropy <- function(A, B) {
  idx <- A > 0
  A <- A[idx]
  B <- B[idx]
  -sum(A * log(B))
}
