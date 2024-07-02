#' Entropy measures
#'
#' Kullback-Leibler and Jensen-Shannon divergence
#'
#' @param A a vector of probability densities
#' @param B a vector of probability densities
#' @return The divergence measure
#' @aliases KLdiv KLDiv_bayes JSDiv JSDiv_bayes
#' @name divergence_measures
NULL

#' @rdname divergence_measures
#' @export
KLDiv <- function(A, B) {
  sum(A * log(A / B))
}

#' @rdname divergence_measures
#' @export
KLDiv_bayes <- function(A, B) {
  length(A) / 2 * log(mean(A / B, na.rm = TRUE) * mean(B / A, na.rm = TRUE))
}

#' @rdname divergence_measures
#' @export
JSDiv <- function(A, B) {
  M <- (A + B) / 2
  jsd <- 0.5 * KLDiv(A, M) + 0.5 * KLDiv(B, M)
  return(jsd)
}

#' @rdname divergence_measures
#' @export
JSDiv_bayes <- function(A, B) {
  M <- (A + B) / 2
  jsd <- 0.5 * KLDiv_bayes(A, M) + 0.5 * KLDiv_bayes(B, M)
  return(jsd)
}
