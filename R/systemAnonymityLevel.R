#' System Anonymity Level
#'
#' Computes the system anonymity level based on the entropy of sender-receiver combinations
#' in anonymous communication networks, normalized by the number of users.
#'
#' @param A A square binary matrix representing sender-receiver adjacency, where
#'        \code{A[i, j] = 1} indicates that sender \code{i} can send to receiver \code{j}.
#' @return A numeric value between 0 and 1 indicating the system anonymity level.
#'         A value of 1 indicates maximum anonymity.
#'
#' @details
#' The system anonymity level quantifies the uncertainty an adversary faces in
#' identifying the correct sender-receiver mapping. It is calculated as:
#'
#' \deqn{privSAL = \frac{1}{\log_2(|U|!)} H(p(x))},
#' where \eqn{H(p(x))} is the Shannon entropy of the adversary's estimated
#' distribution over sender-receiver combinations and \eqn{|U|} is the number of users.
#'
#' The probability \eqn{p(x)} is computed as \eqn{|E| / per(A)}, where
#' \eqn{per(A)} is the matrix permanent of the adjacency matrix \eqn{A},
#' and \eqn{|E|} is the size of the equivalence class.
#'
#' @note
#' Computing the permanent of a matrix is computationally intensive (NP-hard).
#' This implementation is feasible only for small matrices (typically \eqn{<= 6x6}).
#'
#' @references
#' I. Wagner and D. Eckhoff (2018). Technical Privacy Metrics: A Systematic Survey.
#' \emph{ACM Computing Surveys}, 51(3), Article 57. Section 5.2.8.
#'
#' @examples
#' # Example: Anonymity in a simple 3x3 anonymous messaging network
#' A <- matrix(c(
#'   1, 1, 0,
#'   1, 0, 1,
#'   0, 1, 1
#' ), nrow = 3, byrow = TRUE)
#' systemAnonymityLevel(A)
#'
#' # Example for privacy evaluation: Higher values indicate less certainty
#' # about communication patterns in the system, and thus greater privacy.
#' A2 <- matrix(1, nrow = 3, ncol = 3) # Fully connected
#' systemAnonymityLevel(A2) # Should return 1 (maximum uncertainty)
#'
#' @export
systemAnonymityLevel <- function(A) {
  if (!is.matrix(A) || !is.numeric(A)) stop("A must be a numeric matrix.")
  if (nrow(A) != ncol(A)) stop("A must be a square matrix.")
  if (any(A != 0 & A != 1)) stop("A must contain only 0s and 1s.")
  U <- nrow(A)
  if (U == 1) return(0)

  # Approximate permanent using Ryser’s formula
  permanent <- function(mat) {
    n <- nrow(mat)
    total <- 0
    for (s in 0:(2^n - 1)) {
      b <- as.integer(intToBits(s))[1:n]
      row_sum <- rep(0, n)
      for (j in 1:n) {
        if (b[j] == 1) {
          row_sum <- row_sum + mat[, j]
        }
      }
      prod_sum <- prod(row_sum)
      total <- total + (-1)^(sum(b)) * prod_sum
    }
    return((-1)^n * total)
  }

  permA <- permanent(A)
  if (permA <= 0) return(0)  # No valid permutations

  p <- rep(1 / permA, permA)
  entropy <- -sum(p * log2(p))

  return(entropy / log2(factorial(U)))
}
