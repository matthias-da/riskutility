#' @title Mutual Information Privacy Metrics
#'
#' @description
#' Computes various information-theoretic privacy metrics based on entropy, mutual information, and related measures. Includes:
#' - Mutual Information (MI)
#' - Normalized Mutual Information (NMI)
#' - Conditional Privacy Loss (CPL)
#' - Conditional Mutual Information (CMI)
#' - Relative Loss of Anonymity (RLA)
#'
#' @details
#' Mutual Information (MI) quantifies how much information is shared between two random variables X and Y. It is computed as:
#' \deqn{MI(X, Y) = H(X) - H(X | Y)}
#'
#' Conditional Mutual Information (CMI) measures the information shared between X and Y, given prior knowledge Z:
#' \deqn{CMI(X, Y | Z) = H(X | Z) - H(X | Y, Z)}
#'
#' Conditional Privacy Loss (CPL) is a normalized measure based on MI:
#' \deqn{CPL = 1 - 2^{-MI(X, Y)}}
#'
#' Relative Loss of Anonymity (RLA) can be approximated using conditional mutual information. It reflects the maximum leakage given prior knowledge.
#'
#' The normalized mutual information (`normalized = TRUE`) follows:
#' \deqn{NMI(X, Y) = \frac{MI(X, Y)}{H(X)}}
#' or can be normalized by the number of records (rows in `X`) for privacy interpretation.
#'
#' **Interpreting "greater" MI**: Higher MI means more information is leaked. But “greater” only implies higher *relative* risk. If the absolute MI is low (e.g., < 0.1 bits), the disclosure risk might still be negligible, even if it's greater than another case.
#'
#' @param X A matrix, data frame or probability vector representing the true data distribution.
#' @param Y A matrix, data frame or probability vector representing the observed or obfuscated data.
#' @param joint A joint probability matrix for X and Y.
#' @param Z Optional third variable (matrix or data.frame) representing prior knowledge for conditional metrics.
#' @param normalized Logical; whether to normalize mutual information (default: FALSE).
#' @return A list of privacy metrics: MI, NMI, CPL, CMI, and RLA (if Z is provided).
#' @examples
#' # Privacy example with synthetic probability distributions
#' X <- c(0.2, 0.5, 0.3)
#' Y <- c(0.25, 0.45, 0.3)
#' jointXY <- matrix(c(0.1, 0.05, 0.05,
#'                     0.1, 0.2,  0.15,
#'                     0.05, 0.2, 0.1), nrow = 3, byrow = TRUE)
#' mutualInformation(X, Y, jointXY, normalized = TRUE)
#'
#' @export
mutualInformation <- function(X, Y, joint, Z = NULL, normalized = FALSE) {
  # Helper: entropy
  H <- function(p) {
    p <- p[p > 0]
    -sum(p * log2(p))
  }

  # MI = H(X) + H(Y) - H(X,Y)
  Hx <- H(X)
  Hy <- H(Y)
  Hxy <- H(as.vector(joint))
  MI <- Hx + Hy - Hxy

  CPL <- 1 - 2^(-MI)

  # Normalize MI
  NMI <- if (normalized) MI / Hx else MI

  result <- list(
    MI = MI,
    NMI = NMI,
    CPL = CPL
  )

  # Conditional Mutual Information (CMI)
  if (!is.null(Z)) {
    # Estimate joint probabilities: p(x, y, z)
    jointXYZ <- table(X, Y, Z) / length(X)
    pXZ <- apply(jointXYZ, c(1,3), sum)
    pYZ <- apply(jointXYZ, c(2,3), sum)
    pZ  <- apply(jointXYZ, 3, sum)

    H_XZ <- H(as.vector(pXZ))
    H_XYZ <- H(as.vector(jointXYZ))

    CMI <- H_XZ - H_XYZ
    RLA <- CMI  # Interpreted as "max CMI over X" approximation
    result$CMI <- CMI
    result$RLA <- RLA
  }

  return(result)
}
