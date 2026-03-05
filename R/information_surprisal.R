#' Information Surprisal
#'
#' Computes the self-information (surprisal) of an outcome \code{x}, given its probability \code{p(x)}.
#' Surprisal quantifies how unexpected or informative an event is, and can be used to evaluate privacy risks.
#'
#' @details
#' In privacy evaluation, information surprisal measures how surprising an adversary would find
#' a specific combination of attributes. The lower the frequency of this combination in the population,
#' the higher the surprisal and the more identifiable the user becomes.
#'
#' Mathematically, information surprisal is defined as:
#' \deqn{I(x) = -\log_2(p(x))}
#'
#' A low probability \code{p(x)} results in a high surprisal value, indicating uniqueness or rarity,
#' which may correspond to a higher risk of re-identification.
#'
#' @param p A numeric vector of probabilities (or frequencies normalized to sum to 1).
#' @param logbase The logarithm base used for surprisal (default is 2).
#' @return A numeric vector of surprisal values.
#' @examples
#' # Example 1: Basic surprisal values
#' probs <- c(0.5, 0.25, 0.125, 0.125)
#' information_surprisal(probs)
#'
#' # Example 2: Using with a psychological dataset (e.g., personality profiles)
#' # Suppose these are relative frequencies of specific personality profiles
#' profile_freqs <- c("TypeA" = 0.4, "TypeB" = 0.35, "TypeC" = 0.2, "TypeD" = 0.05)
#' information_surprisal(profile_freqs)
#'
#' # Example 3: Privacy evaluation - rare attribute combination
#' # Suppose an individual has a rare attribute combination with p(x) = 0.0001
#' information_surprisal(0.0001)
#' # Output: 13.29 bits - high risk due to uniqueness
#'
#' @references
#' Wagner, I., & Eckhoff, D. (2018). Technical Privacy Metrics: A Systematic Survey.
#' ACM Computing Surveys, 51(3), Article 57. \doi{10.1145/3168389}
#'
#' @family information-theory
#' @export
information_surprisal <- function(p, logbase = 2) {
  if (any(p <= 0)) {
    stop("All probability values must be greater than 0.")
  }
  -log(p, base = logbase)
}
