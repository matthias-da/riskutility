#' Privacy Score
#'
#' Computes the privacy score of a user based on the sensitivity of information items and their visibility.
#'
#' @description
#' The privacy score quantifies a user's potential privacy risk by considering the sensitivity of disclosed information and its visibility.
#' It is defined as:
#'
#' \deqn{privPS(u) = \sum_{x^* \in X^*} \omega_{x^*} \cdot Vis(x^*, u)}
#'
#' where \eqn{\omega_{x^*}} is the sensitivity of item \eqn{x^*} (e.g., how private or revealing it is), and \eqn{Vis(x^*, u)} is the visibility
#' of that item for user \eqn{u} (e.g., the number of people who have access to that item).
#'
#' @param sensitivity Named numeric vector of sensitivities for each information item (e.g., c(gender = 0.3, location = 0.8)).
#' @param visibility Named numeric vector of visibility values (e.g., number of users who can see each item).
#' @param normalize Logical. If TRUE, returns the normalized privacy score (between 0 and 1).
#' @return A single numeric value representing the privacy score for the user.
#'
#' @details
#' The privacy score helps quantify how exposed a user's data is, depending on both what is shared and how widely.
#'
#' If normalization is enabled, the score is divided by the maximum possible score (sum of sensitivities times max visibility).
#'
#' **Privacy interpretation**: A higher privacy score indicates higher exposure and, therefore, higher risk. The score increases with:
#' - More sensitive items being disclosed
#' - Items being visible to more people
#'
#' @references
#' Wagner, I., & Eckhoff, D. (2018). Technical Privacy Metrics: A Systematic Survey. *ACM Computing Surveys*, 51(3), 1–38. https://doi.org/10.1145/3168389
#'
#' @examples
#' # Example 1: Basic usage
#' sensitivity <- c(gender = 0.2, birthday = 0.8, location = 0.6)
#' visibility  <- c(gender = 100, birthday = 20, location = 50)
#' privacy_score(sensitivity, visibility)
#'
#' # Example 2: Normalized privacy score
#' privacy_score(sensitivity, visibility, normalize = TRUE)
#'
#' # Example 3: Privacy evaluation example
#' # Assume visibility corresponds to number of users
#' # who can see each attribute (e.g., social network)
#' psych_profile <- c(stress_score = 0.9, diagnosis = 1.0)
#' visibility <- c(stress_score = 5, diagnosis = 2)
#' privacy_score(psych_profile, visibility)
#'
#' @family information-theory
#' @export
privacy_score <- function(sensitivity, visibility, normalize = FALSE) {
  if (!all(names(sensitivity) %in% names(visibility))) {
    stop("All sensitivity items must also be present in the visibility vector (with matching names).")
  }

  # Align both vectors
  common_items <- intersect(names(sensitivity), names(visibility))
  s <- sensitivity[common_items]
  v <- visibility[common_items]

  score <- sum(s * v)

  if (normalize) {
    max_vis <- max(v)
    max_score <- sum(s * max_vis)
    return(score / max_score)
  }

  return(score)
}
