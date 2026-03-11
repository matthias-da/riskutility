# =============================================================================
# Shared internal helper functions
# =============================================================================

#' Wilson score confidence interval
#'
#' @param n_success integer number of successes
#' @param n_total integer number of trials
#' @param confidence_level numeric confidence level
#' @return list with estimate, ci_lower, ci_upper
#' @keywords internal
.wilson_score <- function(n_success, n_total, confidence_level = 0.95) {
  if (n_total == 0) {
    return(list(estimate = 0, ci_lower = 0, ci_upper = 0))
  }

  p_hat <- n_success / n_total
  z <- qnorm(1 - (1 - confidence_level) / 2)
  z2 <- z^2

  denom <- 1 + z2 / n_total
  center <- (p_hat + z2 / (2 * n_total)) / denom
  margin <- z * sqrt((p_hat * (1 - p_hat) + z2 / (4 * n_total)) / n_total) / denom

  list(
    estimate = p_hat,
    ci_lower = max(0, center - margin),
    ci_upper = min(1, center + margin)
  )
}


#' Compute residual risk from attack and control Wilson scores
#'
#' Shared helper for singling_out and linkability.
#' Returns list with risk, risk_ci, risk_attack, risk_attack_ci,
#' risk_control, risk_control_ci, privacy_pass.
#'
#' @param n_success integer, attack successes
#' @param n_control_success integer, control successes
#' @param n_attacks integer, total attacks
#' @param confidence_level numeric, confidence level
#' @param privacy_threshold numeric, threshold for privacy_pass
#' @return list with risk components
#' @keywords internal
.residual_risk <- function(n_success, n_control_success, n_attacks,
                           confidence_level = 0.95, privacy_threshold = 0.1) {
  attack_wilson <- .wilson_score(n_success, n_attacks, confidence_level)
  control_wilson <- .wilson_score(n_control_success, n_attacks, confidence_level)

  risk_attack <- attack_wilson$estimate
  risk_control <- control_wilson$estimate

  # Residual risk: (attack - control) / (1 - control), bounded [0, 1]
  if (risk_control < 1) {
    risk <- max(0, (risk_attack - risk_control) / (1 - risk_control))
  } else {
    risk <- 0
  }

  # Propagate CI for residual risk
  risk_ci_lower <- max(0, (attack_wilson$ci_lower - control_wilson$ci_upper) /
                          max(1e-10, 1 - control_wilson$ci_upper))
  risk_ci_upper <- min(1, (attack_wilson$ci_upper - control_wilson$ci_lower) /
                          max(1e-10, 1 - control_wilson$ci_lower))
  risk_ci <- c(lower = risk_ci_lower, upper = risk_ci_upper)

  list(
    risk = risk,
    risk_ci = risk_ci,
    risk_attack = risk_attack,
    risk_attack_ci = c(lower = attack_wilson$ci_lower,
                       upper = attack_wilson$ci_upper),
    risk_control = risk_control,
    risk_control_ci = c(lower = control_wilson$ci_lower,
                        upper = control_wilson$ci_upper),
    privacy_pass = risk <= privacy_threshold
  )
}


#' Min-max normalization
#'
#' @param x numeric vector
#' @return normalized vector in \[0, 1\]
#' @keywords internal
.normalize_minmax <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (rng[2] - rng[1] == 0) return(rep(0, length(x)))
  (x - rng[1]) / (rng[2] - rng[1])
}


#' Euclidean distance matrix
#'
#' Computes pairwise Euclidean distances between rows of A and B.
#'
#' @param A numeric matrix
#' @param B numeric matrix
#' @return nrow(A) x nrow(B) distance matrix
#' @keywords internal
.euclidean_dist <- function(A, B) {
  A2 <- rowSums(A^2)
  B2 <- rowSums(B^2)
  cross <- tcrossprod(A, B)
  d2 <- outer(A2, B2, "+") - 2 * cross
  d2[d2 < 0] <- 0
  sqrt(d2)
}


#' Shared holdout preparation for distance-risk functions
#'
#' Validates inputs, intersects variables, handles NAs, and splits holdout.
#' Used by dcr(), nndr(), and rf_privacy().
#'
#' @param X data.frame of original data
#' @param Y data.frame of synthetic data
#' @param holdout data.frame or NULL
#' @param holdout_fraction numeric in (0, 1)
#' @param vars character vector or NULL
#' @param na.rm logical
#' @param seed integer or NULL
#' @param min_holdout integer, minimum holdout size (nndr needs 2)
#' @return list with train, synthetic, holdout, vars, was_split
#' @keywords internal
.distance_risk_prepare <- function(X, Y, holdout = NULL,
                                    holdout_fraction = 0.5,
                                    vars = NULL, na.rm = TRUE,
                                    seed = NULL, min_holdout = 1L) {
  if (!is.data.frame(X)) stop("X must be a data frame", call. = FALSE)
  if (!is.data.frame(Y)) stop("Y must be a data frame", call. = FALSE)
  if (!is.null(holdout) && !is.data.frame(holdout)) {
    stop("holdout must be a data frame or NULL", call. = FALSE)
  }
  if (holdout_fraction <= 0 || holdout_fraction >= 1) {
    stop("holdout_fraction must be in (0, 1)", call. = FALSE)
  }

  # Determine common vars
  common <- intersect(names(X), names(Y))
  if (!is.null(holdout)) {
    common <- intersect(common, names(holdout))
  }
  if (!is.null(vars)) {
    common <- intersect(common, vars)
  }
  if (length(common) == 0) {
    stop("No common variables found between X, Y, and holdout", call. = FALSE)
  }

  X <- X[, common, drop = FALSE]
  Y <- Y[, common, drop = FALSE]

  # NA handling
  if (na.rm) {
    X <- X[stats::complete.cases(X), , drop = FALSE]
    Y <- Y[stats::complete.cases(Y), , drop = FALSE]
  }

  was_split <- FALSE
  if (is.null(holdout)) {
    # Split X into train + holdout
    if (!is.null(seed)) set.seed(seed)
    n <- nrow(X)
    n_holdout <- max(min_holdout, floor(n * holdout_fraction))
    n_holdout <- min(n_holdout, n - 1L)  # keep at least 1 train record
    holdout_idx <- sample(n, n_holdout)
    holdout_df <- X[holdout_idx, , drop = FALSE]
    train_df <- X[-holdout_idx, , drop = FALSE]
    was_split <- TRUE
  } else {
    holdout_df <- holdout[, common, drop = FALSE]
    if (na.rm) {
      holdout_df <- holdout_df[stats::complete.cases(holdout_df), ,
                               drop = FALSE]
    }
    train_df <- X
  }

  list(
    train = train_df,
    synthetic = Y,
    holdout = holdout_df,
    vars = common,
    was_split = was_split
  )
}
