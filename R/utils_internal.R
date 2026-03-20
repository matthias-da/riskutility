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
  # The algebraic expansion ||a-b||^2 = ||a||^2 + ||b||^2 - 2*a.b
  # can produce small negative or small positive residuals for identical
  # points due to floating-point cancellation.  Clamp near-zero values
  # so that identical rows yield distance exactly 0.
  tol <- max(max(abs(A2)), max(abs(B2))) * .Machine$double.eps * 4
  d2[d2 < tol] <- 0
  sqrt(d2)
}


#' Joint min-max normalization for multiple data frames
#'
#' Combines multiple numeric data frames, applies per-column min-max
#' normalization over the union, and returns each subset as a numeric matrix.
#' Used by euclidean-distance branches in nnaa(), hitting_rate(), etc.
#'
#' @param ... data frames to normalize jointly. All must share the same columns
#'   and contain only numeric data.
#' @return A list of numeric matrices in the same order as the inputs,
#'   each min-max normalized using the joint column ranges.
#' @keywords internal
.normalize_and_split <- function(...) {
  dfs <- list(...)
  ns <- vapply(dfs, nrow, integer(1))

  # Stack all data frames
  all_data <- do.call(rbind, dfs)

  # Check all numeric
  if (!all(vapply(all_data, is.numeric, logical(1)))) {
    stop("method='euclidean' requires all variables to be numeric. ",
         "Use method='gower' for mixed types.", call. = FALSE)
  }

  # Normalize each column using joint ranges
  all_norm <- as.data.frame(lapply(all_data, .normalize_minmax))

  # Split back into matrices
  cum_ns <- cumsum(c(0L, ns))
  lapply(seq_along(dfs), function(k) {
    as.matrix(all_norm[(cum_ns[k] + 1L):cum_ns[k + 1L], , drop = FALSE])
  })
}


#' Shared input validation for two-dataset (X vs Y) risk measures
#'
#' Validates data frames, intersects variables, checks types, removes NAs,
#' and optionally splits holdout.  Used by nnaa(), singling_out(),
#' linkability(), domias(), hitting_rate(), and epsilon_identifiability().
#'
#' @param X data.frame of original data
#' @param Y data.frame of synthetic data
#' @param holdout data.frame or NULL; if NULL and \code{holdout_fraction} is
#'   not NULL, X is split into train + holdout
#' @param holdout_fraction numeric in (0, 1) or NULL (no holdout splitting)
#' @param vars character vector of variable names or NULL (auto-intersect)
#' @param na.rm logical, remove incomplete cases
#' @param seed integer or NULL, random seed for holdout sampling
#' @param check_types logical, check that variable classes match in X and Y
#'   (default TRUE)
#' @param min_vars integer, minimum number of variables required (default 1)
#' @return list with X, Y, vars (character), and optionally train, holdout
#'   (if holdout splitting was performed)
#' @keywords internal
.validate_pair_inputs <- function(X, Y,
                                  holdout = NULL,
                                  holdout_fraction = NULL,
                                  vars = NULL,
                                  na.rm = TRUE,
                                  seed = NULL,
                                  check_types = TRUE,
                                  min_vars = 1L) {
  # --- Basic type checks ---
  if (!is.data.frame(X)) stop("X must be a data frame.")
  if (!is.data.frame(Y)) stop("Y must be a data frame.")
  if (!is.null(holdout) && !is.data.frame(holdout)) {
    stop("holdout must be a data frame or NULL.")
  }

  # --- Determine variables to use ---
  if (is.null(vars)) {
    vars <- intersect(names(X), names(Y))
    if (!is.null(holdout)) {
      vars <- intersect(vars, names(holdout))
    }
  }

  if (length(vars) < min_vars) {
    if (min_vars == 1L) {
      stop("No common variables found between datasets.")
    } else {
      stop(sprintf("At least %d variables are required.", min_vars))
    }
  }

  # --- Check variables exist ---
  missing_X <- setdiff(vars, names(X))
  missing_Y <- setdiff(vars, names(Y))
  if (length(missing_X) > 0) {
    stop(paste("Variables missing in X:", paste(missing_X, collapse = ", ")))
  }
  if (length(missing_Y) > 0) {
    stop(paste("Variables missing in Y:", paste(missing_Y, collapse = ", ")))
  }

  # --- Check variable types match ---
  if (check_types) {
    for (var in vars) {
      if (!identical(class(X[[var]]), class(Y[[var]]))) {
        stop(paste("Variable", var, "has different class in X and Y."))
      }
    }
  }

  # --- Subset to selected variables ---
  X <- X[, vars, drop = FALSE]
  Y <- Y[, vars, drop = FALSE]

  # --- Handle missing values ---
  if (na.rm) {
    complete_X <- stats::complete.cases(X)
    complete_Y <- stats::complete.cases(Y)
    X <- X[complete_X, , drop = FALSE]
    Y <- Y[complete_Y, , drop = FALSE]
  }

  if (nrow(X) == 0) stop("No complete cases in X after removing NAs.")
  if (nrow(Y) == 0) stop("No complete cases in Y after removing NAs.")

  # --- Holdout splitting ---
  result <- list(X = X, Y = Y, vars = vars)

  if (!is.null(holdout_fraction)) {
    if (is.null(holdout)) {
      if (!is.null(seed)) set.seed(seed)
      n_holdout <- max(1, floor(nrow(X) * holdout_fraction))
      holdout_idx <- sample(nrow(X), n_holdout)
      result$holdout <- X[holdout_idx, , drop = FALSE]
      result$train   <- X[-holdout_idx, , drop = FALSE]
    } else {
      holdout <- holdout[, vars, drop = FALSE]
      if (na.rm) {
        complete_H <- stats::complete.cases(holdout)
        holdout <- holdout[complete_H, , drop = FALSE]
      }
      if (nrow(holdout) == 0) stop("No complete cases in holdout after removing NAs.")
      result$holdout <- holdout
      result$train   <- X
    }
  }

  result
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
