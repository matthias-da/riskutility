#' Prepare Mahalanobis Distance Components
#'
#' Estimates a (robust) covariance matrix from the data's numeric/ordinal
#' variables and pre-computes its inverse. Also identifies the split between
#' numeric and nominal key variables for the combined distance.
#'
#' @param data data.frame. Training data (typically the original/search data).
#' @param key character. Key variable names.
#' @param type named character. Per-key type ("numeric", "ordinal", "nominal").
#' @param robust logical. If TRUE (default), use MCD via
#'   \code{robustbase::covMcd()} or \code{MASS::cov.rob(method = "mcd")}.
#'   Falls back to classical \code{cov()} if MCD fails.
#' @param alpha numeric. Significance level for chi-squared threshold
#'   (default 0.025, giving the 97.5 percent quantile).
#' @return list with cov_inv, center (diagnostic only; not used in distance
#'   computation), numeric_keys, nominal_keys, alpha (proportion numeric),
#'   robust, chi_sq_threshold, p.
#' @keywords internal
.mahal_prepare <- function(data, key, type, robust = TRUE, alpha = 0.025) {

  numeric_keys <- key[type[key] %in% c("numeric", "ordinal")]
  nominal_keys <- key[type[key] == "nominal"]

  if (length(numeric_keys) == 0)
    stop("method = 'mahalanobis' requires at least one numeric or ordinal ",
         "key variable.", call. = FALSE)

  p <- length(numeric_keys)

  # Build numeric matrix
  num_mat <- data[, numeric_keys, drop = FALSE]
  for (v in numeric_keys) {
    if (is.factor(num_mat[[v]])) num_mat[[v]] <- as.integer(num_mat[[v]])
    num_mat[[v]] <- as.numeric(num_mat[[v]])
  }
  num_mat <- as.matrix(num_mat)

  # Remove rows with NA for covariance estimation
  complete <- stats::complete.cases(num_mat)
  num_complete <- num_mat[complete, , drop = FALSE]
  n <- nrow(num_complete)

  if (n < p + 1) {
    warning("Too few complete observations (", n, ") for ", p,
            " variables. Falling back to classical covariance.",
            call. = FALSE)
    robust <- FALSE
  }

  # Estimate covariance
  cov_est <- NULL

  if (robust && p >= 2) {
    cov_est <- tryCatch({
      if (requireNamespace("robustbase", quietly = TRUE)) {
        mcd <- robustbase::covMcd(num_complete)
        list(cov = mcd$cov, center = mcd$center)
      } else {
        rob <- MASS::cov.rob(num_complete, method = "mcd")
        list(cov = rob$cov, center = rob$center)
      }
    }, error = function(e) {
      warning("MCD estimation failed: ", conditionMessage(e),
              ". Falling back to classical covariance.", call. = FALSE)
      NULL
    })
  }

  if (is.null(cov_est)) {
    if (p == 1) {
      v <- stats::var(num_complete[, 1], na.rm = TRUE)
      if (!is.finite(v) || v <= 0) v <- 1
      cov_est <- list(
        cov = matrix(v, 1, 1),
        center = mean(num_complete[, 1], na.rm = TRUE)
      )
    } else {
      cov_est <- list(
        cov = stats::cov(num_complete),
        center = colMeans(num_complete)
      )
    }
    if (robust) robust <- FALSE
  }

  rob_cov <- cov_est$cov
  center <- cov_est$center

  # Check for singularity and regularize (relative to scale)
  if (rcond(rob_cov) < .Machine$double.eps) {
    warning("Covariance matrix is (near-)singular. ",
            "Regularizing with ridge penalty.", call. = FALSE)
    ridge <- 1e-6 * max(abs(diag(rob_cov)), 1)
    rob_cov <- rob_cov + diag(ridge, p)
  }

  # Invert
  cov_inv <- tryCatch(
    solve(rob_cov),
    error = function(e) {
      warning("Could not invert covariance matrix. ",
              "Using regularized version.", call. = FALSE)
      ridge <- 1e-4 * max(abs(diag(rob_cov)), 1)
      solve(rob_cov + diag(ridge, p))
    }
  )

  # Alpha: proportion of numeric vars among all key vars
  alpha_mix <- length(numeric_keys) / length(key)


  # Chi-squared threshold for pairwise distances:
  # For two independent draws x_i, x_j ~ N(mu, Sigma), the difference
  # has covariance 2*Sigma, so d^2(x_i, x_j)/2 ~ chi^2_p.
  chi_sq_thr <- sqrt(2 * stats::qchisq(1 - alpha, df = p))

  list(
    cov_inv = cov_inv,
    center = center,
    numeric_keys = numeric_keys,
    nominal_keys = nominal_keys,
    alpha = alpha_mix,
    robust = robust,
    chi_sq_threshold = chi_sq_thr,
    p = p
  )
}


#' Compute Mahalanobis-Based Distance from One Record to Candidates
#'
#' For numeric/ordinal variables, computes Mahalanobis distance using the
#' pre-computed inverse covariance from \code{.mahal_prepare()}, normalized
#' by the chi-squared threshold to produce values in approximately \eqn{[0, 1]}.
#' For nominal variables, computes Gower-style exact match (0/1).
#' Combined as weighted average based on proportion of each type.
#'
#' @param x_row data.frame with 1 row. Query record.
#' @param candidates data.frame. Search records (multiple rows).
#' @param prep list. Output from \code{.mahal_prepare()}.
#' @param type named character. Per-key type.
#' @return numeric vector of distances (one per candidate row).
#' @keywords internal
.mahal_dist <- function(x_row, candidates, prep, type) {
  nc <- nrow(candidates)
  if (nc == 0L) return(numeric(0L))

  # --- Numeric part: Mahalanobis distance ---
  d_maha <- rep(0, nc)
  if (length(prep$numeric_keys) > 0) {
    x_num <- numeric(prep$p)
    for (j in seq_along(prep$numeric_keys)) {
      v <- prep$numeric_keys[j]
      val <- x_row[[v]]
      if (is.factor(val)) val <- as.integer(val)
      x_num[j] <- as.numeric(val)
    }

    cand_mat <- candidates[, prep$numeric_keys, drop = FALSE]
    for (v in prep$numeric_keys) {
      if (is.factor(cand_mat[[v]])) cand_mat[[v]] <- as.integer(cand_mat[[v]])
      cand_mat[[v]] <- as.numeric(cand_mat[[v]])
    }
    cand_mat <- as.matrix(cand_mat)

    # diff_mat: nc x p
    diff_mat <- sweep(cand_mat, 2, x_num)
    # Mahalanobis: sqrt(diff %*% Sigma^{-1} %*% diff')
    rmd_sq <- rowSums((diff_mat %*% prep$cov_inv) * diff_mat)
    rmd_sq <- pmax(rmd_sq, 0)
    d_maha <- sqrt(rmd_sq)

    # Normalize by chi-squared threshold
    d_maha <- d_maha / prep$chi_sq_threshold
  }

  # --- Nominal part: exact match (Gower) ---
  d_nom <- rep(0, nc)
  if (length(prep$nominal_keys) > 0) {
    n_nom <- length(prep$nominal_keys)
    for (v in prep$nominal_keys) {
      d_nom <- d_nom + as.numeric(x_row[[v]] != candidates[[v]])
    }
    d_nom <- d_nom / n_nom
  }

  # --- Combine ---
  if (length(prep$nominal_keys) == 0) {
    d_maha
  } else {
    # Clamp numeric part to [0,1] before combining with bounded nominal part
    d_maha_clamped <- pmin(d_maha, 1)
    prep$alpha * d_maha_clamped + (1 - prep$alpha) * d_nom
  }
}
