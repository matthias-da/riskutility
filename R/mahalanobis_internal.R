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
#'   (default 0.025, giving a 97.5\% quantile).
#' @return list with cov_inv, center, numeric_keys, nominal_keys,
#'   alpha (proportion numeric), robust, chi_sq_threshold, p.
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

  # Check for singularity and regularize
  cov_det <- det(rob_cov)
  if (is.na(cov_det) || cov_det < .Machine$double.eps) {
    warning("Covariance matrix is (near-)singular. ",
            "Regularizing with ridge penalty.", call. = FALSE)
    rob_cov <- rob_cov + diag(1e-6, p)
  }

  # Invert
  cov_inv <- tryCatch(
    solve(rob_cov),
    error = function(e) {
      warning("Could not invert covariance matrix. ",
              "Using regularized version.", call. = FALSE)
      solve(rob_cov + diag(1e-4, p))
    }
  )

  # Alpha: proportion of numeric vars among all key vars
  alpha_mix <- length(numeric_keys) / length(key)

  # Chi-squared threshold for p numeric dimensions
  chi_sq_thr <- sqrt(stats::qchisq(1 - alpha, df = p))

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
