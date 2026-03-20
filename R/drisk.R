#' Disclosure Risk for Continuous Variables (dRisk / dRiskRMD)
#'
#' Computes disclosure risk measures for continuous key variables, based on
#' interval overlap (dRisk) and Robust Mahalanobis Distance (dRiskRMD). These
#' methods assess whether original records can be re-identified from synthetic
#' data using numeric quasi-identifiers.
#'
#' @param X data frame of original/training data
#' @param Y data frame of synthetic data
#' @param vars character vector of numeric variable names to use. If NULL
#'   (default), all common numeric variables between X and Y are used.
#' @param method character, which risk measure(s) to compute: "interval"
#'   (dRisk), "rmd" (dRiskRMD), or "both" (default).
#' @param outlier_par numeric, outlier parameter controlling interval width for
#'   the interval method (default: 0.01). Larger values produce wider intervals
#'   and higher risk estimates.
#' @param alpha numeric, significance level for the chi-squared threshold in the
#'   RMD method (default: 0.05). Smaller alpha means a stricter threshold
#'   (fewer flagged records).
#' @param na.rm logical, remove records with NA values (default: TRUE)
#' @param ... additional arguments passed to methods (currently unused)
#'
#' @return An object of class "drisk" containing:
#' \itemize{
#'   \item drisk_interval: fraction of at-risk records (interval method), NA if
#'     not computed
#'   \item drisk_rmd: fraction of at-risk records (RMD method), NA if not
#'     computed
#'   \item at_risk_interval: logical vector per original record (interval
#'     method), NULL if not computed
#'   \item at_risk_rmd: logical vector per original record (RMD method), NULL if
#'     not computed
#'   \item min_rmd: numeric vector of minimum Robust Mahalanobis Distance per
#'     original record, NULL if method does not include "rmd"
#'   \item interval_widths: named numeric vector of interval half-widths d_v per
#'     variable, NULL if method does not include "interval"
#'   \item method: method(s) used
#'   \item privacy_pass: logical, max(drisk_interval, drisk_rmd) <= 0.1
#'   \item n_original: number of original records used
#'   \item n_synthetic: number of synthetic records used
#'   \item vars: variables used
#'   \item outlier_par: outlier parameter used
#'   \item alpha: significance level used
#' }
#'
#' @details
#' Both methods work ONLY on numeric variables. Non-numeric columns are
#' automatically excluded (or cause an error if explicitly requested via
#' \code{vars}).
#'
#' \strong{Interval method (dRisk):}
#' For each variable v, an interval half-width is computed as:
#' \deqn{d_v = \mathrm{outlier\_par} \times \mathrm{IQR}(X_v) / (2n)^{1/3}}
#' where n = nrow(X). For each original record i, the record is flagged as
#' at-risk if there exists at least one synthetic record j where ALL variables
#' of record i fall within the intervals around record j:
#' \deqn{|x_{iv} - y_{jv}| \le d_v \quad \forall v}
#' The dRisk score is the fraction of original records that are at-risk.
#'
#' \strong{RMD method (dRiskRMD):}
#' Uses the Minimum Covariance Determinant (MCD) estimator from
#' \code{\link[MASS]{cov.rob}} to obtain a robust covariance matrix from the
#' original data. For each original record, the minimum Robust Mahalanobis
#' Distance to any synthetic record is computed. Records with minimum RMD below
#' the chi-squared threshold \eqn{\chi^2_{1-\alpha, p}} (where p = number of
#' variables) are flagged as at-risk.
#'
#' @section Interpretation:
#' \itemize{
#'   \item \strong{dRisk close to 0}: Low disclosure risk - few original records
#'     can be matched by synthetic data within the tolerance intervals.
#'   \item \strong{dRisk close to 1}: High disclosure risk - most original
#'     records have close matches in the synthetic data.
#'   \item \strong{dRiskRMD close to 0}: Low risk - original records are far
#'     from synthetic records in robust Mahalanobis distance.
#'   \item \strong{dRiskRMD close to 1}: High risk - many original records have
#'     suspiciously close synthetic counterparts.
#' }
#'
#' @seealso \code{\link{dcr}} for distance-based privacy with holdout comparison,
#'   \code{\link{nndr}} for nearest neighbor distance ratio,
#'   \code{\link{ims}} for identical match detection
#'
#' @references
#' Templ, M. (2017). \emph{Statistical Disclosure Control for Microdata:
#' Methods and Applications in R}. Springer.
#' \doi{10.1007/978-3-319-50272-4}
#'
#' Templ, M., Kowarik, A. & Meindl, B. (2024). sdcMicro: Statistical
#' Disclosure Control Methods for Anonymization of Data and Risk Estimation.
#' R package. \url{https://CRAN.R-project.org/package=sdcMicro}
#'
#' Rousseeuw, P. J. & Van Driessen, K. (1999). A fast algorithm for the
#' Minimum Covariance Determinant estimator.
#' \emph{Technometrics}, 41(3), 212--223.
#' \doi{10.1080/00401706.1999.10485670}
#'
#' @author Matthias Templ
#' @family privacy-models
#' @export
#' @importFrom MASS cov.rob
#' @importFrom stats complete.cases IQR mahalanobis qchisq
#' @importFrom graphics hist abline legend barplot par text stripchart
#' @examples
#' # Create example data
#' set.seed(123)
#' n <- 100
#' X <- data.frame(
#'   age = rnorm(n, 40, 10),
#'   income = rnorm(n, 50000, 15000),
#'   hours = rnorm(n, 40, 8)
#' )
#'
#' # Good synthetic data (independent generation)
#' Y <- data.frame(
#'   age = rnorm(n, 40, 10),
#'   income = rnorm(n, 50000, 15000),
#'   hours = rnorm(n, 40, 8)
#' )
#'
#' result <- drisk(X, Y)
#' print(result)
#' summary(result)
#'
#' \donttest{
#' # Memorized data (Y is copy of X) - should show high risk
#' Y_copy <- X + rnorm(n * 3, sd = 0.01)
#' result_bad <- drisk(X, Y_copy)
#' print(result_bad)
#' }
drisk <- function(X, ...) {
  UseMethod("drisk")
}

#' @rdname drisk
#' @export
drisk.synth_pair <- function(X, ...) {
  # Intentionally X$num_vars (not X$vars): dRisk interval widths and Robust
  # Mahalanobis Distances are defined only for numeric variables; categorical
  # columns would cause errors in IQR() and cov.rob().
  drisk.default(
    X = X$original,
    Y = X$synthetic,
    vars = X$num_vars,
    ...
  )
}

#' @rdname drisk
#' @export
drisk.default <- function(X, Y,
                          vars = NULL,
                          method = c("both", "interval", "rmd"),
                          outlier_par = 0.01,
                          alpha = 0.05,
                          na.rm = TRUE,
                          ...) {

  # Input validation
  if (!is.data.frame(X)) stop("X must be a data frame.")
  if (!is.data.frame(Y)) stop("Y must be a data frame.")

  method <- match.arg(method)

  if (!is.numeric(outlier_par) || length(outlier_par) != 1 || outlier_par <= 0) {
    stop("outlier_par must be a positive number.")
  }
  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 1) {
    stop("alpha must be a number between 0 and 1.")
  }

  # Determine variables to use
  if (is.null(vars)) {
    common <- intersect(names(X), names(Y))
    if (length(common) == 0) {
      stop("No common variables found between datasets.")
    }
    # Select only numeric variables
    vars <- common[vapply(X[, common, drop = FALSE], is.numeric, logical(1)) &
                   vapply(Y[, common, drop = FALSE], is.numeric, logical(1))]
    if (length(vars) == 0) {
      stop("No common numeric variables found between datasets.")
    }
  } else {
    # Check variables exist
    missing_X <- setdiff(vars, names(X))
    missing_Y <- setdiff(vars, names(Y))
    if (length(missing_X) > 0) {
      stop(paste("Variables missing in X:", paste(missing_X, collapse = ", ")))
    }
    if (length(missing_Y) > 0) {
      stop(paste("Variables missing in Y:", paste(missing_Y, collapse = ", ")))
    }
    # Verify all are numeric
    non_num_X <- vars[!vapply(X[, vars, drop = FALSE], is.numeric, logical(1))]
    non_num_Y <- vars[!vapply(Y[, vars, drop = FALSE], is.numeric, logical(1))]
    non_num <- unique(c(non_num_X, non_num_Y))
    if (length(non_num) > 0) {
      stop(paste("All variables must be numeric. Non-numeric:",
                 paste(non_num, collapse = ", ")))
    }
  }

  if (length(vars) < 1) {
    stop("At least one numeric variable is required.")
  }

  # Subset to selected variables
  X <- X[, vars, drop = FALSE]
  Y <- Y[, vars, drop = FALSE]

  # Handle missing values
  if (na.rm) {
    complete_X <- complete.cases(X)
    complete_Y <- complete.cases(Y)
    X <- X[complete_X, , drop = FALSE]
    Y <- Y[complete_Y, , drop = FALSE]
  }

  if (nrow(X) == 0) stop("No complete cases in X after removing NAs.")
  if (nrow(Y) == 0) stop("No complete cases in Y after removing NAs.")

  n_original <- nrow(X)
  n_synthetic <- nrow(Y)
  p <- length(vars)

  do_interval <- method %in% c("interval", "both")
  do_rmd <- method %in% c("rmd", "both")

  # --- Interval method (dRisk) ---
  drisk_interval <- NA_real_
  at_risk_interval <- NULL
  interval_widths <- NULL

  if (do_interval) {
    # Compute interval half-widths
    denom <- (2 * n_original)^(1 / 3)
    d_v <- vapply(vars, function(v) {
      outlier_par * IQR(X[[v]], na.rm = TRUE) / denom
    }, numeric(1))
    names(d_v) <- vars
    interval_widths <- d_v

    # Convert to matrices for vectorized computation
    X_mat <- as.matrix(X)
    Y_mat <- as.matrix(Y)

    # For each original record, check if any synthetic record matches within
    # intervals on ALL variables
    at_risk_interval <- logical(n_original)

    for (i in seq_len(n_original)) {
      # Compute absolute differences: n_synthetic x p matrix
      abs_diff <- abs(sweep(Y_mat, 2, X_mat[i, ]))
      # For each synthetic record j, check if all variables are within d_v
      within_interval <- sweep(abs_diff, 2, d_v, "<=")
      # Record i is at risk if any row (synthetic record) has all TRUE
      at_risk_interval[i] <- any(rowSums(within_interval) == p)
    }

    drisk_interval <- mean(at_risk_interval)
  }

  # --- RMD method (dRiskRMD) ---
  drisk_rmd <- NA_real_
  at_risk_rmd <- NULL
  min_rmd <- NULL

  if (do_rmd) {
    if (p < 2) {
      stop("RMD method requires at least 2 numeric variables.")
    }

    # Robust covariance estimation using MCD on original data
    rob <- MASS::cov.rob(X, method = "mcd")
    rob_center <- rob$center
    rob_cov <- rob$cov

    # Check for singular covariance
    cov_det <- det(rob_cov)
    if (is.na(cov_det) || cov_det < .Machine$double.eps) {
      warning("Robust covariance matrix is (near-)singular. ",
              "RMD results may be unreliable. Consider using fewer or ",
              "less collinear variables.")
      # Fall back to regularized covariance
      rob_cov <- rob_cov + diag(1e-6, p)
    }

    # Compute inverse covariance once
    cov_inv <- tryCatch(
      solve(rob_cov),
      error = function(e) {
        warning("Could not invert robust covariance matrix. ",
                "Using regularized version.")
        solve(rob_cov + diag(1e-6, p))
      }
    )

    # For each original record, find minimum RMD to any synthetic record
    X_mat <- as.matrix(X)
    Y_mat <- as.matrix(Y)
    min_rmd <- numeric(n_original)

    for (i in seq_len(n_original)) {
      # Compute RMD from original record i to all synthetic records
      diff_mat <- sweep(Y_mat, 2, X_mat[i, ])
      # Mahalanobis distance: sqrt(diff %*% cov_inv %*% t(diff))
      # For each row j of diff_mat: d_j = sqrt(diff_j %*% cov_inv %*% diff_j)
      rmd_sq <- rowSums((diff_mat %*% cov_inv) * diff_mat)
      # Ensure non-negative (numerical issues)
      rmd_sq <- pmax(rmd_sq, 0)
      min_rmd[i] <- min(sqrt(rmd_sq))
    }

    # Chi-squared threshold
    threshold <- sqrt(qchisq(1 - alpha, df = p))
    at_risk_rmd <- min_rmd < threshold
    drisk_rmd <- mean(at_risk_rmd)
  }

  # Privacy assessment
  risk_vals <- c(drisk_interval, drisk_rmd)
  risk_vals <- risk_vals[!is.na(risk_vals)]
  privacy_pass <- max(risk_vals) <= 0.1

  results <- list(
    drisk_interval = drisk_interval,
    drisk_rmd = drisk_rmd,
    at_risk_interval = at_risk_interval,
    at_risk_rmd = at_risk_rmd,
    min_rmd = min_rmd,
    interval_widths = interval_widths,
    method = method,
    privacy_pass = privacy_pass,
    n_original = n_original,
    n_synthetic = n_synthetic,
    vars = vars,
    outlier_par = outlier_par,
    alpha = alpha
  )

  class(results) <- "drisk"
  return(results)
}

# --- S3 methods ---

#' Print method for drisk objects
#'
#' @param x an object of class "drisk"
#' @param ... additional arguments (ignored)
#' @export
print.drisk <- function(x, ...) {
  cat("Disclosure Risk for Continuous Variables (dRisk)\n")
  cat("================================================\n")
  cat("Method:", x$method, "\n")
  cat("Variables used:", length(x$vars),
      paste0("(", paste(x$vars, collapse = ", "), ")"), "\n\n")

  cat("Dataset Sizes:\n")
  cat("  Original records:", x$n_original, "\n")
  cat("  Synthetic records:", x$n_synthetic, "\n\n")

  if (!is.na(x$drisk_interval)) {
    n_at_risk <- sum(x$at_risk_interval)
    cat("Interval Method (dRisk):\n")
    cat("  dRisk:", sprintf("%.4f", x$drisk_interval),
        sprintf("(%d of %d records at risk)", n_at_risk, x$n_original), "\n")
    cat("  Outlier parameter:", x$outlier_par, "\n")
    cat("  Interval widths:\n")
    for (v in names(x$interval_widths)) {
      cat("   ", v, ":", sprintf("%.4f", x$interval_widths[v]), "\n")
    }
    cat("\n")
  }

  if (!is.na(x$drisk_rmd)) {
    n_at_risk <- sum(x$at_risk_rmd)
    cat("RMD Method (dRiskRMD):\n")
    cat("  dRiskRMD:", sprintf("%.4f", x$drisk_rmd),
        sprintf("(%d of %d records at risk)", n_at_risk, x$n_original), "\n")
    cat("  Alpha:", x$alpha, "\n")
    cat("  Chi-squared threshold:", sprintf("%.4f",
        sqrt(qchisq(1 - x$alpha, df = length(x$vars)))), "\n")
    cat("  Min RMD range: [",
        sprintf("%.4f", min(x$min_rmd)), ",",
        sprintf("%.4f", max(x$min_rmd)), "]\n\n")
  }

  cat("Privacy Assessment:")
  if (x$privacy_pass) {
    cat(" PASS\n")
    cat("  Disclosure risk is within acceptable bounds (<= 0.1).\n")
  } else {
    cat(" WARNING\n")
    cat("  Elevated disclosure risk detected (> 0.1).\n")
    cat("  Synthetic data may allow re-identification of original records.\n")
  }

  invisible(x)
}

#' Summary method for drisk objects
#'
#' @param object an object of class "drisk"
#' @param ... additional arguments (ignored)
#' @export
summary.drisk <- function(object, ...) {
  summ <- list(
    drisk_interval = object$drisk_interval,
    drisk_rmd = object$drisk_rmd,
    privacy_pass = object$privacy_pass,
    method = object$method,
    n_original = object$n_original,
    n_synthetic = object$n_synthetic,
    vars = object$vars,
    outlier_par = object$outlier_par,
    alpha = object$alpha
  )

  # Interval method details
  if (!is.na(object$drisk_interval)) {
    summ$n_at_risk_interval <- sum(object$at_risk_interval)
    summ$pct_at_risk_interval <- 100 * object$drisk_interval
    summ$interval_widths <- object$interval_widths
  }

  # RMD method details
  if (!is.na(object$drisk_rmd)) {
    summ$n_at_risk_rmd <- sum(object$at_risk_rmd)
    summ$pct_at_risk_rmd <- 100 * object$drisk_rmd
    summ$threshold_rmd <- sqrt(qchisq(1 - object$alpha, df = length(object$vars)))
    summ$min_rmd_quantiles <- quantile(object$min_rmd,
                                        probs = c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1),
                                        na.rm = TRUE)
    summ$mean_min_rmd <- mean(object$min_rmd, na.rm = TRUE)
    summ$sd_min_rmd <- sd(object$min_rmd, na.rm = TRUE)
  }

  class(summ) <- "summary.drisk"
  return(summ)
}

#' Print method for summary.drisk objects
#'
#' @param x an object of class "summary.drisk"
#' @param ... additional arguments (ignored)
#' @export
print.summary.drisk <- function(x, ...) {
  cat("Summary: Disclosure Risk for Continuous Variables (dRisk)\n")
  cat("========================================================\n")
  cat("Method:", x$method, "\n")
  cat("Variables:", length(x$vars),
      paste0("(", paste(x$vars, collapse = ", "), ")"), "\n\n")

  cat("Dataset Sizes:\n")
  cat("  Original:", x$n_original, "| Synthetic:", x$n_synthetic, "\n\n")

  if (!is.na(x$drisk_interval)) {
    cat("Interval Method (dRisk):\n")
    cat("  dRisk:", sprintf("%.4f", x$drisk_interval), "\n")
    cat("  At-risk records:", x$n_at_risk_interval,
        sprintf("(%.1f%%)", x$pct_at_risk_interval), "\n")
    cat("  Outlier parameter:", x$outlier_par, "\n")
    cat("  Interval widths:\n")
    for (v in names(x$interval_widths)) {
      cat("   ", v, ":", sprintf("%.6f", x$interval_widths[v]), "\n")
    }
    cat("\n")
  }

  if (!is.na(x$drisk_rmd)) {
    cat("RMD Method (dRiskRMD):\n")
    cat("  dRiskRMD:", sprintf("%.4f", x$drisk_rmd), "\n")
    cat("  At-risk records:", x$n_at_risk_rmd,
        sprintf("(%.1f%%)", x$pct_at_risk_rmd), "\n")
    cat("  Alpha:", x$alpha, "\n")
    cat("  Chi-squared threshold:", sprintf("%.4f", x$threshold_rmd), "\n")
    cat("  Min RMD: mean =", sprintf("%.4f", x$mean_min_rmd),
        "| sd =", sprintf("%.4f", x$sd_min_rmd), "\n")
    cat("  Min RMD quantiles:\n")
    print(round(x$min_rmd_quantiles, 4))
    cat("\n")
  }

  cat("Privacy:", ifelse(x$privacy_pass, "PASS", "WARNING"), "\n")

  invisible(x)
}

#' Plot method for drisk objects
#'
#' @param x an object of class "drisk"
#' @param y not used
#' @param ... additional arguments passed to plotting functions
#' @param which which plot(s) to show: 1 = per-record risk indicator strip chart,
#'   2 = minimum RMD distance distribution (only available when method includes
#'   "rmd")
#' @export
plot.drisk <- function(x, y = NULL, ..., which = 1) {
  show <- rep(FALSE, 2)
  show[which] <- TRUE

  n_plots <- sum(show)
  if (n_plots > 1) {
    op <- par(mfrow = c(1, 2))
    on.exit(par(op))
  }

  if (show[1]) {
    # Per-record risk indicator strip chart
    has_interval <- !is.na(x$drisk_interval)
    has_rmd <- !is.na(x$drisk_rmd)

    if (has_interval && has_rmd) {
      # Both methods: side by side barplot of risk fractions
      vals <- c(x$drisk_interval, x$drisk_rmd)
      names_vals <- c("Interval\n(dRisk)", "RMD\n(dRiskRMD)")
      cols <- c("coral", "steelblue")

      bp <- barplot(vals, names.arg = names_vals,
                    main = "Disclosure Risk for Continuous Variables",
                    ylab = "Fraction at risk",
                    col = cols,
                    ylim = c(0, max(0.2, max(vals) * 1.3)), ...)
      abline(h = 0.1, col = "red", lwd = 2, lty = 2)
      text(bp, vals + max(vals) * 0.05,
           labels = sprintf("%.3f", vals), cex = 0.9)
      legend("topright", "threshold = 0.1", col = "red", lty = 2, lwd = 2,
             cex = 0.8)

    } else if (has_interval) {
      # Interval only: stripchart of at-risk vs not-at-risk
      risk_status <- ifelse(x$at_risk_interval, "At risk", "Safe")
      risk_colors <- ifelse(x$at_risk_interval, "coral", "steelblue")
      stripchart(seq_len(x$n_original) ~ factor(risk_status,
                                                 levels = c("Safe", "At risk")),
                 main = paste("dRisk: Record Risk Status\n",
                              sprintf("%.1f%% at risk",
                                      100 * x$drisk_interval)),
                 xlab = "Record index",
                 ylab = "",
                 pch = 16, col = c("steelblue", "coral"),
                 method = "jitter", jitter = 0.2, ...)

    } else if (has_rmd) {
      # RMD only: stripchart of at-risk vs not-at-risk
      risk_status <- ifelse(x$at_risk_rmd, "At risk", "Safe")
      stripchart(seq_len(x$n_original) ~ factor(risk_status,
                                                 levels = c("Safe", "At risk")),
                 main = paste("dRiskRMD: Record Risk Status\n",
                              sprintf("%.1f%% at risk",
                                      100 * x$drisk_rmd)),
                 xlab = "Record index",
                 ylab = "",
                 pch = 16, col = c("steelblue", "coral"),
                 method = "jitter", jitter = 0.2, ...)
    }
  }

  if (show[2]) {
    if (is.null(x$min_rmd)) {
      message("Plot 2 (min RMD distribution) requires method = 'rmd' or 'both'.")
    } else {
      # Minimum RMD distribution
      threshold <- sqrt(qchisq(1 - x$alpha, df = length(x$vars)))
      hist(x$min_rmd,
           breaks = 30,
           main = paste("Minimum Robust Mahalanobis Distance\n",
                        sprintf("%.1f%% below threshold",
                                100 * x$drisk_rmd)),
           xlab = "Min RMD to closest synthetic record",
           col = "lightblue",
           border = "white", ...)
      abline(v = threshold, col = "red", lwd = 2, lty = 2)
      legend("topright",
             legend = c(sprintf("threshold = %.2f", threshold),
                        sprintf("alpha = %.2f", x$alpha)),
             col = c("red", NA), lty = c(2, NA), lwd = c(2, NA),
             cex = 0.8)
    }
  }
}
