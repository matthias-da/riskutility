#' Confidence Interval Proximity
#'
#' Computes the proximity of confidence intervals for numeric variables between
#' original and synthetic datasets. This utility measure assesses how well the
#' synthetic data preserves point estimates and their uncertainty.
#'
#' @param X A data.frame or data.table containing the original dataset.
#' @param Y A data.frame or data.table containing the synthetic/anonymized dataset.
#' @param vars Character vector of numeric variable names to compare.
#'   If NULL (default), all common numeric variables are used.
#' @param conf.level Numeric, confidence level for intervals. Default 0.95.
#' @param weight_X Optional character string specifying sampling weight variable in X.
#' @param weight_Y Optional character string specifying sampling weight variable in Y.
#' @param na.rm Logical, whether to remove NA values. Default TRUE.
#' @param ... additional arguments passed to methods (currently unused)
#'
#' @return An object of class "ci_proximity" containing:
#' \itemize{
#'   \item per_variable: data.frame with CI comparison for each variable
#'   \item proximity_mean: mean proximity score across all variables
#'   \item overlap_mean: mean overlap coefficient across all variables
#'   \item relative_error_mean: mean relative error of synthetic means
#'   \item n_vars: number of variables compared
#'   \item vars: variable names used
#'   \item conf.level: confidence level used
#' }
#'
#' @details
#' For each numeric variable, this function computes:
#' \itemize{
#'   \item \strong{CI Overlap}: Proportion of the CIs that overlap, ranging from
#'     0 (no overlap) to 1 (complete overlap or one contains the other)
#'   \item \strong{Relative Error}: |mean_Y - mean_X| / |mean_X|, measuring
#'     how close the synthetic mean is to the original mean
#'   \item \strong{Proximity Score}: Combined measure accounting for both
#'     mean accuracy and CI overlap
#' }
#'
#' The CI overlap coefficient is computed as:
#' \deqn{Overlap = \frac{\max(0, \min(U_X, U_Y) - \max(L_X, L_Y))}{\max(U_X, U_Y) - \min(L_X, L_Y)}}
#'
#' where L and U are lower and upper CI bounds.
#'
#' For utility assessment, higher proximity scores indicate better preservation
#' of statistical properties in synthetic data.
#'
#' @seealso \code{\link{ci_overlap}} for basic CI overlap calculation,
#'   \code{\link{compare_ks_test}} for distribution comparison
#'
#' @author Matthias Templ
#' @family utility
#' @export
#' @importFrom stats t.test weighted.mean complete.cases qt sd
#'
#' @examples
#' set.seed(123)
#' # Original data
#' X <- data.frame(
#'   income = rnorm(500, mean = 50000, sd = 10000),
#'   age = rnorm(500, mean = 40, sd = 10),
#'   score = rnorm(500, mean = 100, sd = 15)
#' )
#'
#' # Good synthetic data (similar means)
#' Y_good <- data.frame(
#'   income = rnorm(500, mean = 50000, sd = 10000),
#'   age = rnorm(500, mean = 40, sd = 10),
#'   score = rnorm(500, mean = 100, sd = 15)
#' )
#'
#' # Poor synthetic data (shifted means)
#' Y_poor <- data.frame(
#'   income = rnorm(500, mean = 55000, sd = 10000),
#'   age = rnorm(500, mean = 45, sd = 10),
#'   score = rnorm(500, mean = 90, sd = 15)
#' )
#'
#' result_good <- ci_proximity(X, Y_good)
#' print(result_good)
#' summary(result_good)
#'
#' result_poor <- ci_proximity(X, Y_poor)
#' print(result_poor)
ci_proximity <- function(X, ...) {
  UseMethod("ci_proximity")
}

#' @rdname ci_proximity
#' @export
ci_proximity.synth_pair <- function(X, ...) {
  ci_proximity.default(
    X = X$original,
    Y = X$synthetic,
    vars = X$num_vars,  # Use auto-detected numeric variables
    weight_X = X$weight_original,
    weight_Y = X$weight_synthetic,
    ...
  )
}

#' @rdname ci_proximity
#' @export
ci_proximity.default <- function(X, Y,
                                 vars = NULL,
                                 conf.level = 0.95,
                                 weight_X = NULL,
                                 weight_Y = NULL,
                                 na.rm = TRUE,
                                 ...) {

  # Convert to data.frame
  X <- as.data.frame(X)
  Y <- as.data.frame(Y)

  # Auto-detect numeric variables if not specified
  if (is.null(vars)) {
    num_vars_X <- names(X)[sapply(X, is.numeric)]
    num_vars_Y <- names(Y)[sapply(Y, is.numeric)]
    vars <- intersect(num_vars_X, num_vars_Y)
    # Exclude weight variables
    if (!is.null(weight_X)) vars <- setdiff(vars, weight_X)
    if (!is.null(weight_Y)) vars <- setdiff(vars, weight_Y)
  }

  if (length(vars) == 0) {
    stop("No numeric variables found or specified for comparison.")
  }

  # Check variables exist
  missing_X <- setdiff(vars, names(X))
  missing_Y <- setdiff(vars, names(Y))
  if (length(missing_X) > 0) {
    stop(paste("Variables missing in X:", paste(missing_X, collapse = ", ")))
  }
  if (length(missing_Y) > 0) {
    stop(paste("Variables missing in Y:", paste(missing_Y, collapse = ", ")))
  }

  # Helper function to compute weighted CI
  compute_ci <- function(x, w = NULL, conf.level) {
    if (na.rm) {
      if (!is.null(w)) {
        valid <- !is.na(x) & !is.na(w)
        x <- x[valid]
        w <- w[valid]
      } else {
        x <- x[!is.na(x)]
      }
    }

    n <- length(x)
    if (n < 2) {
      return(list(mean = NA, lower = NA, upper = NA, se = NA, n = n))
    }

    if (!is.null(w) && length(w) == n) {
      # Weighted mean and SE
      m <- weighted.mean(x, w)
      # Weighted variance (using frequency weights interpretation)
      v <- sum(w * (x - m)^2) / sum(w)
      se <- sqrt(v / n)
    } else {
      m <- mean(x)
      se <- sd(x) / sqrt(n)
    }

    alpha <- 1 - conf.level
    t_crit <- qt(1 - alpha/2, df = n - 1)
    lower <- m - t_crit * se
    upper <- m + t_crit * se

    list(mean = m, lower = lower, upper = upper, se = se, n = n)
  }

  # Helper function to compute CI overlap coefficient
  compute_overlap <- function(ci_X, ci_Y) {
    # Handle missing values
    if (is.na(ci_X$lower) || is.na(ci_Y$lower)) {
      return(NA_real_)
    }

    # Overlap region
    overlap_lower <- max(ci_X$lower, ci_Y$lower)
    overlap_upper <- min(ci_X$upper, ci_Y$upper)

    # Total span
    total_lower <- min(ci_X$lower, ci_Y$lower)
    total_upper <- max(ci_X$upper, ci_Y$upper)

    # Overlap coefficient
    if (overlap_upper > overlap_lower) {
      overlap <- (overlap_upper - overlap_lower) / (total_upper - total_lower)
    } else {
      overlap <- 0  # No overlap
    }

    return(overlap)
  }

  # Compute for each variable
  results_list <- vector("list", length(vars))
  names(results_list) <- vars

  for (v in vars) {
    wx <- if (!is.null(weight_X) && weight_X %in% names(X)) X[[weight_X]] else NULL
    wy <- if (!is.null(weight_Y) && weight_Y %in% names(Y)) Y[[weight_Y]] else NULL

    ci_X <- compute_ci(X[[v]], wx, conf.level)
    ci_Y <- compute_ci(Y[[v]], wy, conf.level)

    # CI overlap
    overlap <- compute_overlap(ci_X, ci_Y)

    # Relative error of mean
    if (!is.na(ci_X$mean) && abs(ci_X$mean) > 1e-10) {
      rel_error <- abs(ci_Y$mean - ci_X$mean) / abs(ci_X$mean)
    } else if (!is.na(ci_X$mean) && !is.na(ci_Y$mean)) {
      # If original mean is near zero, use absolute difference scaled by SD
      pooled_sd <- sqrt((ci_X$se^2 * ci_X$n + ci_Y$se^2 * ci_Y$n) / (ci_X$n + ci_Y$n))
      if (pooled_sd > 0) {
        rel_error <- abs(ci_Y$mean - ci_X$mean) / pooled_sd
      } else {
        rel_error <- 0
      }
    } else {
      rel_error <- NA_real_
    }

    # Proximity score: combination of overlap and accuracy
    # Higher is better; 1 = perfect match
    if (!is.na(overlap) && !is.na(rel_error)) {
      # Transform relative error to [0,1] scale (1 = no error)
      accuracy <- 1 / (1 + rel_error)
      # Combine overlap and accuracy
      proximity <- (overlap + accuracy) / 2
    } else {
      proximity <- NA_real_
      accuracy <- NA_real_
    }

    # Does synthetic CI contain original mean?
    contains_original_mean <- !is.na(ci_Y$lower) && !is.na(ci_X$mean) &&
      (ci_Y$lower <= ci_X$mean && ci_X$mean <= ci_Y$upper)

    results_list[[v]] <- data.frame(
      variable = v,
      mean_X = ci_X$mean,
      mean_Y = ci_Y$mean,
      ci_lower_X = ci_X$lower,
      ci_upper_X = ci_X$upper,
      ci_lower_Y = ci_Y$lower,
      ci_upper_Y = ci_Y$upper,
      overlap = overlap,
      relative_error = rel_error,
      accuracy = if (exists("accuracy")) accuracy else NA_real_,
      proximity = proximity,
      contains_orig_mean = contains_original_mean,
      n_X = ci_X$n,
      n_Y = ci_Y$n
    )
  }

  per_variable <- do.call(rbind, results_list)
  rownames(per_variable) <- NULL

  # Summary statistics
  valid_prox <- per_variable$proximity[!is.na(per_variable$proximity)]
  valid_overlap <- per_variable$overlap[!is.na(per_variable$overlap)]
  valid_error <- per_variable$relative_error[!is.na(per_variable$relative_error)]

  result <- list(
    per_variable = per_variable,
    proximity_mean = if (length(valid_prox) > 0) mean(valid_prox) else NA_real_,
    overlap_mean = if (length(valid_overlap) > 0) mean(valid_overlap) else NA_real_,
    relative_error_mean = if (length(valid_error) > 0) mean(valid_error) else NA_real_,
    pct_containing_orig_mean = 100 * mean(per_variable$contains_orig_mean, na.rm = TRUE),
    n_vars = length(vars),
    vars = vars,
    conf.level = conf.level,
    n_X = nrow(X),
    n_Y = nrow(Y)
  )

  class(result) <- "ci_proximity"
  return(result)
}


#' Print method for ci_proximity objects
#'
#' @param x an object of class "ci_proximity"
#' @param ... additional arguments (ignored)
#' @return The input object, invisibly.
#' @export
print.ci_proximity <- function(x, ...) {
  cat("Confidence Interval Proximity - Statistical Inference Preservation\n")
  cat("===================================================================\n\n")

  cat("Configuration:\n")
  cat("  Confidence level:", x$conf.level * 100, "%\n")
  cat("  Variables compared:", x$n_vars, "\n")
  cat("  Sample sizes: X =", x$n_X, ", Y =", x$n_Y, "\n\n")

  cat("Summary Metrics:\n")
  cat("  Mean proximity score:    ", sprintf("%.4f", x$proximity_mean),
      "(1 = perfect)\n")
  cat("  Mean CI overlap:         ", sprintf("%.4f", x$overlap_mean),
      "(1 = complete)\n")
  cat("  Mean relative error:     ", sprintf("%.4f", x$relative_error_mean),
      "(0 = perfect)\n")
  cat("  CIs containing orig mean:", sprintf("%.1f%%", x$pct_containing_orig_mean), "\n\n")

  cat("Interpretation:\n")
  if (x$proximity_mean >= 0.9) {
    cat("  EXCELLENT: Statistical properties are very well preserved.\n")
  } else if (x$proximity_mean >= 0.8) {
    cat("  GOOD: Statistical properties are reasonably preserved.\n")
  } else if (x$proximity_mean >= 0.7) {
    cat("  MODERATE: Some degradation of statistical properties.\n")
  } else {
    cat("  POOR: Significant differences in statistical properties.\n")
  }

  invisible(x)
}


#' Summary method for ci_proximity objects
#'
#' @param object an object of class "ci_proximity"
#' @param ... additional arguments (ignored)
#' @return A list of summary statistics for the corresponding object.
#' @export
summary.ci_proximity <- function(object, ...) {
  df <- object$per_variable

  summ <- list(
    per_variable = df,
    proximity_mean = object$proximity_mean,
    proximity_sd = sd(df$proximity, na.rm = TRUE),
    proximity_min = min(df$proximity, na.rm = TRUE),
    proximity_max = max(df$proximity, na.rm = TRUE),
    overlap_mean = object$overlap_mean,
    relative_error_mean = object$relative_error_mean,
    relative_error_max = max(df$relative_error, na.rm = TRUE),
    pct_containing_orig_mean = object$pct_containing_orig_mean,
    worst_variable = df$variable[which.min(df$proximity)],
    best_variable = df$variable[which.max(df$proximity)],
    n_vars = object$n_vars,
    conf.level = object$conf.level
  )

  class(summ) <- "summary.ci_proximity"
  return(summ)
}


#' Print method for summary.ci_proximity objects
#'
#' @param x an object of class "summary.ci_proximity"
#' @param ... additional arguments (ignored)
#' @return The input object, invisibly.
#' @export
print.summary.ci_proximity <- function(x, ...) {
  cat("Summary: Confidence Interval Proximity\n")
  cat("======================================\n\n")

  cat("Proximity Score Statistics:\n")
  cat("  Mean:", sprintf("%.4f", x$proximity_mean), "\n")
  cat("  SD:  ", sprintf("%.4f", x$proximity_sd), "\n")
  cat("  Min: ", sprintf("%.4f", x$proximity_min), "\n")
  cat("  Max: ", sprintf("%.4f", x$proximity_max), "\n\n")

  cat("Best preserved variable: ", x$best_variable, "\n")
  cat("Worst preserved variable:", x$worst_variable, "\n\n")

  cat("Per-Variable Results:\n")
  # Print subset of columns for readability
  print_df <- x$per_variable[, c("variable", "mean_X", "mean_Y",
                                  "overlap", "relative_error", "proximity")]
  print_df[, 2:6] <- round(print_df[, 2:6], 4)
  print(print_df, row.names = FALSE)

  invisible(x)
}


#' Plot method for ci_proximity objects
#'
#' @param x an object of class "ci_proximity"
#' @param y not used
#' @param ... additional arguments passed to plotting functions
#' @param which integer, which plot: 1 = CI comparison plot,
#'   2 = proximity scores bar chart, 3 = error vs overlap scatter
#' @importFrom graphics arrows barplot par points abline legend
#' @return No return value; called for the side effect of producing a plot.
#' @export
plot.ci_proximity <- function(x, y = NULL, ..., which = 1) {
  show <- rep(FALSE, 3)
  show[which] <- TRUE

  n_plots <- sum(show)
  if (n_plots > 1) {
    if (n_plots == 2) {
      op <- par(mfrow = c(1, 2))
    } else {
      op <- par(mfrow = c(2, 2))
    }
    on.exit(par(op))
  }

  df <- x$per_variable
  n_vars <- nrow(df)

  if (show[1]) {
    # CI comparison plot
    # Determine y-range
    all_vals <- c(df$ci_lower_X, df$ci_upper_X, df$ci_lower_Y, df$ci_upper_Y)
    y_range <- range(all_vals, na.rm = TRUE)

    plot(1:n_vars, df$mean_X, type = "n",
         xlim = c(0.5, n_vars + 0.5),
         ylim = y_range,
         xlab = "", ylab = "Value",
         main = "Confidence Interval Comparison",
         xaxt = "n", ...)
    axis(1, at = 1:n_vars, labels = df$variable, las = 2)

    # Original CIs (blue)
    arrows(1:n_vars - 0.1, df$ci_lower_X, 1:n_vars - 0.1, df$ci_upper_X,
           angle = 90, code = 3, length = 0.05, col = "steelblue", lwd = 2)
    points(1:n_vars - 0.1, df$mean_X, pch = 19, col = "steelblue")

    # Synthetic CIs (red)
    arrows(1:n_vars + 0.1, df$ci_lower_Y, 1:n_vars + 0.1, df$ci_upper_Y,
           angle = 90, code = 3, length = 0.05, col = "coral", lwd = 2)
    points(1:n_vars + 0.1, df$mean_Y, pch = 17, col = "coral")

    legend("topright", legend = c("Original", "Synthetic"),
           col = c("steelblue", "coral"), pch = c(19, 17), lwd = 2)
  }

  if (show[2]) {
    # Proximity scores bar chart
    colors <- ifelse(df$proximity >= 0.9, "forestgreen",
                     ifelse(df$proximity >= 0.8, "steelblue",
                            ifelse(df$proximity >= 0.7, "orange", "firebrick")))

    barplot(df$proximity,
            names.arg = df$variable,
            main = "Proximity Scores by Variable",
            ylab = "Proximity Score",
            col = colors,
            las = 2,
            ylim = c(0, 1),
            ...)
    abline(h = c(0.7, 0.8, 0.9), lty = 2, col = "gray50")
    abline(h = x$proximity_mean, lty = 1, col = "red", lwd = 2)
  }

  if (show[3]) {
    # Error vs overlap scatter
    plot(df$relative_error, df$overlap,
         xlab = "Relative Error (lower = better)",
         ylab = "CI Overlap (higher = better)",
         main = "Error vs Overlap Trade-off",
         pch = 19, col = "steelblue", cex = 1.5, ...)

    # Add variable labels
    text(df$relative_error, df$overlap, labels = df$variable,
         pos = 3, cex = 0.7)

    # Reference lines
    abline(h = 0.5, lty = 2, col = "gray50")
    abline(v = 0.1, lty = 2, col = "gray50")
  }
}
