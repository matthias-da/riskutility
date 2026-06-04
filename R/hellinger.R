#' Hellinger Distance for Categorical Distributions
#'
#' Computes the Hellinger distance between categorical variable distributions
#' in two datasets. The Hellinger distance measures the similarity between
#' probability distributions and is bounded between 0 (identical) and 1 (no overlap).
#'
#' @param X A data.frame or data.table containing the original dataset.
#' @param Y A data.frame or data.table containing the synthetic/anonymized dataset.
#' @param vars Character vector of categorical variable names to compare.
#'   If NULL (default), all common factor/character variables are used.
#' @param weight_X Optional character string specifying sampling weight variable in X.
#' @param weight_Y Optional character string specifying sampling weight variable in Y.
#' @param na.rm Logical, whether to remove NA values before computation. Default TRUE.
#' @param ... additional arguments passed to methods (currently unused)
#'
#' @return An object of class "hellinger" containing:
#' \itemize{
#'   \item per_variable: data.frame with Hellinger distance for each variable
#'   \item hellinger_mean: mean Hellinger distance across all variables
#'   \item hellinger_max: maximum Hellinger distance across variables
#'   \item n_vars: number of variables compared
#'   \item vars: variable names used
#'   \item utility_score: 1 - hellinger_mean (higher = better utility)
#' }
#'
#' @details
#' The Hellinger distance between two discrete probability distributions P and Q is:
#' \deqn{H(P, Q) = \frac{1}{\sqrt{2}} \sqrt{\sum_{i} (\sqrt{p_i} - \sqrt{q_i})^2}}
#'
#' Properties:
#' \itemize{
#'   \item \strong{H = 0}: Distributions are identical
#'   \item \strong{H = 1}: Distributions have no overlap (disjoint support)
#'   \item Symmetric: H(P,Q) = H(Q,P)
#'   \item Bounded: 0 <= H <= 1
#' }
#'
#' For utility assessment, lower Hellinger distance indicates better preservation
#' of categorical distributions in synthetic data.
#'
#' @seealso \code{\link{compare_chisq_gof}} for chi-squared goodness of fit,
#'   \code{\link{energy_distance}} for multivariate numeric comparison
#'
#' @references
#' Hellinger, E. (1909). Neue Begruendung der Theorie quadratischer Formen von
#' unendlichvielen Veraenderlichen. Journal fuer die reine und angewandte
#' Mathematik, 136, 210-271.
#'
#' @author Matthias Templ
#' @family utility
#' @export
#' @importFrom data.table as.data.table
#' @importFrom stats xtabs complete.cases
#'
#' @examples
#' set.seed(123)
#' # Original data
#' X <- data.frame(
#'   gender = sample(c("Male", "Female"), 500, replace = TRUE),
#'   region = sample(c("North", "South", "East", "West"), 500, replace = TRUE),
#'   education = sample(c("High", "Medium", "Low"), 500, replace = TRUE)
#' )
#'
#' # Good synthetic data (similar distributions)
#' Y_good <- data.frame(
#'   gender = sample(c("Male", "Female"), 500, replace = TRUE),
#'   region = sample(c("North", "South", "East", "West"), 500, replace = TRUE),
#'   education = sample(c("High", "Medium", "Low"), 500, replace = TRUE)
#' )
#'
#' # Poor synthetic data (different distributions)
#' Y_poor <- data.frame(
#'   gender = sample(c("Male", "Female"), 500, replace = TRUE, prob = c(0.9, 0.1)),
#'   region = sample(c("North", "South", "East", "West"), 500, replace = TRUE,
#'                   prob = c(0.7, 0.1, 0.1, 0.1)),
#'   education = sample(c("High", "Medium", "Low"), 500, replace = TRUE,
#'                      prob = c(0.1, 0.1, 0.8))
#' )
#'
#' result_good <- hellinger(X, Y_good)
#' print(result_good)
#'
#' result_poor <- hellinger(X, Y_poor)
#' print(result_poor)
hellinger <- function(X, ...) {
  UseMethod("hellinger")
}

#' @rdname hellinger
#' @export
hellinger.synth_pair <- function(X, ...) {
  hellinger.default(
    X = X$original,
    Y = X$synthetic,
    vars = X$cat_vars,  # Use auto-detected categorical variables
    weight_X = X$weight_original,
    weight_Y = X$weight_synthetic,
    ...
  )
}

#' @rdname hellinger
#' @export
hellinger.default <- function(X, Y,
                              vars = NULL,
                              weight_X = NULL,
                              weight_Y = NULL,
                              na.rm = TRUE,
                              ...) {

  # Convert to data.table for efficiency
  X <- as.data.table(X)
  Y <- as.data.table(Y)

  # Auto-detect categorical variables if not specified
  if (is.null(vars)) {
    cat_vars_X <- names(X)[sapply(X, function(col) is.factor(col) || is.character(col))]
    cat_vars_Y <- names(Y)[sapply(Y, function(col) is.factor(col) || is.character(col))]
    vars <- intersect(cat_vars_X, cat_vars_Y)
    # Exclude weight variables
    if (!is.null(weight_X)) vars <- setdiff(vars, weight_X)
    if (!is.null(weight_Y)) vars <- setdiff(vars, weight_Y)
  }

  if (length(vars) == 0) {
    stop("No categorical variables found or specified for comparison.")
  }

  # Check variables exist in both datasets
  missing_X <- setdiff(vars, names(X))
  missing_Y <- setdiff(vars, names(Y))
  if (length(missing_X) > 0) {
    stop(paste("Variables missing in X:", paste(missing_X, collapse = ", ")))
  }
  if (length(missing_Y) > 0) {
    stop(paste("Variables missing in Y:", paste(missing_Y, collapse = ", ")))
  }

  # Helper function to compute Hellinger distance for a single variable
  compute_hellinger <- function(x, y, wx = NULL, wy = NULL) {
    # Remove NAs if requested
    if (na.rm) {
      if (!is.null(wx)) {
        valid_x <- !is.na(x) & !is.na(wx)
        x <- x[valid_x]
        wx <- wx[valid_x]
      } else {
        x <- x[!is.na(x)]
      }
      if (!is.null(wy)) {
        valid_y <- !is.na(y) & !is.na(wy)
        y <- y[valid_y]
        wy <- wy[valid_y]
      } else {
        y <- y[!is.na(y)]
      }
    }

    if (length(x) == 0 || length(y) == 0) {
      return(NA_real_)
    }

    # Compute weighted frequency tables
    if (!is.null(wx)) {
      tab_x <- tapply(wx, x, sum, default = 0)
    } else {
      tab_x <- table(x)
    }

    if (!is.null(wy)) {
      tab_y <- tapply(wy, y, sum, default = 0)
    } else {
      tab_y <- table(y)
    }

    # Convert to proportions
    p_x <- tab_x / sum(tab_x)
    p_y <- tab_y / sum(tab_y)

    # Align categories (union of all levels)
    all_levels <- union(names(p_x), names(p_y))

    p_x_aligned <- rep(0, length(all_levels))
    names(p_x_aligned) <- all_levels
    p_x_aligned[names(p_x)] <- p_x

    p_y_aligned <- rep(0, length(all_levels))
    names(p_y_aligned) <- all_levels
    p_y_aligned[names(p_y)] <- p_y

    # Hellinger distance formula
    # H(P,Q) = (1/sqrt(2)) * sqrt(sum((sqrt(p) - sqrt(q))^2))
    H <- (1 / sqrt(2)) * sqrt(sum((sqrt(p_x_aligned) - sqrt(p_y_aligned))^2))

    return(H)
  }

  # Compute Hellinger distance for each variable
  results_list <- vector("list", length(vars))
  names(results_list) <- vars

  for (v in vars) {
    wx <- if (!is.null(weight_X) && weight_X %in% names(X)) X[[weight_X]] else NULL
    wy <- if (!is.null(weight_Y) && weight_Y %in% names(Y)) Y[[weight_Y]] else NULL

    H <- compute_hellinger(X[[v]], Y[[v]], wx, wy)

    # Also compute number of categories
    n_cats_x <- length(unique(na.omit(X[[v]])))
    n_cats_y <- length(unique(na.omit(Y[[v]])))
    n_cats_union <- length(union(unique(na.omit(X[[v]])), unique(na.omit(Y[[v]]))))

    results_list[[v]] <- data.frame(
      variable = v,
      hellinger = H,
      n_categories_X = n_cats_x,
      n_categories_Y = n_cats_y,
      n_categories_union = n_cats_union
    )
  }

  per_variable <- do.call(rbind, results_list)
  rownames(per_variable) <- NULL

  # Compute summary statistics
  valid_H <- per_variable$hellinger[!is.na(per_variable$hellinger)]
  hellinger_mean <- if (length(valid_H) > 0) mean(valid_H) else NA_real_
  hellinger_max <- if (length(valid_H) > 0) max(valid_H) else NA_real_
  hellinger_min <- if (length(valid_H) > 0) min(valid_H) else NA_real_

  result <- list(
    per_variable = per_variable,
    hellinger_mean = hellinger_mean,
    hellinger_max = hellinger_max,
    hellinger_min = hellinger_min,
    n_vars = length(vars),
    vars = vars,
    utility_score = 1 - hellinger_mean,
    n_X = nrow(X),
    n_Y = nrow(Y)
  )

  class(result) <- "hellinger"
  return(result)
}


#' Print method for hellinger objects
#'
#' @param x an object of class "hellinger"
#' @param ... additional arguments (ignored)
#' @return The input object, invisibly.
#' @export
print.hellinger <- function(x, ...) {
  cat("Hellinger Distance - Categorical Distribution Comparison\n")
  cat("=========================================================\n\n")

  cat("Dataset Sizes:\n")
  cat("  Original (X):", x$n_X, "| Synthetic (Y):", x$n_Y, "\n")
  cat("  Variables compared:", x$n_vars, "\n\n")

  cat("Summary:\n")
  cat("  Mean Hellinger distance:", sprintf("%.4f", x$hellinger_mean), "\n")
  cat("  Max Hellinger distance: ", sprintf("%.4f", x$hellinger_max), "\n")
  cat("  Min Hellinger distance: ", sprintf("%.4f", x$hellinger_min), "\n")
  cat("  Utility score (1-mean):", sprintf("%.4f", x$utility_score), "\n\n")

  cat("Interpretation:\n")
  if (x$hellinger_mean < 0.1) {
    cat("  EXCELLENT: Categorical distributions are very similar.\n")
  } else if (x$hellinger_mean < 0.2) {
    cat("  GOOD: Categorical distributions are reasonably preserved.\n")
  } else if (x$hellinger_mean < 0.3) {
    cat("  MODERATE: Some differences in categorical distributions.\n")
  } else {
    cat("  POOR: Significant differences in categorical distributions.\n")
  }

  invisible(x)
}


#' Summary method for hellinger objects
#'
#' @param object an object of class "hellinger"
#' @param ... additional arguments (ignored)
#' @return A list of summary statistics for the corresponding object.
#' @export
summary.hellinger <- function(object, ...) {
  summ <- list(
    per_variable = object$per_variable,
    hellinger_mean = object$hellinger_mean,
    hellinger_max = object$hellinger_max,
    hellinger_min = object$hellinger_min,
    hellinger_sd = sd(object$per_variable$hellinger, na.rm = TRUE),
    n_vars = object$n_vars,
    utility_score = object$utility_score,
    worst_variable = object$per_variable$variable[which.max(object$per_variable$hellinger)],
    best_variable = object$per_variable$variable[which.min(object$per_variable$hellinger)]
  )

  class(summ) <- "summary.hellinger"
  return(summ)
}


#' Print method for summary.hellinger objects
#'
#' @param x an object of class "summary.hellinger"
#' @param ... additional arguments (ignored)
#' @return The input object, invisibly.
#' @export
print.summary.hellinger <- function(x, ...) {
  cat("Summary: Hellinger Distance\n")
  cat("===========================\n\n")

  cat("Overall Statistics:\n")
  cat("  Mean:", sprintf("%.4f", x$hellinger_mean), "\n")
  cat("  SD:  ", sprintf("%.4f", x$hellinger_sd), "\n")
  cat("  Min: ", sprintf("%.4f", x$hellinger_min), "\n")
  cat("  Max: ", sprintf("%.4f", x$hellinger_max), "\n\n")

  cat("Best preserved variable: ", x$best_variable, "\n")
  cat("Worst preserved variable:", x$worst_variable, "\n\n")

  cat("Per-Variable Results:\n")
  print(x$per_variable, row.names = FALSE)

  invisible(x)
}


#' Plot method for hellinger objects
#'
#' @param x an object of class "hellinger"
#' @param y not used
#' @param ... additional arguments passed to plotting functions
#' @param which integer, which plot: 1 = bar chart of Hellinger distances,
#'   2 = dotchart sorted by distance
#' @importFrom graphics barplot dotchart abline par
#' @return No return value; called for the side effect of producing a plot.
#' @export
plot.hellinger <- function(x, y = NULL, ..., which = 1) {
  show <- rep(FALSE, 2)
  show[which] <- TRUE

  n_plots <- sum(show)
  if (n_plots > 1) {
    op <- par(mfrow = c(1, 2))
    on.exit(par(op))
  }

  df <- x$per_variable

  if (show[1]) {
    # Bar chart
    colors <- ifelse(df$hellinger < 0.1, "forestgreen",
                     ifelse(df$hellinger < 0.2, "steelblue",
                            ifelse(df$hellinger < 0.3, "orange", "firebrick")))

    barplot(df$hellinger,
            names.arg = df$variable,
            main = "Hellinger Distance by Variable",
            ylab = "Hellinger Distance",
            col = colors,
            las = 2,
            ylim = c(0, max(1, max(df$hellinger, na.rm = TRUE) * 1.1)),
            ...)
    abline(h = c(0.1, 0.2, 0.3), lty = 2, col = "gray50")
    abline(h = x$hellinger_mean, lty = 1, col = "red", lwd = 2)
  }

  if (show[2]) {
    # Sorted dotchart
    ord <- order(df$hellinger)
    dotchart(df$hellinger[ord],
             labels = df$variable[ord],
             main = "Hellinger Distance (sorted)",
             xlab = "Hellinger Distance",
             pch = 19,
             color = "steelblue",
             ...)
    abline(v = x$hellinger_mean, lty = 2, col = "red", lwd = 2)
  }
}
