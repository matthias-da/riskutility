#' Energy Distance for Multivariate Numeric Data
#'
#' Computes the energy distance between multivariate numeric distributions
#' in two datasets. Energy distance is a statistical distance that characterizes
#' equality of distributions and is zero if and only if the distributions are identical.
#'
#' @param X A data.frame or data.table containing the original dataset.
#' @param Y A data.frame or data.table containing the synthetic/anonymized dataset.
#' @param vars Character vector of numeric variable names to compare.
#'   If NULL (default), all common numeric variables are used.
#' @param standardize Logical, whether to standardize variables before computing
#'   distances. Default TRUE (recommended for variables on different scales).
#' @param n_sample Integer, maximum sample size for computation. If datasets are
#'   larger, random sampling is used to reduce computation time. Default 1000.
#'   Set to NULL for no sampling (may be slow for large datasets).
#' @param na.rm Logical, whether to remove rows with NA values. Default TRUE.
#' @param seed Integer, random seed for reproducible sampling. Default NULL.
#' @param ... additional arguments passed to methods (currently unused)
#'
#' @return An object of class "energy_distance" containing:
#' \itemize{
#'   \item energy_distance: the computed energy distance
#'   \item energy_distance_normalized: energy distance divided by a reference value
#'   \item mean_dist_XY: mean distance between X and Y samples
#'   \item mean_dist_XX: mean distance within X samples
#'   \item mean_dist_YY: mean distance within Y samples
#'   \item n_X, n_Y: sample sizes used in computation
#'   \item n_vars: number of variables
#'   \item vars: variable names used
#'   \item standardized: whether standardization was applied
#'   \item utility_score: transformed score (higher = better utility)
#' }
#'
#' @details
#' The energy distance between distributions F and G is defined as:
#' \deqn{D_E(F, G) = 2E||X - Y|| - E||X - X'|| - E||Y - Y'||}
#'
#' where X, X' are independent samples from F, and Y, Y' are independent
#' samples from G. The energy distance:
#' \itemize{
#'   \item Is always non-negative
#'   \item Equals zero if and only if F = G
#'   \item Is sensitive to differences in both location and scale
#'   \item Does not require density estimation
#' }
#'
#' For synthetic data evaluation, lower energy distance indicates better
#' preservation of the multivariate numeric distribution.
#'
#' The computation uses Euclidean distances. For large datasets, random sampling
#' is applied to keep computation tractable (O(n^2) complexity).
#'
#' @seealso \code{\link{hellinger}} for categorical distributions,
#'   \code{\link{compare_wasserstein}} for univariate Wasserstein distance,
#'   \code{\link{gower}} for mixed-type data
#'
#' @references
#' Szekely, G. J. and Rizzo, M. L. (2013). Energy statistics: A class of statistics
#' based on distances. Journal of Statistical Planning and Inference, 143(8), 1249-1272.
#'
#' Rizzo, M. L. and Szekely, G. J. (2016). Energy Distance.
#' WIREs Computational Statistics, 8(1), 27-38.
#'
#' @author Matthias Templ
#' @family utility
#' @export
#' @importFrom stats complete.cases sd
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
#' # Good synthetic data (similar distribution)
#' Y_good <- data.frame(
#'   income = rnorm(500, mean = 50000, sd = 10000),
#'   age = rnorm(500, mean = 40, sd = 10),
#'   score = rnorm(500, mean = 100, sd = 15)
#' )
#'
#' # Poor synthetic data (shifted distribution)
#' Y_poor <- data.frame(
#'   income = rnorm(500, mean = 60000, sd = 15000),
#'   age = rnorm(500, mean = 45, sd = 15),
#'   score = rnorm(500, mean = 90, sd = 20)
#' )
#'
#' result_good <- energy_distance(X, Y_good, seed = 42)
#' print(result_good)
#'
#' result_poor <- energy_distance(X, Y_poor, seed = 42)
#' print(result_poor)
energy_distance <- function(X, ...) {
  UseMethod("energy_distance")
}

#' @rdname energy_distance
#' @export
energy_distance.synth_pair <- function(X, ...) {
  energy_distance.default(
    X = X$original,
    Y = X$synthetic,
    vars = X$num_vars,  # Use auto-detected numeric variables
    ...
  )
}

#' @rdname energy_distance
#' @export
energy_distance.default <- function(X, Y,
                                    vars = NULL,
                                    standardize = TRUE,
                                    n_sample = 1000,
                                    na.rm = TRUE,
                                    seed = NULL,
                                    ...) {

  # Convert to data.frame if needed
  X <- as.data.frame(X)
  Y <- as.data.frame(Y)

  # Auto-detect numeric variables if not specified
  if (is.null(vars)) {
    num_vars_X <- names(X)[sapply(X, is.numeric)]
    num_vars_Y <- names(Y)[sapply(Y, is.numeric)]
    vars <- intersect(num_vars_X, num_vars_Y)
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

  # Subset to selected variables
  X <- X[, vars, drop = FALSE]
  Y <- Y[, vars, drop = FALSE]

  # Handle missing values
  if (na.rm) {
    X <- X[complete.cases(X), , drop = FALSE]
    Y <- Y[complete.cases(Y), , drop = FALSE]
  }

  if (nrow(X) == 0) stop("No complete cases in X after removing NAs.")
  if (nrow(Y) == 0) stop("No complete cases in Y after removing NAs.")

  # Set seed for reproducibility
  if (!is.null(seed)) set.seed(seed)

  # Sample if datasets are large
  n_X_orig <- nrow(X)
  n_Y_orig <- nrow(Y)

  if (!is.null(n_sample)) {
    if (nrow(X) > n_sample) {
      X <- X[sample(nrow(X), n_sample), , drop = FALSE]
    }
    if (nrow(Y) > n_sample) {
      Y <- Y[sample(nrow(Y), n_sample), , drop = FALSE]
    }
  }

  n_X <- nrow(X)
  n_Y <- nrow(Y)

  # Standardize if requested (using pooled mean/sd)
  if (standardize && length(vars) > 0) {
    combined <- rbind(X, Y)
    for (v in vars) {
      m <- mean(combined[[v]], na.rm = TRUE)
      s <- sd(combined[[v]], na.rm = TRUE)
      if (s > 0) {
        X[[v]] <- (X[[v]] - m) / s
        Y[[v]] <- (Y[[v]] - m) / s
      }
    }
  }

  # Convert to matrices for distance computation
  X_mat <- as.matrix(X)
  Y_mat <- as.matrix(Y)

  # Helper function to compute mean pairwise Euclidean distance
  mean_euclidean_dist <- function(A, B) {
    # Compute squared Euclidean distances efficiently
    # ||a - b||^2 = ||a||^2 + ||b||^2 - 2 * a'b
    A_sq <- rowSums(A^2)
    B_sq <- rowSums(B^2)
    AB <- A %*% t(B)
    D_sq <- outer(A_sq, B_sq, "+") - 2 * AB
    D_sq[D_sq < 0] <- 0  # Handle numerical issues
    D <- sqrt(D_sq)
    mean(D)
  }

  # Compute mean distances
  # For within-sample distances, we need to exclude self-distances
  mean_dist_XY <- mean_euclidean_dist(X_mat, Y_mat)

  # Within X (excluding diagonal)
  if (n_X > 1) {
    X_sq <- rowSums(X_mat^2)
    XX <- X_mat %*% t(X_mat)
    D_XX_sq <- outer(X_sq, X_sq, "+") - 2 * XX
    D_XX_sq[D_XX_sq < 0] <- 0
    D_XX <- sqrt(D_XX_sq)
    diag(D_XX) <- NA
    mean_dist_XX <- mean(D_XX, na.rm = TRUE)
  } else {
    mean_dist_XX <- 0
  }

  # Within Y (excluding diagonal)
  if (n_Y > 1) {
    Y_sq <- rowSums(Y_mat^2)
    YY <- Y_mat %*% t(Y_mat)
    D_YY_sq <- outer(Y_sq, Y_sq, "+") - 2 * YY
    D_YY_sq[D_YY_sq < 0] <- 0
    D_YY <- sqrt(D_YY_sq)
    diag(D_YY) <- NA
    mean_dist_YY <- mean(D_YY, na.rm = TRUE)
  } else {
    mean_dist_YY <- 0
  }

  # Energy distance formula
  # E(X,Y) = 2*E||X-Y|| - E||X-X'|| - E||Y-Y'||
  energy_dist <- 2 * mean_dist_XY - mean_dist_XX - mean_dist_YY

  # Ensure non-negative (numerical precision)
  energy_dist <- max(0, energy_dist)

  # Normalized version: divide by sum of within-sample distances
  # This gives a scale-free measure
  reference <- (mean_dist_XX + mean_dist_YY) / 2
  if (reference > 0) {
    energy_normalized <- energy_dist / reference
  } else {
    energy_normalized <- energy_dist
  }

  # Utility score: transform so higher = better
  # Using exponential decay: exp(-energy_dist)
  utility_score <- exp(-energy_dist)

  result <- list(
    energy_distance = energy_dist,
    energy_distance_normalized = energy_normalized,
    mean_dist_XY = mean_dist_XY,
    mean_dist_XX = mean_dist_XX,
    mean_dist_YY = mean_dist_YY,
    n_X = n_X,
    n_Y = n_Y,
    n_X_original = n_X_orig,
    n_Y_original = n_Y_orig,
    n_vars = length(vars),
    vars = vars,
    standardized = standardize,
    utility_score = utility_score,
    sampled = !is.null(n_sample) && (n_X_orig > n_sample || n_Y_orig > n_sample)
  )

  class(result) <- "energy_distance"
  return(result)
}


#' Print method for energy_distance objects
#'
#' @param x an object of class "energy_distance"
#' @param ... additional arguments (ignored)
#' @export
print.energy_distance <- function(x, ...) {
  cat("Energy Distance - Multivariate Numeric Distribution Comparison\n")
  cat("===============================================================\n\n")

  cat("Dataset Sizes:\n")
  if (x$sampled) {
    cat("  Original (X):", x$n_X_original, "-> sampled to", x$n_X, "\n")
    cat("  Synthetic (Y):", x$n_Y_original, "-> sampled to", x$n_Y, "\n")
  } else {
    cat("  Original (X):", x$n_X, "\n")
    cat("  Synthetic (Y):", x$n_Y, "\n")
  }
  cat("  Variables:", x$n_vars, "\n")
  cat("  Standardized:", x$standardized, "\n\n")

  cat("Energy Distance:\n")
  cat("  Raw:       ", sprintf("%.4f", x$energy_distance), "\n")
  cat("  Normalized:", sprintf("%.4f", x$energy_distance_normalized), "\n")
  cat("  Utility:   ", sprintf("%.4f", x$utility_score), "(exp(-E), higher=better)\n\n")

  cat("Distance Components:\n")
  cat("  Mean dist(X,Y):", sprintf("%.4f", x$mean_dist_XY), "\n")
  cat("  Mean dist(X,X):", sprintf("%.4f", x$mean_dist_XX), "\n")
  cat("  Mean dist(Y,Y):", sprintf("%.4f", x$mean_dist_YY), "\n\n")

  cat("Interpretation:\n")
  if (x$energy_distance < 0.1) {
    cat("  EXCELLENT: Multivariate distributions are very similar.\n")
  } else if (x$energy_distance < 0.3) {
    cat("  GOOD: Multivariate distributions are reasonably similar.\n")
  } else if (x$energy_distance < 0.5) {
    cat("  MODERATE: Some differences in multivariate structure.\n")
  } else {
    cat("  POOR: Significant differences in multivariate distributions.\n")
  }

  invisible(x)
}


#' Summary method for energy_distance objects
#'
#' @param object an object of class "energy_distance"
#' @param ... additional arguments (ignored)
#' @export
summary.energy_distance <- function(object, ...) {
  summ <- list(
    energy_distance = object$energy_distance,
    energy_distance_normalized = object$energy_distance_normalized,
    mean_dist_XY = object$mean_dist_XY,
    mean_dist_XX = object$mean_dist_XX,
    mean_dist_YY = object$mean_dist_YY,
    utility_score = object$utility_score,
    n_X = object$n_X,
    n_Y = object$n_Y,
    n_vars = object$n_vars,
    vars = object$vars,
    standardized = object$standardized,
    # Ratio of between to within distances
    between_within_ratio = object$mean_dist_XY / ((object$mean_dist_XX + object$mean_dist_YY) / 2)
  )

  class(summ) <- "summary.energy_distance"
  return(summ)
}


#' Print method for summary.energy_distance objects
#'
#' @param x an object of class "summary.energy_distance"
#' @param ... additional arguments (ignored)
#' @export
print.summary.energy_distance <- function(x, ...) {
  cat("Summary: Energy Distance\n")
  cat("========================\n\n")

  cat("Variables (", x$n_vars, "):", paste(x$vars, collapse = ", "), "\n\n")

  cat("Distance Metrics:\n")
  cat("  Energy distance:       ", sprintf("%.4f", x$energy_distance), "\n")
  cat("  Normalized:            ", sprintf("%.4f", x$energy_distance_normalized), "\n")
  cat("  Between/Within ratio:  ", sprintf("%.4f", x$between_within_ratio), "\n")
  cat("  Utility score:         ", sprintf("%.4f", x$utility_score), "\n\n")

  cat("Sample Sizes: X =", x$n_X, ", Y =", x$n_Y, "\n")
  cat("Standardized:", x$standardized, "\n")

  invisible(x)
}


#' Plot method for energy_distance objects
#'
#' @param x an object of class "energy_distance"
#' @param y not used
#' @param ... additional arguments passed to plotting functions
#' @param which integer, which plot: 1 = distance comparison bar chart,
#'   2 = utility gauge
#' @importFrom graphics barplot par text rect
#' @export
plot.energy_distance <- function(x, y = NULL, ..., which = 1) {
  show <- rep(FALSE, 2)
  show[which] <- TRUE

  n_plots <- sum(show)
  if (n_plots > 1) {
    op <- par(mfrow = c(1, 2))
    on.exit(par(op))
  }

  if (show[1]) {
    # Bar chart comparing distance components
    vals <- c(x$mean_dist_XX, x$mean_dist_YY, x$mean_dist_XY)
    names(vals) <- c("Within X", "Within Y", "Between X-Y")

    barplot(vals,
            main = "Distance Components",
            ylab = "Mean Euclidean Distance",
            col = c("steelblue", "steelblue", "coral"),
            las = 1,
            ...)

    # Add energy distance annotation
    mtext(sprintf("Energy Distance = %.4f", x$energy_distance),
          side = 3, line = 0, cex = 0.9)
  }

  if (show[2]) {
    # Simple utility visualization
    plot(1, type = "n", xlim = c(0, 1), ylim = c(0, 1),
         xlab = "", ylab = "", axes = FALSE,
         main = "Utility Score")

    # Background gradient (poor to good)
    for (i in 1:100) {
      col <- rgb(1 - i/100, i/100, 0, alpha = 0.3)
      rect((i-1)/100, 0.3, i/100, 0.7, col = col, border = NA)
    }

    # Mark current score
    rect(0, 0.3, 1, 0.7, border = "black", lwd = 2)
    points(x$utility_score, 0.5, pch = 18, cex = 3, col = "black")

    # Labels
    text(0.5, 0.1, sprintf("Utility = %.3f", x$utility_score), cex = 1.2)
    text(0, 0.85, "Poor", adj = 0)
    text(1, 0.85, "Good", adj = 1)
  }
}
