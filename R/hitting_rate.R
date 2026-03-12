#' Hitting Rate
#'
#' Computes the Hitting Rate privacy metric for synthetic data. The Hitting Rate
#' measures the fraction of synthetic records that fall within a configurable
#' distance threshold of any original record. It bridges the Identical Match
#' Share (IMS, threshold = 0) and the full Distance to Closest Record (DCR)
#' distribution into a single, interpretable scalar.
#'
#' @param X data frame of original/training data
#' @param Y data frame of synthetic data
#' @param threshold numeric, distance threshold for a "hit" (default: 0.05).
#'   A synthetic record with minimum distance to any original record at or below
#'   this threshold is counted as a hit.
#' @param vars character vector of variable names to use for distance calculation.
#'   If NULL (default), all common variables between X and Y are used.
#' @param method character, distance method: "gower" (default, handles mixed types)
#'   or "euclidean" (numerical variables only)
#' @param na.rm logical, remove records with NA values (default: TRUE)
#' @param ... additional arguments passed to methods (currently unused)
#'
#' @return An object of class "hitting_rate" containing:
#' \itemize{
#'   \item rate: fraction of synthetic records within threshold distance of
#'     any original record
#'   \item min_distances: numeric vector of minimum distances per synthetic record
#'   \item hits: logical vector indicating which synthetic records are hits
#'     (min_distance <= threshold)
#'   \item n_hits: count of hits
#'   \item threshold: the threshold used
#'   \item privacy_pass: logical, TRUE if rate <= 0.1
#'   \item n_original: number of original records
#'   \item n_synthetic: number of synthetic records
#'   \item method: distance method used
#'   \item vars: variables used
#'   \item rate_at_zero: fraction of synthetic records with exact matches
#'     (min_distance == 0), equivalent to IMS
#' }
#'
#' @details
#' The Hitting Rate provides a threshold-based view of synthetic data privacy.
#' For each synthetic record, the minimum distance to any original record is
#' computed. Records with a minimum distance at or below the threshold are
#' counted as "hits" -- near copies that may pose a disclosure risk.
#'
#' The metric is defined as:
#' \deqn{HR(\tau) = \frac{|\{y \in Y : \min_{x \in X} d(x, y) \le \tau\}|}{|Y|}}
#'
#' The hitting rate generalises two existing metrics:
#' \itemize{
#'   \item At threshold = 0, the hitting rate equals the Identical Match Share
#'     (\code{\link{ims}}).
#'   \item The full distribution of minimum distances is the same as the DCR
#'     training distribution in \code{\link{dcr}}, but without requiring a
#'     holdout set.
#' }
#'
#' Interpretation:
#' \itemize{
#'   \item \strong{rate ~ 0}: Good privacy -- few synthetic records are close to
#'     any original record.
#'   \item \strong{rate > 0.1}: Privacy concern -- too many near copies.
#'   \item \strong{rate ~ 1}: Severe memorization -- nearly all synthetic records
#'     are close copies of original records.
#' }
#'
#' The default threshold of 0.05 is appropriate for Gower distances (which are
#' bounded between 0 and 1). For Euclidean distances on normalized data, adjust the
#' threshold according to the dimensionality and scale of the data.
#'
#' @seealso \code{\link{dcr}} for distance to closest record (with holdout comparison),
#'   \code{\link{ims}} for exact match detection (equivalent to hitting rate at threshold 0),
#'   \code{\link{nndr}} for nearest neighbor distance ratio
#'
#' @references
#' Platzer, M. & Reutterer, T. (2021). Holdout-Based Empirical Assessment of
#' Mixed-Type Synthetic Data. \emph{Frontiers in Big Data}, 4, 679939.
#' \doi{10.3389/fdata.2021.679939}
#'
#' @author Matthias Templ
#' @family distance-risk
#' @export
#' @importFrom VIM gowerD
#' @importFrom stats complete.cases quantile sd median
#' @importFrom graphics hist abline legend par
#' @examples
#' # Create example data
#' set.seed(123)
#' n <- 200
#' X <- data.frame(
#'   age = rnorm(n, 40, 10),
#'   income = rnorm(n, 50000, 15000),
#'   gender = sample(c("M", "F"), n, replace = TRUE),
#'   region = sample(c("N", "S", "E", "W"), n, replace = TRUE)
#' )
#'
#' # Good synthetic data (random, no memorization)
#' Y_good <- data.frame(
#'   age = rnorm(n, 40, 10),
#'   income = rnorm(n, 50000, 15000),
#'   gender = sample(c("M", "F"), n, replace = TRUE),
#'   region = sample(c("N", "S", "E", "W"), n, replace = TRUE)
#' )
#'
#' result <- hitting_rate(X, Y_good)
#' print(result)
#' summary(result)
#'
#' \donttest{
#' # Memorized data (Y is copy of X) -- should show high rate
#' Y_bad <- X[sample(nrow(X), n, replace = TRUE), ]
#' result_bad <- hitting_rate(X, Y_bad)
#' print(result_bad)
#'
#' # Sweep thresholds visually
#' plot(result, which = 2)
#' }
hitting_rate <- function(X, ...) {
  UseMethod("hitting_rate")
}

#' @rdname hitting_rate
#' @export
hitting_rate.synth_pair <- function(X, ...) {
  hitting_rate.default(
    X = X$original,
    Y = X$synthetic,
    vars = X$vars,
    ...
  )
}

#' @rdname hitting_rate
#' @export
hitting_rate.default <- function(X, Y,
                                  threshold = 0.05,
                                  vars = NULL,
                                  method = c("gower", "euclidean"),
                                  na.rm = TRUE,
                                  ...) {

  # Input validation
  if (!is.data.frame(X)) stop("X must be a data frame.")
  if (!is.data.frame(Y)) stop("Y must be a data frame.")

  if (!is.numeric(threshold) || length(threshold) != 1 || threshold < 0) {
    stop("threshold must be a single non-negative number.")
  }

  method <- match.arg(method)

  # Determine variables to use
  if (is.null(vars)) {
    vars <- intersect(names(X), names(Y))
  }

  if (length(vars) == 0) {
    stop("No common variables found between datasets.")
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

  # Check variable types match
  for (var in vars) {
    if (!identical(class(X[[var]]), class(Y[[var]]))) {
      stop(paste("Variable", var, "has different class in X and Y."))
    }
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

  # .normalize_minmax is defined in R/utils_internal.R

  # Compute minimum distance from each synthetic record to closest original
  if (method == "gower") {
    # Gower distance (handles mixed types)
    dist_mat <- VIM::gowerD(Y, X)
    min_distances <- apply(dist_mat, 1, min, na.rm = TRUE)

  } else if (method == "euclidean") {
    # Check all variables are numeric
    if (!all(sapply(Y, is.numeric))) {
      stop("method='euclidean' requires all variables to be numeric. Use method='gower' for mixed types.")
    }

    # Combine all data for consistent normalization
    all_data <- rbind(X, Y)
    all_data_norm <- as.data.frame(lapply(all_data, .normalize_minmax))

    X_norm <- all_data_norm[seq_len(n_original), , drop = FALSE]
    Y_norm <- all_data_norm[(n_original + 1):nrow(all_data_norm), , drop = FALSE]

    # Compute distances
    min_distances <- numeric(n_synthetic)
    for (i in seq_len(n_synthetic)) {
      diffs <- sweep(as.matrix(X_norm), 2, as.numeric(Y_norm[i, ]))
      dists <- sqrt(rowSums(diffs^2))
      min_distances[i] <- min(dists)
    }
  }

  # Determine hits
  hits <- min_distances <= threshold
  n_hits <- sum(hits)
  rate <- n_hits / n_synthetic

  # Rate at zero (equivalent to IMS)
  rate_at_zero <- sum(min_distances == 0) / n_synthetic

  # Privacy check
  privacy_pass <- rate <= 0.1

  results <- list(
    rate = rate,
    min_distances = min_distances,
    hits = hits,
    n_hits = n_hits,
    threshold = threshold,
    privacy_pass = privacy_pass,
    n_original = n_original,
    n_synthetic = n_synthetic,
    method = method,
    vars = vars,
    rate_at_zero = rate_at_zero
  )

  class(results) <- "hitting_rate"
  return(results)
}

# --- S3 methods ---

#' Print method for hitting_rate objects
#'
#' @param x an object of class "hitting_rate"
#' @param ... additional arguments (ignored)
#' @export
print.hitting_rate <- function(x, ...) {
  cat("Hitting Rate Privacy Metric\n")
  cat("===========================\n")
  cat("Method:", x$method, "\n")
  cat("Variables used:", length(x$vars), "\n")
  cat("Threshold:", x$threshold, "\n\n")

  cat("Dataset Sizes:\n")
  cat("  Original records:", x$n_original, "\n")
  cat("  Synthetic records:", x$n_synthetic, "\n\n")

  cat("Hitting Rate Results:\n")
  cat("  Hits (min distance <= threshold):", x$n_hits,
      sprintf("(%.1f%%)", 100 * x$rate), "\n")
  cat("  Hitting rate:", round(x$rate, 4), "\n")
  cat("  Rate at zero (exact matches):", round(x$rate_at_zero, 4),
      sprintf("(%.1f%%)", 100 * x$rate_at_zero), "\n\n")

  cat("Min Distance Distribution:\n")
  cat("  Mean:", round(mean(x$min_distances), 4), "\n")
  cat("  Median:", round(median(x$min_distances), 4), "\n")
  cat("  Min:", round(min(x$min_distances), 4),
      "| Max:", round(max(x$min_distances), 4), "\n\n")

  cat("Privacy Assessment:")
  if (x$privacy_pass) {
    cat(" PASS\n")
    cat("  Hitting rate is within acceptable bounds (<= 0.1).\n")
  } else {
    cat(" WARNING\n")
    cat("  Elevated hitting rate detected (> 0.1).\n")
    cat("  Synthetic data may contain too many near copies of original records.\n")
  }

  invisible(x)
}

#' Summary method for hitting_rate objects
#'
#' @param object an object of class "hitting_rate"
#' @param ... additional arguments (ignored)
#' @export
summary.hitting_rate <- function(object, ...) {
  summ <- list(
    rate = object$rate,
    n_hits = object$n_hits,
    threshold = object$threshold,
    privacy_pass = object$privacy_pass,
    rate_at_zero = object$rate_at_zero,
    n_exact = sum(object$min_distances == 0),
    mean_min_distance = mean(object$min_distances, na.rm = TRUE),
    median_min_distance = median(object$min_distances, na.rm = TRUE),
    sd_min_distance = sd(object$min_distances, na.rm = TRUE),
    quantiles = quantile(object$min_distances,
                          probs = c(0, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 1),
                          na.rm = TRUE),
    rates_at_thresholds = .compute_rates_at_thresholds(object$min_distances),
    n_original = object$n_original,
    n_synthetic = object$n_synthetic,
    method = object$method,
    vars = object$vars
  )

  class(summ) <- "summary.hitting_rate"
  return(summ)
}

#' Compute hitting rates at several standard thresholds
#'
#' @param min_distances numeric vector of minimum distances
#' @return named numeric vector of rates at standard thresholds
#' @keywords internal
.compute_rates_at_thresholds <- function(min_distances) {
  thresholds <- c(0, 0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.3)
  n <- length(min_distances)
  rates <- vapply(thresholds, function(tau) {
    sum(min_distances <= tau) / n
  }, numeric(1))
  names(rates) <- paste0("tau=", thresholds)
  rates
}

#' Print method for summary.hitting_rate objects
#'
#' @param x an object of class "summary.hitting_rate"
#' @param ... additional arguments (ignored)
#' @export
print.summary.hitting_rate <- function(x, ...) {
  cat("Summary: Hitting Rate Privacy Metric\n")
  cat("====================================\n")
  cat("Method:", x$method, "\n")
  cat("Variables:", length(x$vars), "\n\n")

  cat("Dataset Sizes:\n")
  cat("  Original:", x$n_original, "| Synthetic:", x$n_synthetic, "\n\n")

  cat("Key Metrics:\n")
  cat("  Hitting rate (at threshold", x$threshold, "):",
      round(x$rate, 4), sprintf("(%d hits)", x$n_hits), "\n")
  cat("  Exact matches (threshold = 0):", x$n_exact,
      sprintf("(%.4f)", x$rate_at_zero), "\n")
  cat("  Privacy:", ifelse(x$privacy_pass, "PASS", "WARNING"), "\n\n")

  cat("Min Distance Distribution:\n")
  cat("  Mean:", round(x$mean_min_distance, 4),
      "| Median:", round(x$median_min_distance, 4),
      "| SD:", round(x$sd_min_distance, 4), "\n")
  cat("  Quantiles:\n")
  print(round(x$quantiles, 4))
  cat("\n")

  cat("Hitting Rate at Standard Thresholds:\n")
  for (i in seq_along(x$rates_at_thresholds)) {
    cat(sprintf("  %-10s  %.4f\n",
                names(x$rates_at_thresholds)[i],
                x$rates_at_thresholds[i]))
  }

  invisible(x)
}

#' Plot method for hitting_rate objects
#'
#' @param x an object of class "hitting_rate"
#' @param y not used
#' @param ... additional arguments passed to plotting functions
#' @param which which plot(s) to show: 1 = histogram of min distances with
#'   threshold line, 2 = hitting rate vs threshold curve (sweeps thresholds
#'   from 0 to max distance)
#' @export
plot.hitting_rate <- function(x, y = NULL, ..., which = 1) {
  show <- rep(FALSE, 2)
  show[which] <- TRUE

  n_plots <- sum(show)
  if (n_plots > 1) {
    op <- par(mfrow = c(1, 2))
    on.exit(par(op))
  }

  if (show[1]) {
    # Histogram of min distances with threshold line
    max_dist <- max(x$min_distances, na.rm = TRUE)
    breaks <- seq(0, max_dist * 1.05, length.out = 31)

    hist(x$min_distances, breaks = breaks, col = rgb(0.3, 0.5, 0.8, 0.6),
         main = paste("Min Distance Distribution\nHitting rate:",
                      sprintf("%.1f%%", 100 * x$rate)),
         xlab = "Min distance to closest original record",
         ylab = "Frequency", ...)
    abline(v = x$threshold, col = "red", lwd = 2, lty = 2)
    legend("topright",
           legend = c(paste("threshold =", x$threshold),
                      paste("hits =", x$n_hits, "/", x$n_synthetic)),
           col = c("red", NA), lty = c(2, NA), lwd = c(2, NA),
           cex = 0.8)
  }

  if (show[2]) {
    # Hitting rate vs threshold curve
    max_dist <- max(x$min_distances, na.rm = TRUE)
    thresholds <- seq(0, max_dist, length.out = 200)
    n <- length(x$min_distances)

    rates <- vapply(thresholds, function(tau) {
      sum(x$min_distances <= tau) / n
    }, numeric(1))

    plot(thresholds, rates, type = "l", lwd = 2,
         main = "Hitting Rate vs Threshold",
         xlab = "Threshold",
         ylab = "Hitting Rate",
         ylim = c(0, 1), ...)
    abline(v = x$threshold, col = "red", lwd = 1.5, lty = 2)
    abline(h = x$rate, col = "red", lwd = 1, lty = 3)
    abline(h = 0.1, col = "grey50", lwd = 1, lty = 2)
    points(x$threshold, x$rate, pch = 16, col = "red", cex = 1.5)
    legend("bottomright",
           legend = c(paste("current threshold =", x$threshold),
                      paste("rate =", round(x$rate, 3)),
                      "privacy threshold (0.1)"),
           col = c("red", "red", "grey50"),
           lty = c(2, 3, 2), lwd = c(1.5, 1, 1),
           pch = c(NA, NA, NA),
           cex = 0.8)
  }
}
