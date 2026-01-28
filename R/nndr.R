#' Nearest Neighbor Distance Ratio (NNDR)
#'
#' Computes the Nearest Neighbor Distance Ratio privacy metric for synthetic data.
#' NNDR detects potential memorization by comparing the distance to the closest
#' record versus the second closest record. Low NNDR values indicate suspicious
#' proximity to specific training records.
#'
#' @param X data frame of original/training data
#' @param Y data frame of synthetic data
#' @param holdout data frame of holdout data (optional). If NULL, a holdout set
#'   is automatically created by splitting X.
#' @param holdout_fraction numeric, fraction of X to use as holdout if holdout
#'   is NULL (default: 0.5)
#' @param vars character vector of variable names to use for distance calculation.
#'   If NULL (default), all common variables are used.
#' @param method character, distance method: "gower" (default, handles mixed types)
#'   or "euclidean" (numerical variables only)
#' @param na.rm logical, remove records with NA values (default: TRUE)
#' @param seed integer, random seed for holdout sampling (default: NULL)
#' @param progress logical, show progress bar for long computations (default: FALSE)
#' @param ... additional arguments passed to methods (currently unused)
#'
#' @return An object of class "nndr" containing:
#' \itemize{
#'   \item nndr_train: NNDR values for synthetic vs training data
#'   \item nndr_holdout: NNDR values for synthetic vs holdout data (reference)
#'   \item nndr_ratio: ratio of 5th percentile NNDR (train/holdout)
#'   \item mean_nndr_train, mean_nndr_holdout: mean NNDR values
#'   \item n_suspicious: count of very low NNDR values (< 0.1)
#'   \item privacy_pass: logical, TRUE if NNDR distribution is acceptable
#' }
#'
#' @details
#' The Nearest Neighbor Distance Ratio (NNDR) is defined for each synthetic
#' record as:
#' \deqn{NNDR = \frac{d_1}{d_2}}
#' where \eqn{d_1} is the distance to the nearest neighbor and \eqn{d_2} is the
#' distance to the second nearest neighbor.
#'
#' Interpretation:
#' \itemize{
#'   \item \strong{NNDR close to 1}: The two nearest neighbors are similarly distant,
#'     suggesting the synthetic record is not unusually close to any specific
#'     training record (good privacy).
#'   \item \strong{NNDR close to 0}: The nearest neighbor is much closer than the
#'     second nearest, suggesting potential memorization or copying of a specific
#'     training record (privacy concern).
#' }
#'
#' The metric compares NNDR distributions for synthetic-to-training versus
#' synthetic-to-holdout. If the training NNDR distribution has significantly
#' more low values, this indicates privacy leakage.
#'
#' @section Holdout splitting (important):
#' Like \code{\link{dcr}}, NNDR internally splits the original data when no
#' external holdout is provided. With default \code{holdout_fraction = 0.5},
#' only half of the original data is used for comparison. This reduces
#' effective sample size and may miss disclosure if the synthetic data
#' memorized records that ended up in the holdout. Providing a separate
#' holdout from the synthesis process is recommended when available.
#' See \code{\link{dcr}} for detailed guidance.
#'
#' @seealso \code{\link{dcr}} for distance to closest record,
#'   \code{\link{ims}} for exact match detection
#'
#' @references
#' MOSTLY AI (2024). Evaluate generator quality.
#' \url{https://docs.mostly.ai/generators/evaluate-quality}
#'
#' Lowe, D.G. (2004). Distinctive Image Features from Scale-Invariant Keypoints.
#' \emph{International Journal of Computer Vision}, 60(2), 91-110.
#'
#' @author Matthias Templ
#' @export
#' @importFrom VIM gowerD
#' @importFrom stats complete.cases quantile sd
#' @importFrom graphics hist abline legend boxplot par
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @examples
#' # Create example data
#' set.seed(123)
#' n <- 300
#' X <- data.frame(
#'   age = rnorm(n, 40, 10),
#'   income = rnorm(n, 50000, 15000),
#'   gender = sample(c("M", "F"), n, replace = TRUE),
#'   region = sample(c("N", "S", "E", "W"), n, replace = TRUE)
#' )
#'
#' # Good synthetic data
#' Y_good <- data.frame(
#'   age = rnorm(200, 40, 10),
#'   income = rnorm(200, 50000, 15000),
#'   gender = sample(c("M", "F"), 200, replace = TRUE),
#'   region = sample(c("N", "S", "E", "W"), 200, replace = TRUE)
#' )
#'
#' result <- nndr(X, Y_good, seed = 42)
#' print(result)
#' summary(result)
#' plot(result)
nndr <- function(X, ...) {
  UseMethod("nndr")
}

#' @rdname nndr
#' @export
nndr.synth_pair <- function(X, ...) {
  nndr.default(
    X = X$original,
    Y = X$synthetic,
    holdout = X$holdout,
    vars = X$vars,
    ...
  )
}

#' @rdname nndr
#' @export
nndr.default <- function(X, Y,
                         holdout = NULL,
                         holdout_fraction = 0.5,
                         vars = NULL,
                         method = c("gower", "euclidean"),
                         na.rm = TRUE,
                         seed = NULL,
                         progress = FALSE,
                         ...) {

  # Input validation
  if (!is.data.frame(X)) stop("X must be a data frame.")
  if (!is.data.frame(Y)) stop("Y must be a data frame.")
  if (!is.null(holdout) && !is.data.frame(holdout)) {
    stop("holdout must be a data frame or NULL.")
  }

  method <- match.arg(method)

  # Determine variables to use
  if (is.null(vars)) {
    vars <- intersect(names(X), names(Y))
    if (!is.null(holdout)) {
      vars <- intersect(vars, names(holdout))
    }
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

  # Create or validate holdout
  if (is.null(holdout)) {
    if (!is.null(seed)) set.seed(seed)
    n_holdout <- max(2, floor(nrow(X) * holdout_fraction))  # Need at least 2 for NNDR
    holdout_idx <- sample(nrow(X), n_holdout)
    holdout <- X[holdout_idx, , drop = FALSE]
    train <- X[-holdout_idx, , drop = FALSE]
  } else {
    holdout <- holdout[, vars, drop = FALSE]
    if (na.rm) {
      complete_H <- complete.cases(holdout)
      holdout <- holdout[complete_H, , drop = FALSE]
    }
    if (nrow(holdout) < 2) stop("Holdout must have at least 2 records for NNDR.")
    train <- X
  }

  if (nrow(train) < 2) stop("Training data must have at least 2 records for NNDR.")

  n_synthetic <- nrow(Y)
  n_train <- nrow(train)
  n_holdout <- nrow(holdout)

  # Helper function to compute NNDR
  compute_nndr <- function(dist_matrix) {
    # For each row (synthetic record), find distances to 1st and 2nd nearest
    nndr_vals <- numeric(nrow(dist_matrix))

    for (i in seq_len(nrow(dist_matrix))) {
      dists <- sort(dist_matrix[i, ], na.last = TRUE)
      d1 <- dists[1]  # Distance to nearest
      d2 <- dists[2]  # Distance to 2nd nearest

      if (is.na(d2) || d2 == 0) {
        nndr_vals[i] <- NA
      } else {
        nndr_vals[i] <- d1 / d2
      }
    }
    return(nndr_vals)
  }

  # Compute distances and NNDR
  if (method == "gower") {
    dist_to_train <- VIM::gowerD(Y, train)
    dist_to_holdout <- VIM::gowerD(Y, holdout)

    nndr_train <- compute_nndr(dist_to_train)
    nndr_holdout <- compute_nndr(dist_to_holdout)

  } else if (method == "euclidean") {
    if (!all(sapply(Y, is.numeric))) {
      stop("method='euclidean' requires all variables to be numeric. Use method='gower' for mixed types.")
    }

    # Normalize for Euclidean distance
    normalize <- function(x) {
      rng <- range(x, na.rm = TRUE)
      if (rng[2] - rng[1] == 0) return(rep(0, length(x)))
      (x - rng[1]) / (rng[2] - rng[1])
    }

    all_data <- rbind(train, holdout, Y)
    all_data_norm <- as.data.frame(lapply(all_data, normalize))

    train_norm <- all_data_norm[1:n_train, , drop = FALSE]
    holdout_norm <- all_data_norm[(n_train + 1):(n_train + n_holdout), , drop = FALSE]
    Y_norm <- all_data_norm[(n_train + n_holdout + 1):nrow(all_data_norm), , drop = FALSE]

    # Compute distance matrices
    dist_to_train <- matrix(0, n_synthetic, n_train)
    dist_to_holdout <- matrix(0, n_synthetic, n_holdout)

    # Setup progress bar if requested
    if (progress) {
      pb <- txtProgressBar(min = 0, max = n_synthetic, style = 3)
    }

    for (i in seq_len(n_synthetic)) {
      diffs_train <- sweep(as.matrix(train_norm), 2, as.numeric(Y_norm[i, ]))
      dist_to_train[i, ] <- sqrt(rowSums(diffs_train^2))

      diffs_holdout <- sweep(as.matrix(holdout_norm), 2, as.numeric(Y_norm[i, ]))
      dist_to_holdout[i, ] <- sqrt(rowSums(diffs_holdout^2))

      if (progress) setTxtProgressBar(pb, i)
    }

    if (progress) close(pb)

    nndr_train <- compute_nndr(dist_to_train)
    nndr_holdout <- compute_nndr(dist_to_holdout)
  }

  # Remove NAs for statistics
  nndr_train_clean <- nndr_train[!is.na(nndr_train)]
  nndr_holdout_clean <- nndr_holdout[!is.na(nndr_holdout)]

  # Compute summary statistics
  mean_nndr_train <- mean(nndr_train_clean)
  mean_nndr_holdout <- mean(nndr_holdout_clean)

  # Compare lower percentiles (where memorization would show)
  p5_train <- quantile(nndr_train_clean, 0.05)
  p5_holdout <- quantile(nndr_holdout_clean, 0.05)
  p10_train <- quantile(nndr_train_clean, 0.10)
  p10_holdout <- quantile(nndr_holdout_clean, 0.10)

  # NNDR ratio: compare 5th percentiles
  nndr_ratio <- as.numeric(p5_train / p5_holdout)

  # Count suspicious records (very low NNDR)
  n_suspicious <- sum(nndr_train_clean < 0.1)
  pct_suspicious <- 100 * n_suspicious / length(nndr_train_clean)

  # Privacy check: NNDR ratio should be close to 1, and few suspicious records
  # If train NNDR is much lower than holdout NNDR, it indicates memorization
  privacy_pass <- nndr_ratio >= 0.8 && pct_suspicious <= 5

  results <- list(
    nndr_train = nndr_train,
    nndr_holdout = nndr_holdout,
    nndr_ratio = nndr_ratio,
    mean_nndr_train = mean_nndr_train,
    mean_nndr_holdout = mean_nndr_holdout,
    median_nndr_train = median(nndr_train_clean),
    median_nndr_holdout = median(nndr_holdout_clean),
    p5_train = as.numeric(p5_train),
    p5_holdout = as.numeric(p5_holdout),
    p10_train = as.numeric(p10_train),
    p10_holdout = as.numeric(p10_holdout),
    n_suspicious = n_suspicious,
    pct_suspicious = pct_suspicious,
    privacy_pass = privacy_pass,
    n_synthetic = n_synthetic,
    n_train = n_train,
    n_holdout = n_holdout,
    method = method,
    vars = vars
  )

  class(results) <- "nndr"
  return(results)
}

#' Print method for nndr objects
#'
#' @param x an object of class "nndr"
#' @param ... additional arguments (ignored)
#' @export
print.nndr <- function(x, ...) {
 cat("Nearest Neighbor Distance Ratio (NNDR) Privacy Metric\n")
  cat("=====================================================\n")
  cat("Method:", x$method, "\n")
  cat("Variables used:", length(x$vars), "\n\n")

  cat("Dataset Sizes:\n")
  cat("  Training records:", x$n_train, "\n")
  cat("  Holdout records:", x$n_holdout, "\n")
  cat("  Synthetic records:", x$n_synthetic, "\n\n")

  cat("NNDR Results:\n")
  cat("  Mean NNDR (to training):", round(x$mean_nndr_train, 4), "\n")
  cat("  Mean NNDR (to holdout):", round(x$mean_nndr_holdout, 4), "\n")
  cat("  5th percentile (training):", round(x$p5_train, 4), "\n")
  cat("  5th percentile (holdout):", round(x$p5_holdout, 4), "\n")
  cat("  NNDR ratio (p5 train/holdout):", round(x$nndr_ratio, 4),
      ifelse(x$nndr_ratio < 0.8, " (privacy concern)", ""), "\n")
  cat("  Suspicious records (NNDR < 0.1):", x$n_suspicious,
      sprintf("(%.1f%%)", x$pct_suspicious),
      ifelse(x$pct_suspicious > 5, " (privacy concern)", ""), "\n\n")

  cat("Privacy Assessment:")
  if (x$privacy_pass) {
    cat(" PASS\n")
    cat("  No evidence of training data memorization.\n")
  } else {
    cat(" WARNING\n")
    cat("  Potential memorization detected - synthetic records are\n")
    cat("  suspiciously close to specific training records.\n")
  }

  invisible(x)
}

#' Summary method for nndr objects
#'
#' @param object an object of class "nndr"
#' @param ... additional arguments (ignored)
#' @export
summary.nndr <- function(object, ...) {
  nndr_train <- object$nndr_train[!is.na(object$nndr_train)]
  nndr_holdout <- object$nndr_holdout[!is.na(object$nndr_holdout)]

  summ <- list(
    nndr_ratio = object$nndr_ratio,
    privacy_pass = object$privacy_pass,
    mean_nndr_train = object$mean_nndr_train,
    mean_nndr_holdout = object$mean_nndr_holdout,
    median_nndr_train = object$median_nndr_train,
    median_nndr_holdout = object$median_nndr_holdout,
    sd_nndr_train = sd(nndr_train),
    sd_nndr_holdout = sd(nndr_holdout),
    quantiles_train = quantile(nndr_train,
                               probs = c(0, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 1)),
    quantiles_holdout = quantile(nndr_holdout,
                                 probs = c(0, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 1)),
    n_suspicious = object$n_suspicious,
    pct_suspicious = object$pct_suspicious,
    n_very_suspicious = sum(nndr_train < 0.05),
    n_synthetic = object$n_synthetic,
    n_train = object$n_train,
    n_holdout = object$n_holdout,
    method = object$method,
    vars = object$vars
  )

  class(summ) <- "summary.nndr"
  return(summ)
}

#' Print method for summary.nndr objects
#'
#' @param x an object of class "summary.nndr"
#' @param ... additional arguments (ignored)
#' @export
print.summary.nndr <- function(x, ...) {
  cat("Summary: Nearest Neighbor Distance Ratio (NNDR)\n")
  cat("================================================\n")
  cat("Method:", x$method, "\n")
  cat("Variables:", length(x$vars), "\n\n")

  cat("Dataset Sizes:\n")
  cat("  Training:", x$n_train, "| Holdout:", x$n_holdout,
      "| Synthetic:", x$n_synthetic, "\n\n")

  cat("Key Metrics:\n")
  cat("  NNDR ratio (5th percentile):", round(x$nndr_ratio, 4), "(ideal: ~1.0)\n")
  cat("  Privacy:", ifelse(x$privacy_pass, "PASS", "WARNING"), "\n\n")

  cat("NNDR to Training:\n")
  cat("  Mean:", round(x$mean_nndr_train, 4),
      "| Median:", round(x$median_nndr_train, 4),
      "| SD:", round(x$sd_nndr_train, 4), "\n")
  cat("  Quantiles:\n")
  print(round(x$quantiles_train, 4))
  cat("\n")

  cat("NNDR to Holdout (reference):\n")
  cat("  Mean:", round(x$mean_nndr_holdout, 4),
      "| Median:", round(x$median_nndr_holdout, 4),
      "| SD:", round(x$sd_nndr_holdout, 4), "\n")
  cat("  Quantiles:\n")
  print(round(x$quantiles_holdout, 4))
  cat("\n")

  cat("Suspicious Records:\n")
  cat("  NNDR < 0.1:", x$n_suspicious, sprintf("(%.1f%%)", x$pct_suspicious), "\n")
  cat("  NNDR < 0.05:", x$n_very_suspicious, "(very suspicious)\n")

  invisible(x)
}

#' Plot method for nndr objects
#'
#' @param x an object of class "nndr"
#' @param y not used
#' @param ... additional arguments passed to plotting functions
#' @param which which plot(s) to show: 1 = overlaid histograms, 2 = boxplot
#'   comparison, 3 = cumulative distribution comparison
#' @export
plot.nndr <- function(x, y = NULL, ..., which = 1) {
  show <- rep(FALSE, 3)
  show[which] <- TRUE

  nndr_train <- x$nndr_train[!is.na(x$nndr_train)]
  nndr_holdout <- x$nndr_holdout[!is.na(x$nndr_holdout)]

  n_plots <- sum(show)
  if (n_plots > 1) {
    if (n_plots == 2) {
      op <- par(mfrow = c(1, 2))
    } else {
      op <- par(mfrow = c(1, 3))
    }
    on.exit(par(op))
  }

  if (show[1]) {
    # Overlaid histograms
    breaks <- seq(0, 1, by = 0.05)

    hist(nndr_train, breaks = breaks, col = rgb(1, 0, 0, 0.5),
         main = "NNDR Distribution",
         xlab = "Nearest Neighbor Distance Ratio",
         xlim = c(0, 1), ...)
    hist(nndr_holdout, breaks = breaks, col = rgb(0, 0, 1, 0.5), add = TRUE)
    abline(v = 0.1, col = "darkred", lwd = 2, lty = 2)
    legend("topleft",
           legend = c("To training", "To holdout", "Suspicious threshold"),
           fill = c(rgb(1, 0, 0, 0.5), rgb(0, 0, 1, 0.5), NA),
           border = c("black", "black", NA),
           lty = c(NA, NA, 2),
           lwd = c(NA, NA, 2),
           col = c(NA, NA, "darkred"),
           cex = 0.8)
  }

  if (show[2]) {
    # Boxplot comparison
    boxplot(list("To Training" = nndr_train, "To Holdout" = nndr_holdout),
            main = "NNDR Comparison",
            ylab = "Nearest Neighbor Distance Ratio",
            col = c("coral", "steelblue"),
            ylim = c(0, 1), ...)
    abline(h = 0.1, col = "darkred", lwd = 2, lty = 2)
  }

  if (show[3]) {
    # Cumulative distribution
    plot(ecdf(nndr_train), col = "red", main = "Cumulative Distribution of NNDR",
         xlab = "NNDR", ylab = "Cumulative Proportion",
         xlim = c(0, 1), ...)
    plot(ecdf(nndr_holdout), col = "blue", add = TRUE)
    abline(v = 0.1, col = "darkred", lwd = 2, lty = 2)
    legend("bottomright",
           legend = c("To training", "To holdout"),
           col = c("red", "blue"), lwd = 2, cex = 0.8)
  }
}
