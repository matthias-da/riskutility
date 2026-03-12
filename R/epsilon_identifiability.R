#' Epsilon Identifiability
#'
#' Computes the Epsilon Identifiability risk metric for synthetic data, a
#' distance-entropy hybrid from the SynthEval framework. For each synthetic
#' record, the metric finds its weighted Gower distance to the closest original
#' record, where variable weights are the inverse of Shannon entropy (rare
#' attributes are penalized more heavily). Records with minimum distance below
#' the epsilon threshold are flagged as identifiable.
#'
#' @param X data frame of original/training data
#' @param Y data frame of synthetic data
#' @param epsilon numeric, distance threshold below which a synthetic record is
#'   considered identifiable (default: 0.05). Must be between 0 and 1.
#' @param vars character vector of variable names to use. If NULL (default),
#'   all common variables between X and Y are used.
#' @param na.rm logical, remove records with NA values (default: TRUE)
#' @param seed integer, random seed for reproducibility (default: NULL).
#'   Currently unused but reserved for future extensions.
#' @param n_bins integer, number of bins for discretizing numeric variables
#'   before computing entropy (default: 20)
#' @param ... additional arguments passed to methods (currently unused)
#'
#' @return An object of class "epsilon_identifiability" containing:
#' \itemize{
#'   \item distances: numeric vector of minimum weighted distances per synthetic record
#'   \item flagged: logical vector, TRUE where distance < epsilon
#'   \item identifiability_rate: fraction of synthetic records flagged
#'   \item entropy_weights: named numeric vector of per-variable weights (normalized, sum to 1)
#'   \item entropies: named numeric vector of per-variable Shannon entropies
#'   \item epsilon: threshold used
#'   \item n_flagged: count of flagged records
#'   \item privacy_pass: logical, TRUE if identifiability_rate <= 0.1
#'   \item n_original: number of original records
#'   \item n_synthetic: number of synthetic records
#'   \item vars: variables used
#' }
#'
#' @details
#' The Epsilon Identifiability metric combines distance-based and
#' information-theoretic approaches:
#'
#' \enumerate{
#'   \item Compute Shannon entropy \eqn{H_v} for each variable \eqn{v} in the
#'     original data. For numeric variables, values are discretized into
#'     \code{n_bins} equal-width bins. For categorical variables, entropy is
#'     computed directly from category frequencies.
#'   \item Compute inverse-entropy weights \eqn{w_v = 1 / H_v}. Variables with
#'     zero entropy (constant) receive zero weight. Weights are normalized to
#'     sum to 1.
#'   \item Compute weighted Gower distance from each synthetic record to all
#'     original records using \code{VIM::gowerD} with the entropy-based weights.
#'   \item For each synthetic record, find the minimum weighted distance.
#'   \item Flag records where distance < epsilon.
#'   \item Compute identifiability rate = fraction of flagged records.
#' }
#'
#' The intuition is that variables with low entropy (few distinct values or
#' skewed distributions) carry more identifying power. A close match on a
#' rare-valued variable is more concerning than a close match on a
#' high-entropy variable.
#'
#' @section Choosing epsilon:
#' The default \code{epsilon = 0.05} is a conservative threshold. Smaller
#' values are more strict (fewer flagged records). The appropriate threshold
#' depends on the number of variables and the data domain. Use
#' \code{plot(result, which = 1)} to visualize the distance distribution and
#' assess whether the threshold is appropriate.
#'
#' @seealso \code{\link{dcr}} for distance to closest record,
#'   \code{\link{nndr}} for nearest neighbor distance ratio,
#'   \code{\link{ims}} for identical match share
#'
#' @references
#' Lautrup, A., Hyrup, T., Zimek, A. & Blockeel, H. (2025).
#' SynthEval: A Framework for Detailed Utility and Privacy Evaluation of
#' Tabular Synthetic Data.
#' \emph{Data Mining and Knowledge Discovery}, 39(1).
#' \doi{10.1007/s10618-024-01081-4}
#'
#' @author Matthias Templ
#' @family privacy-models
#' @export
#' @importFrom VIM gowerD
#' @importFrom stats complete.cases
#' @importFrom graphics hist abline legend barplot par
#' @examples
#' # Create example data
#' set.seed(123)
#' n <- 200
#' X <- data.frame(
#'   age = sample(20:70, n, replace = TRUE),
#'   income = rnorm(n, 50000, 15000),
#'   gender = sample(c("M", "F"), n, replace = TRUE),
#'   region = sample(c("N", "S", "E", "W"), n, replace = TRUE)
#' )
#'
#' # Good synthetic data (random, no memorization)
#' Y_good <- data.frame(
#'   age = sample(20:70, n, replace = TRUE),
#'   income = rnorm(n, 50000, 15000),
#'   gender = sample(c("M", "F"), n, replace = TRUE),
#'   region = sample(c("N", "S", "E", "W"), n, replace = TRUE)
#' )
#'
#' result <- epsilon_identifiability(X, Y_good)
#' print(result)
#' summary(result)
#'
#' \donttest{
#' # Memorized data (Y is a copy of X) - should show high risk
#' Y_copy <- X[sample(nrow(X), n, replace = TRUE), ]
#' result_bad <- epsilon_identifiability(X, Y_copy)
#' print(result_bad)
#' }
epsilon_identifiability <- function(X, ...) {
  UseMethod("epsilon_identifiability")
}

#' @rdname epsilon_identifiability
#' @export
epsilon_identifiability.synth_pair <- function(X, ...) {
  epsilon_identifiability.default(
    X = X$original,
    Y = X$synthetic,
    vars = X$vars,
    ...
  )
}

#' @rdname epsilon_identifiability
#' @export
epsilon_identifiability.default <- function(X, Y,
                                             epsilon = 0.05,
                                             vars = NULL,
                                             na.rm = TRUE,
                                             seed = NULL,
                                             n_bins = 20L,
                                             ...) {

  # Input validation
  if (!is.data.frame(X)) stop("X must be a data frame.")
  if (!is.data.frame(Y)) stop("Y must be a data frame.")

  if (!is.numeric(epsilon) || length(epsilon) != 1 || epsilon <= 0 || epsilon >= 1) {
    stop("epsilon must be a single number between 0 and 1 (exclusive).")
  }

  if (!is.numeric(n_bins) || length(n_bins) != 1 || n_bins < 2) {
    stop("n_bins must be an integer >= 2.")
  }
  n_bins <- as.integer(n_bins)

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

  # Set seed if provided
  if (!is.null(seed)) set.seed(seed)

  # --- Step 1: Compute Shannon entropy for each variable in original data ---
  entropies <- vapply(vars, function(v) {
    col <- X[[v]]
    .compute_entropy(col, n_bins = n_bins)
  }, numeric(1))
  names(entropies) <- vars

  # --- Step 2: Compute inverse-entropy weights ---
  raw_weights <- ifelse(entropies > 0, 1 / entropies, 0)

  # Normalize weights to sum to 1
  w_sum <- sum(raw_weights)
  if (w_sum == 0) {
    # All variables are constant; use equal weights
    entropy_weights <- rep(1 / length(vars), length(vars))
    warning("All variables have zero entropy (constant). Using equal weights.")
  } else {
    entropy_weights <- raw_weights / w_sum
  }
  names(entropy_weights) <- vars

  # --- Step 3: Compute weighted Gower distances ---
  # VIM::gowerD accepts a weights parameter
  dist_mat <- VIM::gowerD(Y, X, weights = entropy_weights)

  # --- Step 4: For each synthetic record, find the minimum distance ---
  distances <- apply(dist_mat, 1, min, na.rm = TRUE)

  # --- Step 5: Flag records where distance < epsilon ---
  flagged <- distances < epsilon
  n_flagged <- sum(flagged)

  # --- Step 6: Compute identifiability rate ---
  identifiability_rate <- n_flagged / nrow(Y)
  privacy_pass <- identifiability_rate <= 0.1

  results <- list(
    distances = distances,
    flagged = flagged,
    identifiability_rate = identifiability_rate,
    entropy_weights = entropy_weights,
    entropies = entropies,
    epsilon = epsilon,
    n_flagged = n_flagged,
    privacy_pass = privacy_pass,
    n_original = nrow(X),
    n_synthetic = nrow(Y),
    vars = vars
  )

  class(results) <- "epsilon_identifiability"
  return(results)
}

# --- Internal helper ---

#' Compute Shannon entropy for a variable
#'
#' For numeric variables, values are discretized into bins before computing
#' entropy. For categorical/factor variables, entropy is computed directly
#' from category frequencies.
#'
#' @param x a vector (numeric, character, or factor)
#' @param n_bins integer, number of bins for discretizing numeric variables
#' @return numeric scalar, the Shannon entropy (natural log)
#' @keywords internal
# Note: handles raw data vectors (discretization + entropy), unlike the
# probability-vector entropy functions in divergence.R/divergence2.R.
.compute_entropy <- function(x, n_bins = 20L) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(0)

  if (is.numeric(x)) {
    # Discretize numeric variables into equal-width bins
    rng <- range(x)
    if (rng[1] == rng[2]) return(0)  # Constant variable
    breaks <- seq(rng[1], rng[2], length.out = n_bins + 1)
    # Adjust endpoints slightly to include all values
    breaks[1] <- breaks[1] - 1e-10
    breaks[n_bins + 1] <- breaks[n_bins + 1] + 1e-10
    binned <- cut(x, breaks = breaks, include.lowest = TRUE)
    freq <- table(binned)
  } else {
    # Categorical or factor: compute directly from frequencies
    freq <- table(x)
  }

  # Remove zero-count bins
  freq <- freq[freq > 0]
  if (length(freq) <= 1) return(0)

  # Shannon entropy: H = -sum(p * log(p))
  p <- as.numeric(freq) / sum(freq)
  entropy <- -sum(p * log(p))
  return(entropy)
}

# --- S3 methods ---

#' Print method for epsilon_identifiability objects
#'
#' @param x an object of class "epsilon_identifiability"
#' @param ... additional arguments (ignored)
#' @export
print.epsilon_identifiability <- function(x, ...) {
  cat("Epsilon Identifiability Risk Assessment\n")
  cat("=======================================\n")
  cat("Variables used:", length(x$vars), "\n")
  cat("Epsilon threshold:", x$epsilon, "\n\n")

  cat("Dataset Sizes:\n")
  cat("  Original records:", x$n_original, "\n")
  cat("  Synthetic records:", x$n_synthetic, "\n\n")

  cat("Results:\n")
  cat("  Identifiability rate:", sprintf("%.4f", x$identifiability_rate),
      sprintf("(%.1f%%)", 100 * x$identifiability_rate), "\n")
  cat("  Flagged records:", x$n_flagged, "of", x$n_synthetic, "\n")
  cat("  Mean min distance:", round(mean(x$distances), 4), "\n")
  cat("  Median min distance:", round(median(x$distances), 4), "\n\n")

  cat("Entropy Weights (top 5):\n")
  sorted_w <- sort(x$entropy_weights, decreasing = TRUE)
  n_show <- min(5, length(sorted_w))
  for (i in seq_len(n_show)) {
    cat(sprintf("  %-20s %.4f (H = %.4f)\n",
                names(sorted_w)[i], sorted_w[i],
                x$entropies[names(sorted_w)[i]]))
  }
  if (length(sorted_w) > 5) {
    cat("  ...\n")
  }
  cat("\n")

  cat("Privacy Assessment:")
  if (x$privacy_pass) {
    cat(" PASS\n")
    cat("  Identifiability rate is within acceptable bounds (<= 0.1).\n")
  } else {
    cat(" WARNING\n")
    cat("  Elevated identifiability risk detected (> 0.1).\n")
    cat("  Synthetic records may be too close to original records.\n")
  }

  invisible(x)
}

#' Summary method for epsilon_identifiability objects
#'
#' @param object an object of class "epsilon_identifiability"
#' @param ... additional arguments (ignored)
#' @export
summary.epsilon_identifiability <- function(object, ...) {
  summ <- list(
    identifiability_rate = object$identifiability_rate,
    n_flagged = object$n_flagged,
    privacy_pass = object$privacy_pass,
    epsilon = object$epsilon,
    mean_distance = mean(object$distances),
    median_distance = median(object$distances),
    sd_distance = sd(object$distances),
    quantiles = quantile(object$distances,
                         probs = c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1),
                         na.rm = TRUE),
    entropy_weights = object$entropy_weights,
    entropies = object$entropies,
    n_original = object$n_original,
    n_synthetic = object$n_synthetic,
    vars = object$vars
  )

  class(summ) <- "summary.epsilon_identifiability"
  return(summ)
}

#' Print method for summary.epsilon_identifiability objects
#'
#' @param x an object of class "summary.epsilon_identifiability"
#' @param ... additional arguments (ignored)
#' @export
print.summary.epsilon_identifiability <- function(x, ...) {
  cat("Summary: Epsilon Identifiability Risk Assessment\n")
  cat("================================================\n")
  cat("Variables:", length(x$vars), "\n")
  cat("Epsilon threshold:", x$epsilon, "\n\n")

  cat("Dataset Sizes:\n")
  cat("  Original:", x$n_original, "| Synthetic:", x$n_synthetic, "\n\n")

  cat("Key Metrics:\n")
  cat("  Identifiability rate:", sprintf("%.4f", x$identifiability_rate),
      sprintf("(%.1f%%)", 100 * x$identifiability_rate), "\n")
  cat("  Flagged records:", x$n_flagged, "of", x$n_synthetic, "\n")
  cat("  Privacy:", ifelse(x$privacy_pass, "PASS", "WARNING"), "\n\n")

  cat("Distance Distribution:\n")
  cat("  Mean:", round(x$mean_distance, 4),
      "| Median:", round(x$median_distance, 4),
      "| SD:", round(x$sd_distance, 4), "\n")
  cat("  Quantiles:\n")
  print(round(x$quantiles, 4))
  cat("\n")

  cat("Entropy Weights per Variable:\n")
  weight_df <- data.frame(
    Variable = names(x$entropy_weights),
    Entropy = round(x$entropies[names(x$entropy_weights)], 4),
    Weight = round(x$entropy_weights, 4)
  )
  weight_df <- weight_df[order(-weight_df$Weight), ]
  rownames(weight_df) <- NULL
  print(weight_df, row.names = FALSE)

  invisible(x)
}

#' Plot method for epsilon_identifiability objects
#'
#' @param x an object of class "epsilon_identifiability"
#' @param y not used
#' @param ... additional arguments passed to plotting functions
#' @param which which plot(s) to show: 1 = histogram of minimum distances with
#'   epsilon threshold line, 2 = entropy weights barplot
#' @export
plot.epsilon_identifiability <- function(x, y = NULL, ..., which = 1) {
  show <- rep(FALSE, 2)
  show[which] <- TRUE

  n_plots <- sum(show)
  if (n_plots > 1) {
    op <- par(mfrow = c(1, 2))
    on.exit(par(op))
  }

  if (show[1]) {
    # Histogram of minimum distances with epsilon threshold
    hist(x$distances,
         breaks = 30,
         col = "steelblue",
         main = paste0("Minimum Weighted Distances\n(epsilon = ", x$epsilon, ")"),
         xlab = "Min Weighted Gower Distance to Original",
         ylab = "Frequency", ...)
    abline(v = x$epsilon, col = "red", lwd = 2, lty = 2)
    legend("topright",
           legend = c(paste("epsilon =", x$epsilon),
                      paste("Flagged:", x$n_flagged, "of", x$n_synthetic,
                            sprintf("(%.1f%%)", 100 * x$identifiability_rate))),
           col = c("red", NA),
           lty = c(2, NA), lwd = c(2, NA),
           cex = 0.8)
  }

  if (show[2]) {
    # Entropy weights barplot
    sorted_w <- sort(x$entropy_weights, decreasing = TRUE)
    barplot(sorted_w,
            main = "Entropy-Based Variable Weights",
            ylab = "Weight (1/H, normalized)",
            col = "coral",
            las = 2,
            cex.names = 0.8, ...)
  }
}
