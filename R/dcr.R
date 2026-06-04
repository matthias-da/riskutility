#' Distance to Closest Record (DCR)
#'
#' Computes the Distance to Closest Record privacy metric for synthetic data.
#' DCR compares the distances from synthetic records to their nearest neighbors
#' in the training data versus a holdout set to detect potential privacy leaks.
#'
#' @param X data frame of original/training data
#' @param Y data frame of synthetic data
#' @param holdout data frame of holdout data (optional). If NULL, a holdout set
#'   is automatically created by splitting X.
#' @param holdout_fraction numeric, fraction of X to use as holdout if holdout
#'   is NULL (default: 0.5). Larger holdouts give more reliable results,
#'   especially for smaller datasets. Use 0.1 only for very large datasets.
#' @param vars character vector of variable names to use for distance calculation.
#'   If NULL (default), all common variables between X, Y, and holdout are used.
#' @param method character, distance method: "gower" (default, handles mixed types)
#'   or "euclidean" (numerical variables only)
#' @param na.rm logical, remove records with NA values (default: TRUE)
#' @param seed integer, random seed for holdout sampling (default: NULL)
#' @param progress logical, show progress bar for long computations (default: FALSE)
#' @param null_test logical, perform permutation test comparing observed DCR share
#'   against a null distribution (default: TRUE). Permutes training/holdout
#'   assignment to estimate expected share under no memorization.
#' @param n_null integer, number of permutation samples for null distribution
#'   estimation (default: 100). Only used when null_test = TRUE.
#' @param ... additional arguments passed to methods (currently unused)
#'
#' @return An object of class "dcr" containing:
#' \itemize{
#'   \item dcr_train: distances from synthetic to closest training record
#'   \item dcr_holdout: distances from synthetic to closest holdout record
#'   \item dcr_ratio: ratio of mean distances (train/holdout), ideally ~1
#'   \item dcr_share: proportion of synthetic closer to training than holdout
#'   \item privacy_pass: logical, TRUE if dcr_share <= 0.55 and Wilcoxon p > 0.05
#'   \item wilcox_test: Wilcoxon test result comparing train vs holdout distances
#'   \item null_distribution: null distribution statistics (if null_test = TRUE)
#'   \item n_synthetic, n_train, n_holdout: dataset sizes
#'   \item method, vars: parameters used
#' }
#'
#' @details
#' Distance to Closest Record (DCR) is a privacy metric that detects whether
#' synthetic data has memorized or copied training records. The key insight is
#' that synthetic records should be equally close to training and holdout data
#' if no information leakage occurred.
#'
#' The metric works by:
#' \enumerate{
#'   \item Splitting original data into training and holdout (or using provided holdout)
#'   \item For each synthetic record, computing distance to nearest training neighbor
#'   \item For each synthetic record, computing distance to nearest holdout neighbor
#'   \item Comparing these distance distributions using statistical tests
#' }
#'
#' Interpretation:
#' \itemize{
#'   \item \strong{DCR ratio ~ 1}: Good privacy - synthetic equally distant from both
#'   \item \strong{DCR ratio < 1}: Privacy concern - synthetic closer to training
#'   \item \strong{DCR share ~ 0.5}: Ideal - 50% closer to training, 50% to holdout
#'   \item \strong{DCR share > 0.5}: Privacy concern - too many closer to training
#' }
#'
#' @section Important limitations (The DCR Delusion):
#' Yao et al. (2025) demonstrate that DCR and related distance-based metrics
#' can fail to detect privacy leakage. Key findings:
#' \itemize{
#'   \item \strong{False sense of security}: Datasets deemed "private" by DCR can
#'     still be vulnerable to Membership Inference Attacks (MIAs).
#'   \item \strong{Null distribution matters}: DCR values must be compared against
#'     a proper null distribution, not interpreted in absolute terms.
#'   \item \strong{Not sufficient alone}: DCR should be used alongside other
#'     privacy metrics and, ideally, actual MIA evaluations.
#'   \item \strong{Metric gaming}: Generators can produce low DCR while still
#'     leaking membership information through other channels.
#' }
#'
#' This implementation includes statistical tests and null distribution comparison
#' to provide more rigorous privacy assessment, but users should be aware that
#' passing DCR tests does not guarantee privacy protection.
#'
#' @section Holdout splitting (important):
#' When no external holdout is provided, DCR internally splits the original data X
#' into training and holdout portions. This has important implications:
#' \itemize{
#'   \item \strong{Reduced training set}: With default \code{holdout_fraction = 0.5},
#'     only half of X is used for distance comparison. This may miss some
#'     disclosure risks if the synthetic data memorized records in the holdout.
#'   \item \strong{Variability}: Results depend on the random split (use \code{seed}
#'     for reproducibility).
#'   \item \strong{Best practice}: Provide a separate holdout set from the synthesis
#'     process if available. Many synthesis frameworks support train/test splits.
#' }
#'
#' For small datasets (< 500 records), consider using a smaller holdout fraction
#' (e.g., 0.2-0.3) to maintain sufficient training data for comparison.
#'
#' @seealso \code{\link{nndr}} for nearest neighbor distance ratio,
#'   \code{\link{ims}} for exact match detection,
#'   \code{\link{nnaa}} for nearest-neighbor adversarial accuracy
#'
#' @references
#' Platzer, M. & Reutterer, T. (2021). Holdout-Based Empirical Assessment of
#' Mixed-Type Synthetic Data. \emph{Frontiers in Big Data}, 4, 679939.
#' \doi{10.3389/fdata.2021.679939}
#'
#' Park, N., et al. (2018). Data Synthesis based on Generative Adversarial
#' Networks. \emph{Proceedings of the VLDB Endowment}, 11(10), 1071--1083.
#' \doi{10.14778/3231751.3231757}
#'
#' Yao, Z., Krco, N., Ganev, G. & de Montjoye, Y.-A. (2025).
#' The DCR Delusion: Measuring the Privacy Risk of Synthetic Data.
#' \url{https://arxiv.org/abs/2505.01524}
#'
#' Zhao, Z., et al. (2021). CTAB-GAN: Effective Table Data Synthesizing.
#' \emph{Asian Conference on Machine Learning}.
#'
#' @author Matthias Templ
#' @family distance-risk
#' @export
#' @importFrom VIM gowerD
#' @importFrom stats complete.cases dist quantile sd wilcox.test median
#' @importFrom graphics hist abline legend boxplot par
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @examples
#' # Create example data
#' set.seed(123)
#' X <- data.frame(
#'   age = rnorm(200, 40, 10),
#'   income = rnorm(200, 50000, 15000),
#'   gender = sample(c("M", "F"), 200, replace = TRUE),
#'   region = sample(c("N", "S", "E", "W"), 200, replace = TRUE)
#' )
#'
#' # Good synthetic data (random, no memorization)
#' Y_good <- data.frame(
#'   age = rnorm(200, 40, 10),
#'   income = rnorm(200, 50000, 15000),
#'   gender = sample(c("M", "F"), 200, replace = TRUE),
#'   region = sample(c("N", "S", "E", "W"), 200, replace = TRUE)
#' )
#'
#' # Compute DCR
#' result <- dcr(X, Y_good, seed = 42)
#' print(result)
#' summary(result)
dcr <- function(X, ...) {
  UseMethod("dcr")
}

#' @rdname dcr
#' @export
dcr.synth_pair <- function(X, ...) {
  dcr.default(
    X = X$original,
    Y = X$synthetic,
    holdout = X$holdout,
    vars = X$vars,
    ...
  )
}

#' @rdname dcr
#' @export
dcr.default <- function(X, Y,
                        holdout = NULL,
                        holdout_fraction = 0.5,
                        vars = NULL,
                        method = c("gower", "euclidean"),
                        na.rm = TRUE,
                        seed = NULL,
                        progress = FALSE,
                        null_test = TRUE,
                        n_null = 100,
                        ...) {

  method <- match.arg(method)

  # Shared input validation, variable intersection, NA handling, holdout split
  prep <- .distance_risk_prepare(X, Y, holdout = holdout,
                                  holdout_fraction = holdout_fraction,
                                  vars = vars, na.rm = na.rm,
                                  seed = seed, min_holdout = 1L)
  train <- prep$train
  Y <- prep$synthetic
  holdout <- prep$holdout
  vars <- prep$vars

  if (nrow(Y) == 0) stop("No complete cases in Y after removing NAs.")
  if (nrow(train) == 0) stop("No complete cases in X after removing NAs.")
  if (nrow(holdout) == 0) stop("No complete cases in holdout after removing NAs.")

  n_synthetic <- nrow(Y)
  n_train <- nrow(train)
  n_holdout <- nrow(holdout)

  # Min-max normalization helper (used by euclidean method and null test)
  normalize <- function(x) {
    rng <- range(x, na.rm = TRUE)
    if (rng[2] - rng[1] == 0) return(rep(0, length(x)))
    (x - rng[1]) / (rng[2] - rng[1])
  }

  # Compute distances
  if (method == "gower") {
    # Gower distance (handles mixed types)
    # Use single gowerD call with combined reference data so that
    # range normalization is consistent for train and holdout distances
    combined_real <- rbind(train, holdout)
    dist_all <- VIM::gowerD(Y, combined_real)

    # Split into train and holdout columns
    dist_to_train <- dist_all[, seq_len(n_train), drop = FALSE]
    dist_to_holdout <- dist_all[, n_train + seq_len(n_holdout), drop = FALSE]

    # Get minimum distance for each synthetic record
    dcr_train <- apply(dist_to_train, 1, min, na.rm = TRUE)
    dcr_holdout <- apply(dist_to_holdout, 1, min, na.rm = TRUE)

  } else if (method == "euclidean") {
    # Check all variables are numeric
    if (!all(sapply(Y, is.numeric))) {
      stop("method='euclidean' requires all variables to be numeric. Use method='gower' for mixed types.")
    }

    # Combine all data for consistent normalization
    all_data <- rbind(train, holdout, Y)
    all_data_norm <- as.data.frame(lapply(all_data, normalize))

    train_norm <- all_data_norm[1:n_train, , drop = FALSE]
    holdout_norm <- all_data_norm[(n_train + 1):(n_train + n_holdout), , drop = FALSE]
    Y_norm <- all_data_norm[(n_train + n_holdout + 1):nrow(all_data_norm), , drop = FALSE]

    # Compute distances
    dcr_train <- numeric(n_synthetic)
    dcr_holdout <- numeric(n_synthetic)

    # Setup progress bar if requested
    if (progress) {
      pb <- txtProgressBar(min = 0, max = n_synthetic, style = 3)
    }

    for (i in seq_len(n_synthetic)) {
      # Distance to training
      diffs_train <- sweep(as.matrix(train_norm), 2, as.numeric(Y_norm[i, ]))
      dists_train <- sqrt(rowSums(diffs_train^2))
      dcr_train[i] <- min(dists_train)

      # Distance to holdout
      diffs_holdout <- sweep(as.matrix(holdout_norm), 2, as.numeric(Y_norm[i, ]))
      dists_holdout <- sqrt(rowSums(diffs_holdout^2))
      dcr_holdout[i] <- min(dists_holdout)

      if (progress) setTxtProgressBar(pb, i)
    }

    if (progress) close(pb)
  }

  # Compute summary statistics
  dcr_ratio <- mean(dcr_train, na.rm = TRUE) / mean(dcr_holdout, na.rm = TRUE)
  dcr_share <- mean(dcr_train < dcr_holdout, na.rm = TRUE)

  # Statistical test: Wilcoxon signed-rank test
  # H0: synthetic records are equally close to training and holdout
  # H1: synthetic records are systematically closer to training (privacy concern)
  wilcox_result <- tryCatch({
    wilcox.test(dcr_train, dcr_holdout, paired = TRUE, alternative = "less")
  }, error = function(e) {
    list(statistic = NA, p.value = NA, method = "Wilcoxon signed-rank test (failed)")
  })

  # Null distribution comparison (following Yao et al. 2025 recommendations)
  null_distribution <- NULL
  null_share_pvalue <- NA

  if (null_test && n_null > 0) {
    # Generate null distribution by computing DCR share for random subsets
    # Under null hypothesis (no memorization), share should be ~0.5
    null_shares <- numeric(n_null)

    for (b in seq_len(n_null)) {
      # Randomly permute training/holdout assignment
      all_real <- rbind(train, holdout)
      perm_idx <- sample(nrow(all_real))
      perm_train <- all_real[perm_idx[1:n_train], , drop = FALSE]
      perm_holdout <- all_real[perm_idx[(n_train + 1):nrow(all_real)], , drop = FALSE]

      # Compute distances for a random subset of synthetic records (for speed)
      sample_size <- min(50, n_synthetic)
      sample_idx <- sample(n_synthetic, sample_size)

      if (method == "gower") {
        # Single gowerD call for consistent range normalization
        Y_sub <- Y[sample_idx, , drop = FALSE]
        perm_combined <- rbind(perm_train, perm_holdout)
        d_all_perm <- VIM::gowerD(Y_sub, perm_combined)
        d_train_perm <- d_all_perm[, seq_len(n_train), drop = FALSE]
        d_holdout_perm <- d_all_perm[, n_train + seq_len(n_holdout), drop = FALSE]
        dcr_t <- apply(d_train_perm, 1, min, na.rm = TRUE)
        dcr_h <- apply(d_holdout_perm, 1, min, na.rm = TRUE)
      } else {
        # Re-normalize using permuted data for consistent scaling
        perm_all <- rbind(perm_train, perm_holdout, Y)
        perm_all_norm <- as.data.frame(lapply(perm_all, normalize))
        perm_train_norm <- perm_all_norm[seq_len(n_train), , drop = FALSE]
        perm_holdout_norm <- perm_all_norm[n_train + seq_len(n_holdout), , drop = FALSE]
        perm_Y_norm <- perm_all_norm[(n_train + n_holdout + 1):nrow(perm_all_norm), , drop = FALSE]

        dcr_t <- numeric(sample_size)
        dcr_h <- numeric(sample_size)
        for (j in seq_len(sample_size)) {
          idx <- sample_idx[j]
          diffs_t <- sweep(as.matrix(perm_train_norm), 2, as.numeric(perm_Y_norm[idx, ]))
          dcr_t[j] <- min(sqrt(rowSums(diffs_t^2)))
          diffs_h <- sweep(as.matrix(perm_holdout_norm), 2, as.numeric(perm_Y_norm[idx, ]))
          dcr_h[j] <- min(sqrt(rowSums(diffs_h^2)))
        }
      }

      null_shares[b] <- mean(dcr_t < dcr_h, na.rm = TRUE)
    }

    null_distribution <- list(
      shares = null_shares,
      mean = mean(null_shares),
      sd = sd(null_shares),
      quantiles = quantile(null_shares, probs = c(0.025, 0.5, 0.975))
    )

    # P-value: proportion of null shares >= observed share
    # (one-sided test: is observed share unusually high?)
    null_share_pvalue <- mean(null_shares >= dcr_share)
  }

  # Basic privacy check: share should be around 0.5
  # Also consider statistical significance
  privacy_pass <- dcr_share <= 0.55 &&
    (is.na(wilcox_result$p.value) || wilcox_result$p.value > 0.05)

  results <- list(
    dcr_train = dcr_train,
    dcr_holdout = dcr_holdout,
    dcr_ratio = dcr_ratio,
    dcr_share = dcr_share,
    privacy_pass = privacy_pass,
    wilcox_test = wilcox_result,
    null_distribution = null_distribution,
    null_share_pvalue = null_share_pvalue,
    mean_dcr_train = mean(dcr_train, na.rm = TRUE),
    mean_dcr_holdout = mean(dcr_holdout, na.rm = TRUE),
    median_dcr_train = median(dcr_train, na.rm = TRUE),
    median_dcr_holdout = median(dcr_holdout, na.rm = TRUE),
    n_synthetic = n_synthetic,
    n_train = n_train,
    n_holdout = n_holdout,
    method = method,
    vars = vars,
    holdout_fraction = if (prep$was_split) holdout_fraction else NA
  )

  class(results) <- "dcr"
  return(results)
}

#' Print method for dcr objects
#'
#' @param x an object of class "dcr"
#' @param ... additional arguments (ignored)
#' @return The input object, invisibly.
#' @export
print.dcr <- function(x, ...) {
  cat("Distance to Closest Record (DCR) Privacy Metric\n")
  cat("================================================\n")

  # Warning about DCR limitations
  cat("NOTE: DCR has known limitations - see Yao et al. (2025)\n")
  cat("      'The DCR Delusion' (arXiv:2505.01524)\n\n")

  cat("Method:", x$method, "\n")
  cat("Variables used:", length(x$vars), "\n\n")

  cat("Dataset Sizes:\n")
  cat("  Training records:", x$n_train, "\n")
  cat("  Holdout records:", x$n_holdout, "\n")
  cat("  Synthetic records:", x$n_synthetic, "\n\n")

  cat("DCR Results:\n")
  cat("  Mean DCR to training:", round(x$mean_dcr_train, 4), "\n")
  cat("  Mean DCR to holdout:", round(x$mean_dcr_holdout, 4), "\n")
  cat("  DCR ratio (train/holdout):", round(x$dcr_ratio, 4),
      ifelse(x$dcr_ratio < 0.9, " (privacy concern)", ""), "\n")
  cat("  DCR share (closer to train):", sprintf("%.1f%%", 100 * x$dcr_share),
      ifelse(x$dcr_share > 0.55, " (privacy concern)", ""), "\n\n")

  # Statistical test results
  cat("Statistical Tests:\n")
  if (!is.null(x$wilcox_test) && !is.na(x$wilcox_test$p.value)) {
    cat("  Wilcoxon signed-rank test p-value:", format.pval(x$wilcox_test$p.value), "\n")
    cat("  (H0: equal distance to train/holdout; H1: closer to train)\n")
  }
  if (!is.null(x$null_distribution)) {
    cat("  Null distribution share: mean =", round(x$null_distribution$mean, 3),
        ", sd =", round(x$null_distribution$sd, 3), "\n")
    cat("  Observed vs null p-value:", format.pval(x$null_share_pvalue), "\n")
  }
  cat("\n")

  cat("Privacy Assessment:")
  if (x$privacy_pass) {
    cat(" PASS\n")
    cat("  DCR test passed, but this does NOT guarantee privacy.\n")
    cat("  Consider additional metrics (DCAP, TCAP) and MIA evaluation.\n")
  } else {
    cat(" WARNING\n")
    cat("  Synthetic data may have memorized some training records.\n")
  }

  invisible(x)
}

#' Summary method for dcr objects
#'
#' @param object an object of class "dcr"
#' @param ... additional arguments (ignored)
#' @return A list of summary statistics for the corresponding object.
#' @export
summary.dcr <- function(object, ...) {
  summ <- list(
    dcr_ratio = object$dcr_ratio,
    dcr_share = object$dcr_share,
    privacy_pass = object$privacy_pass,
    wilcox_pvalue = if (!is.null(object$wilcox_test)) object$wilcox_test$p.value else NA,
    null_share_pvalue = object$null_share_pvalue,
    null_distribution = object$null_distribution,
    mean_dcr_train = object$mean_dcr_train,
    mean_dcr_holdout = object$mean_dcr_holdout,
    median_dcr_train = object$median_dcr_train,
    median_dcr_holdout = object$median_dcr_holdout,
    sd_dcr_train = sd(object$dcr_train, na.rm = TRUE),
    sd_dcr_holdout = sd(object$dcr_holdout, na.rm = TRUE),
    quantiles_train = quantile(object$dcr_train,
                               probs = c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1),
                               na.rm = TRUE),
    quantiles_holdout = quantile(object$dcr_holdout,
                                 probs = c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1),
                                 na.rm = TRUE),
    n_closer_to_train = sum(object$dcr_train < object$dcr_holdout, na.rm = TRUE),
    n_closer_to_holdout = sum(object$dcr_train >= object$dcr_holdout, na.rm = TRUE),
    n_identical_train = sum(object$dcr_train == 0, na.rm = TRUE),
    n_synthetic = object$n_synthetic,
    n_train = object$n_train,
    n_holdout = object$n_holdout,
    method = object$method,
    vars = object$vars
  )

  class(summ) <- "summary.dcr"
  return(summ)
}

#' Print method for summary.dcr objects
#'
#' @param x an object of class "summary.dcr"
#' @param ... additional arguments (ignored)
#' @return The input object, invisibly.
#' @export
print.summary.dcr <- function(x, ...) {
  cat("Summary: Distance to Closest Record (DCR)\n")
  cat("==========================================\n")
  cat("WARNING: DCR has known limitations - see arXiv:2505.01524\n\n")

  cat("Method:", x$method, "\n")
  cat("Variables:", length(x$vars), "\n\n")

  cat("Dataset Sizes:\n")
  cat("  Training:", x$n_train, "| Holdout:", x$n_holdout,
      "| Synthetic:", x$n_synthetic, "\n\n")

  cat("Key Metrics:\n")
  cat("  DCR ratio:", round(x$dcr_ratio, 4), "(ideal: ~1.0)\n")
  cat("  DCR share:", sprintf("%.1f%%", 100 * x$dcr_share), "(ideal: ~50%)\n")
  cat("  Privacy:", ifelse(x$privacy_pass, "PASS", "WARNING"), "\n\n")

  cat("Statistical Tests:\n")
  if (!is.na(x$wilcox_pvalue)) {
    cat("  Wilcoxon p-value:", format.pval(x$wilcox_pvalue), "\n")
  }
  if (!is.na(x$null_share_pvalue)) {
    cat("  Null comparison p-value:", format.pval(x$null_share_pvalue), "\n")
  }
  if (!is.null(x$null_distribution)) {
    cat("  Null share distribution: mean =", round(x$null_distribution$mean, 3),
        ", 95% CI = [", round(x$null_distribution$quantiles[1], 3), ",",
        round(x$null_distribution$quantiles[3], 3), "]\n")
  }
  cat("\n")

  cat("DCR to Training:\n")
  cat("  Mean:", round(x$mean_dcr_train, 4),
      "| Median:", round(x$median_dcr_train, 4),
      "| SD:", round(x$sd_dcr_train, 4), "\n")
  cat("  Quantiles:\n")
  print(round(x$quantiles_train, 4))
  cat("\n")

  cat("DCR to Holdout:\n")
  cat("  Mean:", round(x$mean_dcr_holdout, 4),
      "| Median:", round(x$median_dcr_holdout, 4),
      "| SD:", round(x$sd_dcr_holdout, 4), "\n")
  cat("  Quantiles:\n")
  print(round(x$quantiles_holdout, 4))
  cat("\n")

  cat("Proximity Analysis:\n")
  cat("  Closer to training:", x$n_closer_to_train,
      sprintf("(%.1f%%)", 100 * x$n_closer_to_train / x$n_synthetic), "\n")
  cat("  Closer to holdout:", x$n_closer_to_holdout,
      sprintf("(%.1f%%)", 100 * x$n_closer_to_holdout / x$n_synthetic), "\n")
  cat("  Identical to training (DCR=0):", x$n_identical_train, "\n")

  invisible(x)
}

#' Plot method for dcr objects
#'
#' @param x an object of class "dcr"
#' @param y not used
#' @param ... additional arguments passed to plotting functions
#' @param which which plot(s) to show: 1 = overlaid histograms, 2 = boxplot
#'   comparison, 3 = scatter plot of train vs holdout DCR
#' @return No return value; called for the side effect of producing a plot.
#' @export
plot.dcr <- function(x, y = NULL, ..., which = 1) {
  show <- rep(FALSE, 3)
  show[which] <- TRUE

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
    breaks <- seq(0, max(c(x$dcr_train, x$dcr_holdout), na.rm = TRUE) * 1.1,
                  length.out = 31)

    hist(x$dcr_train, breaks = breaks, col = rgb(1, 0, 0, 0.5),
         main = "DCR Distribution",
         xlab = "Distance to Closest Record", ...)
    hist(x$dcr_holdout, breaks = breaks, col = rgb(0, 0, 1, 0.5), add = TRUE)
    legend("topright",
           legend = c(paste("To training (mean:", round(x$mean_dcr_train, 3), ")"),
                      paste("To holdout (mean:", round(x$mean_dcr_holdout, 3), ")")),
           fill = c(rgb(1, 0, 0, 0.5), rgb(0, 0, 1, 0.5)),
           cex = 0.8)
  }

  if (show[2]) {
    # Boxplot comparison
    boxplot(list("To Training" = x$dcr_train, "To Holdout" = x$dcr_holdout),
            main = "DCR Comparison",
            ylab = "Distance to Closest Record",
            col = c("coral", "steelblue"), ...)
  }

  if (show[3]) {
    # Scatter plot
    plot(x$dcr_holdout, x$dcr_train,
         main = paste("DCR: Training vs Holdout\nShare closer to train:",
                      sprintf("%.1f%%", 100 * x$dcr_share)),
         xlab = "DCR to Holdout",
         ylab = "DCR to Training",
         pch = 16, col = rgb(0, 0, 0, 0.3), ...)
    abline(0, 1, col = "red", lwd = 2, lty = 2)
    legend("topleft", "y = x (equal distance)", col = "red", lty = 2, lwd = 2,
           cex = 0.8)
  }
}
