#' Density-based Membership Inference Attack (DOMIAS)
#'
#' Detects local overfitting in synthetic data by estimating density ratios
#' between the synthetic data distribution and a reference distribution at each
#' test point. Records where the synthetic density is much higher than the
#' reference density are likely memorized. The AUC of classifying training
#' versus holdout records by their density ratios serves as the attack metric.
#'
#' @param X data frame of original/training data
#' @param Y data frame of synthetic data
#' @param holdout data frame of holdout data (optional). If NULL, a holdout set
#'   is automatically created by splitting X.
#' @param holdout_fraction numeric, fraction of X to use as holdout if holdout
#'   is NULL (default: 0.5)
#' @param radius numeric, Gower distance radius for neighbor counting
#'   (default: 0.1). Smaller values increase sensitivity but require more data.
#' @param vars character vector of variable names to use. If NULL (default),
#'   all common variables between X and Y are used.
#' @param na.rm logical, remove records with NA values (default: TRUE)
#' @param seed integer, random seed for holdout sampling (default: NULL)
#' @param ... additional arguments passed to methods (currently unused)
#'
#' @return An object of class "domias" containing:
#' \itemize{
#'   \item density_ratios_train: numeric vector of density ratios for training records
#'   \item density_ratios_holdout: numeric vector of density ratios for holdout records
#'   \item auc: AUC of train-vs-holdout classification by density ratio
#'   \item mean_ratio_train: mean density ratio for training records
#'   \item mean_ratio_holdout: mean density ratio for holdout records
#'   \item memorization_score: mean_ratio_train / mean_ratio_holdout (ideally ~1)
#'   \item privacy_pass: logical, TRUE if auc <= 0.6
#'   \item n_train, n_holdout, n_synthetic: dataset sizes
#'   \item radius: radius used for density estimation
#'   \item vars: variables used
#' }
#'
#' @details
#' DOMIAS (Density-based Overfitting Membership Inference Attack on Synthetic data)
#' detects whether a synthetic data generator has memorized specific training records.
#' The key insight is that if a generator overfits to training data, the synthetic
#' distribution will have higher density around training points compared to unseen
#' holdout points.
#'
#' The algorithm works by:
#' \enumerate{
#'   \item Splitting original data into training (used for synthesis) and holdout
#'   \item For each test record, estimating synthetic density by counting neighbors
#'     within radius \code{r} in the synthetic data
#'   \item Estimating reference density by counting neighbors within radius \code{r}
#'     in the full original data (training + holdout)
#'   \item Computing the density ratio: (density_synth + 1) / (density_ref + 1),
#'     with Laplace smoothing to avoid division by zero
#'   \item Comparing density ratios for training versus holdout records
#' }
#'
#' Interpretation:
#' \itemize{
#'   \item \strong{AUC ~ 0.5}: Good privacy - density ratios cannot distinguish
#'     training from holdout (no memorization detected)
#'   \item \strong{AUC > 0.5}: Privacy concern - training records have higher
#'     density ratios, suggesting local overfitting
#'   \item \strong{AUC > 0.6}: Likely memorization detected
#'   \item \strong{memorization_score ~ 1}: No differential memorization
#'   \item \strong{memorization_score >> 1}: Synthetic data concentrated around
#'     training records
#' }
#'
#' This implementation uses Gower distance for neighbor counting, which handles
#' mixed data types (numeric, categorical, ordinal). For purely numeric data,
#' the radius parameter corresponds to the Gower distance (normalized between 0 and 1).
#'
#' @section Holdout splitting:
#' When no external holdout is provided, the original data is split internally.
#' The holdout serves as a control group: if the generator has not memorized
#' training data, density ratios should be similar for training and holdout records.
#' For best results, provide a separate holdout set from the synthesis process.
#'
#' @section Choosing the radius:
#' The \code{radius} parameter controls the bandwidth of the kernel density
#' estimate. Guidelines:
#' \itemize{
#'   \item Smaller radius (0.01-0.05): More sensitive to local memorization
#'     but noisier estimates, needs more data
#'   \item Default radius (0.1): Good balance for most datasets
#'   \item Larger radius (0.2-0.5): Smoother estimates but may miss localized
#'     overfitting
#' }
#'
#' @seealso \code{\link{dcr}} for distance to closest record,
#'   \code{\link{nndr}} for nearest neighbor distance ratio,
#'   \code{\link{nnaa}} for nearest-neighbor adversarial accuracy
#'
#' @references
#' Van Breugel, B., Sun, H., Qian, Z. & van der Schaar, M. (2023).
#' Membership Inference Attacks against Synthetic Data through
#' Overfitting Detection.
#' \emph{Proceedings of the 26th International Conference on Artificial
#' Intelligence and Statistics (AISTATS)}, PMLR 206, 3493--3514.
#' \url{https://proceedings.mlr.press/v206/breugel23a.html}
#'
#' @author Matthias Templ
#' @family privacy-models
#' @export
#' @importFrom VIM gowerD
#' @importFrom stats complete.cases wilcox.test median quantile sd
#' @importFrom graphics hist abline legend barplot par boxplot
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
#' result <- domias(X, Y_good, seed = 42)
#' print(result)
#' summary(result)
#'
#' \donttest{
#' # Memorized data (Y is copy of X) - should show high AUC
#' Y_copy <- X[sample(nrow(X), n, replace = TRUE), ]
#' result_bad <- domias(X, Y_copy, seed = 42)
#' print(result_bad)
#' }
domias <- function(X, ...) {
  UseMethod("domias")
}

#' @rdname domias
#' @export
domias.synth_pair <- function(X, ...) {
  if (!is.null(X$source) && X$source == "sdcMicro") {
    stop("domias is designed for synthetic data evaluation and is not applicable to ",
         "traditionally anonymized data (sdcMicro objects). ",
         "Use dcr, nndr, or ims for distance-based privacy evaluation instead.",
         call. = FALSE)
  }
  domias.default(
    X = X$original,
    Y = X$synthetic,
    holdout = X$holdout,
    vars = X$vars,
    ...
  )
}

#' @rdname domias
#' @export
domias.default <- function(X, Y,
                            holdout = NULL,
                            holdout_fraction = 0.5,
                            radius = 0.1,
                            vars = NULL,
                            na.rm = TRUE,
                            seed = NULL,
                            ...) {

  # Input validation
  if (!is.data.frame(X)) stop("X must be a data frame.")
  if (!is.data.frame(Y)) stop("Y must be a data frame.")
  if (!is.null(holdout) && !is.data.frame(holdout)) {
    stop("holdout must be a data frame or NULL.")
  }

  if (!is.numeric(radius) || length(radius) != 1 || radius <= 0 || radius >= 1) {
    stop("radius must be a number between 0 and 1 (exclusive).")
  }

  if (!is.numeric(holdout_fraction) || length(holdout_fraction) != 1 ||
      holdout_fraction <= 0 || holdout_fraction >= 1) {
    stop("holdout_fraction must be a number between 0 and 1 (exclusive).")
  }

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

  # Create or validate holdout
  if (is.null(holdout)) {
    if (!is.null(seed)) set.seed(seed)
    n_holdout <- max(1, floor(nrow(X) * holdout_fraction))
    holdout_idx <- sample(nrow(X), n_holdout)
    holdout <- X[holdout_idx, , drop = FALSE]
    train <- X[-holdout_idx, , drop = FALSE]
  } else {
    holdout <- holdout[, vars, drop = FALSE]
    if (na.rm) {
      complete_H <- complete.cases(holdout)
      holdout <- holdout[complete_H, , drop = FALSE]
    }
    if (nrow(holdout) == 0) stop("No complete cases in holdout after removing NAs.")
    train <- X
  }

  n_train <- nrow(train)
  n_holdout <- nrow(holdout)
  n_synthetic <- nrow(Y)

  # Reference data = full original (training + holdout)
  reference <- rbind(train, holdout)
  n_ref <- nrow(reference)

  # Use single gowerD calls with all test points combined against each target
  # so that range normalization is consistent across train and holdout
  test_points <- rbind(train, holdout)  # rows 1:n_train = train, rest = holdout

  # Distances from all test points to synthetic data (consistent normalization)
  dist_to_synth <- VIM::gowerD(test_points, Y)

  # Distances from all test points to reference data (consistent normalization)
  dist_to_ref <- VIM::gowerD(test_points, reference)

  # Split into train and holdout portions
  dist_train_to_synth <- dist_to_synth[seq_len(n_train), , drop = FALSE]
  dist_holdout_to_synth <- dist_to_synth[n_train + seq_len(n_holdout), , drop = FALSE]
  dist_train_to_ref <- dist_to_ref[seq_len(n_train), , drop = FALSE]
  dist_holdout_to_ref <- dist_to_ref[n_train + seq_len(n_holdout), , drop = FALSE]

  # Count neighbors within radius for training records
  density_synth_train <- rowSums(dist_train_to_synth <= radius)
  # Self-exclude: training records appear in the reference at positions 1:n_train.
  # Record i in train maps to column i in reference. Set self-distance to > radius.
  for (i in seq_len(n_train)) {
    dist_train_to_ref[i, i] <- radius + 1
  }
  density_ref_train <- rowSums(dist_train_to_ref <= radius)

  # Count neighbors within radius for holdout records
  density_synth_holdout <- rowSums(dist_holdout_to_synth <= radius)
  # Self-exclude: holdout records appear in reference at positions
  # (n_train + 1) : (n_train + n_holdout). Record j in holdout maps to
  # column (n_train + j) in reference.
  for (j in seq_len(n_holdout)) {
    dist_holdout_to_ref[j, n_train + j] <- radius + 1
  }
  density_ref_holdout <- rowSums(dist_holdout_to_ref <= radius)

  # Compute density ratios with Laplace smoothing
  density_ratios_train <- (density_synth_train + 1) / (density_ref_train + 1)
  density_ratios_holdout <- (density_synth_holdout + 1) / (density_ref_holdout + 1)

  # Compute AUC via Wilcoxon-Mann-Whitney statistic
  # AUC = P(ratio_train > ratio_holdout)
  auc <- tryCatch({
    w <- wilcox.test(density_ratios_train, density_ratios_holdout,
                     exact = FALSE)$statistic
    as.numeric(w / (n_train * n_holdout))
  }, error = function(e) NA_real_)

  # Summary statistics
  mean_ratio_train <- mean(density_ratios_train, na.rm = TRUE)
  mean_ratio_holdout <- mean(density_ratios_holdout, na.rm = TRUE)
  memorization_score <- mean_ratio_train / mean_ratio_holdout

  # Privacy assessment: AUC close to 0.5 = good

  privacy_pass <- !is.na(auc) && auc <= 0.6

  results <- list(
    density_ratios_train = density_ratios_train,
    density_ratios_holdout = density_ratios_holdout,
    auc = auc,
    mean_ratio_train = mean_ratio_train,
    mean_ratio_holdout = mean_ratio_holdout,
    memorization_score = memorization_score,
    privacy_pass = privacy_pass,
    n_train = n_train,
    n_holdout = n_holdout,
    n_synthetic = n_synthetic,
    radius = radius,
    vars = vars
  )

  class(results) <- "domias"
  return(results)
}

# --- S3 methods ---

#' Print method for domias objects
#'
#' @param x an object of class "domias"
#' @param ... additional arguments (ignored)
#' @export
print.domias <- function(x, ...) {
  cat("DOMIAS: Density-based Membership Inference Attack\n")
  cat("==================================================\n")
  cat("Variables used:", length(x$vars), "\n")
  cat("Density estimation radius:", x$radius, "\n\n")

  cat("Dataset Sizes:\n")
  cat("  Training records:", x$n_train, "\n")
  cat("  Holdout records:", x$n_holdout, "\n")
  cat("  Synthetic records:", x$n_synthetic, "\n\n")

  cat("Density Ratio Results:\n")
  cat("  Mean ratio (training):", round(x$mean_ratio_train, 4), "\n")
  cat("  Mean ratio (holdout):", round(x$mean_ratio_holdout, 4), "\n")
  cat("  Memorization score:", round(x$memorization_score, 4),
      ifelse(x$memorization_score > 1.5, " (privacy concern)", ""), "\n\n")

  cat("Attack Performance:\n")
  if (!is.na(x$auc)) {
    cat("  AUC (train vs holdout):", round(x$auc, 4), "\n")
    cat("  (0.5 = no memorization, 1.0 = complete memorization)\n\n")
  } else {
    cat("  AUC: could not be computed\n\n")
  }

  cat("Privacy Assessment:")
  if (x$privacy_pass) {
    cat(" PASS\n")
    cat("  No significant memorization detected (AUC <= 0.6).\n")
  } else {
    cat(" WARNING\n")
    cat("  Elevated memorization risk detected (AUC > 0.6).\n")
    cat("  Synthetic data may have overfitted to training records.\n")
  }

  invisible(x)
}

#' Summary method for domias objects
#'
#' @param object an object of class "domias"
#' @param ... additional arguments (ignored)
#' @export
summary.domias <- function(object, ...) {
  summ <- list(
    auc = object$auc,
    mean_ratio_train = object$mean_ratio_train,
    mean_ratio_holdout = object$mean_ratio_holdout,
    memorization_score = object$memorization_score,
    privacy_pass = object$privacy_pass,
    median_ratio_train = median(object$density_ratios_train, na.rm = TRUE),
    median_ratio_holdout = median(object$density_ratios_holdout, na.rm = TRUE),
    sd_ratio_train = sd(object$density_ratios_train, na.rm = TRUE),
    sd_ratio_holdout = sd(object$density_ratios_holdout, na.rm = TRUE),
    quantiles_train = quantile(object$density_ratios_train,
                               probs = c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1),
                               na.rm = TRUE),
    quantiles_holdout = quantile(object$density_ratios_holdout,
                                 probs = c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1),
                                 na.rm = TRUE),
    n_high_ratio_train = sum(object$density_ratios_train > 2, na.rm = TRUE),
    n_high_ratio_holdout = sum(object$density_ratios_holdout > 2, na.rm = TRUE),
    n_train = object$n_train,
    n_holdout = object$n_holdout,
    n_synthetic = object$n_synthetic,
    radius = object$radius,
    vars = object$vars
  )

  class(summ) <- "summary.domias"
  return(summ)
}

#' Print method for summary.domias objects
#'
#' @param x an object of class "summary.domias"
#' @param ... additional arguments (ignored)
#' @export
print.summary.domias <- function(x, ...) {
  cat("Summary: DOMIAS Membership Inference Attack\n")
  cat("============================================\n")
  cat("Variables:", length(x$vars), "\n")
  cat("Radius:", x$radius, "\n\n")

  cat("Dataset Sizes:\n")
  cat("  Training:", x$n_train, "| Holdout:", x$n_holdout,
      "| Synthetic:", x$n_synthetic, "\n\n")

  cat("Key Metrics:\n")
  cat("  AUC:", round(x$auc, 4), "(ideal: ~0.5)\n")
  cat("  Memorization score:", round(x$memorization_score, 4), "(ideal: ~1.0)\n")
  cat("  Privacy:", ifelse(x$privacy_pass, "PASS", "WARNING"), "\n\n")

  cat("Density Ratios (Training):\n")
  cat("  Mean:", round(x$mean_ratio_train, 4),
      "| Median:", round(x$median_ratio_train, 4),
      "| SD:", round(x$sd_ratio_train, 4), "\n")
  cat("  High ratio (> 2):", x$n_high_ratio_train,
      sprintf("(%.1f%%)", 100 * x$n_high_ratio_train / x$n_train), "\n")
  cat("  Quantiles:\n")
  print(round(x$quantiles_train, 4))
  cat("\n")

  cat("Density Ratios (Holdout):\n")
  cat("  Mean:", round(x$mean_ratio_holdout, 4),
      "| Median:", round(x$median_ratio_holdout, 4),
      "| SD:", round(x$sd_ratio_holdout, 4), "\n")
  cat("  High ratio (> 2):", x$n_high_ratio_holdout,
      sprintf("(%.1f%%)", 100 * x$n_high_ratio_holdout / x$n_holdout), "\n")
  cat("  Quantiles:\n")
  print(round(x$quantiles_holdout, 4))

  invisible(x)
}

#' Plot method for domias objects
#'
#' @param x an object of class "domias"
#' @param y not used
#' @param ... additional arguments passed to plotting functions
#' @param which which plot(s) to show: 1 = density ratio distributions
#'   (train vs holdout overlaid histograms), 2 = boxplot comparison of
#'   density ratios
#' @export
plot.domias <- function(x, y = NULL, ..., which = 1) {
  show <- rep(FALSE, 2)
  show[which] <- TRUE

  n_plots <- sum(show)
  if (n_plots > 1) {
    op <- par(mfrow = c(1, 2))
    on.exit(par(op))
  }

  if (show[1]) {
    # Overlaid histograms of density ratios
    all_ratios <- c(x$density_ratios_train, x$density_ratios_holdout)
    breaks <- seq(0, max(all_ratios, na.rm = TRUE) * 1.1, length.out = 31)

    hist(x$density_ratios_train, breaks = breaks,
         col = rgb(1, 0, 0, 0.5),
         main = paste("DOMIAS Density Ratios\nAUC:",
                      round(x$auc, 3)),
         xlab = "Density Ratio (synthetic / reference)",
         ylab = "Frequency", ...)
    hist(x$density_ratios_holdout, breaks = breaks,
         col = rgb(0, 0, 1, 0.5), add = TRUE)
    abline(v = 1, col = "grey40", lwd = 2, lty = 2)
    legend("topright",
           legend = c(
             paste("Training (mean:", round(x$mean_ratio_train, 3), ")"),
             paste("Holdout (mean:", round(x$mean_ratio_holdout, 3), ")"),
             "Ratio = 1 (no overfitting)"
           ),
           fill = c(rgb(1, 0, 0, 0.5), rgb(0, 0, 1, 0.5), NA),
           border = c("black", "black", NA),
           lty = c(NA, NA, 2), lwd = c(NA, NA, 2),
           col = c(NA, NA, "grey40"),
           cex = 0.8)
  }

  if (show[2]) {
    # Boxplot comparison
    boxplot(list("Training" = x$density_ratios_train,
                 "Holdout" = x$density_ratios_holdout),
            main = paste("DOMIAS Density Ratios\nMemorization Score:",
                         round(x$memorization_score, 3)),
            ylab = "Density Ratio",
            col = c("coral", "steelblue"), ...)
    abline(h = 1, col = "grey40", lwd = 2, lty = 2)
  }
}
