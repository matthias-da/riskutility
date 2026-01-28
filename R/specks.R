#' SPECKS - Propensity Score Comparison via Kolmogorov-Smirnov Test
#'
#' Computes the SPECKS (Synthetic data generation; Propensity score matching;
#' Empirical Comparison via the Kolmogorov-Smirnov distance) utility measure
#' for synthetic data. This metric assesses how distinguishable synthetic data
#' is from original data using propensity scores.
#'
#' @param X data frame of original data, or a \code{\link{synth_pair}} object
#' @param Y data frame of synthetic data (not needed if X is a synth_pair)
#' @param vars character vector of variable names to use. If NULL, all common
#'   variables are used.
#' @param method character, method for propensity score estimation:
#'   "cart" (classification tree, default), "logit" (logistic regression),
#'   or "rf" (random forest)
#' @param maxorder integer, maximum order of interactions for logit method (default: 1)
#' @param k integer, number of percentiles to use for KS calculation (default: NULL
#'   uses all unique values)
#' @param na.rm logical, remove records with NA values (default: TRUE)
#' @param seed integer, random seed for reproducibility (default: NULL)
#' @param ... additional arguments passed to methods (currently unused)
#'
#' @return An object of class "specks" containing:
#' \itemize{
#'   \item specks: the SPECKS statistic (KS test statistic on propensity scores)
#'   \item pMSE: propensity score Mean Squared Error
#'   \item pMSE_ratio: ratio of pMSE to null expectation
#'   \item ks_pvalue: p-value from the KS test
#'   \item propensity_original: propensity scores for original records
#'   \item propensity_synthetic: propensity scores for synthetic records
#'   \item utility_score: 1 - specks (higher = better utility)
#'   \item method: the method used
#'   \item n_original, n_synthetic: dataset sizes
#' }
#'
#' @details
#' SPECKS (Woo et al. 2009) is a utility measure that assesses how well synthetic
#' data preserves the joint distribution of variables. It works by:
#'
#' 1. Combining original (X) and synthetic (Y) data with labels
#' 2. Fitting a propensity score model to predict the probability of being synthetic
#' 3. Computing the Kolmogorov-Smirnov statistic between the propensity score
#'    distributions of original and synthetic records
#'
#' Interpretation:
#' \itemize{
#'   \item \strong{SPECKS near 0}: Excellent utility - synthetic indistinguishable from original
#'   \item \strong{SPECKS < 0.1}: Good utility - minor distributional differences
#'   \item \strong{SPECKS < 0.2}: Moderate utility - noticeable differences
#'   \item \strong{SPECKS > 0.3}: Poor utility - substantial differences
#' }
#'
#' The pMSE (propensity Mean Squared Error) measures the average squared deviation
#' of propensity scores from the expected value under random assignment:
#' \deqn{pMSE = \frac{1}{N} \sum_{i=1}^{N} (\hat{p}_i - c)^2}
#' where \eqn{c = n_Y / (n_X + n_Y)} is the proportion of synthetic records.
#'
#' The pMSE ratio compares pMSE to its null expectation, with values near 1
#' indicating good utility.
#'
#' @seealso \code{\link{propscore}} for related propensity score utility,
#'   \code{\link{compare_distributions_cont}} for distributional comparisons
#'
#' @references
#' Woo, M. J., Reiter, J. P., Oganian, A., & Karr, A. F. (2009). Global measures
#' of data utility for microdata masked for disclosure limitation.
#' \emph{Journal of Privacy and Confidentiality}, 1(1), 111-124.
#'
#' Snoke, J., Raab, G. M., Nowok, B., Dibben, C., & Slavkovic, A. (2018).
#' General and specific utility measures for synthetic data.
#' \emph{Journal of the Royal Statistical Society: Series A}, 181(3), 663-688.
#'
#' @author Matthias Templ
#' @export
#' @importFrom stats ks.test glm binomial predict complete.cases as.formula
#' @examples
#' # Create example data
#' set.seed(123)
#' n <- 500
#' original <- data.frame(
#'   age = sample(18:80, n, replace = TRUE),
#'   gender = sample(c("M", "F"), n, replace = TRUE),
#'   income = exp(rnorm(n, 10, 1)),
#'   education = sample(c("low", "medium", "high"), n, replace = TRUE)
#' )
#'
#' # Good synthetic data (similar distributions)
#' synthetic_good <- data.frame(
#'   age = sample(18:80, n, replace = TRUE),
#'   gender = sample(c("M", "F"), n, replace = TRUE),
#'   income = exp(rnorm(n, 10, 1)),
#'   education = sample(c("low", "medium", "high"), n, replace = TRUE)
#' )
#'
#' result_good <- specks(original, synthetic_good)
#' print(result_good)
#' summary(result_good)
#'
#' # Poor synthetic data (different distributions)
#' synthetic_poor <- data.frame(
#'   age = sample(18:80, n, replace = TRUE, prob = c(rep(0.01, 30), rep(0.02, 33))),
#'   gender = sample(c("M", "F"), n, replace = TRUE, prob = c(0.8, 0.2)),
#'   income = exp(rnorm(n, 11, 0.5)),
#'   education = sample(c("low", "medium", "high"), n, replace = TRUE, prob = c(0.1, 0.2, 0.7))
#' )
#'
#' result_poor <- specks(original, synthetic_poor)
#' print(result_poor)
#'
#' # Compare results
#' cat("Good synthetic SPECKS:", result_good$specks, "\n")
#' cat("Poor synthetic SPECKS:", result_poor$specks, "\n")
specks <- function(X, ...) {
  UseMethod("specks")
}

#' @rdname specks
#' @export
specks.synth_pair <- function(X, ...) {
  specks.default(
    X = X$original,
    Y = X$synthetic,
    vars = X$vars,
    ...
  )
}

#' @rdname specks
#' @export
specks.default <- function(X, Y,
                           vars = NULL,
                           method = c("cart", "logit", "rf"),
                           maxorder = 1,
                           k = NULL,
                           na.rm = TRUE,
                           seed = NULL,
                           ...) {

  # Input validation
  if (!is.data.frame(X)) stop("X must be a data frame.")
  if (!is.data.frame(Y)) stop("Y must be a data frame.")

  method <- match.arg(method)

  # Set seed if provided
  if (!is.null(seed)) set.seed(seed)

  # Determine variables to use
  if (is.null(vars)) {
    vars <- intersect(names(X), names(Y))
  }

  if (length(vars) == 0) {
    stop("No common variables found between X and Y.")
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
  X_sub <- X[, vars, drop = FALSE]
  Y_sub <- Y[, vars, drop = FALSE]

  # Handle missing values
  if (na.rm) {
    complete_X <- complete.cases(X_sub)
    complete_Y <- complete.cases(Y_sub)
    X_sub <- X_sub[complete_X, , drop = FALSE]
    Y_sub <- Y_sub[complete_Y, , drop = FALSE]
  }

  if (nrow(X_sub) == 0) stop("No complete cases in X after removing NAs.")
  if (nrow(Y_sub) == 0) stop("No complete cases in Y after removing NAs.")

  n_X <- nrow(X_sub)
  n_Y <- nrow(Y_sub)
  N <- n_X + n_Y

  # Expected proportion of synthetic records
  c_expected <- n_Y / N

  # Combine datasets with indicator
  combined <- rbind(X_sub, Y_sub)
  combined$..synthetic.. <- factor(c(rep(0, n_X), rep(1, n_Y)))

  # Fit propensity score model
  if (method == "logit") {
    # Logistic regression with interactions up to maxorder
    if (maxorder == 1) {
      form <- as.formula(paste("..synthetic.. ~", paste(vars, collapse = " + ")))
    } else {
      # Create interaction terms
      form <- as.formula(paste("..synthetic.. ~ (", paste(vars, collapse = " + "), ")^", maxorder))
    }

    model <- tryCatch({
      glm(form, data = combined, family = binomial())
    }, error = function(e) {
      # Fall back to simple model if interactions fail
      warning("Full interaction model failed, using additive model")
      glm(as.formula(paste("..synthetic.. ~", paste(vars, collapse = " + "))),
          data = combined, family = binomial())
    })

    propensity <- predict(model, type = "response")

  } else if (method == "cart") {
    # Classification tree
    if (!requireNamespace("rpart", quietly = TRUE)) {
      stop("Package 'rpart' required for method='cart'. Please install it.")
    }

    form <- as.formula(paste("..synthetic.. ~", paste(vars, collapse = " + ")))
    model <- rpart::rpart(form, data = combined, method = "class")
    propensity <- predict(model, type = "prob")[, 2]

  } else if (method == "rf") {
    # Random forest
    form <- as.formula(paste("..synthetic.. ~", paste(vars, collapse = " + ")))
    model <- randomForest::randomForest(form, data = combined)
    propensity <- predict(model, type = "prob")[, 2]
  }

  # Split propensity scores
  prop_original <- propensity[1:n_X]
  prop_synthetic <- propensity[(n_X + 1):N]

  # Compute SPECKS (KS statistic)
  ks_result <- ks.test(prop_original, prop_synthetic)
  specks_stat <- as.numeric(ks_result$statistic)
  ks_pvalue <- ks_result$p.value

  # Compute pMSE
  pMSE <- mean((propensity - c_expected)^2)

  # Compute null expectation of pMSE
  # Under null (no difference), expected pMSE = c * (1 - c) / N approximately
  # More precisely, for a good fit: E[pMSE] ~ c(1-c) * k / N where k = number of parameters
  # For simplicity, use c(1-c)/N as baseline
  pMSE_null <- c_expected * (1 - c_expected) / N
  pMSE_ratio <- pMSE / pMSE_null

  # Alternative standardized pMSE (Snoke et al. 2018)
  # This divides by c(1-c) to get a scale-free measure
  pMSE_standardized <- pMSE / (c_expected * (1 - c_expected))

  # Utility score (1 - specks, higher is better)
  utility_score <- 1 - specks_stat

  # Build results
  results <- list(
    specks = specks_stat,
    pMSE = pMSE,
    pMSE_ratio = pMSE_ratio,
    pMSE_standardized = pMSE_standardized,
    ks_pvalue = ks_pvalue,
    propensity_original = prop_original,
    propensity_synthetic = prop_synthetic,
    c_expected = c_expected,
    utility_score = utility_score,
    method = method,
    vars = vars,
    n_original = n_X,
    n_synthetic = n_Y
  )

  class(results) <- "specks"
  return(results)
}


#' Print method for specks objects
#'
#' @param x an object of class "specks"
#' @param ... additional arguments (ignored)
#' @export
print.specks <- function(x, ...) {
  cat("SPECKS - Propensity Score Utility Measure\n")
  cat("==========================================\n\n")

  cat("Dataset Sizes:\n")
  cat("  Original (X):", x$n_original, "| Synthetic (Y):", x$n_synthetic, "\n")
  cat("  Variables used:", length(x$vars), "\n")
  cat("  Method:", x$method, "\n\n")

  cat("Results:\n")
  cat("  SPECKS (KS statistic):", sprintf("%.4f", x$specks), "\n")
  cat("  KS test p-value:", sprintf("%.4f", x$ks_pvalue), "\n")
  cat("  pMSE:", sprintf("%.6f", x$pMSE), "\n")
  cat("  pMSE ratio:", sprintf("%.2f", x$pMSE_ratio), "\n")
  cat("  Utility score (1-SPECKS):", sprintf("%.4f", x$utility_score), "\n\n")

  cat("Interpretation:\n")
  if (x$specks < 0.05) {
    cat("  EXCELLENT: Synthetic data is virtually indistinguishable from original.\n")
  } else if (x$specks < 0.1) {
    cat("  GOOD: Minor differences; synthetic data preserves distributions well.\n")
  } else if (x$specks < 0.2) {
    cat("  MODERATE: Noticeable differences in joint distributions.\n")
  } else if (x$specks < 0.3) {
    cat("  FAIR: Substantial differences; utility may be limited for some analyses.\n")
  } else {
    cat("  POOR: Synthetic data differs significantly from original.\n")
  }

  invisible(x)
}


#' Summary method for specks objects
#'
#' @param object an object of class "specks"
#' @param ... additional arguments (ignored)
#' @export
summary.specks <- function(object, ...) {
  prop_orig <- object$propensity_original
  prop_synth <- object$propensity_synthetic

  summ <- list(
    specks = object$specks,
    pMSE = object$pMSE,
    pMSE_ratio = object$pMSE_ratio,
    pMSE_standardized = object$pMSE_standardized,
    ks_pvalue = object$ks_pvalue,
    utility_score = object$utility_score,
    prop_original_summary = c(
      mean = mean(prop_orig),
      sd = sd(prop_orig),
      min = min(prop_orig),
      q25 = quantile(prop_orig, 0.25),
      median = median(prop_orig),
      q75 = quantile(prop_orig, 0.75),
      max = max(prop_orig)
    ),
    prop_synthetic_summary = c(
      mean = mean(prop_synth),
      sd = sd(prop_synth),
      min = min(prop_synth),
      q25 = quantile(prop_synth, 0.25),
      median = median(prop_synth),
      q75 = quantile(prop_synth, 0.75),
      max = max(prop_synth)
    ),
    c_expected = object$c_expected,
    method = object$method,
    vars = object$vars,
    n_original = object$n_original,
    n_synthetic = object$n_synthetic
  )

  class(summ) <- "summary.specks"
  return(summ)
}


#' Print method for summary.specks objects
#'
#' @param x an object of class "summary.specks"
#' @param ... additional arguments (ignored)
#' @export
print.summary.specks <- function(x, ...) {
  cat("Summary: SPECKS Utility Measure\n")
  cat("================================\n\n")

  cat("Method:", x$method, "\n")
  cat("Variables:", paste(x$vars, collapse = ", "), "\n")
  cat("Sample sizes: Original =", x$n_original, ", Synthetic =", x$n_synthetic, "\n\n")

  cat("Main Results:\n")
  cat("  SPECKS:", sprintf("%.4f", x$specks), "\n")
  cat("  pMSE:", sprintf("%.6f", x$pMSE), "\n")
  cat("  pMSE ratio:", sprintf("%.2f", x$pMSE_ratio), "\n")
  cat("  pMSE standardized:", sprintf("%.4f", x$pMSE_standardized), "\n")
  cat("  KS p-value:", sprintf("%.4f", x$ks_pvalue), "\n")
  cat("  Utility score:", sprintf("%.4f", x$utility_score), "\n\n")

  cat("Propensity Score Distribution (Original):\n")
  cat("  Mean:", sprintf("%.4f", x$prop_original_summary["mean"]),
      " SD:", sprintf("%.4f", x$prop_original_summary["sd"]), "\n")
  cat("  Range: [", sprintf("%.4f", x$prop_original_summary["min"]),
      ",", sprintf("%.4f", x$prop_original_summary["max"]), "]\n")
  cat("  Quartiles:", sprintf("%.4f", x$prop_original_summary["q25"]),
      sprintf("%.4f", x$prop_original_summary["median"]),
      sprintf("%.4f", x$prop_original_summary["q75"]), "\n\n")

  cat("Propensity Score Distribution (Synthetic):\n")
  cat("  Mean:", sprintf("%.4f", x$prop_synthetic_summary["mean"]),
      " SD:", sprintf("%.4f", x$prop_synthetic_summary["sd"]), "\n")
  cat("  Range: [", sprintf("%.4f", x$prop_synthetic_summary["min"]),
      ",", sprintf("%.4f", x$prop_synthetic_summary["max"]), "]\n")
  cat("  Quartiles:", sprintf("%.4f", x$prop_synthetic_summary["q25"]),
      sprintf("%.4f", x$prop_synthetic_summary["median"]),
      sprintf("%.4f", x$prop_synthetic_summary["q75"]), "\n\n")

  cat("Expected proportion (c):", sprintf("%.4f", x$c_expected), "\n")
  cat("  (Ideal propensity scores cluster around this value)\n")

  invisible(x)
}


#' Plot method for specks objects
#'
#' @param x an object of class "specks"
#' @param y not used
#' @param ... additional arguments passed to plotting functions
#' @param which integer, which plot: 1 = density comparison (default),
#'   2 = ECDF comparison, 3 = histogram comparison
#' @importFrom graphics hist lines legend abline par
#' @importFrom stats density ecdf
#' @export
plot.specks <- function(x, y = NULL, ..., which = 1) {
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

  prop_orig <- x$propensity_original
  prop_synth <- x$propensity_synthetic

  if (show[1]) {
    # Density comparison
    dens_orig <- density(prop_orig, from = 0, to = 1)
    dens_synth <- density(prop_synth, from = 0, to = 1)

    plot(dens_orig, main = "Propensity Score Densities",
         xlab = "Propensity Score", ylab = "Density",
         col = "steelblue", lwd = 2,
         xlim = c(0, 1),
         ylim = c(0, max(c(dens_orig$y, dens_synth$y)) * 1.1), ...)
    lines(dens_synth, col = "coral", lwd = 2)
    abline(v = x$c_expected, col = "gray40", lty = 2, lwd = 1.5)
    legend("topright",
           legend = c("Original", "Synthetic", paste0("Expected (", round(x$c_expected, 2), ")")),
           col = c("steelblue", "coral", "gray40"),
           lty = c(1, 1, 2), lwd = c(2, 2, 1.5))
  }

  if (show[2]) {
    # ECDF comparison
    ecdf_orig <- ecdf(prop_orig)
    ecdf_synth <- ecdf(prop_synth)

    plot(ecdf_orig, main = paste0("ECDF Comparison (SPECKS = ", round(x$specks, 3), ")"),
         xlab = "Propensity Score", ylab = "Cumulative Probability",
         col = "steelblue", lwd = 2,
         xlim = c(0, 1), ...)
    lines(ecdf_synth, col = "coral", lwd = 2)
    abline(v = x$c_expected, col = "gray40", lty = 2)
    legend("bottomright",
           legend = c("Original", "Synthetic"),
           col = c("steelblue", "coral"),
           lty = 1, lwd = 2)
  }

  if (show[3]) {
    # Histogram comparison
    breaks <- seq(0, 1, by = 0.05)

    hist(prop_orig, breaks = breaks, col = adjustcolor("steelblue", 0.5),
         main = "Propensity Score Distributions",
         xlab = "Propensity Score", ylab = "Frequency",
         xlim = c(0, 1), ...)
    hist(prop_synth, breaks = breaks, col = adjustcolor("coral", 0.5), add = TRUE)
    abline(v = x$c_expected, col = "gray40", lty = 2, lwd = 2)
    legend("topright",
           legend = c("Original", "Synthetic"),
           fill = adjustcolor(c("steelblue", "coral"), 0.5))
  }
}
