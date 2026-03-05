#' Chi-Square Utility Measures
#'
#' Computes chi-square based utility measures comparing original and synthetic data
#' frequency tables. Implements VW (Voas-Williamson), FT (Freeman-Tukey), G (log-likelihood),
#' and JSD (Jensen-Shannon Divergence) statistics.
#'
#' @param X original data frame, or a \code{\link{synth_pair}} object
#' @param Y synthetic data frame (not required if X is a synth_pair)
#' @param vars character vector of variable names to include in tables
#' @param max_cells integer, maximum number of cells in cross-tabulation (default: 10000)
#' @param weight_X optional, name of weight variable in X or numeric vector
#' @param weight_Y optional, name of weight variable in Y or numeric vector
#' @param ... additional arguments passed to methods
#'
#' @return An object of class "chisq_utility" containing:
#' \itemize{
#'   \item chi2: standard chi-square statistic
#'   \item df: degrees of freedom
#'   \item p_value: p-value for chi-square test
#'   \item VW: Voas-Williamson statistic (normalized chi-square)
#'   \item FT: Freeman-Tukey statistic
#'   \item G: G-test (log-likelihood ratio) statistic
#'   \item JSD: Jensen-Shannon Divergence
#'   \item n_cells: number of cells in the table
#'   \item n_empty_orig: empty cells in original
#'   \item n_empty_synth: empty cells in synthetic
#'   \item pct_utility: overall utility percentage (based on VW)
#'   \item table_orig: frequency table for original data
#'   \item table_synth: frequency table for synthetic data
#' }
#'
#' @details
#' These measures compare frequency tables from original and synthetic data.
#'
#' \strong{Chi-Square (\eqn{\chi^2}):}
#' The standard Pearson chi-square statistic:
#' \deqn{\chi^2 = \sum \frac{(O - E)^2}{E}}
#' where O is observed (synthetic) and E is expected (original proportions scaled).
#'
#' \strong{Voas-Williamson (VW):}
#' A normalized chi-square measure (Voas & Williamson, 2001):
#' \deqn{VW = \frac{\chi^2 - df}{N}}
#' where df is degrees of freedom and N is sample size. Lower is better; 0 indicates
#' perfect replication.
#'
#' \strong{Freeman-Tukey (FT):}
#' Uses a variance-stabilizing transformation:
#' \deqn{FT = 4 \sum (\sqrt{O} - \sqrt{E})^2}
#'
#' \strong{G-Test (Log-Likelihood Ratio):}
#' \deqn{G = 2 \sum O \log(O/E)}
#'
#' \strong{Jensen-Shannon Divergence (JSD):}
#' A symmetric, bounded divergence measure:
#' \deqn{JSD = \frac{1}{2}[KL(P||M) + KL(Q||M)]}
#' where M = (P + Q)/2, and KL is Kullback-Leibler divergence.
#'
#' \strong{Interpretation:}
#' \itemize{
#'   \item VW close to 0: excellent utility
#'   \item VW < 0.01: good utility
#'   \item VW > 0.05: potential utility concerns
#'   \item JSD ranges 0-1: 0 = identical, 1 = completely different
#' }
#'
#' @seealso \code{\link{propscore}} for propensity score utility,
#'   \code{\link{specks}} for SPECKS utility measure
#'
#' @references
#' Voas, D., & Williamson, P. (2001). Evaluating goodness-of-fit measures for
#' synthetic microdata. \emph{Geographical and Environmental Modelling}, 5(2), 177-200.
#'
#' Freeman, M. F., & Tukey, J. W. (1950). Transformations related to the angular
#' and the square root. \emph{Annals of Mathematical Statistics}, 21, 607-611.
#'
#' @author Matthias Templ
#' @importFrom stats pchisq xtabs
#' @family utility
#' @export
#' @examples
#' # Create example data
#' set.seed(123)
#' orig <- data.frame(
#'   age = sample(c("young", "middle", "old"), 500, replace = TRUE),
#'   gender = sample(c("M", "F"), 500, replace = TRUE),
#'   region = sample(c("N", "S", "E", "W"), 500, replace = TRUE)
#' )
#'
#' # Synthetic data with some differences
#' synth <- data.frame(
#'   age = sample(c("young", "middle", "old"), 500, replace = TRUE,
#'                prob = c(0.35, 0.4, 0.25)),
#'   gender = sample(c("M", "F"), 500, replace = TRUE),
#'   region = sample(c("N", "S", "E", "W"), 500, replace = TRUE)
#' )
#'
#' # Compute chi-square utility
#' result <- chisq_utility(orig, synth, vars = c("age", "gender"))
#' print(result)
#' summary(result)
#' plot(result)
chisq_utility <- function(X, ...) {
  UseMethod("chisq_utility")
}

#' @rdname chisq_utility
#' @export
chisq_utility.synth_pair <- function(X, vars = NULL, max_cells = 10000, ...) {
  if (is.null(vars)) {
    # Use all common categorical variables
    common_vars <- intersect(names(X$original), names(X$synthetic))
    vars <- common_vars[sapply(X$original[, common_vars, drop = FALSE],
                               function(x) is.factor(x) || is.character(x))]
    if (length(vars) == 0) {
      stop("No categorical variables found. Specify 'vars' explicitly.")
    }
  }

  chisq_utility.default(
    X = X$original,
    Y = X$synthetic,
    vars = vars,
    max_cells = max_cells,
    weight_X = X$weight_original,
    weight_Y = X$weight_synthetic,
    ...
  )
}

#' @rdname chisq_utility
#' @export
chisq_utility.default <- function(X,
                                  Y,
                                  vars,
                                  max_cells = 10000,
                                  weight_X = NULL,
                                  weight_Y = NULL,
                                  ...) {

  # Input validation
  if (!is.data.frame(X)) stop("X must be a data frame.")
  if (!is.data.frame(Y)) stop("Y must be a data frame.")

  # Check variables exist
  missing_X <- setdiff(vars, names(X))
  missing_Y <- setdiff(vars, names(Y))
  if (length(missing_X) > 0) {
    stop(paste("Variables missing in X:", paste(missing_X, collapse = ", ")))
  }
  if (length(missing_Y) > 0) {
    stop(paste("Variables missing in Y:", paste(missing_Y, collapse = ", ")))
  }

  # Extract weights if specified
  wt_X <- NULL
  wt_Y <- NULL

  if (!is.null(weight_X)) {
    if (is.character(weight_X) && length(weight_X) == 1) {
      if (!weight_X %in% names(X)) stop(paste("Weight variable", weight_X, "not found in X"))
      wt_X <- X[[weight_X]]
    } else if (is.numeric(weight_X)) {
      if (length(weight_X) != nrow(X)) stop("weight_X must have same length as nrow(X)")
      wt_X <- weight_X
    }
  }

  if (!is.null(weight_Y)) {
    if (is.character(weight_Y) && length(weight_Y) == 1) {
      if (!weight_Y %in% names(Y)) stop(paste("Weight variable", weight_Y, "not found in Y"))
      wt_Y <- Y[[weight_Y]]
    } else if (is.numeric(weight_Y)) {
      if (length(weight_Y) != nrow(Y)) stop("weight_Y must have same length as nrow(Y)")
      wt_Y <- weight_Y
    }
  }

  # Subset to vars
  X_sub <- X[, vars, drop = FALSE]
  Y_sub <- Y[, vars, drop = FALSE]

  # Convert to factors with consistent levels
  for (v in vars) {
    all_levels <- union(unique(as.character(X_sub[[v]])),
                        unique(as.character(Y_sub[[v]])))
    X_sub[[v]] <- factor(X_sub[[v]], levels = all_levels)
    Y_sub[[v]] <- factor(Y_sub[[v]], levels = all_levels)
  }

  # Check cell count
  n_cells_possible <- prod(sapply(X_sub, nlevels))
  if (n_cells_possible > max_cells) {
    warning(sprintf("Cross-tabulation would have %d cells (max: %d). Consider fewer variables.",
                    n_cells_possible, max_cells))
  }

  # Create frequency tables
  if (is.null(wt_X)) {
    tab_X <- table(X_sub)
  } else {
    # Weighted table using xtabs
    X_sub$..wt.. <- wt_X
    formula_str <- paste("..wt.. ~", paste(vars, collapse = " + "))
    tab_X <- xtabs(as.formula(formula_str), data = X_sub)
  }

  if (is.null(wt_Y)) {
    tab_Y <- table(Y_sub)
  } else {
    Y_sub$..wt.. <- wt_Y
    formula_str <- paste("..wt.. ~", paste(vars, collapse = " + "))
    tab_Y <- xtabs(as.formula(formula_str), data = Y_sub)
  }

  # Flatten to vectors
  obs_X <- as.vector(tab_X)
  obs_Y <- as.vector(tab_Y)

  # Convert to proportions
  prop_X <- obs_X / sum(obs_X)
  prop_Y <- obs_Y / sum(obs_Y)

  # Expected counts (original proportions scaled to synthetic N)
  N_Y <- sum(obs_Y)
  expected <- prop_X * N_Y

  # Avoid division by zero
  eps <- 1e-10
  expected[expected < eps] <- eps

  # Number of cells
  n_cells <- length(obs_X)
  n_empty_orig <- sum(obs_X == 0)
  n_empty_synth <- sum(obs_Y == 0)

  # Degrees of freedom (cells with non-zero expected - 1)
  df <- sum(prop_X > 0) - 1

  # Chi-square statistic
  chi2 <- sum((obs_Y - expected)^2 / expected)

  # p-value
  p_value <- pchisq(chi2, df = df, lower.tail = FALSE)

  # Voas-Williamson: (chi2 - df) / N
  VW <- (chi2 - df) / N_Y

  # Freeman-Tukey
  FT <- 4 * sum((sqrt(obs_Y) - sqrt(expected))^2)

  # G-test (log-likelihood ratio)
  # Handle zeros: only include cells where obs_Y > 0
  mask_G <- obs_Y > 0
  G <- 2 * sum(obs_Y[mask_G] * log(obs_Y[mask_G] / expected[mask_G]))

  # Jensen-Shannon Divergence
  # Add small epsilon for numerical stability
  p <- prop_X + eps
  q <- prop_Y + eps
  p <- p / sum(p)
  q <- q / sum(q)
  m <- (p + q) / 2

  kl_pm <- sum(p * log(p / m))
  kl_qm <- sum(q * log(q / m))
  JSD <- (kl_pm + kl_qm) / 2

  # Utility percentage (inverse of VW, capped)
  # VW of 0 = 100% utility; use exponential decay
  pct_utility <- 100 * exp(-abs(VW) * 10)

  results <- list(
    chi2 = chi2,
    df = df,
    p_value = p_value,
    VW = VW,
    FT = FT,
    G = G,
    JSD = JSD,
    n_cells = n_cells,
    n_empty_orig = n_empty_orig,
    n_empty_synth = n_empty_synth,
    pct_utility = pct_utility,
    n_orig = nrow(X),
    n_synth = nrow(Y),
    vars = vars,
    table_orig = tab_X,
    table_synth = tab_Y,
    prop_orig = prop_X,
    prop_synth = prop_Y
  )

  class(results) <- "chisq_utility"
  return(results)
}


#' Print method for chisq_utility objects
#' @param x an object of class "chisq_utility"
#' @param ... additional arguments (ignored)
#' @export
print.chisq_utility <- function(x, ...) {
  cat("Chi-Square Utility Assessment\n")
  cat("=============================\n\n")

  cat("Variables:", paste(x$vars, collapse = ", "), "\n")
  cat("Cells:", x$n_cells, "(", x$n_empty_orig, "empty in original,",
      x$n_empty_synth, "empty in synthetic )\n\n")

  cat("Utility Measures:\n")
  cat("  Chi-square:", sprintf("%.2f", x$chi2), "(df =", x$df,
      ", p =", format.pval(x$p_value), ")\n")
  cat("  Voas-Williamson (VW):", sprintf("%.4f", x$VW), "\n")
  cat("  Freeman-Tukey (FT):", sprintf("%.2f", x$FT), "\n")
  cat("  G-test:", sprintf("%.2f", x$G), "\n")
  cat("  Jensen-Shannon Divergence:", sprintf("%.4f", x$JSD), "\n\n")

  cat("Interpretation:\n")
  if (x$VW < 0.005) {
    cat("  Excellent utility - synthetic closely matches original.\n")
  } else if (x$VW < 0.01) {
    cat("  Good utility - minor distributional differences.\n")
  } else if (x$VW < 0.05) {
    cat("  Moderate utility - noticeable differences exist.\n")
  } else {
    cat("  Poor utility - significant distributional differences.\n")
  }

  invisible(x)
}


#' Summary method for chisq_utility objects
#' @param object an object of class "chisq_utility"
#' @param ... additional arguments (ignored)
#' @export
summary.chisq_utility <- function(object, ...) {
  # Calculate cell-level differences
  diff_abs <- abs(object$prop_orig - object$prop_synth)
  diff_rel <- ifelse(object$prop_orig > 0,
                     abs(object$prop_orig - object$prop_synth) / object$prop_orig,
                     NA)

  summ <- list(
    chi2 = object$chi2,
    df = object$df,
    p_value = object$p_value,
    VW = object$VW,
    FT = object$FT,
    G = object$G,
    JSD = object$JSD,
    n_cells = object$n_cells,
    vars = object$vars,
    diff_summary = data.frame(
      metric = c("Mean abs diff", "Max abs diff", "Mean rel diff", "Max rel diff"),
      value = c(mean(diff_abs), max(diff_abs),
                mean(diff_rel, na.rm = TRUE), max(diff_rel, na.rm = TRUE))
    ),
    pct_utility = object$pct_utility
  )

  class(summ) <- "summary.chisq_utility"
  return(summ)
}


#' Print method for summary.chisq_utility objects
#' @param x an object of class "summary.chisq_utility"
#' @param ... additional arguments (ignored)
#' @export
print.summary.chisq_utility <- function(x, ...) {
  cat("Summary: Chi-Square Utility Assessment\n")
  cat("======================================\n\n")

  cat("Variables:", paste(x$vars, collapse = ", "), "\n")
  cat("Number of cells:", x$n_cells, "\n\n")

  cat("Test Statistics:\n")
  cat(sprintf("  Chi-square: %.2f (df=%d, p=%s)\n",
              x$chi2, x$df, format.pval(x$p_value)))
  cat(sprintf("  Voas-Williamson: %.5f\n", x$VW))
  cat(sprintf("  Freeman-Tukey: %.2f\n", x$FT))
  cat(sprintf("  G-test: %.2f\n", x$G))
  cat(sprintf("  JSD: %.5f\n\n", x$JSD))

  cat("Cell-Level Differences:\n")
  print(x$diff_summary, row.names = FALSE)
  cat("\n")

  cat(sprintf("Estimated Utility: %.1f%%\n", x$pct_utility))

  invisible(x)
}


#' Plot method for chisq_utility objects
#' @param x an object of class "chisq_utility"
#' @param y not used
#' @param ... additional arguments passed to plotting functions
#' @param which integer, which plot: 1 = proportion comparison, 2 = residuals
#' @importFrom graphics barplot plot abline par legend text
#' @export
plot.chisq_utility <- function(x, y = NULL, ..., which = 1) {
  show <- rep(FALSE, 2)
  show[which] <- TRUE

  if (sum(show) > 1) {
    op <- par(mfrow = c(1, 2))
    on.exit(par(op))
  }

  if (show[1]) {
    # Scatter plot of proportions
    max_prop <- max(c(x$prop_orig, x$prop_synth))

    plot(x$prop_orig, x$prop_synth,
         main = "Cell Proportions: Original vs Synthetic",
         xlab = "Original Proportion",
         ylab = "Synthetic Proportion",
         pch = 19, col = rgb(0.2, 0.4, 0.8, 0.6),
         xlim = c(0, max_prop), ylim = c(0, max_prop), ...)
    abline(0, 1, col = "red", lty = 2, lwd = 2)
    legend("topleft",
           legend = c("Perfect match", sprintf("VW = %.4f", x$VW)),
           lty = c(2, NA), col = c("red", NA), lwd = 2, cex = 0.8)
  }

  if (show[2]) {
    # Standardized residuals
    expected <- x$prop_orig * x$n_synth
    expected[expected < 0.5] <- 0.5
    observed <- x$prop_synth * x$n_synth
    std_resid <- (observed - expected) / sqrt(expected)

    # Sort by absolute value for display
    ord <- order(abs(std_resid), decreasing = TRUE)
    top_n <- min(20, length(std_resid))

    barplot(std_resid[ord[1:top_n]],
            main = "Standardized Residuals (Top 20 Cells)",
            ylab = "Standardized Residual",
            col = ifelse(std_resid[ord[1:top_n]] > 0, "steelblue", "firebrick"),
            las = 2, ...)
    abline(h = c(-2, 2), col = "orange", lty = 2)
  }
}
