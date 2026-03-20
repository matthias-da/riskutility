#' Regression Fidelity — Coefficient Comparison
#'
#' Fits the same regression model on original and synthetic data, then compares
#' coefficients to assess specific-analysis fidelity. The utility score is the
#' mean confidence interval overlap across all coefficients.
#'
#' @param X A data.frame containing the original dataset.
#' @param Y A data.frame containing the synthetic/anonymized dataset.
#' @param formula An R formula specifying the regression model. Required;
#'   an error is raised if not provided. For example, \code{y ~ x1 + x2}.
#' @param model Character, the type of model to fit. One of \code{"lm"}
#'   (default) or \code{"glm"}.
#' @param family A family object for GLM fitting (e.g., \code{binomial()}).
#'   Only used when \code{model = "glm"}. Default \code{binomial()}.
#' @param conf_level Numeric, the confidence level for confidence intervals.
#'   Default 0.95.
#' @param na.rm Logical, whether to remove rows with NA values in the
#'   variables used by the formula. Default TRUE.
#' @param ... Additional arguments passed to the model fitting function
#'   (\code{\link{lm}} or \code{\link{glm}}).
#'
#' @return An object of class \code{"regression_fidelity"} containing:
#' \itemize{
#'   \item \code{utility_score}: mean CI overlap across coefficients, in \eqn{[0, 1]}
#'   \item \code{coefficients}: data.frame with columns: term, estimate_orig,
#'     estimate_synth, se_orig, se_synth, bias, std_bias, ci_overlap, sig_orig,
#'     sig_synth, sig_agreement
#'   \item \code{mean_ci_overlap}: same as utility_score
#'   \item \code{mean_abs_std_bias}: mean absolute standardized bias
#'   \item \code{sig_agreement_rate}: proportion of coefficients with matching
#'     significance
#'   \item \code{formula}: the formula used
#'   \item \code{model}: the model type used
#'   \item \code{conf_level}: the confidence level used
#'   \item \code{n_X}: number of observations in original data
#'   \item \code{n_Y}: number of observations in synthetic data
#'   \item \code{n_coef}: number of coefficients compared
#' }
#'
#' @details
#' For each coefficient, the following quantities are computed:
#' \itemize{
#'   \item \strong{Bias}: difference between synthetic and original estimates
#'     (\code{synth - orig})
#'   \item \strong{Standardized bias}: bias divided by the original standard
#'     error (\code{bias / SE_orig})
#'   \item \strong{CI overlap}: overlap length divided by the average of the
#'     two interval widths, clamped to \eqn{[0, 1]}
#'   \item \strong{Significance agreement}: whether both models agree on
#'     statistical significance at the given confidence level
#' }
#'
#' The CI overlap for two intervals \eqn{[lo_1, hi_1]} and \eqn{[lo_2, hi_2]}
#' is:
#' \deqn{overlap = \frac{\max(0, \min(hi_1, hi_2) - \max(lo_1, lo_2))}
#'   {0.5 \cdot (w_1 + w_2)}}
#' where \eqn{w_k = hi_k - lo_k}. The result is clamped to \eqn{[0, 1]}.
#'
#' For \code{model = "lm"}, confidence intervals use the t-distribution with
#' the appropriate residual degrees of freedom. For \code{model = "glm"},
#' Wald intervals based on the normal distribution are used.
#'
#' A utility score of 1 indicates perfect overlap (coefficients are identical),
#' while 0 indicates no overlap at all. The \code{utility_score} field provides
#' a standard interface compatible with \code{\link{subgroup_utility}} and
#' \code{\link{rumap}}.
#'
#' \strong{Interpretation (heuristic thresholds):}
#' \itemize{
#'   \item utility_score > 0.9: EXCELLENT -- regression results very well preserved
#'   \item utility_score > 0.7: GOOD -- reasonably preserved
#'   \item utility_score > 0.4: MODERATE -- some differences
#'   \item utility_score <= 0.4: POOR -- significant differences
#' }
#'
#' @seealso \code{\link{propscore}} for propensity score utility,
#'   \code{\link{compare_model_performance}} for predictive performance comparison,
#'   \code{\link{contingency_fidelity}} for categorical dependence comparison,
#'   \code{\link{subgroup_utility}} for stratified utility assessment
#'
#' @references
#' Karr, A. F., Kohnen, C. N., Oganian, A., Reiter, J. P., and Sanil, A. P.
#' (2006). A Framework for Evaluating the Utility of Data Altered to Protect
#' Confidentiality. The American Statistician, 60(3), 224-232.
#'
#' Snoke, J., Raab, G. M., Nowok, B., Dibben, C., and Slavkovic, A. (2018).
#' General and Specific Utility Measures for Synthetic Data. Journal of the
#' Royal Statistical Society: Series A, 181(3), 663-688.
#'
#' @author Matthias Templ
#' @family utility
#' @export
#' @importFrom stats lm glm confint coef summary.lm vcov binomial qnorm qt
#'
#' @examples
#' set.seed(123)
#' n <- 200
#' x1 <- rnorm(n)
#' x2 <- rnorm(n)
#' y <- 2 + 3 * x1 - 1.5 * x2 + rnorm(n)
#' X <- data.frame(y = y, x1 = x1, x2 = x2)
#'
#' # Good synthetic data (similar DGP)
#' x1s <- rnorm(n)
#' x2s <- rnorm(n)
#' ys <- 2 + 3 * x1s - 1.5 * x2s + rnorm(n)
#' Y_good <- data.frame(y = ys, x1 = x1s, x2 = x2s)
#'
#' result <- regression_fidelity(X, Y_good, formula = y ~ x1 + x2)
#' print(result)
#' summary(result)
#'
#' # Using synth_pair
#' pair <- synth_pair(X, Y_good)
#' result2 <- regression_fidelity(pair, formula = y ~ x1 + x2)
#'
#' \donttest{
#' # Poor synthetic data (wrong coefficients)
#' ys_bad <- 10 + 0.1 * x1s + 5 * x2s + rnorm(n)
#' Y_bad <- data.frame(y = ys_bad, x1 = x1s, x2 = x2s)
#' result_bad <- regression_fidelity(X, Y_bad, formula = y ~ x1 + x2)
#' print(result_bad)
#'
#' # GLM example
#' y_bin <- rbinom(n, 1, plogis(0.5 + 1.5 * x1))
#' X_bin <- data.frame(y = y_bin, x1 = x1, x2 = x2)
#' y_bin_s <- rbinom(n, 1, plogis(0.5 + 1.5 * x1s))
#' Y_bin <- data.frame(y = y_bin_s, x1 = x1s, x2 = x2s)
#' result_glm <- regression_fidelity(X_bin, Y_bin,
#'   formula = y ~ x1 + x2, model = "glm")
#' print(result_glm)
#' }
regression_fidelity <- function(X, ...) {
  UseMethod("regression_fidelity")
}

#' @rdname regression_fidelity
#' @export
regression_fidelity.synth_pair <- function(X, ...) {
  regression_fidelity.default(
    X = X$original,
    Y = X$synthetic,
    ...
  )
}

#' @rdname regression_fidelity
#' @export
regression_fidelity.default <- function(X, Y,
                                        formula = NULL,
                                        model = c("lm", "glm"),
                                        family = binomial(),
                                        conf_level = 0.95,
                                        na.rm = TRUE,
                                        ...) {

  model <- match.arg(model)

  # Validate formula

  if (is.null(formula)) {
    stop("'formula' is required. Specify the regression model, e.g., y ~ x1 + x2.")
  }
  if (!inherits(formula, "formula")) {
    stop("'formula' must be an R formula object.")
  }

  # Convert to data.frame if needed
  X <- as.data.frame(X)
  Y <- as.data.frame(Y)

  # Extract variable names from formula
  formula_vars <- all.vars(formula)
  missing_X <- setdiff(formula_vars, names(X))
  missing_Y <- setdiff(formula_vars, names(Y))
  if (length(missing_X) > 0) {
    stop(paste("Variables missing in X:", paste(missing_X, collapse = ", ")))
  }
  if (length(missing_Y) > 0) {
    stop(paste("Variables missing in Y:", paste(missing_Y, collapse = ", ")))
  }

  # Handle missing values
  if (na.rm) {
    X <- X[complete.cases(X[, formula_vars, drop = FALSE]), , drop = FALSE]
    Y <- Y[complete.cases(Y[, formula_vars, drop = FALSE]), , drop = FALSE]
  }

  if (nrow(X) == 0) stop("No complete cases in X after removing NAs.")
  if (nrow(Y) == 0) stop("No complete cases in Y after removing NAs.")

  n_X <- nrow(X)
  n_Y <- nrow(Y)

  # Fit models
  if (model == "lm") {
    fit_orig  <- lm(formula, data = X, ...)
    fit_synth <- lm(formula, data = Y, ...)
  } else {
    fit_orig  <- glm(formula, data = X, family = family, ...)
    fit_synth <- glm(formula, data = Y, family = family, ...)
  }

  # Extract coefficients and standard errors
  summ_orig  <- summary(fit_orig)
  summ_synth <- summary(fit_synth)

  coef_orig  <- coef(summ_orig)
  coef_synth <- coef(summ_synth)

  # Align terms — use the intersection of terms present in both models
  common_terms <- intersect(rownames(coef_orig), rownames(coef_synth))
  if (length(common_terms) == 0) {
    stop("No common coefficient terms between the two models.")
  }

  coef_orig  <- coef_orig[common_terms, , drop = FALSE]
  coef_synth <- coef_synth[common_terms, , drop = FALSE]

  # Extract estimates and SEs
  est_orig  <- coef_orig[, 1]
  est_synth <- coef_synth[, 1]
  se_orig   <- coef_orig[, 2]
  se_synth  <- coef_synth[, 2]

  # Compute CIs (t-distribution for lm, normal for glm)
  alpha <- 1 - conf_level
  if (model == "lm") {
    q_orig  <- qt(1 - alpha / 2, df = fit_orig$df.residual)
    q_synth <- qt(1 - alpha / 2, df = fit_synth$df.residual)
  } else {
    q_orig  <- qnorm(1 - alpha / 2)
    q_synth <- q_orig
  }

  lo_orig  <- est_orig  - q_orig  * se_orig
  hi_orig  <- est_orig  + q_orig  * se_orig
  lo_synth <- est_synth - q_synth * se_synth
  hi_synth <- est_synth + q_synth * se_synth

  # CI overlap: overlap length / average width, clamped to [0, 1]
  width_orig  <- hi_orig  - lo_orig
  width_synth <- hi_synth - lo_synth
  avg_width   <- 0.5 * (width_orig + width_synth)
  overlap_num <- pmax(0, pmin(hi_orig, hi_synth) - pmax(lo_orig, lo_synth))
  ci_overlap  <- ifelse(avg_width > 0, pmin(overlap_num / avg_width, 1), 1)

  # Bias and standardized bias
  bias <- est_synth - est_orig
  std_bias <- ifelse(se_orig > 0, bias / se_orig, NA_real_)

  # Significance: p-value < alpha
  # p-values are in column 4 of the summary coefficient table
  pval_orig  <- coef_orig[, 4]
  pval_synth <- coef_synth[, 4]
  sig_orig  <- pval_orig  < alpha
  sig_synth <- pval_synth < alpha
  sig_agreement <- sig_orig == sig_synth

  # Build coefficients data.frame
  coef_df <- data.frame(
    term           = common_terms,
    estimate_orig  = est_orig,
    estimate_synth = est_synth,
    se_orig        = se_orig,
    se_synth       = se_synth,
    bias           = bias,
    std_bias       = std_bias,
    ci_overlap     = ci_overlap,
    sig_orig       = sig_orig,
    sig_synth      = sig_synth,
    sig_agreement  = sig_agreement,
    row.names = NULL
  )

  # Summary statistics
  mean_ci_overlap    <- mean(ci_overlap)
  mean_abs_std_bias  <- mean(abs(std_bias), na.rm = TRUE)
  sig_agreement_rate <- mean(sig_agreement)

  result <- list(
    utility_score      = mean_ci_overlap,
    coefficients       = coef_df,
    mean_ci_overlap    = mean_ci_overlap,
    mean_abs_std_bias  = mean_abs_std_bias,
    sig_agreement_rate = sig_agreement_rate,
    formula            = formula,
    model              = model,
    conf_level         = conf_level,
    n_X                = n_X,
    n_Y                = n_Y,
    n_coef             = length(common_terms)
  )

  class(result) <- "regression_fidelity"
  return(result)
}


#' Print method for regression_fidelity objects
#'
#' @param x an object of class \code{"regression_fidelity"}
#' @param ... additional arguments (ignored)
#' @export
print.regression_fidelity <- function(x, ...) {
  cat("Regression Fidelity - Coefficient Comparison\n")
  cat("=============================================\n\n")

  cat("Formula:", deparse(x$formula), "\n")
  cat("Model:  ", x$model, "\n")
  cat("Conf. level:", x$conf_level, "\n\n")

  cat("Dataset Sizes:\n")
  cat("  Original (X):", x$n_X, "\n")
  cat("  Synthetic (Y):", x$n_Y, "\n")
  cat("  Coefficients:", x$n_coef, "\n\n")

  cat("Results:\n")
  cat("  Utility score (mean CI overlap):", sprintf("%.4f", x$utility_score), "\n")
  cat("  Mean |std. bias|:               ", sprintf("%.4f", x$mean_abs_std_bias), "\n")
  cat("  Significance agreement rate:    ", sprintf("%.4f", x$sig_agreement_rate),
      sprintf(" (%d/%d)", sum(x$coefficients$sig_agreement), x$n_coef), "\n\n")

  # Show disagreements
  disagree <- x$coefficients[!x$coefficients$sig_agreement, ]
  if (nrow(disagree) > 0) {
    cat("Significance Disagreements:\n")
    for (i in seq_len(nrow(disagree))) {
      cat("  ", disagree$term[i], ": orig",
          ifelse(disagree$sig_orig[i], "significant", "not significant"),
          ", synth",
          ifelse(disagree$sig_synth[i], "significant", "not significant"),
          "\n")
    }
    cat("\n")
  }

  cat("Interpretation:\n")
  if (x$utility_score > 0.9) {
    cat("  EXCELLENT: Regression results are very well preserved.\n")
  } else if (x$utility_score > 0.7) {
    cat("  GOOD: Regression results are reasonably preserved.\n")
  } else if (x$utility_score > 0.4) {
    cat("  MODERATE: Some differences in regression coefficients.\n")
  } else {
    cat("  POOR: Significant differences in regression coefficients.\n")
  }

  invisible(x)
}


#' Summary method for regression_fidelity objects
#'
#' @param object an object of class \code{"regression_fidelity"}
#' @param ... additional arguments (ignored)
#' @export
summary.regression_fidelity <- function(object, ...) {
  summ <- list(
    utility_score      = object$utility_score,
    coefficients       = object$coefficients,
    mean_ci_overlap    = object$mean_ci_overlap,
    mean_abs_std_bias  = object$mean_abs_std_bias,
    sig_agreement_rate = object$sig_agreement_rate,
    formula            = object$formula,
    model              = object$model,
    conf_level         = object$conf_level,
    n_X                = object$n_X,
    n_Y                = object$n_Y,
    n_coef             = object$n_coef
  )

  class(summ) <- "summary.regression_fidelity"
  return(summ)
}


#' Print method for summary.regression_fidelity objects
#'
#' @param x an object of class \code{"summary.regression_fidelity"}
#' @param ... additional arguments (ignored)
#' @export
print.summary.regression_fidelity <- function(x, ...) {
  cat("Summary: Regression Fidelity\n")
  cat("============================\n\n")

  cat("Formula:", deparse(x$formula), "\n")
  cat("Model:  ", x$model, " | Conf. level:", x$conf_level, "\n")
  cat("Samples: X =", x$n_X, ", Y =", x$n_Y, "\n\n")

  cat("Coefficient Comparison:\n")
  # Print formatted table
  df <- x$coefficients
  fmt <- data.frame(
    Term       = df$term,
    Est.Orig   = sprintf("%.4f", df$estimate_orig),
    Est.Synth  = sprintf("%.4f", df$estimate_synth),
    Bias       = sprintf("%.4f", df$bias),
    Std.Bias   = sprintf("%.4f", df$std_bias),
    CI.Overlap = sprintf("%.4f", df$ci_overlap),
    Sig.Agree  = ifelse(df$sig_agreement, "yes", "NO")
  )
  print(fmt, row.names = FALSE, right = FALSE)

  cat("\nSummary Statistics:\n")
  cat("  Utility score (mean CI overlap):", sprintf("%.4f", x$utility_score), "\n")
  cat("  Mean |standardized bias|:       ", sprintf("%.4f", x$mean_abs_std_bias), "\n")
  cat("  Significance agreement rate:    ", sprintf("%.1f%%", x$sig_agreement_rate * 100), "\n")

  invisible(x)
}


#' Plot method for regression_fidelity objects
#'
#' @param x an object of class \code{"regression_fidelity"}
#' @param y not used
#' @param ... additional arguments passed to plotting functions
#' @param which integer, which plot(s) to produce:
#'   1 = forest plot comparing original and synthetic coefficients with CIs,
#'   2 = CI overlap bar chart per coefficient
#' @importFrom graphics abline arrows axis barplot legend mtext par plot points segments
#' @export
plot.regression_fidelity <- function(x, y = NULL, ..., which = 1) {
  show <- rep(FALSE, 2)
  show[which] <- TRUE

  n_plots <- sum(show)
  if (n_plots > 1) {
    op <- par(mfrow = c(1, 2))
    on.exit(par(op))
  }

  df <- x$coefficients
  n <- nrow(df)
  z <- qnorm(1 - (1 - x$conf_level) / 2)

  if (show[1]) {
    # Forest plot: original and synthetic coefficients with CIs
    lo_orig  <- df$estimate_orig  - z * df$se_orig
    hi_orig  <- df$estimate_orig  + z * df$se_orig
    lo_synth <- df$estimate_synth - z * df$se_synth
    hi_synth <- df$estimate_synth + z * df$se_synth

    all_vals <- c(lo_orig, hi_orig, lo_synth, hi_synth)
    xlim <- range(all_vals, na.rm = TRUE)
    xlim <- xlim + diff(xlim) * c(-0.1, 0.1)

    # Y positions: two rows per coefficient (offset for orig vs synth)
    offset <- 0.15
    ypos_orig  <- seq_len(n) + offset
    ypos_synth <- seq_len(n) - offset

    plot(NA, xlim = xlim, ylim = c(0.5, n + 0.5),
         xlab = "Estimate", ylab = "",
         yaxt = "n", main = "Regression Coefficients",
         ...)
    axis(2, at = seq_len(n), labels = df$term, las = 1)
    abline(v = 0, lty = 2, col = "grey50")

    # Original CIs
    segments(lo_orig, ypos_orig, hi_orig, ypos_orig,
             col = "steelblue", lwd = 2)
    points(df$estimate_orig, ypos_orig, pch = 19, col = "steelblue", cex = 1.2)

    # Synthetic CIs
    segments(lo_synth, ypos_synth, hi_synth, ypos_synth,
             col = "coral", lwd = 2)
    points(df$estimate_synth, ypos_synth, pch = 17, col = "coral", cex = 1.2)

    legend("topright",
           legend = c("Original", "Synthetic"),
           col = c("steelblue", "coral"),
           pch = c(19, 17),
           lty = 1, lwd = 2,
           bty = "n")
  }

  if (show[2]) {
    # CI overlap bar chart
    cols <- ifelse(df$ci_overlap > 0.5, "steelblue", "coral")
    bp <- barplot(df$ci_overlap,
                  names.arg = df$term,
                  main = "CI Overlap per Coefficient",
                  ylab = "CI Overlap",
                  ylim = c(0, 1),
                  col = cols,
                  las = 2,
                  ...)
    abline(h = 0.5, lty = 2, col = "grey50")
    mtext(sprintf("Mean CI overlap = %.3f", x$mean_ci_overlap),
          side = 3, line = 0, cex = 0.9)
  }

  invisible(x)
}
