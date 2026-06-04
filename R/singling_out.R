#' Singling Out Risk
#'
#' Computes the Singling Out Risk for synthetic data, one of three explicit
#' GDPR anonymization failure criteria (Article 29 Working Party). An attacker
#' uses synthetic data to craft logical predicates (conjunctions of column
#' conditions) and tests whether these predicates uniquely identify exactly one
#' record in the original data.
#'
#' @param X data frame of original/training data
#' @param Y data frame of synthetic data
#' @param holdout data frame of holdout data (optional). If NULL, a holdout set
#'   is automatically created by splitting X.
#' @param holdout_fraction numeric, fraction of X to use as holdout if holdout
#'   is NULL (default: 0.5)
#' @param n_attacks integer, number of predicate attacks to generate (default: 2000)
#' @param n_cols integer, number of columns to use per predicate in multivariate
#'   mode (default: 3)
#' @param mode character, attack mode: "multivariate" (default) builds predicates
#'   from random column subsets; "univariate" targets single columns
#' @param vars character vector of variable names to use. If NULL (default),
#'   all common variables between X and Y are used.
#' @param na.rm logical, remove records with NA values (default: TRUE)
#' @param seed integer, random seed for reproducibility (default: NULL)
#' @param confidence_level numeric, confidence level for Wilson score intervals
#'   (default: 0.95)
#' @param ... additional arguments passed to methods (currently unused)
#'
#' @return An object of class "singling_out" containing:
#' \itemize{
#'   \item risk: residual risk score (or absolute if no holdout), bounded between 0 and 1
#'   \item risk_ci: confidence interval for risk
#'   \item risk_attack: attack success rate (fraction singling out in original)
#'   \item risk_attack_ci: Wilson CI for attack success rate
#'   \item risk_control: control success rate (fraction singling out in holdout, NA if none)
#'   \item risk_control_ci: Wilson CI for control success rate (NA if none)
#'   \item n_attacks: number of predicates generated
#'   \item n_success: number singling out in original
#'   \item n_control_success: number singling out in holdout
#'   \item match_counts: integer vector of match counts per predicate in original
#'   \item match_counts_control: integer vector of match counts per predicate in holdout
#'   \item privacy_pass: logical, risk <= 0.1
#'   \item n_original: number of original records (training portion)
#'   \item n_synthetic: number of synthetic records
#'   \item n_holdout: number of holdout records
#'   \item mode: attack mode used
#'   \item n_cols: number of columns per predicate
#'   \item vars: variables used
#'   \item confidence_level: confidence level used
#' }
#'
#' @details
#' The singling out attack works by sampling random records from the synthetic
#' data and constructing logical predicates from their column values. For each
#' predicate, the algorithm checks how many records in the original data match.
#' A predicate that matches exactly one record constitutes a successful singling
#' out attack.
#'
#' In multivariate mode (default), each predicate is a conjunction of conditions
#' on \code{n_cols} randomly selected columns:
#' \itemize{
#'   \item Categorical/factor columns: exact match (\code{col == val})
#'   \item Numeric columns: inequality based on median (\code{col >= val} if
#'     val > median, \code{col <= val} otherwise), targeting distribution tails
#' }
#'
#' In univariate mode, each predicate tests a single column condition.
#'
#' The risk score uses Wilson score intervals for robust estimation:
#' \deqn{r_{attack} = n_{success} / n_{attacks}}
#' \deqn{r_{control} = n_{control\_success} / n_{attacks}}
#' \deqn{risk = (r_{attack} - r_{control}) / (1 - r_{control})}
#'
#' Without holdout data, the absolute attack rate is reported as the risk.
#'
#' @section Holdout splitting:
#' When no external holdout is provided, the original data is split internally.
#' The holdout serves as a control: predicates that single out records in the
#' holdout (which were not used for synthesis) represent baseline singling out
#' risk due to data structure rather than information leakage. The residual
#' risk subtracts this baseline.
#'
#' @seealso \code{\link{dcap}} for differential correct attribution probability,
#'   \code{\link{tcap}} for targeted CAP,
#'   \code{\link{nnaa}} for nearest-neighbor adversarial accuracy
#'
#' @references
#' Giomi, M., Boenisch, F., Wehmeyer, C. & Tasnadi, B. (2023).
#' A Unified Framework for Quantifying Privacy Risk in Synthetic Data.
#' \emph{Proceedings on Privacy Enhancing Technologies (PoPETs)}, 2023(2), 312--328.
#' \doi{10.56553/popets-2023-0055}
#'
#' Cohen, A. & Nissim, K. (2020).
#' Towards Formalizing the GDPR's Notion of Singling Out.
#' \emph{Proceedings of the National Academy of Sciences}, 117(15), 8344--8352.
#' \doi{10.1073/pnas.1914598117}
#'
#' Article 29 Data Protection Working Party (2014).
#' Opinion 05/2014 on Anonymisation Techniques. WP216.
#'
#' @author Matthias Templ
#' @family privacy-models
#' @export
#' @importFrom stats complete.cases median qnorm
#' @importFrom graphics hist abline legend barplot par text
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
#' # Good synthetic data (independent, no memorization)
#' Y <- data.frame(
#'   age = sample(20:70, n, replace = TRUE),
#'   income = rnorm(n, 50000, 15000),
#'   gender = sample(c("M", "F"), n, replace = TRUE),
#'   region = sample(c("N", "S", "E", "W"), n, replace = TRUE)
#' )
#'
#' result <- singling_out(X, Y, n_attacks = 500, seed = 42)
#' print(result)
#' summary(result)
#'
#' \donttest{
#' # Memorized data (Y is copy of X) - should show high risk
#' Y_copy <- X[sample(nrow(X), n, replace = TRUE), ]
#' result_bad <- singling_out(X, Y_copy, n_attacks = 500, seed = 42)
#' print(result_bad)
#' }
singling_out <- function(X, ...) {
  UseMethod("singling_out")
}

#' @rdname singling_out
#' @export
singling_out.synth_pair <- function(X, ...) {
  singling_out.default(
    X = X$original,
    Y = X$synthetic,
    holdout = X$holdout,
    vars = X$vars,
    ...
  )
}

#' @rdname singling_out
#' @export
singling_out.default <- function(X, Y,
                                  holdout = NULL,
                                  holdout_fraction = 0.5,
                                  n_attacks = 2000,
                                  n_cols = 3,
                                  mode = c("multivariate", "univariate"),
                                  vars = NULL,
                                  na.rm = TRUE,
                                  seed = NULL,
                                  confidence_level = 0.95,
                                  ...) {

  mode <- match.arg(mode)

  if (!is.numeric(n_attacks) || length(n_attacks) != 1 || n_attacks < 1) {
    stop("n_attacks must be a positive integer.")
  }
  n_attacks <- as.integer(n_attacks)

  if (!is.numeric(n_cols) || length(n_cols) != 1 || n_cols < 1) {
    stop("n_cols must be a positive integer.")
  }
  n_cols <- as.integer(n_cols)

  if (!is.numeric(confidence_level) || length(confidence_level) != 1 ||
      confidence_level <= 0 || confidence_level >= 1) {
    stop("confidence_level must be a number between 0 and 1.")
  }

  # Shared input validation, variable intersection, NA removal, holdout split
  v <- .validate_pair_inputs(X, Y, holdout = holdout,
                             holdout_fraction = holdout_fraction,
                             vars = vars, na.rm = na.rm, seed = seed)
  Y       <- v$Y
  vars    <- v$vars
  train   <- v$train
  holdout <- v$holdout

  # Cap n_cols at number of available variables
  if (n_cols > length(vars)) {
    n_cols <- length(vars)
  }

  n_train <- nrow(train)
  n_synthetic <- nrow(Y)
  n_holdout_final <- nrow(holdout)

  # Determine column types and precompute medians for numeric columns
  col_types <- vapply(Y, function(col) {
    if (is.numeric(col)) "numeric"
    else "categorical"
  }, character(1))

  medians <- setNames(
    vapply(vars, function(v) {
      if (col_types[v] == "numeric") median(Y[[v]], na.rm = TRUE)
      else NA_real_
    }, numeric(1)),
    vars
  )

  # Set seed for reproducible attacks (after holdout splitting)
  if (!is.null(seed)) set.seed(seed + 1L)

  # Generate and evaluate predicates
  match_counts <- integer(n_attacks)
  match_counts_control <- integer(n_attacks)

  if (mode == "multivariate") {
    for (i in seq_len(n_attacks)) {
      # Sample a random record from Y
      rec_idx <- sample(n_synthetic, 1)
      # Sample random columns
      cols <- sample(vars, n_cols)

      # Build and evaluate predicate
      pred <- .build_predicate(Y[rec_idx, , drop = FALSE], cols, medians, col_types)
      match_counts[i] <- .evaluate_predicate(pred, train)
      match_counts_control[i] <- .evaluate_predicate(pred, holdout)
    }
  } else {
    # Univariate mode: one column per predicate
    for (i in seq_len(n_attacks)) {
      rec_idx <- sample(n_synthetic, 1)
      col <- sample(vars, 1)

      pred <- .build_predicate(Y[rec_idx, , drop = FALSE], col, medians, col_types)
      match_counts[i] <- .evaluate_predicate(pred, train)
      match_counts_control[i] <- .evaluate_predicate(pred, holdout)
    }
  }

  # Count successes (exactly 1 match = singling out)
  n_success <- sum(match_counts == 1L)
  n_control_success <- sum(match_counts_control == 1L)

  # Compute residual risk with Wilson CIs
  rr <- .residual_risk(n_success, n_control_success, n_attacks, confidence_level)

  results <- list(
    risk = rr$risk,
    risk_ci = rr$risk_ci,
    risk_attack = rr$risk_attack,
    risk_attack_ci = rr$risk_attack_ci,
    risk_control = rr$risk_control,
    risk_control_ci = rr$risk_control_ci,
    n_attacks = n_attacks,
    n_success = n_success,
    n_control_success = n_control_success,
    match_counts = match_counts,
    match_counts_control = match_counts_control,
    privacy_pass = rr$privacy_pass,
    n_original = n_train,
    n_synthetic = n_synthetic,
    n_holdout = n_holdout_final,
    mode = mode,
    n_cols = n_cols,
    vars = vars,
    confidence_level = confidence_level
  )

  class(results) <- "singling_out"
  return(results)
}

# --- Internal helpers ---

#' Build a predicate from a synthetic record
#'
#' @param record one-row data frame
#' @param cols character vector of column names to include
#' @param medians named numeric vector of medians for numeric columns
#' @param col_types named character vector of column types
#' @return list of conditions, each with elements: col, op, value
#' @keywords internal
.build_predicate <- function(record, cols, medians, col_types) {
  lapply(cols, function(col) {
    val <- record[[col]]
    if (col_types[col] == "numeric") {
      if (val > medians[col]) {
        list(col = col, op = ">=", value = val)
      } else {
        list(col = col, op = "<=", value = val)
      }
    } else {
      list(col = col, op = "==", value = val)
    }
  })
}

#' Evaluate a predicate against a dataset
#'
#' @param predicate list of conditions from .build_predicate
#' @param data data frame to evaluate against
#' @return integer count of matching rows
#' @keywords internal
.evaluate_predicate <- function(predicate, data) {
  matches <- rep(TRUE, nrow(data))
  for (cond in predicate) {
    col_vals <- data[[cond$col]]
    if (cond$op == "==") {
      matches <- matches & (col_vals == cond$value)
    } else if (cond$op == ">=") {
      matches <- matches & (col_vals >= cond$value)
    } else if (cond$op == "<=") {
      matches <- matches & (col_vals <= cond$value)
    }
  }
  sum(matches, na.rm = TRUE)
}

# .wilson_score is defined in R/utils_internal.R

# --- S3 methods ---

#' Print method for singling_out objects
#'
#' @param x an object of class "singling_out"
#' @param ... additional arguments (ignored)
#' @return The input object, invisibly.
#' @export
print.singling_out <- function(x, ...) {
  cat("Singling Out Risk Assessment\n")
  cat("============================\n")
  cat("Mode:", x$mode, "| Columns per predicate:", x$n_cols, "\n")
  cat("Variables used:", length(x$vars), "\n\n")

  cat("Dataset Sizes:\n")
  cat("  Original (training):", x$n_original, "\n")
  cat("  Holdout:", x$n_holdout, "\n")
  cat("  Synthetic:", x$n_synthetic, "\n\n")

  cat("Attack Results (", x$n_attacks, " predicates):\n", sep = "")
  cat("  Singling out in original:", x$n_success,
      sprintf("(%.1f%%)", 100 * x$risk_attack), "\n")
  cat("  Singling out in holdout:", x$n_control_success,
      sprintf("(%.1f%%)", 100 * x$risk_control), "\n\n")

  cat("Risk Score:\n")
  cat("  Residual risk:", round(x$risk, 4),
      sprintf("[%.4f, %.4f]", x$risk_ci["lower"], x$risk_ci["upper"]),
      sprintf("(%d%% CI)\n", round(100 * x$confidence_level)))

  cat("\nPrivacy Assessment:")
  if (x$privacy_pass) {
    cat(" PASS\n")
    cat("  Singling out risk is within acceptable bounds (<= 0.1).\n")
  } else {
    cat(" WARNING\n")
    cat("  Elevated singling out risk detected (> 0.1).\n")
    cat("  Synthetic data may enable identification of individuals.\n")
  }

  invisible(x)
}

#' Summary method for singling_out objects
#'
#' @param object an object of class "singling_out"
#' @param ... additional arguments (ignored)
#' @return A list of summary statistics for the corresponding object.
#' @export
summary.singling_out <- function(object, ...) {
  # Match count distribution
  mc_table <- table(factor(pmin(object$match_counts, 10),
                           levels = 0:10))
  names(mc_table)[11] <- "10+"

  mc_table_ctrl <- table(factor(pmin(object$match_counts_control, 10),
                                levels = 0:10))
  names(mc_table_ctrl)[11] <- "10+"

  summ <- list(
    risk = object$risk,
    risk_ci = object$risk_ci,
    risk_attack = object$risk_attack,
    risk_attack_ci = object$risk_attack_ci,
    risk_control = object$risk_control,
    risk_control_ci = object$risk_control_ci,
    n_attacks = object$n_attacks,
    n_success = object$n_success,
    n_control_success = object$n_control_success,
    match_count_dist = mc_table,
    match_count_dist_control = mc_table_ctrl,
    mean_matches = mean(object$match_counts),
    median_matches = median(object$match_counts),
    mean_matches_control = mean(object$match_counts_control),
    median_matches_control = median(object$match_counts_control),
    privacy_pass = object$privacy_pass,
    n_original = object$n_original,
    n_synthetic = object$n_synthetic,
    n_holdout = object$n_holdout,
    mode = object$mode,
    n_cols = object$n_cols,
    vars = object$vars,
    confidence_level = object$confidence_level
  )

  class(summ) <- "summary.singling_out"
  return(summ)
}

#' Print method for summary.singling_out objects
#'
#' @param x an object of class "summary.singling_out"
#' @param ... additional arguments (ignored)
#' @return The input object, invisibly.
#' @export
print.summary.singling_out <- function(x, ...) {
  cat("Summary: Singling Out Risk Assessment\n")
  cat("=====================================\n")
  cat("Mode:", x$mode, "| Columns per predicate:", x$n_cols, "\n")
  cat("Variables:", length(x$vars), "\n\n")

  cat("Dataset Sizes:\n")
  cat("  Original:", x$n_original, "| Holdout:", x$n_holdout,
      "| Synthetic:", x$n_synthetic, "\n\n")

  cat("Key Metrics:\n")
  cat("  Residual risk:", round(x$risk, 4),
      sprintf("[%.4f, %.4f]", x$risk_ci["lower"], x$risk_ci["upper"]), "\n")
  cat("  Attack rate:", sprintf("%.4f", x$risk_attack),
      sprintf("[%.4f, %.4f]", x$risk_attack_ci["lower"],
              x$risk_attack_ci["upper"]), "\n")
  cat("  Control rate:", sprintf("%.4f", x$risk_control),
      sprintf("[%.4f, %.4f]", x$risk_control_ci["lower"],
              x$risk_control_ci["upper"]), "\n")
  cat("  Privacy:", ifelse(x$privacy_pass, "PASS", "WARNING"), "\n\n")

  cat("Match Count Distribution (original):\n")
  cat("  Mean matches:", round(x$mean_matches, 2),
      "| Median:", x$median_matches, "\n")
  cat("  Counts (0, 1, 2, ...):\n")
  print(x$match_count_dist)
  cat("\n")

  cat("Match Count Distribution (holdout):\n")
  cat("  Mean matches:", round(x$mean_matches_control, 2),
      "| Median:", x$median_matches_control, "\n")
  cat("  Counts (0, 1, 2, ...):\n")
  print(x$match_count_dist_control)

  invisible(x)
}

#' Plot method for singling_out objects
#'
#' @param x an object of class "singling_out"
#' @param y not used
#' @param ... additional arguments passed to plotting functions
#' @param which which plot(s) to show: 1 = Risk comparison barplot (attack vs control),
#'   2 = Match count distribution histogram
#' @return No return value; called for the side effect of producing a plot.
#' @export
plot.singling_out <- function(x, y = NULL, ..., which = 1) {
  show <- rep(FALSE, 2)
  show[which] <- TRUE

  n_plots <- sum(show)
  if (n_plots > 1) {
    op <- par(mfrow = c(1, 2))
    on.exit(par(op))
  }

  if (show[1]) {
    # Risk comparison barplot
    vals <- c(x$risk_attack, x$risk_control, x$risk)
    names_vals <- c("Attack", "Control", "Residual")
    cols <- c("coral", "steelblue", "darkgreen")

    bp <- barplot(vals, names.arg = names_vals,
                  main = paste("Singling Out Risk\n(",
                               x$n_attacks, "predicates)"),
                  ylab = "Rate",
                  col = cols,
                  ylim = c(0, max(0.2, max(vals) * 1.3)), ...)
    abline(h = 0.1, col = "red", lwd = 2, lty = 2)
    text(bp, vals + max(vals) * 0.05,
         labels = round(vals, 3), cex = 0.9)
    legend("topright", "threshold = 0.1", col = "red", lty = 2, lwd = 2,
           cex = 0.8)
  }

  if (show[2]) {
    # Match count distribution
    max_count <- min(max(c(x$match_counts, x$match_counts_control)), 20)
    breaks <- seq(-0.5, max_count + 0.5, by = 1)

    hist(x$match_counts, breaks = breaks, col = rgb(1, 0, 0, 0.5),
         main = "Match Count Distribution",
         xlab = "Number of matching records per predicate",
         ylab = "Frequency", ...)
    hist(x$match_counts_control, breaks = breaks,
         col = rgb(0, 0, 1, 0.5), add = TRUE)
    abline(v = 1, col = "black", lwd = 2, lty = 2)
    legend("topright",
           legend = c("Original (attack)", "Holdout (control)",
                      "Singling out (count = 1)"),
           fill = c(rgb(1, 0, 0, 0.5), rgb(0, 0, 1, 0.5), NA),
           border = c("black", "black", NA),
           lty = c(NA, NA, 2), lwd = c(NA, NA, 2),
           col = c(NA, NA, "black"),
           cex = 0.8)
  }
}
