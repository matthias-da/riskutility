#' delta-Presence Risk Assessment
#'
#' Computes the delta-presence risk measure for synthetic/released data relative
#' to the original/population data. delta-Presence (Nergiz et al., 2007) bounds
#' the membership probability for any record in the population, measuring how
#' confidently an adversary can determine whether a specific individual contributed
#' to the released dataset.
#'
#' @param X data frame of original/population data, or a \code{\link{synth_pair}} object
#' @param Y data frame of synthetic/released data
#' @param key_vars character vector of quasi-identifier variable names
#' @param delta_min numeric, lower bound on acceptable membership probability
#'   (default: 0.0)
#' @param delta_max numeric, upper bound on acceptable membership probability
#'   (default: 1.0)
#' @param na.rm logical, remove records with NA in key variables (default: TRUE)
#' @param ... additional arguments passed to methods (currently unused)
#'
#' @return An object of class "delta_presence" containing:
#' \itemize{
#'   \item membership_prob: numeric vector of f_k/F_k per record in X
#'   \item per_combination: data.frame with QI combination, f_k, F_k, prob
#'   \item delta_min, delta_max: bounds used
#'   \item n_violations_lower: number of records below delta_min
#'   \item n_violations_upper: number of records above delta_max
#'   \item pct_violations: fraction of records violating bounds
#'   \item satisfies_delta: logical, all records within bounds
#'   \item privacy_pass: same as satisfies_delta
#'   \item n_original, n_synthetic: dataset sizes
#'   \item key_vars: quasi-identifier variables used
#' }
#'
#' @details
#' delta-Presence is a formal privacy model that bounds the membership disclosure
#' probability. Given a population dataset (the original data X) and a published
#' subset (the synthetic/released data Y), delta-presence bounds the probability
#' that any record from X appears in Y.
#'
#' For each quasi-identifier combination k:
#' \itemize{
#'   \item \eqn{f_k} = frequency of combination k in the synthetic data Y
#'   \item \eqn{F_k} = frequency of combination k in the original data X
#'   \item Membership probability: \eqn{Pr(t \in Y | t \in X) = f_k / F_k}
#' }
#'
#' The dataset satisfies \eqn{(\delta_{min}, \delta_{max})}-presence if for
#' all records: \eqn{\delta_{min} \le f_k / F_k \le \delta_{max}}.
#'
#' Membership probabilities are capped at 1.0 (when f_k > F_k, which can occur
#' when the synthetic data overrepresents certain combinations).
#'
#' \strong{Interpretation:}
#' \itemize{
#'   \item Values close to 0: The QI combination is underrepresented in Y
#'     relative to X, meaning membership is unlikely
#'   \item Values close to 1: The QI combination has similar frequency in both
#'     datasets, meaning membership is highly probable
#'   \item Values = 0: The QI combination exists in X but not in Y (no
#'     membership risk for these records)
#' }
#'
#' \strong{QI combinations in Y but not X:}
#' If QI combinations appear in the synthetic data Y that do not exist in the
#' original data X, a warning is issued. These represent fabricated combinations
#' that do not correspond to any real records.
#'
#' @seealso \code{\link{kanonymity}} for k-anonymity assessment,
#'   \code{\link{ldiversity}} for l-diversity assessment,
#'   \code{\link{tcloseness}} for t-closeness assessment,
#'   \code{\link{individual_risk}} for probabilistic risk assessment
#'
#' @references
#' Nergiz, M. E., Atzori, M. & Clifton, C. (2007).
#' Hiding the Presence of Individuals from Shared Databases.
#' \emph{Proceedings of the 2007 ACM SIGMOD International Conference on
#' Management of Data}, 665--676.
#' \doi{10.1145/1247480.1247554}
#'
#' @author Matthias Templ
#' @family privacy-models
#' @export
#' @importFrom stats complete.cases
#' @importFrom graphics hist abline barplot par legend
#' @examples
#' # Create example data
#' set.seed(123)
#' n <- 200
#' X <- data.frame(
#'   age = sample(c("young", "middle", "old"), n, replace = TRUE),
#'   gender = sample(c("M", "F"), n, replace = TRUE),
#'   region = sample(c("N", "S", "E", "W"), n, replace = TRUE)
#' )
#'
#' # Synthetic data with similar distribution
#' Y <- data.frame(
#'   age = sample(c("young", "middle", "old"), n, replace = TRUE),
#'   gender = sample(c("M", "F"), n, replace = TRUE),
#'   region = sample(c("N", "S", "E", "W"), n, replace = TRUE)
#' )
#'
#' result <- delta_presence(X, Y,
#'                          key_vars = c("age", "gender", "region"),
#'                          delta_min = 0.1, delta_max = 0.9)
#' print(result)
#' summary(result)
#' plot(result)
#'
#' \donttest{
#' # Memorized data - all probabilities near 1
#' Y_copy <- X
#' result_bad <- delta_presence(X, Y_copy,
#'                              key_vars = c("age", "gender", "region"),
#'                              delta_max = 0.5)
#' print(result_bad)
#' }
delta_presence <- function(X, ...) {
  UseMethod("delta_presence")
}

#' @rdname delta_presence
#' @export
delta_presence.synth_pair <- function(X, ...) {
  if (is.null(X$key_vars)) {
    stop("synth_pair must have 'key_vars' set for delta_presence()")
  }

  delta_presence.default(
    X = X$original,
    Y = X$synthetic,
    key_vars = X$key_vars,
    ...
  )
}

#' @rdname delta_presence
#' @export
delta_presence.default <- function(X, Y,
                                   key_vars,
                                   delta_min = 0.0,
                                   delta_max = 1.0,
                                   na.rm = TRUE,
                                   ...) {

  # Input validation
  if (!is.data.frame(X)) stop("X must be a data frame.")
  if (!is.data.frame(Y)) stop("Y must be a data frame.")

  if (!is.character(key_vars) || length(key_vars) == 0) {
    stop("key_vars must be a non-empty character vector.")
  }

  if (!is.numeric(delta_min) || length(delta_min) != 1 ||
      delta_min < 0 || delta_min > 1) {
    stop("delta_min must be a single numeric value between 0 and 1.")
  }
  if (!is.numeric(delta_max) || length(delta_max) != 1 ||
      delta_max < 0 || delta_max > 1) {
    stop("delta_max must be a single numeric value between 0 and 1.")
  }
  if (delta_min > delta_max) {
    stop("delta_min must be less than or equal to delta_max.")
  }

  # Check key variables exist in both datasets
  missing_X <- setdiff(key_vars, names(X))
  missing_Y <- setdiff(key_vars, names(Y))
  if (length(missing_X) > 0) {
    stop(paste("Key variables missing in X:", paste(missing_X, collapse = ", ")))
  }
  if (length(missing_Y) > 0) {
    stop(paste("Key variables missing in Y:", paste(missing_Y, collapse = ", ")))
  }

  # Subset to key variables
  X_keys <- X[, key_vars, drop = FALSE]
  Y_keys <- Y[, key_vars, drop = FALSE]

  # Handle missing values
  if (na.rm) {
    complete_X <- complete.cases(X_keys)
    complete_Y <- complete.cases(Y_keys)
    if (sum(!complete_X) > 0) {
      message(sprintf("Removing %d records with NA in key variables from X.",
                      sum(!complete_X)))
    }
    if (sum(!complete_Y) > 0) {
      message(sprintf("Removing %d records with NA in key variables from Y.",
                      sum(!complete_Y)))
    }
    X_keys <- X_keys[complete_X, , drop = FALSE]
    Y_keys <- Y_keys[complete_Y, , drop = FALSE]
  }

  if (nrow(X_keys) == 0) stop("No complete cases remaining in X.")
  if (nrow(Y_keys) == 0) stop("No complete cases remaining in Y.")

  n_original <- nrow(X_keys)
  n_synthetic <- nrow(Y_keys)

  # Create key signatures
  key_sig_X <- apply(X_keys, 1, paste, collapse = "|")
  key_sig_Y <- apply(Y_keys, 1, paste, collapse = "|")

  # Compute frequencies
  F_k_table <- table(key_sig_X)  # frequencies in X (population)
  f_k_table <- table(key_sig_Y)  # frequencies in Y (released)

  # Check for fabricated combinations (in Y but not in X)
  fabricated <- setdiff(names(f_k_table), names(F_k_table))
  if (length(fabricated) > 0) {
    n_fabricated_records <- sum(f_k_table[fabricated])
    warning(sprintf(
      paste0("%d QI combination(s) in synthetic data not found in original ",
             "data (%d records). These are fabricated combinations."),
      length(fabricated), n_fabricated_records
    ))
  }

  # Build per-combination table
  all_combos_X <- names(F_k_table)
  combo_df <- data.frame(
    combination = all_combos_X,
    F_k = as.integer(F_k_table[all_combos_X])
  )

  # Match f_k: combos in X that also appear in Y get their Y count, else 0
  combo_df$f_k <- ifelse(
    combo_df$combination %in% names(f_k_table),
    as.integer(f_k_table[combo_df$combination]),
    0L
  )

  # Compute membership probability, capped at 1
  combo_df$prob <- pmin(combo_df$f_k / combo_df$F_k, 1.0)

  # Sort by probability descending

  combo_df <- combo_df[order(-combo_df$prob, -combo_df$F_k), ]
  rownames(combo_df) <- NULL

  # Compute per-record membership probability
  # Map each record in X to its combination's probability
  combo_prob_lookup <- combo_df$prob
  names(combo_prob_lookup) <- combo_df$combination
  membership_prob <- as.numeric(combo_prob_lookup[key_sig_X])

  # Count violations
  n_violations_lower <- sum(membership_prob < delta_min)
  n_violations_upper <- sum(membership_prob > delta_max)
  n_violations_total <- n_violations_lower + n_violations_upper
  pct_violations <- n_violations_total / n_original

  satisfies_delta <- n_violations_total == 0L

  results <- list(
    membership_prob = membership_prob,
    per_combination = combo_df,
    delta_min = delta_min,
    delta_max = delta_max,
    n_violations_lower = n_violations_lower,
    n_violations_upper = n_violations_upper,
    pct_violations = pct_violations,
    satisfies_delta = satisfies_delta,
    privacy_pass = satisfies_delta,
    n_original = n_original,
    n_synthetic = n_synthetic,
    key_vars = key_vars
  )

  class(results) <- "delta_presence"
  return(results)
}


#' Print method for delta_presence objects
#' @param x an object of class "delta_presence"
#' @param ... additional arguments (ignored)
#' @return The input object, invisibly.
#' @export
print.delta_presence <- function(x, ...) {
  cat("delta-Presence Risk Assessment\n")
  cat("==============================\n\n")

  cat("Key variables:", paste(x$key_vars, collapse = ", "), "\n")
  cat("Original records:", x$n_original,
      "| Synthetic records:", x$n_synthetic, "\n")
  cat("QI combinations:", nrow(x$per_combination), "\n\n")

  cat("Bounds: delta_min =", x$delta_min,
      ", delta_max =", x$delta_max, "\n\n")

  cat("Membership Probability:\n")
  cat("  Mean:", sprintf("%.4f", mean(x$membership_prob)), "\n")
  cat("  Median:", sprintf("%.4f", median(x$membership_prob)), "\n")
  cat("  Range: [", sprintf("%.4f", min(x$membership_prob)), ",",
      sprintf("%.4f", max(x$membership_prob)), "]\n\n")

  cat("Violations:\n")
  cat("  Below delta_min:", x$n_violations_lower, "records\n")
  cat("  Above delta_max:", x$n_violations_upper, "records\n")
  cat("  Total violations:", x$n_violations_lower + x$n_violations_upper,
      sprintf("(%.1f%%)", 100 * x$pct_violations), "\n\n")

  cat("Privacy Assessment:")
  if (x$privacy_pass) {
    cat(" PASS\n")
    cat("  All records satisfy (", x$delta_min, ", ", x$delta_max,
        ")-presence.\n", sep = "")
  } else {
    cat(" WARNING\n")
    cat("  Some records violate the delta-presence bounds.\n")
  }

  invisible(x)
}


#' Summary method for delta_presence objects
#' @param object an object of class "delta_presence"
#' @param ... additional arguments (ignored)
#' @return A list of summary statistics for the corresponding object.
#' @export
summary.delta_presence <- function(object, ...) {
  prob_summary <- data.frame(
    metric = c("Min", "Q1", "Median", "Mean", "Q3", "Max"),
    value = as.numeric(summary(object$membership_prob)[1:6])
  )

  # Combinations with zero membership (in X but not in Y)
  n_zero <- sum(object$per_combination$f_k == 0)

  # Combinations at cap (prob == 1)
  n_capped <- sum(object$per_combination$prob == 1)

  # Worst combinations (highest membership probability)
  worst <- head(object$per_combination[order(-object$per_combination$prob), ], 10)

  summ <- list(
    delta_min = object$delta_min,
    delta_max = object$delta_max,
    satisfies_delta = object$satisfies_delta,
    n_violations_lower = object$n_violations_lower,
    n_violations_upper = object$n_violations_upper,
    pct_violations = object$pct_violations,
    n_original = object$n_original,
    n_synthetic = object$n_synthetic,
    n_combinations = nrow(object$per_combination),
    n_zero_membership = n_zero,
    n_capped = n_capped,
    prob_summary = prob_summary,
    worst_combinations = worst,
    key_vars = object$key_vars
  )

  class(summ) <- "summary.delta_presence"
  return(summ)
}


#' Print method for summary.delta_presence objects
#' @param x an object of class "summary.delta_presence"
#' @param ... additional arguments (ignored)
#' @return The input object, invisibly.
#' @export
print.summary.delta_presence <- function(x, ...) {
  cat("Summary: delta-Presence Risk Assessment\n")
  cat("=======================================\n\n")

  cat("Key variables:", paste(x$key_vars, collapse = ", "), "\n")
  cat("Original records:", x$n_original,
      "| Synthetic records:", x$n_synthetic, "\n")
  cat("Bounds: delta_min =", x$delta_min,
      ", delta_max =", x$delta_max, "\n\n")

  cat("Overall Result:\n")
  cat("  Satisfies delta-presence:",
      ifelse(x$satisfies_delta, "YES", "NO"), "\n")
  cat("  Violations below delta_min:", x$n_violations_lower, "\n")
  cat("  Violations above delta_max:", x$n_violations_upper, "\n")
  cat("  Total violation rate:", sprintf("%.1f%%", 100 * x$pct_violations), "\n\n")

  cat("QI Combinations:\n")
  cat("  Total:", x$n_combinations, "\n")
  cat("  Zero membership (in X only):", x$n_zero_membership, "\n")
  cat("  Capped at 1.0:", x$n_capped, "\n\n")

  cat("Membership Probability Distribution:\n")
  print(x$prob_summary, row.names = FALSE)
  cat("\n")

  cat("Highest-Risk Combinations (by membership probability):\n")
  worst <- x$worst_combinations[, c("combination", "F_k", "f_k", "prob")]
  worst$prob <- sprintf("%.4f", worst$prob)
  print(worst, row.names = FALSE)

  invisible(x)
}


#' Plot method for delta_presence objects
#' @param x an object of class "delta_presence"
#' @param y not used
#' @param ... additional arguments passed to plotting functions
#' @param which integer, which plot(s) to show:
#'   1 = histogram of membership probabilities with delta bounds,
#'   2 = per-QI-combination barplot sorted by probability
#' @importFrom graphics hist abline barplot par legend
#' @return No return value; called for the side effect of producing a plot.
#' @export
plot.delta_presence <- function(x, y = NULL, ..., which = 1) {
  show <- rep(FALSE, 2)
  show[which] <- TRUE

  n_plots <- sum(show)
  if (n_plots > 1) {
    op <- par(mfrow = c(1, 2))
    on.exit(par(op))
  }

  if (show[1]) {
    # Histogram of membership probabilities
    hist(x$membership_prob,
         breaks = 20,
         main = paste0("Membership Probability Distribution\n(",
                        x$delta_min, ", ", x$delta_max, ")-presence"),
         xlab = "Membership Probability (f_k / F_k)",
         ylab = "Number of Records",
         col = "steelblue", border = "white",
         xlim = c(0, max(1, max(x$membership_prob))),
         ...)
    if (x$delta_min > 0) {
      abline(v = x$delta_min, col = "orange", lwd = 2, lty = 2)
    }
    if (x$delta_max < 1) {
      abline(v = x$delta_max, col = "red", lwd = 2, lty = 2)
    }
    leg_labels <- character(0)
    leg_cols <- character(0)
    if (x$delta_min > 0) {
      leg_labels <- c(leg_labels, paste0("delta_min = ", x$delta_min))
      leg_cols <- c(leg_cols, "orange")
    }
    if (x$delta_max < 1) {
      leg_labels <- c(leg_labels, paste0("delta_max = ", x$delta_max))
      leg_cols <- c(leg_cols, "red")
    }
    if (length(leg_labels) > 0) {
      legend("topright", legend = leg_labels,
             col = leg_cols, lty = 2, lwd = 2, cex = 0.8)
    }
  }

  if (show[2]) {
    # Per-combination barplot (top 30, sorted by prob)
    combo <- x$per_combination[order(-x$per_combination$prob), ]
    n_show <- min(30, nrow(combo))
    combo <- combo[seq_len(n_show), ]

    bar_cols <- ifelse(
      combo$prob < x$delta_min | combo$prob > x$delta_max,
      "firebrick", "steelblue"
    )

    bp <- barplot(combo$prob,
                  names.arg = if (n_show <= 15) combo$combination else rep("", n_show),
                  main = "Membership Probability by QI Combination",
                  xlab = if (n_show <= 15) "QI Combination" else paste0("Top ", n_show, " QI Combinations"),
                  ylab = "Membership Probability",
                  col = bar_cols,
                  ylim = c(0, min(1.2, max(combo$prob) * 1.3)),
                  las = if (n_show <= 15) 2 else 0,
                  ...)
    if (x$delta_max < 1) {
      abline(h = x$delta_max, col = "red", lwd = 2, lty = 2)
    }
    if (x$delta_min > 0) {
      abline(h = x$delta_min, col = "orange", lwd = 2, lty = 2)
    }
    legend("topright",
           legend = c("Within bounds", "Violating"),
           fill = c("steelblue", "firebrick"), cex = 0.8)
  }
}
