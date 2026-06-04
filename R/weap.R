#' Within Equivalence Class Attribution Probability (WEAP)
#'
#' Computes the Within Equivalence Class Attribution Probability for synthetic
#' data records. WEAP identifies potentially disclosive synthetic records by
#' measuring the conditional probability of the target value given the key
#' variables within the synthetic dataset.
#'
#' @param X data frame of synthetic data
#' @param key_vars character vector of quasi-identifier variable names
#' @param target_var character, name of the sensitive target variable
#' @param na.rm logical, remove records with NA in key or target (default: TRUE)
#' @param ... additional arguments passed to methods (currently unused)
#'
#' @return An object of class "weap" containing:
#' \itemize{
#'   \item weap_scores: numeric vector of per-record WEAP values in Y
#'   \item weap_mean: mean WEAP across all records
#'   \item equivalence_classes: data frame with equivalence class information
#'   \item n_disclosive: number of records with WEAP = 1 (potentially disclosive)
#'   \item pct_disclosive: percentage of potentially disclosive records
#'   \item n_unique_keys: number of unique key combinations
#'   \item key_vars, target_var: input parameters
#' }
#'
#' @details
#' WEAP (Within Equivalence Class Attribution Probability) is computed for each
#' record in the synthetic dataset. An equivalence class is the set of all
#' records sharing the same key variable values.
#'
#' For each synthetic record j with key values k_j and target value t_j:
#' \deqn{WEAP_j = \frac{|\{i \in Y : keys_i = k_j \land target_i = t_j\}|}{|\{i \in Y : keys_i = k_j\}|}}
#'
#' Records with WEAP = 1 are potentially disclosive because the target value
#' is uniquely determined by the key values within the synthetic data.
#'
#' This measure is used to identify risky synthetic records before release.
#' High WEAP values indicate that if an intruder knows the key values, they
#' can confidently infer the target value.
#'
#' @seealso \code{\link{tcap}} for targeted correct attribution probability,
#'   \code{\link{dcap}} for differential correct attribution probability,
#'   \code{\link{disco}} for identifying disclosive records
#'
#' @references
#' Elliot, M. (2014). Final Report on the Disclosure Risk Associated with the
#' Synthetic Data Produced by the SYLLS Team. Report 2015-2.
#'
#' Taub, J., Elliot, M., Pampaka, M., & Smith, D. (2018). Differential Correct
#' Attribution Probability for Synthetic Data: An Exploration.
#' \emph{Privacy in Statistical Databases}, 122-137.
#'
#' @author Matthias Templ
#' @family attribution-risk
#' @export
#' @importFrom stats complete.cases aggregate median sd
#' @importFrom graphics hist barplot abline legend pie
#' @examples
#' # Create synthetic data
#' set.seed(123)
#' X <- data.frame(
#'   age = sample(c("young", "middle", "old"), 200, replace = TRUE),
#'   gender = sample(c("M", "F"), 200, replace = TRUE),
#'   region = sample(c("N", "S"), 200, replace = TRUE),
#'   income = sample(c("low", "medium", "high"), 200, replace = TRUE)
#' )
#'
#' # Compute WEAP
#' result <- weap(X,
#'                key_vars = c("age", "gender", "region"),
#'                target_var = "income")
#' print(result)
#' summary(result)
#' plot(result)
#'
#' # Identify disclosive records
#' disclosive_idx <- which(result$weap_scores == 1)
#' if (length(disclosive_idx) > 0) {
#'   head(X[disclosive_idx, ])
#' }
weap <- function(X, ...) {
  UseMethod("weap")
}

#' @rdname weap
#' @export
weap.synth_pair <- function(X, ...) {
  if (is.null(X$key_vars)) {
    stop("synth_pair must have 'key_vars' set for weap()")
  }
  if (is.null(X$target_var)) {
    stop("synth_pair must have 'target_var' set for weap()")
  }

  weap.default(
    X = X$synthetic,
    key_vars = X$key_vars,
    target_var = X$target_var,
    ...
  )
}

#' @rdname weap
#' @export
weap.default <- function(X,
                         key_vars,
                         target_var,
                         na.rm = TRUE,
                         ...) {

  # Input validation
  if (!is.data.frame(X)) stop("X must be a data frame.")

  # Check variables exist
  all_vars <- c(key_vars, target_var)
  missing_vars <- setdiff(all_vars, names(X))
  if (length(missing_vars) > 0) {
    stop(paste("Variables missing in X:", paste(missing_vars, collapse = ", ")))
  }

  # Handle missing values
  if (na.rm) {
    complete_X <- complete.cases(X[, all_vars, drop = FALSE])
    X <- X[complete_X, , drop = FALSE]
  }

  if (nrow(X) == 0) stop("No complete cases in X after removing NAs.")

  n <- nrow(X)
  target_X <- X[[target_var]]

  # Create key signature for each record
  make_key <- function(df, vars) {
    apply(df[, vars, drop = FALSE], 1, paste, collapse = "|")
  }
  keys_X <- make_key(X, key_vars)

  # Create combined key-target signature
  key_target <- paste(keys_X, target_X, sep = "||")

  # Count occurrences of each key and each key-target combination
  key_counts <- table(keys_X)
  key_target_counts <- table(key_target)

  # Compute WEAP for each record
  weap_scores <- numeric(n)
  for (i in seq_len(n)) {
    key_i <- keys_X[i]
    key_target_i <- key_target[i]
    weap_scores[i] <- as.numeric(key_target_counts[key_target_i]) /
                      as.numeric(key_counts[key_i])
  }

  # Build equivalence class summary
  eq_classes <- data.frame(
    key = names(key_counts),
    size = as.numeric(key_counts)
  )

  # For each equivalence class, count unique target values
  eq_classes$n_unique_targets <- sapply(eq_classes$key, function(k) {
    length(unique(target_X[keys_X == k]))
  })
  eq_classes$is_disclosive <- eq_classes$n_unique_targets == 1

  # Summary statistics
  n_disclosive <- sum(weap_scores == 1)
  pct_disclosive <- 100 * n_disclosive / n

  results <- list(
    weap_scores = weap_scores,
    weap_mean = mean(weap_scores),
    weap_median = median(weap_scores),
    equivalence_classes = eq_classes,
    n_disclosive = n_disclosive,
    pct_disclosive = pct_disclosive,
    n_total = n,
    n_unique_keys = nrow(eq_classes),
    key_vars = key_vars,
    target_var = target_var
  )

  class(results) <- "weap"
  return(results)
}

#' Print method for weap objects
#'
#' @param x an object of class "weap"
#' @param ... additional arguments (ignored)
#' @return The input object, invisibly.
#' @export
print.weap <- function(x, ...) {
  cat("Within Equivalence Class Attribution Probability (WEAP)\n")
  cat("=======================================================\n")
  cat("Key variables:", paste(x$key_vars, collapse = ", "), "\n")
  cat("Target variable:", x$target_var, "\n\n")

  cat("Results:\n")
  cat("  Mean WEAP:", round(x$weap_mean, 4), "\n")
  cat("  Median WEAP:", round(x$weap_median, 4), "\n")
  cat("  Total records:", x$n_total, "\n")
  cat("  Unique key combinations:", x$n_unique_keys, "\n\n")

  cat("Disclosure Risk:\n")
  cat("  Records with WEAP = 1:", x$n_disclosive,
      sprintf("(%.1f%%)", x$pct_disclosive), "\n")

  if (x$pct_disclosive > 10) {
    cat("  WARNING: High proportion of potentially disclosive records\n")
  } else if (x$pct_disclosive > 0) {
    cat("  Note: Some records may be disclosive\n")
  } else {
    cat("  No records with unique key-target combinations\n")
  }

  invisible(x)
}

#' Summary method for weap objects
#'
#' @param object an object of class "weap"
#' @param ... additional arguments (ignored)
#' @return A list of summary statistics for the corresponding object.
#' @export
summary.weap <- function(object, ...) {
  scores <- object$weap_scores

  # Analyze equivalence class sizes
  ec <- object$equivalence_classes

  summ <- list(
    weap_mean = object$weap_mean,
    weap_median = object$weap_median,
    weap_sd = sd(scores),
    weap_quantiles = quantile(scores, probs = c(0, 0.25, 0.5, 0.75, 1)),
    n_disclosive = object$n_disclosive,
    pct_disclosive = object$pct_disclosive,
    n_total = object$n_total,
    n_unique_keys = object$n_unique_keys,
    ec_size_mean = mean(ec$size),
    ec_size_median = median(ec$size),
    ec_size_range = range(ec$size),
    n_singleton_keys = sum(ec$size == 1),
    n_disclosive_classes = sum(ec$is_disclosive),
    key_vars = object$key_vars,
    target_var = object$target_var
  )

  class(summ) <- "summary.weap"
  return(summ)
}

#' Print method for summary.weap objects
#'
#' @param x an object of class "summary.weap"
#' @param ... additional arguments (ignored)
#' @return The input object, invisibly.
#' @export
print.summary.weap <- function(x, ...) {
  cat("Summary: Within Equivalence Class Attribution Probability (WEAP)\n")
  cat("================================================================\n")
  cat("Key variables:", paste(x$key_vars, collapse = ", "), "\n")
  cat("Target variable:", x$target_var, "\n\n")

  cat("WEAP Distribution:\n")
  cat("  Mean:", round(x$weap_mean, 4), "\n")
  cat("  Median:", round(x$weap_median, 4), "\n")
  cat("  SD:", round(x$weap_sd, 4), "\n")
  cat("  Quantiles:\n")
  print(round(x$weap_quantiles, 4))
  cat("\n")

  cat("Equivalence Class Statistics:\n")
  cat("  Number of unique key combinations:", x$n_unique_keys, "\n")
  cat("  Mean class size:", round(x$ec_size_mean, 1), "\n")
  cat("  Median class size:", x$ec_size_median, "\n")
  cat("  Class size range:", x$ec_size_range[1], "-", x$ec_size_range[2], "\n")
  cat("  Singleton classes (size=1):", x$n_singleton_keys, "\n\n")

  cat("Disclosure Risk:\n")
  cat("  Disclosive records (WEAP=1):", x$n_disclosive,
      sprintf("(%.1f%% of total)", x$pct_disclosive), "\n")
  cat("  Disclosive equivalence classes:", x$n_disclosive_classes,
      sprintf("(%.1f%% of classes)", 100 * x$n_disclosive_classes / x$n_unique_keys), "\n")

  invisible(x)
}

#' Plot method for weap objects
#'
#' @param x an object of class "weap"
#' @param y not used
#' @param ... additional arguments passed to plotting functions
#' @param which which plot(s) to show: 1 = WEAP histogram, 2 = equivalence class
#'   sizes, 3 = disclosive vs non-disclosive
#' @return No return value; called for the side effect of producing a plot.
#' @export
plot.weap <- function(x, y = NULL, ..., which = 1) {
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
    # Histogram of WEAP scores
    hist(x$weap_scores, breaks = 20,
         main = "Distribution of WEAP Scores",
         xlab = "WEAP Score", col = "steelblue", border = "white", ...)
    abline(v = x$weap_mean, col = "red", lwd = 2, lty = 2)
  }

  if (show[2]) {
    # Equivalence class size distribution
    ec_sizes <- x$equivalence_classes$size
    hist(ec_sizes, breaks = 30,
         main = "Equivalence Class Sizes",
         xlab = "Number of Records per Class",
         col = "coral", border = "white", ...)
    abline(v = mean(ec_sizes), col = "red", lwd = 2, lty = 2)
  }

  if (show[3]) {
    # Disclosive vs non-disclosive
    counts <- c(x$n_total - x$n_disclosive, x$n_disclosive)
    names(counts) <- c("Non-disclosive", "Disclosive\n(WEAP=1)")
    barplot(counts,
            main = "Record Disclosure Status",
            ylab = "Number of Records",
            col = c("forestgreen", "firebrick"), ...)
  }
}
