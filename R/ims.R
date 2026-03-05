#' Identical Match Share (IMS)
#'
#' Computes the Identical Match Share privacy metric for synthetic data.
#' IMS measures the proportion of synthetic records that are exact copies
#' of training records, indicating potential memorization or data leakage.
#'
#' For a variant that focuses only on unique (singleton) training records,
#' use \code{\link{repu}} (Replicated Uniques).
#'
#' @param X data frame of original/training data
#' @param Y data frame of synthetic data
#' @param vars character vector of variable names to use for matching.
#'   If NULL (default), all common variables are used.
#' @param na.rm logical, remove records with NA values (default: TRUE)
#' @param ... additional arguments passed to methods (currently unused)
#'
#' @return An object of class "ims" containing:
#' \itemize{
#'   \item ims: identical match share (proportion of synthetic with exact match)
#'   \item ims_pct: IMS as percentage
#'   \item n_identical: number of synthetic records with exact match
#'   \item identical_idx: indices of identical synthetic records in Y
#'   \item identical_records: the actual identical records from Y
#'   \item match_counts: number of training matches per identical synthetic record
#'   \item n_unique_matched: number of unique training records that were copied
#'   \item privacy_pass: logical, TRUE if IMS is acceptably low (< 1%)
#'   \item n_synthetic, n_train: dataset sizes
#' }
#'
#' @details
#' Identical Match Share (IMS) is a straightforward privacy metric that detects
#' exact copies of training records in synthetic data. Even well-designed
#' synthetic data generators can occasionally produce records identical to
#' training data by chance, especially for categorical-heavy datasets.
#'
#' The metric computes:
#' \deqn{IMS = \frac{|\{y \in Y : \exists x \in X, y = x\}|}{|Y|}}
#'
#' Interpretation:
#' \itemize{
#'   \item \strong{IMS = 0%}: No exact copies - ideal for privacy
#'   \item \strong{IMS < 1%}: Acceptable - likely coincidental matches
#'   \item \strong{IMS > 1%}: Concerning - may indicate memorization
#'   \item \strong{IMS > 5%}: Serious privacy issue - significant copying
#' }
#'
#' For datasets with many categorical variables and few unique combinations,
#' some identical matches may occur by chance. Compare IMS against the expected
#' collision rate for your data characteristics.
#'
#' @seealso \code{\link{dcr}} for distance to closest record,
#'   \code{\link{nndr}} for nearest neighbor distance ratio,
#'   \code{\link{disco}} for disclosive records with matching target values
#'
#' @references
#' MOSTLY AI (2024). Evaluate generator quality.
#' \url{https://docs.mostly.ai/generators/evaluate-quality}
#'
#' @author Matthias Templ
#' @family distance-risk
#' @export
#' @importFrom stats complete.cases
#' @importFrom graphics barplot pie
#' @examples
#' # Create example data
#' set.seed(123)
#' X <- data.frame(
#'   age = sample(20:60, 100, replace = TRUE),
#'   gender = sample(c("M", "F"), 100, replace = TRUE),
#'   region = sample(c("N", "S", "E", "W"), 100, replace = TRUE),
#'   income = sample(c("low", "medium", "high"), 100, replace = TRUE)
#' )
#'
#' # Good synthetic data (random)
#' Y_good <- data.frame(
#'   age = sample(20:60, 100, replace = TRUE),
#'   gender = sample(c("M", "F"), 100, replace = TRUE),
#'   region = sample(c("N", "S", "E", "W"), 100, replace = TRUE),
#'   income = sample(c("low", "medium", "high"), 100, replace = TRUE)
#' )
#'
#' result_good <- ims(X, Y_good)
#' print(result_good)
#'
#' # Bad synthetic data (direct copies)
#' Y_bad <- X[sample(nrow(X), 100, replace = TRUE), ]
#' result_bad <- ims(X, Y_bad)
#' print(result_bad)
ims <- function(X, ...) {
  UseMethod("ims")
}

#' @rdname ims
#' @export
ims.synth_pair <- function(X, ...) {
  ims.default(
    X = X$original,
    Y = X$synthetic,
    vars = X$vars,
    ...
  )
}

#' @rdname ims
#' @export
ims.default <- function(X, Y,
                        vars = NULL,
                        na.rm = TRUE,
                        ...) {

  # Input validation
  if (!is.data.frame(X)) stop("X must be a data frame.")
  if (!is.data.frame(Y)) stop("Y must be a data frame.")

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

  n_train <- nrow(X)
  n_synthetic <- nrow(Y)


  # Create string signatures for exact matching using interaction()
  # This is safer than paste() with a delimiter that might appear in data
  make_signature <- function(df) {
    if (ncol(df) == 1) {
      return(as.character(df[[1]]))
    }
    # Use interaction with drop=TRUE to create unique factor levels
    as.character(interaction(df, drop = TRUE, sep = "\x1F"))
  }

  sig_X <- make_signature(X)
  sig_Y <- make_signature(Y)

  # Find identical matches
  is_identical <- sig_Y %in% sig_X
  identical_idx <- which(is_identical)
  n_identical <- length(identical_idx)

  # Count matches per identical synthetic record
  match_counts <- integer(n_identical)
  matched_train_idx <- vector("list", n_identical)

  for (i in seq_along(identical_idx)) {
    matches <- which(sig_X == sig_Y[identical_idx[i]])
    match_counts[i] <- length(matches)
    matched_train_idx[[i]] <- matches
  }

  # Count unique training records that were copied
  all_matched_train <- unique(unlist(matched_train_idx))
  n_unique_matched <- length(all_matched_train)

  # Compute IMS
  ims_value <- n_identical / n_synthetic
  ims_pct <- 100 * ims_value

  # Privacy check: IMS should be very low (< 1%)
  privacy_pass <- ims_pct < 1

  results <- list(
    ims = ims_value,
    ims_pct = ims_pct,
    n_identical = n_identical,
    identical_idx = identical_idx,
    identical_records = if (n_identical > 0) Y[identical_idx, , drop = FALSE] else Y[0, , drop = FALSE],
    match_counts = match_counts,
    matched_train_idx = matched_train_idx,
    n_unique_matched = n_unique_matched,
    pct_train_copied = 100 * n_unique_matched / n_train,
    privacy_pass = privacy_pass,
    n_synthetic = n_synthetic,
    n_train = n_train,
    vars = vars
  )

  class(results) <- "ims"
  return(results)
}

#' Print method for ims objects
#'
#' @param x an object of class "ims"
#' @param ... additional arguments (ignored)
#' @export
print.ims <- function(x, ...) {
  cat("Identical Match Share (IMS) Privacy Metric\n")
  cat("==========================================\n")
  cat("Variables used:", length(x$vars), "\n\n")

  cat("Dataset Sizes:\n")
  cat("  Training records:", x$n_train, "\n")
  cat("  Synthetic records:", x$n_synthetic, "\n\n")

  cat("IMS Results:\n")
  cat("  Identical synthetic records:", x$n_identical,
      sprintf("(%.2f%%)", x$ims_pct), "\n")
  cat("  Unique training records copied:", x$n_unique_matched,
      sprintf("(%.2f%% of training)", x$pct_train_copied), "\n\n")

  cat("Privacy Assessment:")
  if (x$privacy_pass) {
    cat(" PASS\n")
    if (x$n_identical == 0) {
      cat("  No exact copies found in synthetic data.\n")
    } else {
      cat("  Very few exact copies - likely coincidental matches.\n")
    }
  } else {
    cat(" WARNING\n")
    if (x$ims_pct >= 5) {
      cat("  SERIOUS: High proportion of exact copies detected.\n")
      cat("  Synthetic data may be directly copying training records.\n")
    } else {
      cat("  Elevated number of exact copies detected.\n")
      cat("  Review synthetic data generation process.\n")
    }
  }

  invisible(x)
}

#' Summary method for ims objects
#'
#' @param object an object of class "ims"
#' @param ... additional arguments (ignored)
#' @export
summary.ims <- function(object, ...) {
  # Analyze match multiplicities
  if (length(object$match_counts) > 0) {
    match_table <- table(object$match_counts)
  } else {
    match_table <- table(integer(0))
  }

  summ <- list(
    ims = object$ims,
    ims_pct = object$ims_pct,
    n_identical = object$n_identical,
    n_unique_matched = object$n_unique_matched,
    pct_train_copied = object$pct_train_copied,
    privacy_pass = object$privacy_pass,
    match_distribution = match_table,
    mean_matches_per_copy = if (object$n_identical > 0) mean(object$match_counts) else NA,
    max_matches = if (object$n_identical > 0) max(object$match_counts) else 0,
    n_synthetic = object$n_synthetic,
    n_train = object$n_train,
    vars = object$vars
  )

  class(summ) <- "summary.ims"
  return(summ)
}

#' Print method for summary.ims objects
#'
#' @param x an object of class "summary.ims"
#' @param ... additional arguments (ignored)
#' @export
print.summary.ims <- function(x, ...) {
  cat("Summary: Identical Match Share (IMS)\n")
  cat("====================================\n")
  cat("Variables:", length(x$vars), "\n\n")

  cat("Dataset Sizes:\n")
  cat("  Training:", x$n_train, "| Synthetic:", x$n_synthetic, "\n\n")

  cat("Key Metrics:\n")
  cat("  IMS:", sprintf("%.4f", x$ims), sprintf("(%.2f%%)", x$ims_pct), "\n")
  cat("  Privacy:", ifelse(x$privacy_pass, "PASS", "WARNING"), "\n\n")

  cat("Identical Records Analysis:\n")
  cat("  Synthetic records with exact match:", x$n_identical, "\n")
  cat("  Unique training records copied:", x$n_unique_matched,
      sprintf("(%.2f%% of training)", x$pct_train_copied), "\n")

  if (x$n_identical > 0) {
    cat("  Avg training matches per copy:", round(x$mean_matches_per_copy, 2), "\n")
    cat("  Max training matches:", x$max_matches, "\n\n")

    cat("Match Multiplicity Distribution:\n")
    cat("  (How many training records match each identical synthetic record)\n")
    print(x$match_distribution)
  }

  invisible(x)
}

#' Plot method for ims objects
#'
#' @param x an object of class "ims"
#' @param y not used
#' @param ... additional arguments passed to plotting functions
#' @param which which plot to show: 1 = pie chart of identical vs unique,
#'   2 = bar chart of match counts (if any identical records exist)
#' @export
plot.ims <- function(x, y = NULL, ..., which = 1) {
  show <- rep(FALSE, 2)
  show[which] <- TRUE

  n_plots <- sum(show)
  if (n_plots > 1) {
    op <- par(mfrow = c(1, 2))
    on.exit(par(op))
  }

  if (show[1]) {
    # Pie chart or bar chart of identical vs non-identical
    counts <- c(x$n_synthetic - x$n_identical, x$n_identical)
    names(counts) <- c("Unique", "Identical")

    if (x$n_identical == 0) {
      barplot(counts,
              main = "Synthetic Record Status",
              ylab = "Number of Records",
              col = c("forestgreen", "firebrick"), ...)
    } else {
      colors <- c("forestgreen", "firebrick")
      pie(counts,
          main = paste("Identical Match Share:", sprintf("%.2f%%", x$ims_pct)),
          col = colors,
          labels = paste(names(counts), "\n", counts,
                         sprintf("\n(%.1f%%)", 100 * counts / sum(counts))), ...)
    }
  }

  if (show[2]) {
    if (x$n_identical > 0) {
      # Distribution of match counts
      match_table <- table(x$match_counts)
      barplot(match_table,
              main = "Training Matches per Identical Record",
              xlab = "Number of Training Matches",
              ylab = "Count of Synthetic Records",
              col = "coral", ...)
    } else {
      # Empty plot with message
      plot(1, type = "n", axes = FALSE, xlab = "", ylab = "",
           main = "No Identical Records Found")
      text(1, 1, "No identical matches to analyze", cex = 1.2)
    }
  }
}


#' Replicated Uniques (RepU)
#'
#' Alias for \code{\link{ims}}. Computes the share of synthetic records that
#' are identical to unique (singleton) records in the training data. This is
#' particularly concerning as unique records are more identifiable.
#'
#' @inheritParams ims
#' @param uniques_only logical, if TRUE (default), only count matches to
#'   records that are unique in the training data. If FALSE, equivalent to
#'   \code{\link{ims}}.
#'
#' @return An object of class "ims" (same structure as \code{\link{ims}})
#'
#' @details
#' Replicated Uniques (RepU) is a variant of IMS that focuses specifically on
#' synthetic records that match unique (singleton) training records. Unique
#' records are more identifiable and thus pose a higher privacy risk if copied.
#'
#' @seealso \code{\link{ims}} for identical match share
#'
#' @family distance-risk
#' @export
#' @examples
#' set.seed(123)
#' X <- data.frame(
#'   age = sample(20:60, 100, replace = TRUE),
#'   gender = sample(c("M", "F"), 100, replace = TRUE),
#'   income = sample(c("low", "medium", "high"), 100, replace = TRUE)
#' )
#' Y <- X[sample(nrow(X), 50, replace = TRUE), ]
#'
#' # All identical matches
#' result_ims <- ims(X, Y)
#'
#' # Only matches to unique training records
#' result_repu <- repu(X, Y)
repu <- function(X, Y,
                 vars = NULL,
                 uniques_only = TRUE,
                 na.rm = TRUE) {

  if (!uniques_only) {
    return(ims(X, Y, vars = vars, na.rm = na.rm))
  }

  # Input validation
  if (!is.data.frame(X)) stop("X must be a data frame.")
  if (!is.data.frame(Y)) stop("Y must be a data frame.")

  # Determine variables to use
  if (is.null(vars)) {
    vars <- intersect(names(X), names(Y))
  }

  if (length(vars) == 0) {
    stop("No common variables found between X and Y.")
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

  # Create string signatures for exact matching using interaction()
  # This is safer than paste() with a delimiter that might appear in data
  make_signature <- function(df) {
    if (ncol(df) == 1) {
      return(as.character(df[[1]]))
    }
    as.character(interaction(df, drop = TRUE, sep = "\x1F"))
  }

  sig_X <- make_signature(X)
  sig_Y <- make_signature(Y)

  # Find unique (singleton) records in training data
  sig_counts <- table(sig_X)
  unique_sigs <- names(sig_counts[sig_counts == 1])

  n_train <- nrow(X)
  n_synthetic <- nrow(Y)
  n_unique_train <- length(unique_sigs)

  # Find synthetic records matching unique training records
  is_repu <- sig_Y %in% unique_sigs
  identical_idx <- which(is_repu)
  n_identical <- length(identical_idx)

  # Compute metrics
  ims_value <- n_identical / n_synthetic
  ims_pct <- 100 * ims_value

  # Count unique training records that were copied
  matched_sigs <- unique(sig_Y[is_repu])
  n_unique_matched <- length(matched_sigs)

  # Privacy check
  privacy_pass <- ims_pct < 1

  results <- list(
    ims = ims_value,
    ims_pct = ims_pct,
    n_identical = n_identical,
    identical_idx = identical_idx,
    identical_records = if (n_identical > 0) Y[identical_idx, , drop = FALSE] else Y[0, , drop = FALSE],
    match_counts = rep(1L, n_identical),  # All matches are to uniques
    matched_train_idx = as.list(match(sig_Y[identical_idx], sig_X)),
    n_unique_matched = n_unique_matched,
    pct_train_copied = 100 * n_unique_matched / n_unique_train,
    n_unique_train = n_unique_train,
    privacy_pass = privacy_pass,
    n_synthetic = n_synthetic,
    n_train = n_train,
    vars = vars,
    uniques_only = TRUE
  )

  class(results) <- "ims"
  return(results)
}
