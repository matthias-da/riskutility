#' k-Anonymity Assessment
#'
#' Measures k-anonymity of a dataset based on specified quasi-identifier (key)
#' variables. A dataset satisfies k-anonymity if every record is indistinguishable
#' from at least k-1 other records with respect to the key variables.
#'
#' @param X data frame to assess, or a \code{\link{synth_pair}} object
#' @param key_vars character vector of quasi-identifier variable names
#' @param k integer, the k-anonymity threshold to check against (default: 5)
#' @param na.rm logical, remove records with NA in key variables (default: TRUE)
#' @param ... additional arguments passed to methods (currently unused)
#'
#' @return An object of class "kanonymity" containing:
#' \itemize{
#'   \item k_level: the achieved k-anonymity level (minimum equivalence class size)
#'   \item satisfies_k: logical, whether the data satisfies k-anonymity for given k
#'   \item n_violating: number of records in equivalence classes smaller than k
#'   \item pct_violating: percentage of records violating k-anonymity
#'   \item n_ec: number of equivalence classes
#'   \item ec_sizes: table of equivalence class size distribution
#'   \item violating_records: indices of records violating k-anonymity
#'   \item equivalence_classes: data frame with equivalence class details
#'   \item risk_summary: summary of re-identification risk
#' }
#'
#' @details
#' k-Anonymity (Sweeney, 2002) is a privacy model that requires each record to be
#' indistinguishable from at least k-1 other records based on quasi-identifiers.
#'
#' An equivalence class is the set of all records sharing the same values for
#' all key variables. The k-anonymity level is the size of the smallest
#' equivalence class.
#'
#' \strong{Interpretation:}
#' \itemize{
#'   \item \strong{k >= 5}: Generally considered acceptable for public release
#'   \item \strong{k >= 3}: Minimum for most applications
#'   \item \strong{k = 1}: Unique records exist - high re-identification risk
#'   \item \strong{k = 2}: Some records have only one match - elevated risk
#' }
#'
#' Records in small equivalence classes (size < k) are at higher risk of
#' re-identification and may need additional protection.
#'
#' @seealso \code{\link{ldiversity}} for l-diversity assessment,
#'   \code{\link{tcloseness}} for t-closeness assessment,
#'   \code{\link{suda}} for special uniques detection,
#'   \code{\link{individual_risk}} for probabilistic risk assessment
#'
#' @references
#' Sweeney, L. (2002). k-Anonymity: A Model for Protecting Privacy.
#' \emph{International Journal of Uncertainty, Fuzziness and Knowledge-Based Systems},
#' 10(5), 557-570.
#'
#' @author Matthias Templ
#' @family privacy-models
#' @export
#' @examples
#' # Create example data
#' set.seed(123)
#' data <- data.frame(
#'   age = sample(c("young", "middle", "old"), 100, replace = TRUE),
#'   gender = sample(c("M", "F"), 100, replace = TRUE),
#'   region = sample(c("N", "S", "E", "W"), 100, replace = TRUE),
#'   income = sample(c("low", "medium", "high"), 100, replace = TRUE)
#' )
#'
#' # Check k-anonymity
#' result <- kanonymity(data, key_vars = c("age", "gender", "region"), k = 5)
#' print(result)
#' summary(result)
#' plot(result)
#'
#' # Check with more key variables (likely lower k)
#' result2 <- kanonymity(data, key_vars = c("age", "gender", "region", "income"), k = 3)
#' print(result2)
kanonymity <- function(X, ...) {

  UseMethod("kanonymity")
}

#' @rdname kanonymity
#' @param data character, which dataset to assess: "synthetic" (default) or "original".
#'   Only used by the synth_pair method.
#' @export
kanonymity.synth_pair <- function(X, k = 5, data = c("synthetic", "original"), ...) {
  if (is.null(X$key_vars)) {
    stop("synth_pair must have 'key_vars' set for kanonymity()")
  }
  data <- match.arg(data)
  dataset <- if (data == "synthetic") X$synthetic else X$original

  kanonymity.default(
    X = dataset,
    key_vars = X$key_vars,
    k = k,
    ...
  )
}

#' @rdname kanonymity
#' @export
kanonymity.default <- function(X,
                               key_vars,
                               k = 5,
                               na.rm = TRUE,
                               ...) {

  # Input validation
 if (!is.data.frame(X)) stop("X must be a data frame.")

  # Check key variables exist
  missing_vars <- setdiff(key_vars, names(X))
  if (length(missing_vars) > 0) {
    stop(paste("Key variables missing in X:", paste(missing_vars, collapse = ", ")))
  }

  # Handle missing values
  X_keys <- X[, key_vars, drop = FALSE]
  if (na.rm) {
    complete_idx <- complete.cases(X_keys)
    if (sum(!complete_idx) > 0) {
      message(sprintf("Removing %d records with NA in key variables.", sum(!complete_idx)))
    }
    X_keys <- X_keys[complete_idx, , drop = FALSE]
    original_idx <- which(complete_idx)
  } else {
    original_idx <- seq_len(nrow(X))
  }

  if (nrow(X_keys) == 0) stop("No complete cases remaining.")

  n <- nrow(X_keys)

  # Create key signature for each record
  key_signature <- apply(X_keys, 1, paste, collapse = "|")

  # Count equivalence class sizes
  ec_table <- table(key_signature)
  ec_sizes <- as.numeric(ec_table)

  # Compute k-anonymity level
  k_level <- min(ec_sizes)

  # Check if satisfies k-anonymity
  satisfies_k <- k_level >= k

  # Find violating records
  small_ec <- names(ec_table)[ec_table < k]
  violating_mask <- key_signature %in% small_ec
  violating_idx <- original_idx[violating_mask]
  n_violating <- length(violating_idx)
  pct_violating <- 100 * n_violating / n

  # Build equivalence class summary
  ec_df <- data.frame(
    key = names(ec_table),
    size = as.numeric(ec_table),
    stringsAsFactors = FALSE
  )
  ec_df <- ec_df[order(ec_df$size), ]
  ec_df$violates_k <- ec_df$size < k
  rownames(ec_df) <- NULL

  # Size distribution
  size_dist <- table(factor(pmin(ec_sizes, 10), levels = 1:10))
  names(size_dist) <- c(as.character(1:9), "10+")

  # Risk summary
  n_unique <- sum(ec_sizes == 1)
  n_small <- sum(ec_sizes < k)
  risk_summary <- list(
    n_unique = n_unique,
    pct_unique = 100 * sum(ec_sizes == 1) / n,
    n_small_ec = n_small,
    mean_ec_size = mean(ec_sizes),
    median_ec_size = median(ec_sizes)
  )

  results <- list(
    k_level = k_level,
    k_threshold = k,
    satisfies_k = satisfies_k,
    n_violating = n_violating,
    pct_violating = pct_violating,
    n_records = n,
    n_ec = length(ec_sizes),
    ec_size_distribution = size_dist,
    violating_records = violating_idx,
    equivalence_classes = ec_df,
    risk_summary = risk_summary,
    key_vars = key_vars
  )

  class(results) <- "kanonymity"
  return(results)
}


#' Print method for kanonymity objects
#' @param x an object of class "kanonymity"
#' @param ... additional arguments (ignored)
#' @export
print.kanonymity <- function(x, ...) {
 cat("k-Anonymity Assessment\n")
  cat("======================\n\n")

  cat("Key variables:", paste(x$key_vars, collapse = ", "), "\n")
  cat("Records:", x$n_records, "| Equivalence classes:", x$n_ec, "\n\n")

  cat("Results:\n")
  cat("  Achieved k-anonymity level:", x$k_level, "\n")
  cat("  Target k:", x$k_threshold, "\n")
  cat("  Satisfies", x$k_threshold, "-anonymity:",
      ifelse(x$satisfies_k, "YES", "NO"), "\n\n")

  if (x$n_violating > 0) {
    cat("Violations:\n")
    cat("  Records violating k =", x$k_threshold, ":", x$n_violating,
        sprintf("(%.1f%%)", x$pct_violating), "\n")
    cat("  Unique records (k=1):", x$risk_summary$n_unique, "\n")
  } else {
    cat("No violations: All records satisfy", x$k_threshold, "-anonymity.\n")
  }

  cat("\nRisk Assessment:\n")
  if (x$k_level == 1) {
    cat("  HIGH RISK: Unique records exist that can be directly identified.\n")
  } else if (x$k_level < 3) {
    cat("  ELEVATED RISK: Small equivalence classes increase re-identification risk.\n")
  } else if (x$k_level < 5) {
    cat("  MODERATE RISK: Consider increasing anonymization for sensitive data.\n")
  } else {
    cat("  LOW RISK: Data has reasonable k-anonymity protection.\n")
  }

  invisible(x)
}


#' Summary method for kanonymity objects
#' @param object an object of class "kanonymity"
#' @param ... additional arguments (ignored)
#' @export
summary.kanonymity <- function(object, ...) {
  summ <- list(
    k_level = object$k_level,
    k_threshold = object$k_threshold,
    satisfies_k = object$satisfies_k,
    n_violating = object$n_violating,
    pct_violating = object$pct_violating,
    n_records = object$n_records,
    n_ec = object$n_ec,
    ec_size_distribution = object$ec_size_distribution,
    risk_summary = object$risk_summary,
    key_vars = object$key_vars,
    smallest_ec = head(object$equivalence_classes, 10)
  )

  class(summ) <- "summary.kanonymity"
  return(summ)
}


#' Print method for summary.kanonymity objects
#' @param x an object of class "summary.kanonymity"
#' @param ... additional arguments (ignored)
#' @export
print.summary.kanonymity <- function(x, ...) {
  cat("Summary: k-Anonymity Assessment\n")
  cat("================================\n\n")

  cat("Key variables:", paste(x$key_vars, collapse = ", "), "\n\n")

  cat("k-Anonymity:\n")
  cat("  Achieved level:", x$k_level, "\n")
  cat("  Target threshold:", x$k_threshold, "\n")
  cat("  Satisfies target:", ifelse(x$satisfies_k, "YES", "NO"), "\n\n")

  cat("Dataset Statistics:\n")
  cat("  Total records:", x$n_records, "\n")
  cat("  Equivalence classes:", x$n_ec, "\n")
  cat("  Mean EC size:", round(x$risk_summary$mean_ec_size, 1), "\n")
  cat("  Median EC size:", x$risk_summary$median_ec_size, "\n\n")

  cat("Equivalence Class Size Distribution:\n")
  print(x$ec_size_distribution)
  cat("\n")

  cat("Smallest Equivalence Classes:\n")
  print(x$smallest_ec[, c("key", "size", "violates_k")], row.names = FALSE)

  invisible(x)
}


#' Plot method for kanonymity objects
#' @param x an object of class "kanonymity"
#' @param y not used
#' @param ... additional arguments passed to plotting functions
#' @param which integer, which plot: 1 = EC size distribution, 2 = cumulative distribution
#' @importFrom graphics barplot hist abline par legend
#' @export
plot.kanonymity <- function(x, y = NULL, ..., which = 1) {
  show <- rep(FALSE, 2)
  show[which] <- TRUE

  if (sum(show) > 1) {
    op <- par(mfrow = c(1, 2))
    on.exit(par(op))
  }

  if (show[1]) {
    # Bar chart of EC size distribution
    sizes <- x$ec_size_distribution
    colors <- ifelse(as.numeric(names(sizes)) < x$k_threshold | names(sizes) == "10+",
                     "firebrick", "steelblue")
    colors[names(sizes) == "10+"] <- "steelblue"
    colors[1:(x$k_threshold - 1)] <- "firebrick"

    barplot(sizes,
            main = paste0("Equivalence Class Sizes (k = ", x$k_threshold, ")"),
            xlab = "Equivalence Class Size",
            ylab = "Number of Classes",
            col = colors, ...)
    abline(v = x$k_threshold - 0.5, col = "red", lty = 2, lwd = 2)
    legend("topright",
           legend = c(paste0("< k (violating)"), paste0(">= k (ok)")),
           fill = c("firebrick", "steelblue"), cex = 0.8)
  }

  if (show[2]) {
    # Cumulative distribution of EC sizes
    ec_sizes <- x$equivalence_classes$size
    sorted_sizes <- sort(ec_sizes)
    cum_pct <- 100 * seq_along(sorted_sizes) / length(sorted_sizes)

    plot(sorted_sizes, cum_pct, type = "s",
         main = "Cumulative Distribution of EC Sizes",
         xlab = "Equivalence Class Size",
         ylab = "Cumulative % of Classes",
         col = "steelblue", lwd = 2, ...)
    abline(v = x$k_threshold, col = "red", lty = 2, lwd = 2)
    abline(h = 100 * sum(ec_sizes < x$k_threshold) / length(ec_sizes),
           col = "orange", lty = 3)
  }
}
