#' SUDA - Special Uniques Detection Algorithm
#'
#' Implements the Special Uniques Detection Algorithm (SUDA) to identify records
#' that are unique on subsets of key variables. SUDA scores indicate the risk
#' level of individual records based on their uniqueness patterns.
#'
#' @param X data frame to assess, or a \code{\link{synth_pair}} object
#' @param key_vars character vector of quasi-identifier variable names
#' @param max_msu integer, maximum size of Minimal Sample Uniques to consider
#'   (default: NULL, uses all variables up to length(key_vars))
#' @param missing numeric, value to treat as missing (default: NA handling via na.rm)
#' @param na.rm logical, remove records with NA in key variables (default: TRUE)
#' @param ... additional arguments passed to methods (currently unused)
#'
#' @return An object of class "suda" containing:
#' \itemize{
#'   \item suda_scores: numeric vector of SUDA scores for each record
#'   \item dis_scores: numeric vector of Data Intrusion Simulation (DIS) scores
#'   \item msu_counts: number of MSUs each record contributes to
#'   \item n_msu: total number of Minimal Sample Uniques found
#'   \item msu_by_size: count of MSUs by size
#'   \item high_risk_records: indices of records with elevated SUDA scores
#'   \item summary_stats: summary statistics of SUDA scores
#' }
#'
#' @details
#' SUDA (Elliot et al., 2002) identifies records that are unique on subsets of
#' key variables. A Minimal Sample Unique (MSU) is a record that is unique on
#' a set of variables, but not unique on any proper subset of those variables.
#'
#' \strong{SUDA Score:}
#' The SUDA score for a record is the sum of contributions from all MSUs it
#' belongs to, where smaller MSUs contribute more (they are more identifying):
#' \deqn{SUDA_i = \sum_{j \in MSU_i} \frac{1}{2^{|MSU_j| - 1}}}
#'
#' \strong{DIS Score (Data Intrusion Simulation):}
#' A normalized version of SUDA that accounts for the number of possible
#' variable combinations, making it more comparable across different datasets.
#'
#' \strong{Interpretation:}
#' \itemize{
#'   \item \strong{SUDA = 0}: Record is not unique on any variable subset
#'   \item \strong{SUDA > 0}: Record is unique on at least one variable subset
#'   \item \strong{Higher SUDA}: More/smaller uniqueness patterns - higher risk
#' }
#'
#' Records with high SUDA scores are at elevated risk of re-identification
#' because they can be identified using fewer variables.
#'
#' @seealso \code{\link{kanonymity}} for k-anonymity assessment,
#'   \code{\link{individual_risk}} for probabilistic risk measures
#'
#' @references
#' Elliot, M.J., Manning, A.M., Mayes, K., Gurd, J., & Bane, M. (2005).
#' SUDA: A Program for Detecting Special Uniques.
#' Joint UNECE/Eurostat Work Session on Statistical Data Confidentiality.
#'
#' Manning, A.M., Haglin, D.J., & Keane, J.A. (2008).
#' A Recursive Search Algorithm for Statistical Disclosure Assessment.
#' \emph{Data Mining and Knowledge Discovery}, 16(2), 165-196.
#'
#' @author Matthias Templ
#' @importFrom utils combn
#' @family privacy-models
#' @export
#' @examples
#' # Create example data
#' set.seed(123)
#' data <- data.frame(
#'   age = sample(c("young", "middle", "old"), 100, replace = TRUE),
#'   gender = sample(c("M", "F"), 100, replace = TRUE),
#'   region = sample(c("N", "S", "E", "W"), 100, replace = TRUE),
#'   education = sample(c("low", "medium", "high"), 100, replace = TRUE)
#' )
#'
#' # Run SUDA analysis
#' result <- suda(data, key_vars = c("age", "gender", "region", "education"))
#' print(result)
#' summary(result)
#' plot(result)
#'
#' # Identify high-risk records
#' high_risk <- which(result$suda_scores > quantile(result$suda_scores, 0.95))
#' data[high_risk, ]
suda <- function(X, ...) {
  UseMethod("suda")
}

#' @rdname suda
#' @param data character, which dataset to assess: "synthetic" (default) or "original".
#'   Only used by the synth_pair method.
#' @export
suda.synth_pair <- function(X, max_msu = NULL, data = c("synthetic", "original"), ...) {
  if (is.null(X$key_vars)) {
    stop("synth_pair must have 'key_vars' set for suda()")
  }
  data <- match.arg(data)
  dataset <- if (data == "synthetic") X$synthetic else X$original

  suda.default(
    X = dataset,
    key_vars = X$key_vars,
    max_msu = max_msu,
    ...
  )
}

#' @rdname suda
#' @export
suda.default <- function(X,
                         key_vars,
                         max_msu = NULL,
                         missing = NA,
                         na.rm = TRUE,
                         ...) {

  # Input validation
  if (!is.data.frame(X)) stop("X must be a data frame.")

  # Check key variables exist
  missing_vars <- setdiff(key_vars, names(X))
  if (length(missing_vars) > 0) {
    stop(paste("Key variables missing in X:", paste(missing_vars, collapse = ", ")))
  }

  # Subset to key variables
  X_keys <- X[, key_vars, drop = FALSE]

  # Handle missing values
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
  p <- length(key_vars)

  # Set max_msu
  if (is.null(max_msu)) {
    max_msu <- p
  } else {
    max_msu <- min(max_msu, p)
  }

  # Initialize scores
  suda_scores <- numeric(n)
  msu_counts <- integer(n)
  dis_scores <- numeric(n)

  # Track MSUs found
  msu_list <- list()
  msu_by_size <- integer(max_msu)
  names(msu_by_size) <- as.character(1:max_msu)

  # Convert to character for consistent comparison
  X_char <- as.data.frame(lapply(X_keys, as.character), stringsAsFactors = FALSE)

  # For each subset size (starting from smallest = most dangerous)
  for (size in 1:max_msu) {
    # Generate all combinations of 'size' variables
    var_combos <- combn(key_vars, size, simplify = FALSE)

    for (vars in var_combos) {
      # Create key for this variable subset
      if (length(vars) == 1) {
        subset_keys <- X_char[[vars]]
      } else {
        subset_keys <- apply(X_char[, vars, drop = FALSE], 1, paste, collapse = "|")
      }

      # Find unique values (sample uniques)
      key_counts <- table(subset_keys)
      unique_keys <- names(key_counts)[key_counts == 1]

      if (length(unique_keys) > 0) {
        # Find which records are unique on this subset
        unique_mask <- subset_keys %in% unique_keys

        for (idx in which(unique_mask)) {
          # Check if this is a MINIMAL sample unique
          # (i.e., not unique on any proper subset)
          is_minimal <- TRUE

          if (size > 1) {
            # Check all smaller subsets
            for (smaller_size in 1:(size - 1)) {
              smaller_combos <- combn(vars, smaller_size, simplify = FALSE)
              for (smaller_vars in smaller_combos) {
                if (length(smaller_vars) == 1) {
                  smaller_key <- X_char[[smaller_vars]][idx]
                  smaller_keys <- X_char[[smaller_vars]]
                } else {
                  smaller_key <- paste(X_char[idx, smaller_vars], collapse = "|")
                  smaller_keys <- apply(X_char[, smaller_vars, drop = FALSE], 1, paste, collapse = "|")
                }
                if (sum(smaller_keys == smaller_key) == 1) {
                  # Already unique on smaller subset - not minimal
                  is_minimal <- FALSE
                  break
                }
              }
              if (!is_minimal) break
            }
          }

          if (is_minimal) {
            # This is an MSU - add contribution to SUDA score
            # Smaller MSUs contribute more (more identifying)
            contribution <- 1 / (2^(size - 1))
            suda_scores[idx] <- suda_scores[idx] + contribution
            msu_counts[idx] <- msu_counts[idx] + 1
            msu_by_size[size] <- msu_by_size[size] + 1

            # Store MSU info
            msu_list[[length(msu_list) + 1]] <- list(
              record_idx = original_idx[idx],
              variables = vars,
              size = size
            )
          }
        }
      }
    }
  }

  # Compute DIS scores (normalized SUDA)
  # DIS accounts for the number of possible variable combinations
  if (max(suda_scores) > 0) {
    # Normalize by maximum possible score
    max_possible <- sum(sapply(1:max_msu, function(k) choose(p, k) / 2^(k-1)))
    dis_scores <- suda_scores / max_possible
  }

  # Map back to original indices if we removed NAs
  full_suda_scores <- rep(0, nrow(X))
  full_suda_scores[original_idx] <- suda_scores
  full_msu_counts <- rep(0L, nrow(X))
  full_msu_counts[original_idx] <- msu_counts
  full_dis_scores <- rep(0, nrow(X))
  full_dis_scores[original_idx] <- dis_scores

  # Identify high-risk records (top 5% or SUDA > 0)
  if (sum(full_suda_scores > 0) > 0) {
    threshold <- quantile(full_suda_scores[full_suda_scores > 0], 0.75)
    high_risk_idx <- which(full_suda_scores >= threshold)
  } else {
    high_risk_idx <- integer(0)
  }

  # Summary statistics
  summary_stats <- list(
    n_with_msu = sum(full_suda_scores > 0),
    pct_with_msu = 100 * sum(full_suda_scores > 0) / nrow(X),
    mean_suda = mean(full_suda_scores),
    median_suda = median(full_suda_scores),
    max_suda = max(full_suda_scores),
    mean_msu_count = mean(full_msu_counts),
    max_msu_count = max(full_msu_counts)
  )

  results <- list(
    suda_scores = full_suda_scores,
    dis_scores = full_dis_scores,
    msu_counts = full_msu_counts,
    n_msu = sum(msu_by_size),
    msu_by_size = msu_by_size,
    high_risk_records = high_risk_idx,
    n_records = nrow(X),
    n_complete = n,
    summary_stats = summary_stats,
    max_msu = max_msu,
    key_vars = key_vars
  )

  class(results) <- "suda"
  return(results)
}


#' Print method for suda objects
#' @param x an object of class "suda"
#' @param ... additional arguments (ignored)
#' @export
print.suda <- function(x, ...) {
  cat("SUDA - Special Uniques Detection Algorithm\n")
  cat("==========================================\n\n")

  cat("Key variables:", paste(x$key_vars, collapse = ", "), "\n")
  cat("Records:", x$n_records, "| Max MSU size:", x$max_msu, "\n\n")

  cat("Results:\n")
  cat("  Total MSUs found:", x$n_msu, "\n")
  cat("  Records with MSU:", x$summary_stats$n_with_msu,
      sprintf("(%.1f%%)", x$summary_stats$pct_with_msu), "\n")
  cat("  High-risk records:", length(x$high_risk_records), "\n\n")

  cat("MSUs by Size:\n")
  for (i in seq_along(x$msu_by_size)) {
    if (x$msu_by_size[i] > 0) {
      cat(sprintf("  Size %d: %d MSUs\n", i, x$msu_by_size[i]))
    }
  }

  cat("\nSUDA Score Statistics:\n")
  cat("  Mean:", sprintf("%.4f", x$summary_stats$mean_suda), "\n")
  cat("  Median:", sprintf("%.4f", x$summary_stats$median_suda), "\n")
  cat("  Max:", sprintf("%.4f", x$summary_stats$max_suda), "\n")

  cat("\nRisk Assessment:\n")
  if (x$msu_by_size[1] > 0) {
    cat("  HIGH RISK: Records with single-variable uniqueness detected.\n")
  } else if (sum(x$msu_by_size[1:2]) > 0) {
    cat("  ELEVATED RISK: Records unique on small variable subsets.\n")
  } else if (x$n_msu > 0) {
    cat("  MODERATE RISK: Some MSUs found but require many variables.\n")
  } else {
    cat("  LOW RISK: No MSUs detected.\n")
  }

  invisible(x)
}


#' Summary method for suda objects
#' @param object an object of class "suda"
#' @param ... additional arguments (ignored)
#' @export
summary.suda <- function(object, ...) {
  # Get top risky records
  top_n <- min(10, length(object$high_risk_records))
  if (top_n > 0) {
    top_idx <- order(object$suda_scores, decreasing = TRUE)[1:top_n]
    top_records <- data.frame(
      record = top_idx,
      suda_score = object$suda_scores[top_idx],
      dis_score = object$dis_scores[top_idx],
      msu_count = object$msu_counts[top_idx]
    )
  } else {
    top_records <- NULL
  }

  summ <- list(
    n_msu = object$n_msu,
    msu_by_size = object$msu_by_size,
    summary_stats = object$summary_stats,
    suda_quantiles = quantile(object$suda_scores, probs = c(0, 0.25, 0.5, 0.75, 0.9, 0.95, 1)),
    top_records = top_records,
    key_vars = object$key_vars,
    n_records = object$n_records
  )

  class(summ) <- "summary.suda"
  return(summ)
}


#' Print method for summary.suda objects
#' @param x an object of class "summary.suda"
#' @param ... additional arguments (ignored)
#' @export
print.summary.suda <- function(x, ...) {
  cat("Summary: SUDA Analysis\n")
  cat("======================\n\n")

  cat("Key variables:", paste(x$key_vars, collapse = ", "), "\n")
  cat("Records:", x$n_records, "\n\n")

  cat("MSU Summary:\n")
  cat("  Total MSUs:", x$n_msu, "\n")
  cat("  Records with MSU:", x$summary_stats$n_with_msu,
      sprintf("(%.1f%%)", x$summary_stats$pct_with_msu), "\n\n")

  cat("MSUs by Size:\n")
  print(x$msu_by_size[x$msu_by_size > 0])
  cat("\n")

  cat("SUDA Score Distribution:\n")
  print(round(x$suda_quantiles, 4))
  cat("\n")

  if (!is.null(x$top_records)) {
    cat("Highest Risk Records:\n")
    print(x$top_records, row.names = FALSE)
  }

  invisible(x)
}


#' Plot method for suda objects
#' @param x an object of class "suda"
#' @param y not used
#' @param ... additional arguments passed to plotting functions
#' @param which integer, which plot: 1 = SUDA score distribution, 2 = MSU count by size
#' @importFrom graphics barplot hist abline par
#' @export
plot.suda <- function(x, y = NULL, ..., which = 1) {
  show <- rep(FALSE, 2)
  show[which] <- TRUE

  if (sum(show) > 1) {
    op <- par(mfrow = c(1, 2))
    on.exit(par(op))
  }

  if (show[1]) {
    # Histogram of SUDA scores
    scores <- x$suda_scores[x$suda_scores > 0]
    if (length(scores) > 0) {
      hist(scores,
           main = "Distribution of SUDA Scores (non-zero)",
           xlab = "SUDA Score",
           ylab = "Number of Records",
           col = "steelblue", border = "white", ...)
      abline(v = mean(scores), col = "red", lty = 2, lwd = 2)
    } else {
      plot(1, type = "n", main = "No MSUs Found",
           xlab = "", ylab = "", axes = FALSE)
      text(1, 1, "No records with SUDA > 0")
    }
  }

  if (show[2]) {
    # Bar chart of MSUs by size
    msu_counts <- x$msu_by_size
    if (sum(msu_counts) > 0) {
      colors <- c("firebrick", "orange", rep("steelblue", length(msu_counts) - 2))
      colors <- colors[1:length(msu_counts)]

      barplot(msu_counts,
              main = "MSUs by Variable Subset Size",
              xlab = "Number of Variables",
              ylab = "Number of MSUs",
              col = colors, ...)
    } else {
      plot(1, type = "n", main = "No MSUs Found",
           xlab = "", ylab = "", axes = FALSE)
      text(1, 1, "No MSUs detected")
    }
  }
}
