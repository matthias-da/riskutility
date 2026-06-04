#' t-Closeness Assessment
#'
#' Measures t-closeness of a dataset, which extends l-diversity by requiring
#' that the distribution of a sensitive attribute within each equivalence class
#' is close to the overall distribution, using Earth Mover's Distance (EMD).
#'
#' @param X data frame to assess, or a \code{\link{synth_pair}} object
#' @param key_vars character vector of quasi-identifier variable names
#' @param sensitive_var character, name of the sensitive attribute
#' @param t numeric, the t-closeness threshold (default: 0.2)
#' @param na.rm logical, remove records with NA values (default: TRUE)
#' @param ... additional arguments passed to methods (currently unused)
#'
#' @return An object of class "tcloseness" containing:
#' \itemize{
#'   \item t_achieved: maximum EMD across all equivalence classes
#'   \item t_threshold: user-specified t
#'   \item satisfies_t: logical, whether data satisfies t-closeness
#'   \item n_violating: number of records in violating equivalence classes
#'   \item pct_violating: percentage of records violating t-closeness
#'   \item n_records: total number of records
#'   \item n_ec: number of equivalence classes
#'   \item per_ec: data.frame with key, size, emd, violates_t per equivalence class
#'   \item sensitive_type: "numeric" or "categorical"
#'   \item key_vars, sensitive_var: input parameters
#' }
#'
#' @details
#' t-Closeness (Li et al., 2007) addresses limitations of l-diversity by
#' requiring that the distribution of a sensitive attribute within each
#' equivalence class is close to the attribute's overall distribution.
#'
#' The distance is measured using Earth Mover's Distance (EMD), also known
#' as the Wasserstein distance:
#'
#' \strong{Numeric attributes:}
#' EMD is computed as the normalized sum of absolute CDF differences over
#' ordered unique values:
#' \deqn{EMD = \frac{1}{m - 1} \sum_{i=1}^{m} |F_{EC}(x_i) - F(x_i)|}
#' where \eqn{m} is the number of unique values, \eqn{F_{EC}} is the empirical
#' CDF within the equivalence class, and \eqn{F} is the overall empirical CDF.
#' This CDF-difference form corresponds to the ordered-distance EMD of Li et al.
#' (2007) up to the choice of normalisation (here \eqn{1/(m-1)} over the observed
#' support); it is a close approximation rather than the exact published
#' per-rank-weighted quantity.
#'
#' \strong{Categorical attributes:}
#' EMD is the variational distance (half the L1 distance between distributions):
#' \deqn{EMD = \frac{1}{2} \sum_{i} |p_{EC,i} - p_i|}
#'
#' \strong{Interpretation:}
#' \itemize{
#'   \item \strong{t close to 0}: EC distributions nearly identical to overall
#'   \item \strong{t = 0.2}: Common threshold for reasonable protection
#'   \item \strong{t close to 1}: Weak protection, distributions may differ greatly
#' }
#'
#' @seealso \code{\link{kanonymity}} for k-anonymity assessment,
#'   \code{\link{ldiversity}} for l-diversity assessment
#'
#' @references
#' Li, N., Li, T. & Venkatasubramanian, S. (2007).
#' t-Closeness: Privacy Beyond k-Anonymity and l-Diversity.
#' \emph{IEEE 23rd International Conference on Data Engineering}, 106--115.
#' \doi{10.1109/ICDE.2007.367856}
#'
#' @author Matthias Templ
#' @family privacy-models
#' @export
#' @examples
#' # Create example data
#' set.seed(123)
#' data <- data.frame(
#'   age = sample(c("young", "middle", "old"), 200, replace = TRUE),
#'   gender = sample(c("M", "F"), 200, replace = TRUE),
#'   region = sample(c("N", "S"), 200, replace = TRUE),
#'   salary = rnorm(200, mean = 50000, sd = 15000)
#' )
#'
#' # Check t-closeness with numeric sensitive attribute
#' result <- tcloseness(data,
#'                      key_vars = c("age", "gender", "region"),
#'                      sensitive_var = "salary",
#'                      t = 0.2)
#' print(result)
#' summary(result)
#' plot(result)
#'
#' # Check with categorical sensitive attribute
#' data$disease <- sample(c("healthy", "cold", "flu", "covid"), 200,
#'                        replace = TRUE, prob = c(0.5, 0.2, 0.2, 0.1))
#' result2 <- tcloseness(data,
#'                       key_vars = c("age", "gender", "region"),
#'                       sensitive_var = "disease",
#'                       t = 0.3)
#' print(result2)
tcloseness <- function(X, ...) {
  UseMethod("tcloseness")
}

#' @rdname tcloseness
#' @param data character, which dataset to assess: "synthetic" (default) or "original".
#'   Only used by the synth_pair method.
#' @export
tcloseness.synth_pair <- function(X, t = 0.2, data = c("synthetic", "original"), ...) {
  if (is.null(X$key_vars)) {
    stop("synth_pair must have 'key_vars' set for tcloseness()")
  }
  if (is.null(X$target_var)) {
    stop("synth_pair must have 'target_var' set as sensitive variable for tcloseness()")
  }
  data <- match.arg(data)
  dataset <- if (data == "synthetic") X$synthetic else X$original

  tcloseness.default(
    X = dataset,
    key_vars = X$key_vars,
    sensitive_var = X$target_var,
    t = t,
    ...
  )
}

#' @rdname tcloseness
#' @export
tcloseness.default <- function(X,
                               key_vars,
                               sensitive_var,
                               t = 0.2,
                               na.rm = TRUE,
                               ...) {

  # Input validation
  if (!is.data.frame(X)) stop("X must be a data frame.")

  if (!is.numeric(t) || length(t) != 1 || t < 0 || t > 1) {
    stop("t must be a single numeric value between 0 and 1.")
  }

  # Check variables exist
  all_vars <- c(key_vars, sensitive_var)
  missing_vars <- setdiff(all_vars, names(X))
  if (length(missing_vars) > 0) {
    stop(paste("Variables missing in X:", paste(missing_vars, collapse = ", ")))
  }

  # Handle missing values
  if (na.rm) {
    complete_idx <- complete.cases(X[, all_vars, drop = FALSE])
    if (sum(!complete_idx) > 0) {
      message(sprintf("Removing %d records with NA values.", sum(!complete_idx)))
    }
    X <- X[complete_idx, , drop = FALSE]
  }

  if (nrow(X) == 0) stop("No complete cases remaining.")

  n <- nrow(X)
  sensitive <- X[[sensitive_var]]

  # Determine sensitive type
  sensitive_type <- if (is.numeric(sensitive)) "numeric" else "categorical"

  # Pre-compute global distribution
  if (sensitive_type == "numeric") {
    ordered_vals <- sort(unique(sensitive))
    m <- length(ordered_vals)
    # Global CDF at each ordered unique value
    global_ecdf <- ecdf(sensitive)
    global_cdf <- global_ecdf(ordered_vals)
  } else {
    # Global proportions for all categories
    global_tab <- table(sensitive)
    global_props <- as.numeric(global_tab) / sum(global_tab)
    names(global_props) <- names(global_tab)
  }

  # Create key signature
  key_signature <- apply(X[, key_vars, drop = FALSE], 1, paste, collapse = "|")

  # Get unique equivalence classes
  unique_keys <- unique(key_signature)
  n_ec <- length(unique_keys)

  # Compute EMD for each equivalence class
  ec_results <- lapply(unique_keys, function(key) {
    mask <- key_signature == key
    sens_values <- sensitive[mask]
    ec_size <- length(sens_values)

    if (sensitive_type == "numeric") {
      # EMD for numeric: normalized sum of absolute CDF differences
      ec_ecdf <- ecdf(sens_values)
      ec_cdf <- ec_ecdf(ordered_vals)
      if (m > 1) {
        emd <- sum(abs(ec_cdf - global_cdf)) / (m - 1)
      } else {
        emd <- 0
      }
    } else {
      # EMD for categorical: variational distance
      ec_tab <- table(factor(sens_values, levels = names(global_props)))
      ec_props <- as.numeric(ec_tab) / sum(ec_tab)
      emd <- 0.5 * sum(abs(ec_props - global_props))
    }

    data.frame(
      key = key,
      size = ec_size,
      emd = emd,
      violates_t = emd > t
    )
  })

  per_ec <- do.call(rbind, ec_results)
  rownames(per_ec) <- NULL

  # Overall t-closeness
  t_achieved <- max(per_ec$emd)
  satisfies_t <- t_achieved <= t

  # Count violations
  violating_keys <- per_ec$key[per_ec$violates_t]
  n_violating <- sum(key_signature %in% violating_keys)
  pct_violating <- 100 * n_violating / n

  results <- list(
    t_achieved = t_achieved,
    t_threshold = t,
    satisfies_t = satisfies_t,
    n_violating = n_violating,
    pct_violating = pct_violating,
    n_records = n,
    n_ec = n_ec,
    per_ec = per_ec,
    sensitive_type = sensitive_type,
    key_vars = key_vars,
    sensitive_var = sensitive_var
  )

  class(results) <- "tcloseness"
  return(results)
}


#' Print method for tcloseness objects
#' @param x an object of class "tcloseness"
#' @param ... additional arguments (ignored)
#' @return The input object, invisibly.
#' @export
print.tcloseness <- function(x, ...) {
  cat("t-Closeness Assessment\n")
  cat("======================\n\n")

  cat("Key variables:", paste(x$key_vars, collapse = ", "), "\n")
  cat("Sensitive variable:", x$sensitive_var,
      paste0("(", x$sensitive_type, ")"), "\n")
  cat("Records:", x$n_records, "| Equivalence classes:", x$n_ec, "\n\n")

  cat("Target t =", x$t_threshold, "\n\n")

  cat("Results:\n")
  cat("  Maximum EMD:", sprintf("%.4f", x$t_achieved), "\n")
  cat("  Satisfies t-closeness:",
      ifelse(x$satisfies_t, "YES", "NO"), "\n")

  if (x$n_violating > 0) {
    cat("  Violating records:", x$n_violating,
        sprintf("(%.1f%%)", x$pct_violating), "\n")
  }

  cat("\nInterpretation:\n")
  if (x$satisfies_t) {
    cat("  All equivalence class distributions are within t =",
        x$t_threshold, "of the overall distribution.\n")
  } else {
    n_violating_ec <- sum(x$per_ec$violates_t)
    cat("  ", n_violating_ec, "of", x$n_ec,
        "equivalence classes exceed the threshold.\n")
  }

  invisible(x)
}


#' Summary method for tcloseness objects
#' @param object an object of class "tcloseness"
#' @param ... additional arguments (ignored)
#' @return A list of summary statistics for the corresponding object.
#' @export
summary.tcloseness <- function(object, ...) {
  summ <- list(
    t_achieved = object$t_achieved,
    t_threshold = object$t_threshold,
    satisfies_t = object$satisfies_t,
    n_violating = object$n_violating,
    pct_violating = object$pct_violating,
    n_records = object$n_records,
    n_ec = object$n_ec,
    sensitive_type = object$sensitive_type,
    emd_summary = data.frame(
      metric = c("Min", "Q1", "Median", "Mean", "Q3", "Max"),
      value = as.numeric(summary(object$per_ec$emd)[1:6])
    ),
    worst_ec = head(object$per_ec[order(-object$per_ec$emd), ], 10),
    key_vars = object$key_vars,
    sensitive_var = object$sensitive_var
  )

  class(summ) <- "summary.tcloseness"
  return(summ)
}


#' Print method for summary.tcloseness objects
#' @param x an object of class "summary.tcloseness"
#' @param ... additional arguments (ignored)
#' @return The input object, invisibly.
#' @export
print.summary.tcloseness <- function(x, ...) {
  cat("Summary: t-Closeness Assessment\n")
  cat("================================\n\n")

  cat("Key variables:", paste(x$key_vars, collapse = ", "), "\n")
  cat("Sensitive variable:", x$sensitive_var,
      paste0("(", x$sensitive_type, ")"), "\n")
  cat("Target t:", x$t_threshold, "\n\n")

  cat("Overall Result:\n")
  cat("  Maximum EMD:", sprintf("%.4f", x$t_achieved),
      ifelse(x$satisfies_t, "(PASS)", "(FAIL)"), "\n")
  cat("  Violating records:", x$n_violating,
      sprintf("(%.1f%%)", x$pct_violating), "\n\n")

  cat("EMD Distribution Across Equivalence Classes:\n")
  print(x$emd_summary, row.names = FALSE)
  cat("\n")

  cat("Worst Equivalence Classes (by EMD):\n")
  worst <- x$worst_ec[, c("key", "size", "emd", "violates_t")]
  worst$emd <- sprintf("%.4f", worst$emd)
  print(worst, row.names = FALSE)

  invisible(x)
}


#' Plot method for tcloseness objects
#' @param x an object of class "tcloseness"
#' @param y not used
#' @param ... additional arguments passed to plotting functions
#' @param which integer, which plot: 1 = EMD distribution histogram,
#'   2 = EC size vs EMD scatter
#' @importFrom graphics barplot hist abline par legend points
#' @return No return value; called for the side effect of producing a plot.
#' @export
plot.tcloseness <- function(x, y = NULL, ..., which = 1) {
  show <- rep(FALSE, 2)
  show[which] <- TRUE

  if (sum(show) > 1) {
    op <- par(mfrow = c(1, 2))
    on.exit(par(op))
  }

  if (show[1]) {
    # Histogram of EMD values
    hist(x$per_ec$emd,
         main = paste0("EMD Distribution (t = ", x$t_threshold, ")"),
         xlab = "Earth Mover's Distance",
         ylab = "Number of Equivalence Classes",
         col = "steelblue", border = "white", ...)
    abline(v = x$t_threshold, col = "red", lty = 2, lwd = 2)
    legend("topright",
           legend = paste0("t = ", x$t_threshold),
           col = "red", lty = 2, lwd = 2, cex = 0.8)
  }

  if (show[2]) {
    # EC size vs EMD scatter
    cols <- ifelse(x$per_ec$violates_t, "firebrick", "steelblue")
    plot(x$per_ec$size, x$per_ec$emd,
         main = "EC Size vs EMD",
         xlab = "Equivalence Class Size",
         ylab = "Earth Mover's Distance",
         col = cols, pch = 19, ...)
    abline(h = x$t_threshold, col = "red", lty = 2, lwd = 2)
    legend("topright",
           legend = c("Violating", "OK"),
           col = c("firebrick", "steelblue"),
           pch = 19, cex = 0.8)
  }
}
