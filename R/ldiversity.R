#' l-Diversity Assessment
#'
#' Measures l-diversity of a dataset, which extends k-anonymity by requiring
#' diversity in sensitive attribute values within each equivalence class.
#' Implements distinct l-diversity, entropy l-diversity, and recursive (c,l)-diversity.
#'
#' @param X data frame to assess, or a \code{\link{synth_pair}} object
#' @param key_vars character vector of quasi-identifier variable names
#' @param sensitive_var character, name of the sensitive attribute
#' @param l integer, the l-diversity threshold to check against (default: 2)
#' @param c numeric, the c parameter for recursive (c,l)-diversity (default: 2)
#' @param na.rm logical, remove records with NA values (default: TRUE)
#' @param ... additional arguments passed to methods (currently unused)
#'
#' @return An object of class "ldiversityRisk" containing:
#' \itemize{
#'   \item distinct_l: achieved distinct l-diversity level
#'   \item entropy_l: achieved entropy l-diversity level
#'   \item recursive_cl: whether data satisfies recursive (c,l)-diversity
#'   \item satisfies_l: logical, whether data satisfies l-diversity for given l
#'   \item n_violating_distinct: records in EC with fewer than l distinct values
#'   \item n_violating_entropy: records in EC with entropy < log(l)
#'   \item per_ec: data frame with per-equivalence-class diversity measures
#'   \item key_vars, sensitive_var: input parameters
#' }
#'
#' @details
#' l-Diversity (Machanavajjhala et al., 2007) addresses limitations of k-anonymity
#' by requiring diversity in sensitive attributes within equivalence classes.
#'
#' \strong{Distinct l-Diversity:}
#' Each equivalence class must have at least l distinct values of the sensitive
#' attribute. The achieved level is the minimum number of distinct values across
#' all equivalence classes.
#'
#' \strong{Entropy l-Diversity:}
#' Each equivalence class must have entropy of the sensitive attribute >= log(l).
#' Entropy is computed as: \eqn{H = -\sum p_i \log(p_i)}
#'
#' \strong{Recursive (c,l)-Diversity:}
#' The most common sensitive value should not appear too frequently relative
#' to other values. Specifically: \eqn{r_1 < c(r_l + r_{l+1} + ... + r_m)}
#' where r_i is the frequency of the i-th most common value.
#'
#' \strong{Interpretation:}
#' \itemize{
#'   \item \strong{l >= 3}: Good diversity protection
#'   \item \strong{l = 2}: Minimum meaningful diversity
#'   \item \strong{l = 1}: No diversity - sensitive value can be inferred
#' }
#'
#' @seealso \code{\link{kanonymity}} for k-anonymity assessment,
#'   \code{\link{tcloseness}} for t-closeness assessment,
#'   \code{\link{tcap}} for targeted correct attribution probability
#'
#' @references
#' Machanavajjhala, A., Kifer, D., Gehrke, J., & Venkitasubramaniam, M. (2007).
#' l-Diversity: Privacy Beyond k-Anonymity.
#' \emph{ACM Transactions on Knowledge Discovery from Data}, 1(1), Article 3.
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
#'   disease = sample(c("healthy", "cold", "flu", "covid"), 200, replace = TRUE,
#'                    prob = c(0.5, 0.2, 0.2, 0.1))
#' )
#'
#' # Check l-diversity
#' result <- ldiversity(data,
#'                      key_vars = c("age", "gender", "region"),
#'                      sensitive_var = "disease",
#'                      l = 2)
#' print(result)
#' summary(result)
#' plot(result)
ldiversity <- function(X, ...) {
  UseMethod("ldiversity")
}

#' @rdname ldiversity
#' @param data character, which dataset to assess: "synthetic" (default) or "original".
#'   Only used by the synth_pair method.
#' @export
ldiversity.synth_pair <- function(X, l = 2, c = 2, data = c("synthetic", "original"), ...) {
  if (is.null(X$key_vars)) {
    stop("synth_pair must have 'key_vars' set for ldiversity()")
  }
  if (is.null(X$target_var)) {
    stop("synth_pair must have 'target_var' set as sensitive variable for ldiversity()")
  }
  data <- match.arg(data)
  dataset <- if (data == "synthetic") X$synthetic else X$original

  ldiversity.default(
    X = dataset,
    key_vars = X$key_vars,
    sensitive_var = X$target_var,
    l = l,
    c = c,
    ...
  )
}

#' @rdname ldiversity
#' @export
ldiversity.default <- function(X,
                               key_vars,
                               sensitive_var,
                               l = 2,
                               c = 2,
                               na.rm = TRUE,
                               ...) {

  # Input validation
  if (!is.data.frame(X)) stop("X must be a data frame.")

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

  # Create key signature
  key_signature <- apply(X[, key_vars, drop = FALSE], 1, paste, collapse = "|")

  # Get unique equivalence classes
  unique_keys <- unique(key_signature)
  n_ec <- length(unique_keys)

  # Compute diversity for each equivalence class
  ec_results <- lapply(unique_keys, function(key) {
    mask <- key_signature == key
    sens_values <- sensitive[mask]
    ec_size <- length(sens_values)

    # Distinct l-diversity
    distinct_values <- unique(sens_values)
    n_distinct <- length(distinct_values)

    # Entropy l-diversity
    freq_table <- table(sens_values)
    props <- as.numeric(freq_table) / sum(freq_table)
    # Avoid log(0)
    props_positive <- props[props > 0]
    entropy <- -sum(props_positive * log(props_positive))
    # Convert to equivalent l (entropy >= log(l) means l <= exp(entropy))
    entropy_l_equiv <- exp(entropy)

    # Recursive (c,l)-diversity
    sorted_freqs <- sort(as.numeric(freq_table), decreasing = TRUE)
    if (length(sorted_freqs) >= l) {
      r1 <- sorted_freqs[1]
      rest_sum <- sum(sorted_freqs[l:length(sorted_freqs)])
      recursive_satisfied <- r1 < c * rest_sum
    } else {
      recursive_satisfied <- FALSE
    }

    data.frame(
      key = key,
      size = ec_size,
      n_distinct = n_distinct,
      entropy = entropy,
      entropy_l_equiv = entropy_l_equiv,
      recursive_cl = recursive_satisfied,
      violates_distinct = n_distinct < l,
      violates_entropy = entropy < log(l)
    )
  })

  per_ec <- do.call(rbind, ec_results)
  rownames(per_ec) <- NULL

  # Compute overall diversity levels
  distinct_l <- min(per_ec$n_distinct)
  entropy_l <- min(per_ec$entropy_l_equiv)
  recursive_cl_satisfied <- all(per_ec$recursive_cl)

  # Check if satisfies l-diversity
  satisfies_distinct <- distinct_l >= l
  satisfies_entropy <- all(per_ec$entropy >= log(l))

  # Count violations
  n_violating_distinct <- sum(X[key_signature %in% per_ec$key[per_ec$violates_distinct], , drop = FALSE] |> nrow(),
                              na.rm = TRUE)
  # Recalculate properly
  violating_keys_distinct <- per_ec$key[per_ec$violates_distinct]
  n_violating_distinct <- sum(key_signature %in% violating_keys_distinct)

  violating_keys_entropy <- per_ec$key[per_ec$violates_entropy]
  n_violating_entropy <- sum(key_signature %in% violating_keys_entropy)

  results <- list(
    distinct_l = distinct_l,
    entropy_l = entropy_l,
    recursive_cl = recursive_cl_satisfied,
    l_threshold = l,
    c_param = c,
    satisfies_distinct_l = satisfies_distinct,
    satisfies_entropy_l = satisfies_entropy,
    n_violating_distinct = n_violating_distinct,
    n_violating_entropy = n_violating_entropy,
    pct_violating_distinct = 100 * n_violating_distinct / n,
    pct_violating_entropy = 100 * n_violating_entropy / n,
    n_records = n,
    n_ec = n_ec,
    per_ec = per_ec,
    key_vars = key_vars,
    sensitive_var = sensitive_var
  )

  # Class is "ldiversityRisk" (not "ldiversity") to avoid an S3 dispatch
  # collision with sdcMicro::print.ldiversity when both packages are loaded.
  class(results) <- "ldiversityRisk"
  return(results)
}


#' Print method for ldiversityRisk objects
#' @param x an object of class "ldiversityRisk"
#' @param ... additional arguments (ignored)
#' @return The input object, invisibly.
#' @export
print.ldiversityRisk <- function(x, ...) {
  cat("l-Diversity Assessment\n")
  cat("======================\n\n")

  cat("Key variables:", paste(x$key_vars, collapse = ", "), "\n")
  cat("Sensitive variable:", x$sensitive_var, "\n")
  cat("Records:", x$n_records, "| Equivalence classes:", x$n_ec, "\n\n")

  cat("Target l =", x$l_threshold, "\n\n")

  cat("Distinct l-Diversity:\n")
  cat("  Achieved level:", x$distinct_l, "\n")
  cat("  Satisfies", x$l_threshold, "-diversity:",
      ifelse(x$satisfies_distinct_l, "YES", "NO"), "\n")
  if (x$n_violating_distinct > 0) {
    cat("  Violating records:", x$n_violating_distinct,
        sprintf("(%.1f%%)", x$pct_violating_distinct), "\n")
  }

  cat("\nEntropy l-Diversity:\n")
  cat("  Achieved level:", sprintf("%.2f", x$entropy_l), "\n")
  cat("  Satisfies", x$l_threshold, "-diversity:",
      ifelse(x$satisfies_entropy_l, "YES", "NO"), "\n")
  if (x$n_violating_entropy > 0) {
    cat("  Violating records:", x$n_violating_entropy,
        sprintf("(%.1f%%)", x$pct_violating_entropy), "\n")
  }

  cat("\nRecursive (c,l)-Diversity (c =", x$c_param, "):\n")
  cat("  Satisfied:", ifelse(x$recursive_cl, "YES", "NO"), "\n")

  invisible(x)
}


#' Summary method for ldiversityRisk objects
#' @param object an object of class "ldiversityRisk"
#' @param ... additional arguments (ignored)
#' @return A list of summary statistics for the corresponding object.
#' @export
summary.ldiversityRisk <- function(object, ...) {
  summ <- list(
    distinct_l = object$distinct_l,
    entropy_l = object$entropy_l,
    recursive_cl = object$recursive_cl,
    l_threshold = object$l_threshold,
    satisfies_distinct_l = object$satisfies_distinct_l,
    satisfies_entropy_l = object$satisfies_entropy_l,
    n_violating_distinct = object$n_violating_distinct,
    n_violating_entropy = object$n_violating_entropy,
    n_records = object$n_records,
    n_ec = object$n_ec,
    ec_summary = data.frame(
      metric = c("Min distinct", "Max distinct", "Mean distinct",
                 "Min entropy", "Max entropy", "Mean entropy"),
      value = c(min(object$per_ec$n_distinct), max(object$per_ec$n_distinct),
                mean(object$per_ec$n_distinct),
                min(object$per_ec$entropy), max(object$per_ec$entropy),
                mean(object$per_ec$entropy))
    ),
    worst_ec = head(object$per_ec[order(object$per_ec$n_distinct), ], 10),
    key_vars = object$key_vars,
    sensitive_var = object$sensitive_var
  )

  class(summ) <- "summary.ldiversityRisk"
  return(summ)
}


#' Print method for summary.ldiversityRisk objects
#' @param x an object of class "summary.ldiversityRisk"
#' @param ... additional arguments (ignored)
#' @return The input object, invisibly.
#' @export
print.summary.ldiversityRisk <- function(x, ...) {
  cat("Summary: l-Diversity Assessment\n")
  cat("================================\n\n")

  cat("Key variables:", paste(x$key_vars, collapse = ", "), "\n")
  cat("Sensitive variable:", x$sensitive_var, "\n")
  cat("Target l:", x$l_threshold, "\n\n")

  cat("Overall Results:\n")
  cat("  Distinct l-diversity:", x$distinct_l,
      ifelse(x$satisfies_distinct_l, "(PASS)", "(FAIL)"), "\n")
  cat("  Entropy l-diversity:", sprintf("%.2f", x$entropy_l),
      ifelse(x$satisfies_entropy_l, "(PASS)", "(FAIL)"), "\n")
  cat("  Recursive (c,l)-diversity:", ifelse(x$recursive_cl, "PASS", "FAIL"), "\n\n")

  cat("Equivalence Class Statistics:\n")
  print(x$ec_summary, row.names = FALSE)
  cat("\n")

  cat("Worst Equivalence Classes (by distinct values):\n")
  print(x$worst_ec[, c("key", "size", "n_distinct", "entropy")], row.names = FALSE)

  invisible(x)
}


#' Plot method for ldiversityRisk objects
#' @param x an object of class "ldiversityRisk"
#' @param y not used
#' @param ... additional arguments passed to plotting functions
#' @param which integer, which plot: 1 = distinct values distribution,
#'   2 = entropy distribution
#' @importFrom graphics barplot hist abline par legend
#' @return No return value; called for the side effect of producing a plot.
#' @export
plot.ldiversityRisk <- function(x, y = NULL, ..., which = 1) {
  show <- rep(FALSE, 2)
  show[which] <- TRUE

  if (sum(show) > 1) {
    op <- par(mfrow = c(1, 2))
    on.exit(par(op))
  }

  if (show[1]) {
    # Distribution of distinct values
    dist_vals <- table(factor(pmin(x$per_ec$n_distinct, 10), levels = 1:10))
    names(dist_vals) <- c(as.character(1:9), "10+")

    colors <- ifelse(1:10 < x$l_threshold, "firebrick", "steelblue")

    barplot(dist_vals,
            main = paste0("Distinct Sensitive Values per EC (l = ", x$l_threshold, ")"),
            xlab = "Number of Distinct Values",
            ylab = "Number of Equivalence Classes",
            col = colors, ...)
    legend("topright",
           legend = c(paste0("< l (violating)"), paste0(">= l (ok)")),
           fill = c("firebrick", "steelblue"), cex = 0.8)
  }

  if (show[2]) {
    # Histogram of entropy values
    hist(x$per_ec$entropy,
         main = paste0("Entropy Distribution (l = ", x$l_threshold, ")"),
         xlab = "Entropy",
         ylab = "Number of Equivalence Classes",
         col = "steelblue", border = "white", ...)
    abline(v = log(x$l_threshold), col = "red", lty = 2, lwd = 2)
    legend("topright",
           legend = paste0("log(", x$l_threshold, ") = ", round(log(x$l_threshold), 2)),
           col = "red", lty = 2, lwd = 2, cex = 0.8)
  }
}
