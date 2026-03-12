#' Attacker Risk Models (Prosecutor/Journalist/Marketer)
#'
#' Computes re-identification risk under three canonical attacker models
#' from Statistical Disclosure Control (SDC). Each model makes different
#' assumptions about the attacker's background knowledge, leading to
#' different risk quantifications based on quasi-identifier (QI) frequencies.
#'
#' @param X data frame to assess, or a \code{\link{synth_pair}} object
#' @param key_vars character vector of quasi-identifier variable names
#' @param sampling_fraction numeric, the fraction of the population
#'   represented by the data (default: 0.01). Used for the journalist model
#'   to estimate population frequencies.
#' @param model character, which attacker model(s) to compute. One of
#'   "prosecutor", "journalist", "marketer", or "all" (default: "all")
#' @param na.rm logical, remove records with NA in key variables (default: TRUE)
#' @param ... additional arguments passed to methods (currently unused)
#'
#' @return An object of class "attacker_risk" containing:
#' \itemize{
#'   \item risk_per_record: data.frame with columns: prosecutor, journalist (if applicable), freq (f_k)
#'   \item global_risk: named list with prosecutor (mean 1/f_k), journalist (mean 1/F_k), marketer (1/n * sum(1/f_k))
#'   \item n_uniques: count of records with f_k = 1
#'   \item pct_uniques: fraction of unique records
#'   \item freq_table: table of QI combination frequencies
#'   \item model: model(s) used
#'   \item key_vars: quasi-identifier variables used
#'   \item sampling_fraction: sampling fraction used
#'   \item n_records: number of records assessed
#'   \item privacy_pass: logical, prosecutor global risk <= 0.1
#' }
#'
#' @details
#' The three attacker models differ in their assumptions:
#'
#' \strong{Prosecutor Model:}
#' The attacker knows a specific individual is in the dataset and attempts
#' to re-identify them. For each record in equivalence class k with
#' frequency f_k, the individual risk is:
#' \deqn{risk_i = 1 / f_k}
#' The global prosecutor risk is the mean over all records.
#'
#' \strong{Journalist Model:}
#' The attacker does not know whether the target is in the dataset. This
#' requires estimating the population frequency F_k from the sample
#' frequency f_k using the sampling fraction:
#' \deqn{F_k \approx f_k / sampling\_fraction}
#' \deqn{risk_i = 1 / F_k = sampling\_fraction / f_k}
#' This typically yields lower risks than the prosecutor model because
#' population frequencies are larger than sample frequencies.
#'
#' \strong{Marketer Model:}
#' The attacker does not target any specific individual but wants to
#' re-identify as many records as possible. The global risk is the
#' expected number of successful re-identifications divided by the
#' total number of records:
#' \deqn{risk = (1/n) \sum_{i=1}^{n} 1/f_k(i)}
#' Note that for the marketer model, the global risk equals the
#' prosecutor global risk (both are mean 1/f_k).
#'
#' \strong{Interpretation:}
#' \itemize{
#'   \item Risk close to 1: High re-identification risk (many unique records)
#'   \item Risk > 0.1: Elevated risk, privacy protection may be insufficient
#'   \item Risk < 0.05: Generally acceptable
#'   \item Risk close to 0: Low re-identification risk (large equivalence classes)
#' }
#'
#' @seealso \code{\link{kanonymity}} for k-anonymity assessment,
#'   \code{\link{individual_risk}} for model-based individual risk,
#'   \code{\link{suda}} for special uniques detection
#'
#' @references
#' Hundepool, A., Domingo-Ferrer, J., Franconi, L., Giessing, S.,
#' Schulte Nordholt, E., Spicer, K. & de Wolf, P.-P. (2012).
#' \emph{Statistical Disclosure Control}. John Wiley & Sons.
#' \doi{10.1002/9781118348239}
#'
#' Elliot, M. (2006). A Computational Algorithm for Handling the Special
#' Uniques Problem. \emph{International Journal of Uncertainty, Fuzziness
#' and Knowledge-Based Systems}, 14(Supp. 1), 45-59.
#'
#' Templ, M. (2017). \emph{Statistical Disclosure Control for Microdata:
#' Methods and Applications in R}. Springer.
#' \doi{10.1007/978-3-319-50272-4}
#'
#' @author Matthias Templ
#' @family privacy-models
#' @export
#' @importFrom stats complete.cases median quantile
#' @importFrom graphics barplot hist abline par legend text
#' @examples
#' # Create example data
#' set.seed(123)
#' data <- data.frame(
#'   age = sample(c("young", "middle", "old"), 200, replace = TRUE),
#'   gender = sample(c("M", "F"), 200, replace = TRUE),
#'   region = sample(c("N", "S", "E", "W"), 200, replace = TRUE),
#'   income = sample(c("low", "medium", "high"), 200, replace = TRUE)
#' )
#'
#' # Compute all three attacker models
#' result <- attacker_risk(data, key_vars = c("age", "gender", "region"))
#' print(result)
#' summary(result)
#'
#' \donttest{
#' plot(result)
#'
#' # Prosecutor model only
#' result_p <- attacker_risk(data, key_vars = c("age", "gender", "region"),
#'                           model = "prosecutor")
#' print(result_p)
#'
#' # Journalist model with different sampling fraction
#' result_j <- attacker_risk(data, key_vars = c("age", "gender"),
#'                           model = "journalist", sampling_fraction = 0.05)
#' print(result_j)
#' }
attacker_risk <- function(X, ...) {
  UseMethod("attacker_risk")
}

#' @rdname attacker_risk
#' @param data character, which dataset to assess: "synthetic" (default) or "original".
#'   Only used by the synth_pair method.
#' @export
attacker_risk.synth_pair <- function(X, data = c("synthetic", "original"), ...) {
  if (is.null(X$key_vars)) {
    stop("synth_pair must have 'key_vars' set for attacker_risk()")
  }
  data <- match.arg(data)
  dataset <- if (data == "synthetic") X$synthetic else X$original

  attacker_risk.default(
    X = dataset,
    key_vars = X$key_vars,
    ...
  )
}

#' @rdname attacker_risk
#' @export
attacker_risk.default <- function(X,
                                  key_vars,
                                  sampling_fraction = 0.01,
                                  model = c("all", "prosecutor",
                                            "journalist", "marketer"),
                                  na.rm = TRUE,
                                  ...) {

  # Input validation
  if (!is.data.frame(X)) stop("X must be a data frame.")

  model <- match.arg(model)

  # Check key variables exist
  missing_vars <- setdiff(key_vars, names(X))
  if (length(missing_vars) > 0) {
    stop(paste("Key variables missing in X:",
               paste(missing_vars, collapse = ", ")))
  }

  if (!is.numeric(sampling_fraction) || length(sampling_fraction) != 1 ||
      sampling_fraction <= 0 || sampling_fraction > 1) {
    stop("sampling_fraction must be a single numeric value in (0, 1].")
  }

  # Handle missing values
  X_keys <- X[, key_vars, drop = FALSE]
  if (na.rm) {
    complete_idx <- complete.cases(X_keys)
    if (sum(!complete_idx) > 0) {
      message(sprintf("Removing %d records with NA in key variables.",
                      sum(!complete_idx)))
    }
    X_keys <- X_keys[complete_idx, , drop = FALSE]
  }

  if (nrow(X_keys) == 0) stop("No complete cases remaining.")

  n <- nrow(X_keys)

  # Create key signature for each record
  key_signature <- apply(X_keys, 1, paste, collapse = "|")

  # Count equivalence class frequencies
  freq_table <- table(key_signature)
  ec_sizes <- as.numeric(freq_table)

  # Map each record to its frequency f_k
  f_k <- as.numeric(freq_table[key_signature])

  # Determine which models to compute
  compute_prosecutor <- model %in% c("all", "prosecutor")
  compute_journalist <- model %in% c("all", "journalist")
  compute_marketer   <- model %in% c("all", "marketer")

  # --- Per-record risks ---
  risk_per_record <- data.frame(freq = f_k)

  if (compute_prosecutor || compute_marketer) {
    risk_per_record$prosecutor <- 1 / f_k
  }

  if (compute_journalist) {
    # F_k = f_k / sampling_fraction (estimated population frequency)
    F_k <- f_k / sampling_fraction
    risk_per_record$journalist <- 1 / F_k  # = sampling_fraction / f_k
  }

  # --- Global risks ---
  global_risk <- list()

  if (compute_prosecutor) {
    global_risk$prosecutor <- mean(1 / f_k)
  }

  if (compute_journalist) {
    global_risk$journalist <- mean(sampling_fraction / f_k)
  }

  if (compute_marketer) {
    # Marketer: (1/n) * sum(1/f_k)
    global_risk$marketer <- (1 / n) * sum(1 / f_k)
  }

  # --- Unique records ---
  n_uniques <- sum(ec_sizes == 1)
  pct_uniques <- n_uniques / length(ec_sizes)

  # --- Privacy pass: based on prosecutor global risk ---
  if (compute_prosecutor) {
    privacy_pass <- global_risk$prosecutor <= 0.1
  } else if (compute_journalist) {
    privacy_pass <- global_risk$journalist <= 0.1
  } else {
    privacy_pass <- global_risk$marketer <= 0.1
  }

  results <- list(
    risk_per_record = risk_per_record,
    global_risk = global_risk,
    n_uniques = n_uniques,
    pct_uniques = pct_uniques,
    freq_table = freq_table,
    model = model,
    key_vars = key_vars,
    sampling_fraction = sampling_fraction,
    n_records = n,
    privacy_pass = privacy_pass
  )

  class(results) <- "attacker_risk"
  return(results)
}


# --- S3 methods ---

#' Print method for attacker_risk objects
#' @param x an object of class "attacker_risk"
#' @param ... additional arguments (ignored)
#' @export
print.attacker_risk <- function(x, ...) {
  cat("Attacker Risk Models\n")
  cat("====================\n\n")

  cat("Key variables:", paste(x$key_vars, collapse = ", "), "\n")
  cat("Records:", x$n_records, "\n")
  cat("Model:", x$model, "\n\n")

  cat("Global Risk:\n")
  if (!is.null(x$global_risk$prosecutor)) {
    cat("  Prosecutor:", sprintf("%.4f", x$global_risk$prosecutor), "\n")
  }
  if (!is.null(x$global_risk$journalist)) {
    cat("  Journalist:", sprintf("%.4f", x$global_risk$journalist),
        sprintf("(sampling fraction = %.4f)", x$sampling_fraction), "\n")
  }
  if (!is.null(x$global_risk$marketer)) {
    cat("  Marketer:  ", sprintf("%.4f", x$global_risk$marketer), "\n")
  }

  cat("\nUnique Records:", x$n_uniques, "of", length(x$freq_table),
      "equivalence classes",
      sprintf("(%.1f%%)", 100 * x$pct_uniques), "\n")

  cat("\nPrivacy Assessment:")
  if (x$privacy_pass) {
    cat(" PASS\n")
    cat("  Global re-identification risk is within acceptable bounds (<= 0.1).\n")
  } else {
    cat(" WARNING\n")
    cat("  Elevated re-identification risk detected (> 0.1).\n")
    cat("  Consider additional anonymization measures.\n")
  }

  invisible(x)
}


#' Summary method for attacker_risk objects
#' @param object an object of class "attacker_risk"
#' @param ... additional arguments (ignored)
#' @export
summary.attacker_risk <- function(object, ...) {
  # Frequency distribution
  ec_sizes <- as.numeric(object$freq_table)
  size_dist <- table(factor(pmin(ec_sizes, 10), levels = 1:10))
  names(size_dist) <- c(as.character(1:9), "10+")

  # Risk distribution for prosecutor model
  if (!is.null(object$risk_per_record$prosecutor)) {
    risk_vec <- object$risk_per_record$prosecutor
    risk_quantiles <- quantile(risk_vec, probs = c(0, 0.25, 0.5, 0.75, 0.9, 1))
    n_high_risk <- sum(risk_vec > 0.1)
    pct_high_risk <- 100 * n_high_risk / object$n_records
  } else if (!is.null(object$risk_per_record$journalist)) {
    risk_vec <- object$risk_per_record$journalist
    risk_quantiles <- quantile(risk_vec, probs = c(0, 0.25, 0.5, 0.75, 0.9, 1))
    n_high_risk <- sum(risk_vec > 0.1)
    pct_high_risk <- 100 * n_high_risk / object$n_records
  } else {
    risk_vec <- NULL
    risk_quantiles <- NULL
    n_high_risk <- NA
    pct_high_risk <- NA
  }

  summ <- list(
    global_risk = object$global_risk,
    model = object$model,
    key_vars = object$key_vars,
    sampling_fraction = object$sampling_fraction,
    n_records = object$n_records,
    n_ec = length(ec_sizes),
    n_uniques = object$n_uniques,
    pct_uniques = object$pct_uniques,
    ec_size_distribution = size_dist,
    mean_ec_size = mean(ec_sizes),
    median_ec_size = median(ec_sizes),
    risk_quantiles = risk_quantiles,
    n_high_risk = n_high_risk,
    pct_high_risk = pct_high_risk,
    privacy_pass = object$privacy_pass
  )

  class(summ) <- "summary.attacker_risk"
  return(summ)
}


#' Print method for summary.attacker_risk objects
#' @param x an object of class "summary.attacker_risk"
#' @param ... additional arguments (ignored)
#' @export
print.summary.attacker_risk <- function(x, ...) {
  cat("Summary: Attacker Risk Models\n")
  cat("=============================\n\n")

  cat("Key variables:", paste(x$key_vars, collapse = ", "), "\n")
  cat("Model:", x$model, "\n\n")

  cat("Dataset Statistics:\n")
  cat("  Total records:", x$n_records, "\n")
  cat("  Equivalence classes:", x$n_ec, "\n")
  cat("  Unique records (f_k = 1):", x$n_uniques,
      sprintf("(%.1f%% of ECs)", 100 * x$pct_uniques), "\n")
  cat("  Mean EC size:", round(x$mean_ec_size, 1), "\n")
  cat("  Median EC size:", x$median_ec_size, "\n\n")

  cat("Global Risk:\n")
  if (!is.null(x$global_risk$prosecutor)) {
    cat("  Prosecutor:", sprintf("%.4f", x$global_risk$prosecutor), "\n")
  }
  if (!is.null(x$global_risk$journalist)) {
    cat("  Journalist:", sprintf("%.4f", x$global_risk$journalist),
        sprintf("(sampling fraction = %.4f)", x$sampling_fraction), "\n")
  }
  if (!is.null(x$global_risk$marketer)) {
    cat("  Marketer:  ", sprintf("%.4f", x$global_risk$marketer), "\n")
  }

  if (!is.null(x$risk_quantiles)) {
    cat("\nPer-Record Risk Distribution:\n")
    cat("  Min:", sprintf("%.4f", x$risk_quantiles[1]), "\n")
    cat("  Q1: ", sprintf("%.4f", x$risk_quantiles[2]), "\n")
    cat("  Median:", sprintf("%.4f", x$risk_quantiles[3]), "\n")
    cat("  Q3: ", sprintf("%.4f", x$risk_quantiles[4]), "\n")
    cat("  P90:", sprintf("%.4f", x$risk_quantiles[5]), "\n")
    cat("  Max:", sprintf("%.4f", x$risk_quantiles[6]), "\n")
    cat("  Records with risk > 0.1:", x$n_high_risk,
        sprintf("(%.1f%%)", x$pct_high_risk), "\n")
  }

  cat("\nEquivalence Class Size Distribution:\n")
  print(x$ec_size_distribution)

  cat("\nPrivacy:", ifelse(x$privacy_pass, "PASS", "WARNING"), "\n")

  invisible(x)
}


#' Plot method for attacker_risk objects
#' @param x an object of class "attacker_risk"
#' @param y not used
#' @param ... additional arguments passed to plotting functions
#' @param which integer, which plot(s) to show:
#'   1 = Per-record risk distribution (histogram),
#'   2 = Comparison across models (barplot of global risks)
#' @importFrom graphics barplot hist abline par legend text
#' @export
plot.attacker_risk <- function(x, y = NULL, ..., which = 1) {
  show <- rep(FALSE, 2)
  show[which] <- TRUE

  n_plots <- sum(show)
  if (n_plots > 1) {
    op <- par(mfrow = c(1, 2))
    on.exit(par(op))
  }

  if (show[1]) {
    # Per-record risk distribution (histogram)
    # Use prosecutor risk if available, else journalist
    if (!is.null(x$risk_per_record$prosecutor)) {
      risk_vec <- x$risk_per_record$prosecutor
      risk_label <- "Prosecutor Risk"
    } else if (!is.null(x$risk_per_record$journalist)) {
      risk_vec <- x$risk_per_record$journalist
      risk_label <- "Journalist Risk"
    } else {
      # Marketer only: recompute per-record prosecutor for histogram
      risk_vec <- 1 / x$risk_per_record$freq
      risk_label <- "Per-Record Risk (1/f_k)"
    }

    hist(risk_vec,
         main = paste("Per-Record Risk Distribution\n(", risk_label, ")"),
         xlab = "Re-identification Risk",
         ylab = "Number of Records",
         col = "steelblue", border = "white", ...)
    abline(v = 0.1, col = "red", lty = 2, lwd = 2)
    legend("topright",
           legend = "threshold = 0.1",
           col = "red", lty = 2, lwd = 2, cex = 0.8)
  }

  if (show[2]) {
    # Global risk comparison across models
    vals <- numeric(0)
    names_vals <- character(0)

    if (!is.null(x$global_risk$prosecutor)) {
      vals <- c(vals, x$global_risk$prosecutor)
      names_vals <- c(names_vals, "Prosecutor")
    }
    if (!is.null(x$global_risk$journalist)) {
      vals <- c(vals, x$global_risk$journalist)
      names_vals <- c(names_vals, "Journalist")
    }
    if (!is.null(x$global_risk$marketer)) {
      vals <- c(vals, x$global_risk$marketer)
      names_vals <- c(names_vals, "Marketer")
    }

    cols <- c("coral", "steelblue", "seagreen")[seq_along(vals)]

    bp <- barplot(vals, names.arg = names_vals,
                  main = "Global Risk by Attacker Model",
                  ylab = "Global Risk",
                  col = cols,
                  ylim = c(0, max(0.2, max(vals) * 1.3)), ...)
    abline(h = 0.1, col = "red", lwd = 2, lty = 2)
    text(bp, vals + max(vals) * 0.05,
         labels = sprintf("%.4f", vals), cex = 0.9)
    legend("topright", "threshold = 0.1", col = "red", lty = 2, lwd = 2,
           cex = 0.8)
  }
}
