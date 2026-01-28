#' Individual Re-identification Risk
#'
#' Estimates individual re-identification risk based on population frequency models.
#' Computes sample-based and model-based risk estimates using log-linear or Poisson
#' regression models for population frequency estimation.
#'
#' @param X data frame to assess, or a \code{\link{synth_pair}} object
#' @param key_vars character vector of quasi-identifier variable names
#' @param weight optional, name of weight variable or numeric vector of sampling weights
#' @param method character, risk estimation method: "sample" (sample frequencies only),
#'   "model" (log-linear model), or "both" (default)
#' @param threshold numeric, risk threshold for flagging high-risk records (default: 0.1)
#' @param na.rm logical, remove records with NA in key variables (default: TRUE)
#' @param ... additional arguments passed to methods
#'
#' @return An object of class "individual_risk" containing:
#' \itemize{
#'   \item risk: vector of individual risk values
#'   \item method: method used for risk estimation
#'   \item mean_risk: mean risk across all records
#'   \item max_risk: maximum individual risk
#'   \item n_high_risk: number of records exceeding threshold
#'   \item pct_high_risk: percentage of high-risk records
#'   \item threshold: risk threshold used
#'   \item risk_distribution: summary statistics of risk distribution
#'   \item ec_info: equivalence class information
#' }
#'
#' @details
#' Individual risk quantifies the probability that a specific record can be
#' re-identified by an intruder who has access to the quasi-identifiers.
#'
#' \strong{Sample-based Risk:}
#' For a record in an equivalence class of size f_k, the sample-based risk is:
#' \deqn{r_k = \frac{1}{f_k}}
#' This is a simple estimate assuming uniform probability within each class.
#'
#' \strong{Model-based Risk:}
#' When sampling weights are available or for small samples, model-based approaches
#' estimate the population frequency F_k using log-linear models:
#' \deqn{r_k = \frac{1}{F_k}}
#'
#' The function fits a Poisson log-linear model to estimate expected frequencies
#' in the population, accounting for the sampling design.
#'
#' \strong{Risk Interpretation:}
#' \itemize{
#'   \item Risk = 1: Unique record, definitely identifiable
#'   \item Risk > 0.5: High risk, appears in small equivalence class
#'   \item Risk > 0.1: Elevated risk, deserves attention
#'   \item Risk < 0.05: Generally acceptable for most applications
#' }
#'
#' @seealso \code{\link{kanonymity}} for k-anonymity assessment,
#'   \code{\link{suda}} for Special Uniques Detection,
#'   \code{\link{ldiversity}} for l-diversity assessment
#'
#' @references
#' Skinner, C. J., & Elliot, M. J. (2002). A measure of disclosure risk for microdata.
#' \emph{Journal of the Royal Statistical Society: Series B}, 64(4), 855-867.
#'
#' Rinott, Y., & Shlomo, N. (2006). A generalized negative binomial smoothing model
#' for sample disclosure risk estimation. In \emph{Privacy in Statistical Databases}.
#'
#' @author Matthias Templ
#' @importFrom stats poisson glm predict as.formula
#' @export
#' @examples
#' # Create example data
#' set.seed(123)
#' data <- data.frame(
#'   age = sample(c("young", "middle", "old"), 500, replace = TRUE),
#'   gender = sample(c("M", "F"), 500, replace = TRUE),
#'   region = sample(c("N", "S", "E", "W"), 500, replace = TRUE),
#'   income = sample(c("low", "medium", "high"), 500, replace = TRUE)
#' )
#'
#' # Compute individual risk
#' result <- individual_risk(data,
#'                           key_vars = c("age", "gender", "region"),
#'                           threshold = 0.2)
#' print(result)
#' summary(result)
#' plot(result)
#'
#' # With sampling weights
#' data$weight <- runif(500, 0.5, 2)
#' result_weighted <- individual_risk(data,
#'                                    key_vars = c("age", "gender", "region"),
#'                                    weight = "weight")
#' print(result_weighted)
individual_risk <- function(X, ...) {
  UseMethod("individual_risk")
}

#' @rdname individual_risk
#' @export
individual_risk.synth_pair <- function(X, threshold = 0.1, ...) {
  if (is.null(X$key_vars)) {
    stop("synth_pair must have 'key_vars' set for individual_risk()")
  }

  individual_risk.default(
    X = X$synthetic,
    key_vars = X$key_vars,
    weight = X$weight_var,
    threshold = threshold,
    ...
  )
}

#' @rdname individual_risk
#' @export
individual_risk.default <- function(X,
                                    key_vars,
                                    weight = NULL,
                                    method = c("both", "sample", "model"),
                                    threshold = 0.1,
                                    na.rm = TRUE,
                                    ...) {

  method <- match.arg(method)

  # Input validation
  if (!is.data.frame(X)) stop("X must be a data frame.")

  # Check key variables exist
  missing_vars <- setdiff(key_vars, names(X))
  if (length(missing_vars) > 0) {
    stop(paste("Key variables missing in X:", paste(missing_vars, collapse = ", ")))
  }

  # Extract weights
  wt <- NULL
  if (!is.null(weight)) {
    if (is.character(weight) && length(weight) == 1) {
      if (!weight %in% names(X)) stop(paste("Weight variable", weight, "not found in X"))
      wt <- X[[weight]]
    } else if (is.numeric(weight)) {
      if (length(weight) != nrow(X)) stop("weight must have same length as nrow(X)")
      wt <- weight
    }
  }

  # Handle missing values
  X_keys <- X[, key_vars, drop = FALSE]
  if (na.rm) {
    complete_idx <- complete.cases(X_keys)
    if (sum(!complete_idx) > 0) {
      message(sprintf("Removing %d records with NA in key variables.", sum(!complete_idx)))
    }
    X_keys <- X_keys[complete_idx, , drop = FALSE]
    if (!is.null(wt)) wt <- wt[complete_idx]
    original_idx <- which(complete_idx)
  } else {
    original_idx <- seq_len(nrow(X))
  }

  if (nrow(X_keys) == 0) stop("No complete cases remaining.")

  n <- nrow(X_keys)

  # Create key signature for each record
  key_signature <- apply(X_keys, 1, paste, collapse = "|")

  # Sample-based frequencies
  ec_table <- table(key_signature)
  ec_sizes <- as.numeric(ec_table[key_signature])

  # Sample-based risk: 1/f_k
  risk_sample <- 1 / ec_sizes

  # Model-based risk (if requested)
  risk_model <- NULL

  if (method %in% c("model", "both")) {
    tryCatch({
      # Prepare data for Poisson model
      # Aggregate to cell counts
      agg_df <- as.data.frame(table(X_keys), stringsAsFactors = TRUE)
      names(agg_df)[ncol(agg_df)] <- "count"

      if (nrow(agg_df) > 1) {
        # Fit Poisson log-linear model with main effects
        formula_str <- paste("count ~", paste(key_vars, collapse = " + "))
        model <- glm(as.formula(formula_str),
                     data = agg_df,
                     family = poisson())

        # Predict expected counts
        fitted_counts <- predict(model, newdata = agg_df, type = "response")
        fitted_counts[fitted_counts < 1] <- 1  # Minimum expected count

        # Map back to individual records
        # Create lookup
        agg_df$fitted <- fitted_counts
        agg_key <- apply(agg_df[, key_vars, drop = FALSE], 1, paste, collapse = "|")
        names(fitted_counts) <- agg_key

        # Model-based risk
        expected_size <- fitted_counts[key_signature]
        risk_model <- 1 / expected_size

        # If weights available, adjust for sampling fraction
        if (!is.null(wt)) {
          # Estimate population size per cell
          weighted_counts <- tapply(wt, key_signature, sum)
          pop_size <- weighted_counts[key_signature]
          risk_model <- 1 / pmax(pop_size, 1)
        }
      } else {
        risk_model <- risk_sample
        message("Only one cell - using sample-based risk.")
      }
    }, error = function(e) {
      warning(paste("Model fitting failed:", e$message, "- using sample-based risk only."))
      risk_model <<- risk_sample
    })
  }

  # Select final risk based on method
  if (method == "sample") {
    risk_final <- risk_sample
    risk_type <- "sample"
  } else if (method == "model") {
    risk_final <- if (is.null(risk_model)) risk_sample else risk_model
    risk_type <- if (is.null(risk_model)) "sample (model failed)" else "model"
  } else {  # both
    risk_final <- if (is.null(risk_model)) risk_sample else pmin(risk_sample, risk_model)
    risk_type <- "combined"
  }

  # Ensure risk is bounded [0, 1]
  risk_final <- pmin(pmax(risk_final, 0), 1)

  # Count high-risk records
  n_high_risk <- sum(risk_final > threshold)
  pct_high_risk <- 100 * n_high_risk / n

  # Risk distribution summary
  risk_dist <- list(
    min = min(risk_final),
    q25 = quantile(risk_final, 0.25),
    median = median(risk_final),
    mean = mean(risk_final),
    q75 = quantile(risk_final, 0.75),
    q90 = quantile(risk_final, 0.90),
    q95 = quantile(risk_final, 0.95),
    max = max(risk_final)
  )

  # Equivalence class info
  ec_df <- data.frame(
    key = names(ec_table),
    size = as.numeric(ec_table),
    n_records = as.numeric(ec_table),
    stringsAsFactors = FALSE
  )
  ec_df$sample_risk <- 1 / ec_df$size
  ec_df <- ec_df[order(ec_df$size), ]
  rownames(ec_df) <- NULL

  # Build results
  results <- list(
    risk = risk_final,
    risk_sample = risk_sample,
    risk_model = risk_model,
    method = risk_type,
    mean_risk = mean(risk_final),
    max_risk = max(risk_final),
    n_high_risk = n_high_risk,
    pct_high_risk = pct_high_risk,
    threshold = threshold,
    n_records = n,
    n_unique = sum(ec_sizes == 1),
    pct_unique = 100 * sum(ec_sizes == 1) / n,
    n_ec = length(unique(key_signature)),
    risk_distribution = risk_dist,
    ec_info = ec_df,
    key_vars = key_vars,
    record_indices = original_idx
  )

  class(results) <- "individual_risk"
  return(results)
}


#' Print method for individual_risk objects
#' @param x an object of class "individual_risk"
#' @param ... additional arguments (ignored)
#' @export
print.individual_risk <- function(x, ...) {
  cat("Individual Re-identification Risk Assessment\n")
  cat("============================================\n\n")

  cat("Key variables:", paste(x$key_vars, collapse = ", "), "\n")
  cat("Records:", x$n_records, "| Equivalence classes:", x$n_ec, "\n")
  cat("Method:", x$method, "\n\n")

  cat("Risk Summary:\n")
  cat("  Mean risk:", sprintf("%.4f", x$mean_risk), "\n")
  cat("  Max risk:", sprintf("%.4f", x$max_risk), "\n")
  cat("  Median risk:", sprintf("%.4f", x$risk_distribution$median), "\n\n")

  cat("High-Risk Records (threshold =", x$threshold, "):\n")
  cat("  Count:", x$n_high_risk, sprintf("(%.1f%%)\n", x$pct_high_risk))
  cat("  Unique records:", x$n_unique, sprintf("(%.1f%%)\n", x$pct_unique))

  cat("\nRisk Interpretation:\n")
  if (x$mean_risk < 0.05) {
    cat("  LOW RISK: Average risk is acceptable for most applications.\n")
  } else if (x$mean_risk < 0.1) {
    cat("  MODERATE RISK: Some records may need additional protection.\n")
  } else if (x$mean_risk < 0.2) {
    cat("  ELEVATED RISK: Consider additional anonymization measures.\n")
  } else {
    cat("  HIGH RISK: Data requires significant privacy protection.\n")
  }

  invisible(x)
}


#' Summary method for individual_risk objects
#' @param object an object of class "individual_risk"
#' @param ... additional arguments (ignored)
#' @export
summary.individual_risk <- function(object, ...) {
  # Risk bands
  risk_bands <- c(
    "Unique (risk=1)" = sum(object$risk == 1),
    "Very high (>0.5)" = sum(object$risk > 0.5 & object$risk < 1),
    "High (0.2-0.5)" = sum(object$risk > 0.2 & object$risk <= 0.5),
    "Moderate (0.1-0.2)" = sum(object$risk > 0.1 & object$risk <= 0.2),
    "Low (0.05-0.1)" = sum(object$risk > 0.05 & object$risk <= 0.1),
    "Very low (<0.05)" = sum(object$risk <= 0.05)
  )

  summ <- list(
    method = object$method,
    n_records = object$n_records,
    n_ec = object$n_ec,
    mean_risk = object$mean_risk,
    max_risk = object$max_risk,
    risk_distribution = object$risk_distribution,
    risk_bands = risk_bands,
    pct_bands = 100 * risk_bands / object$n_records,
    threshold = object$threshold,
    n_high_risk = object$n_high_risk,
    pct_high_risk = object$pct_high_risk,
    key_vars = object$key_vars,
    smallest_ec = head(object$ec_info, 10)
  )

  class(summ) <- "summary.individual_risk"
  return(summ)
}


#' Print method for summary.individual_risk objects
#' @param x an object of class "summary.individual_risk"
#' @param ... additional arguments (ignored)
#' @export
print.summary.individual_risk <- function(x, ...) {
  cat("Summary: Individual Re-identification Risk\n")
  cat("==========================================\n\n")

  cat("Key variables:", paste(x$key_vars, collapse = ", "), "\n")
  cat("Method:", x$method, "\n\n")

  cat("Risk Distribution:\n")
  cat(sprintf("  Min: %.4f | Q25: %.4f | Median: %.4f | Mean: %.4f\n",
              x$risk_distribution$min, x$risk_distribution$q25,
              x$risk_distribution$median, x$risk_distribution$mean))
  cat(sprintf("  Q75: %.4f | Q90: %.4f | Q95: %.4f | Max: %.4f\n",
              x$risk_distribution$q75, x$risk_distribution$q90,
              x$risk_distribution$q95, x$risk_distribution$max))

  cat("\nRisk Bands:\n")
  for (i in seq_along(x$risk_bands)) {
    cat(sprintf("  %s: %d (%.1f%%)\n",
                names(x$risk_bands)[i], x$risk_bands[i], x$pct_bands[i]))
  }

  cat("\nSmallest Equivalence Classes:\n")
  print(x$smallest_ec[, c("key", "size", "sample_risk")], row.names = FALSE)

  invisible(x)
}


#' Plot method for individual_risk objects
#' @param x an object of class "individual_risk"
#' @param y not used
#' @param ... additional arguments passed to plotting functions
#' @param which integer, which plot: 1 = risk distribution histogram,
#'   2 = risk by equivalence class size
#' @importFrom graphics hist barplot abline par legend
#' @export
plot.individual_risk <- function(x, y = NULL, ..., which = 1) {
  show <- rep(FALSE, 2)
  show[which] <- TRUE

  if (sum(show) > 1) {
    op <- par(mfrow = c(1, 2))
    on.exit(par(op))
  }

  if (show[1]) {
    # Histogram of individual risks
    hist(x$risk,
         main = "Distribution of Individual Risk",
         xlab = "Re-identification Risk",
         ylab = "Number of Records",
         col = "steelblue", border = "white",
         breaks = 20, ...)
    abline(v = x$threshold, col = "red", lty = 2, lwd = 2)
    abline(v = x$mean_risk, col = "orange", lty = 3, lwd = 2)
    legend("topright",
           legend = c(paste("Threshold =", x$threshold),
                      paste("Mean =", round(x$mean_risk, 3))),
           col = c("red", "orange"), lty = c(2, 3), lwd = 2, cex = 0.8)
  }

  if (show[2]) {
    # Risk by equivalence class size
    # Aggregate by EC size
    ec_sizes <- x$ec_info$size
    risk_by_size <- tapply(1/ec_sizes, factor(pmin(ec_sizes, 10),
                                               levels = 1:10), mean)
    names(risk_by_size) <- c(as.character(1:9), "10+")

    barplot(risk_by_size,
            main = "Average Risk by Equivalence Class Size",
            xlab = "Equivalence Class Size",
            ylab = "Average Risk",
            col = ifelse(risk_by_size > x$threshold, "firebrick", "steelblue"),
            ...)
    abline(h = x$threshold, col = "red", lty = 2, lwd = 2)
  }
}


#' Get high-risk records
#'
#' Extract records with risk above a specified threshold.
#'
#' @param x an object of class "individual_risk"
#' @param threshold numeric, risk threshold (default: uses threshold from object)
#' @return data frame with high-risk record indices and their risk values
#' @export
high_risk_records <- function(x, threshold = NULL) {
  if (!inherits(x, "individual_risk")) {
    stop("x must be an 'individual_risk' object")
  }

  if (is.null(threshold)) {
    threshold <- x$threshold
  }

  high_idx <- which(x$risk > threshold)

  data.frame(
    index = x$record_indices[high_idx],
    risk = x$risk[high_idx],
    stringsAsFactors = FALSE
  )
}
