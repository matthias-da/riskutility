#' Population Uniqueness Risk
#'
#' Estimates the probability that records that are unique in the sample are also
#' unique in the population, using super-population models. Three estimation
#' methods are available: Pitman, Zayatz, and Sample-Negative-Binomial (SNB).
#'
#' @param X data frame to assess, or a \code{\link{synth_pair}} object
#' @param key_vars character vector of quasi-identifier variable names
#' @param sampling_fraction numeric, the fraction of the population represented
#'   by the sample (default: 0.01)
#' @param method character, one of "pitman", "zayatz", "snb", or "all" (default:
#'   "all"). If "all", all three methods are computed and compared.
#' @param na.rm logical, remove records with NA in key variables (default: TRUE)
#' @param data character, which dataset to evaluate when using a
#'   \code{\link{synth_pair}} object: "synthetic" (default) or "original".
#'   Ignored for the default method.
#' @param ... additional arguments passed to methods (currently unused)
#'
#' @return An object of class "population_uniqueness" containing:
#' \itemize{
#'   \item risk_per_record: numeric vector of P(F_k=1|f_k=1) per record
#'   \item global_risk: mean risk across all records
#'   \item n_sample_uniques: count of records with f_k=1
#'   \item pct_sample_uniques: percentage of sample uniques
#'   \item n_population_uniques_est: estimated number of population uniques
#'   \item freq_table: table of QI combination frequencies
#'   \item method: method used
#'   \item key_vars: key variables used
#'   \item sampling_fraction: sampling fraction used
#'   \item n_records: number of records assessed
#'   \item privacy_pass: logical, global_risk <= 0.1
#'   \item comparison: if method="all", list of results per method
#' }
#'
#' @details
#' Population uniqueness estimation asks: given that a record has a unique
#' combination of quasi-identifiers in the sample (f_k = 1), what is the
#' probability that this combination is also unique in the entire population
#' (F_k = 1)?
#'
#' Three super-population models are supported:
#'
#' \strong{Zayatz (1991):}
#' The simplest model. For each sample unique (f_k = 1):
#' \deqn{P(F_k = 1 | f_k = 1) = (1 - \pi)^{1/\pi - 1}}
#' where \eqn{\pi} is the sampling fraction. This gives a single constant
#' probability for all sample uniques. Records with f_k > 1 have risk 0.
#'
#' \strong{Pitman (2003):}
#' Uses a Pitman partition model with parameter \eqn{\alpha} estimated from
#' the frequency distribution. For sample uniques:
#' \deqn{P(F_k = 1 | f_k = 1) \approx \frac{\alpha}{\alpha + n \cdot \pi}}
#' where \eqn{n} is the sample size. The parameter \eqn{\alpha} is estimated
#' from the ratio of frequency classes in the sample.
#'
#' \strong{SNB (Bethlehem et al., 1990):}
#' Sample-Negative-Binomial model. Fits a negative binomial distribution to
#' the population frequency distribution using method of moments. Parameters
#' \eqn{r} and \eqn{q} are estimated from the sample frequency distribution.
#' For sample uniques:
#' \deqn{P(F_k = 1 | f_k = 1) = \frac{(1 - \pi q / (1 + q))^r}{\text{norm}}}
#'
#' \strong{Interpretation:}
#' \itemize{
#'   \item \strong{Global risk < 0.05}: Low population uniqueness risk
#'   \item \strong{Global risk 0.05--0.10}: Moderate risk, caution advised
#'   \item \strong{Global risk > 0.10}: Elevated risk, records may be identifiable
#' }
#'
#' @seealso \code{\link{kanonymity}} for k-anonymity assessment,
#'   \code{\link{suda}} for special uniques detection,
#'   \code{\link{individual_risk}} for individual re-identification risk
#'
#' @references
#' Zayatz, L. (1991). Estimation of the Number of Unique Population Elements
#' Using a Sample. \emph{Proceedings of the Section on Survey Research Methods,
#' American Statistical Association}, 369--373.
#'
#' Pitman, J. (2003). Poisson-Kingman Partitions. In \emph{Statistics and
#' Science: A Festschrift for Terry Speed}, IMS Lecture Notes Monograph Series,
#' 40, 1--34.
#'
#' Bethlehem, J.G., Keller, W.J. & Pannekoek, J. (1990).
#' Disclosure Control of Microdata.
#' \emph{Journal of the American Statistical Association}, 85(409), 38--45.
#' \doi{10.1080/01621459.1990.10475304}
#'
#' Skinner, C.J. & Shlomo, N. (2008).
#' Assessing Identification Risk in Survey Microdata Using Log-Linear Models.
#' \emph{Journal of the American Statistical Association}, 103(483), 989--1001.
#' \doi{10.1198/016214507000001328}
#'
#' Templ, M. (2017). \emph{Statistical Disclosure Control for Microdata:
#' Methods and Applications in R}. Springer.
#' \doi{10.1007/978-3-319-50272-4}
#'
#' @author Matthias Templ
#' @family privacy-models
#' @export
#' @importFrom stats complete.cases dnbinom var median
#' @importFrom graphics hist abline legend barplot par text
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
#' # Estimate population uniqueness with all methods
#' result <- population_uniqueness(data,
#'                                  key_vars = c("age", "gender", "region"),
#'                                  sampling_fraction = 0.01)
#' print(result)
#' summary(result)
#'
#' # Use a specific method
#' result_z <- population_uniqueness(data,
#'                                    key_vars = c("age", "gender", "region"),
#'                                    sampling_fraction = 0.05,
#'                                    method = "zayatz")
#' print(result_z)
#'
#' \donttest{
#' # More quasi-identifiers = higher risk of uniqueness
#' result_more <- population_uniqueness(data,
#'                                       key_vars = c("age", "gender", "region", "income"),
#'                                       sampling_fraction = 0.01)
#' print(result_more)
#' plot(result_more)
#' }
population_uniqueness <- function(X, ...) {
  UseMethod("population_uniqueness")
}

#' @rdname population_uniqueness
#' @export
population_uniqueness.synth_pair <- function(X, data = c("synthetic", "original"), ...) {
  if (is.null(X$key_vars)) {
    stop("synth_pair must have 'key_vars' set for population_uniqueness()")
  }

  data <- match.arg(data)
  dataset <- if (data == "synthetic") X$synthetic else X$original

  population_uniqueness.default(
    X = dataset,
    key_vars = X$key_vars,
    ...
  )
}

#' @rdname population_uniqueness
#' @export
population_uniqueness.default <- function(X,
                                           key_vars,
                                           sampling_fraction = 0.01,
                                           method = c("all", "pitman", "zayatz", "snb"),
                                           na.rm = TRUE,
                                           ...) {

  # Input validation
  if (!is.data.frame(X)) stop("X must be a data frame.")

  method <- match.arg(method)

  if (!is.numeric(sampling_fraction) || length(sampling_fraction) != 1 ||
      sampling_fraction <= 0 || sampling_fraction >= 1) {
    stop("sampling_fraction must be a single numeric value between 0 and 1 (exclusive).")
  }

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
  }

  if (nrow(X_keys) == 0) stop("No complete cases remaining.")

  n <- nrow(X_keys)
  pi_frac <- sampling_fraction

  # Compute frequency of each QI combination
  key_signature <- apply(X_keys, 1, paste, collapse = "|")
  freq_tab <- table(key_signature)
  freq_vec <- as.integer(freq_tab)

  # Map each record to its frequency
  record_freq <- as.integer(freq_tab[key_signature])

  # Number of sample uniques
  n_sample_uniques <- sum(record_freq == 1L)
  pct_sample_uniques <- 100 * n_sample_uniques / n

  # Run the requested method(s)
  if (method == "all") {
    res_pitman <- .pu_pitman(record_freq, freq_vec, n, pi_frac)
    res_zayatz <- .pu_zayatz(record_freq, n, pi_frac)
    res_snb <- .pu_snb(record_freq, freq_vec, n, pi_frac)

    # Use Pitman as the primary result (most commonly used)
    risk_per_record <- res_pitman$risk_per_record
    global_risk <- res_pitman$global_risk
    n_population_uniques_est <- res_pitman$n_population_uniques_est

    comparison <- list(
      pitman = res_pitman,
      zayatz = res_zayatz,
      snb = res_snb
    )
  } else if (method == "pitman") {
    res <- .pu_pitman(record_freq, freq_vec, n, pi_frac)
    risk_per_record <- res$risk_per_record
    global_risk <- res$global_risk
    n_population_uniques_est <- res$n_population_uniques_est
    comparison <- NULL
  } else if (method == "zayatz") {
    res <- .pu_zayatz(record_freq, n, pi_frac)
    risk_per_record <- res$risk_per_record
    global_risk <- res$global_risk
    n_population_uniques_est <- res$n_population_uniques_est
    comparison <- NULL
  } else if (method == "snb") {
    res <- .pu_snb(record_freq, freq_vec, n, pi_frac)
    risk_per_record <- res$risk_per_record
    global_risk <- res$global_risk
    n_population_uniques_est <- res$n_population_uniques_est
    comparison <- NULL
  }

  privacy_pass <- global_risk <= 0.1

  results <- list(
    risk_per_record = risk_per_record,
    global_risk = global_risk,
    n_sample_uniques = n_sample_uniques,
    pct_sample_uniques = pct_sample_uniques,
    n_population_uniques_est = n_population_uniques_est,
    freq_table = freq_tab,
    method = method,
    key_vars = key_vars,
    sampling_fraction = pi_frac,
    n_records = n,
    privacy_pass = privacy_pass,
    comparison = comparison
  )

  class(results) <- "population_uniqueness"
  return(results)
}


# --- Internal estimation methods ---

#' Pitman model for population uniqueness
#'
#' @param record_freq integer vector, frequency of each record's QI combination
#' @param freq_vec integer vector, frequencies of all unique QI combinations
#' @param n integer, sample size
#' @param pi_frac numeric, sampling fraction
#' @return list with risk_per_record, global_risk, n_population_uniques_est, alpha
#' @keywords internal
.pu_pitman <- function(record_freq, freq_vec, n, pi_frac) {
  # Estimate alpha from frequency distribution
  # Use the ratio of sample uniques (f1) to pairs (f2)
  f1 <- sum(freq_vec == 1L)
  f2 <- sum(freq_vec == 2L)

  if (f2 > 0) {
    # Method of moments estimator: alpha ~ f1^2 / (2 * f2)
    alpha <- f1^2 / (2 * f2)
  } else if (f1 > 0) {
    # If no pairs, use conservative estimate
    alpha <- f1
  } else {
    # No uniques at all
    alpha <- 1
  }

  # Ensure alpha is positive
  alpha <- max(alpha, 0.01)

  # Risk for sample uniques: P(F_k=1|f_k=1) ~ alpha / (alpha + n * pi)
  p_unique <- alpha / (alpha + n * pi_frac)

  # Assign risk per record
  risk_per_record <- ifelse(record_freq == 1L, p_unique, 0)

  # Global risk: mean across all records
  global_risk <- mean(risk_per_record)

  # Estimated number of population uniques
  n_population_uniques_est <- f1 * p_unique

  list(
    risk_per_record = risk_per_record,
    global_risk = global_risk,
    n_population_uniques_est = n_population_uniques_est,
    alpha = alpha,
    method = "pitman"
  )
}


#' Zayatz model for population uniqueness
#'
#' @param record_freq integer vector, frequency of each record's QI combination
#' @param n integer, sample size
#' @param pi_frac numeric, sampling fraction
#' @return list with risk_per_record, global_risk, n_population_uniques_est
#' @keywords internal
.pu_zayatz <- function(record_freq, n, pi_frac) {
  # P(F_k=1|f_k=1) = (1 - pi)^(1/pi - 1)
  p_unique <- (1 - pi_frac)^(1 / pi_frac - 1)

  # Assign risk per record
  risk_per_record <- ifelse(record_freq == 1L, p_unique, 0)

  # Global risk: mean across all records
  global_risk <- mean(risk_per_record)

  # Estimated number of population uniques
  f1 <- sum(record_freq == 1L)
  n_population_uniques_est <- f1 * p_unique

  list(
    risk_per_record = risk_per_record,
    global_risk = global_risk,
    n_population_uniques_est = n_population_uniques_est,
    p_unique = p_unique,
    method = "zayatz"
  )
}


#' SNB (Sample-Negative-Binomial) model for population uniqueness
#'
#' @param record_freq integer vector, frequency of each record's QI combination
#' @param freq_vec integer vector, frequencies of all unique QI combinations
#' @param n integer, sample size
#' @param pi_frac numeric, sampling fraction
#' @return list with risk_per_record, global_risk, n_population_uniques_est, r, q
#' @keywords internal
.pu_snb <- function(record_freq, freq_vec, n, pi_frac) {
  # Estimate negative binomial parameters from frequency distribution
  # using method of moments
  K <- length(freq_vec)  # number of unique QI combinations
  mean_f <- mean(freq_vec)
  var_f <- var(freq_vec)

  # Method of moments for negative binomial:
  # mean = r * q / (1 - q) ... but we need to relate sample to population
  # Using the approach from Bethlehem et al.:
  # Population frequency F_k ~ NB(r, q)
  # Sample frequency f_k ~ Binomial(F_k, pi) integrated over F_k

  # The marginal distribution of f_k is also negative binomial with
  # modified parameters:
  #   E[f_k] = r * pi * q / (1 - q)
  #   Var[f_k] = r * pi * q * (1 + pi * q) / (1 - q)^2

  # From mean and variance of sample frequencies, estimate r and q:
  # E[f_k|f_k>0] is the mean of observed frequencies
  # But we observe frequencies of unique keys, not all cells.

  # Simpler estimation via the ratio:
  #   var_f / mean_f = (1 + pi * q) / (1 - q)
  # Let R = var_f / mean_f
  # Then: R * (1 - q) = 1 + pi * q
  #       R - R*q = 1 + pi*q
  #       R - 1 = q * (R + pi)
  #       q = (R - 1) / (R + pi)

  if (var_f > mean_f) {
    R <- var_f / mean_f
    q_est <- (R - 1) / (R + pi_frac)
  } else {
    # Variance <= mean: no overdispersion, use small q
    q_est <- 0.01
  }

  # Bound q away from 0 and 1
  q_est <- max(0.001, min(q_est, 0.999))

  # Estimate r from mean:
  # mean_f = r * pi * q / (1 - q)
  # r = mean_f * (1 - q) / (pi * q)
  r_est <- mean_f * (1 - q_est) / (pi_frac * q_est)

  # Ensure r is positive
  r_est <- max(r_est, 0.01)

  # P(F_k = 1 | f_k = 1) for the SNB model
  # Using the posterior probability:
  # P(F_k = 1 | f_k = 1) = P(f_k = 1 | F_k = 1) * P(F_k = 1) / P(f_k = 1)
  #
  # P(F_k = 1) = r * q * (1 - q)^(r-1) / (1 - (1-q)^r)  [truncated at F>0]
  #            simplified for untruncated: r * q / (1 + q)^(r+1) ... NB pmf
  #
  # Using the standard NB pmf: P(F_k = j) = choose(j+r-1, j) * (q/(1+q))^j * (1/(1+q))^r
  #
  # P(f_k = 1 | F_k) = F_k * pi * (1 - pi)^(F_k - 1)
  #
  # The posterior for sample uniques:
  # P(F_k = 1 | f_k = 1) = [pi * P(F_k=1)] / sum_{j>=1}[j*pi*(1-pi)^(j-1) * P(F_k=j)]
  #                       = P(F_k=1) / sum_{j>=1}[j*(1-pi)^(j-1) * P(F_k=j)]

  # NB pmf (parametrized as size=r, prob = 1/(1+q)):
  # P(F_k = j) = choose(j + r - 1, j) * (q/(1+q))^j * (1/(1+q))^r
  # = dnbinom(j, size = r, prob = 1/(1+q))

  prob_nb <- 1 / (1 + q_est)

  # Numerator: P(F_k = 1) in NB
  p_Fk_1 <- dnbinom(1, size = r_est, prob = prob_nb)

  # Denominator: E[F_k * (1-pi)^(F_k-1) | F_k >= 1] * P(F_k >= 1)
  # We compute the sum over a finite range
  max_j <- max(100, 5 * ceiling(r_est * q_est))
  max_j <- min(max_j, 10000)  # cap for safety
  j_vals <- seq_len(max_j)

  p_Fk_j <- dnbinom(j_vals, size = r_est, prob = prob_nb)
  weights_j <- j_vals * (1 - pi_frac)^(j_vals - 1)
  denom <- sum(weights_j * p_Fk_j)

  if (denom > 0) {
    p_unique <- p_Fk_1 / denom
  } else {
    p_unique <- 0
  }

  # Bound to [0, 1]
  p_unique <- max(0, min(p_unique, 1))

  # Assign risk per record
  risk_per_record <- ifelse(record_freq == 1L, p_unique, 0)

  # Global risk: mean across all records
  global_risk <- mean(risk_per_record)

  # Estimated number of population uniques
  f1 <- sum(record_freq == 1L)
  n_population_uniques_est <- f1 * p_unique

  list(
    risk_per_record = risk_per_record,
    global_risk = global_risk,
    n_population_uniques_est = n_population_uniques_est,
    r = r_est,
    q = q_est,
    p_unique = p_unique,
    method = "snb"
  )
}


# --- S3 methods ---

#' Print method for population_uniqueness objects
#'
#' @param x an object of class "population_uniqueness"
#' @param ... additional arguments (ignored)
#' @export
print.population_uniqueness <- function(x, ...) {
  cat("Population Uniqueness Risk Assessment\n")
  cat("======================================\n")
  cat("Method:", x$method, "\n")
  cat("Key variables:", paste(x$key_vars, collapse = ", "), "\n")
  cat("Sampling fraction:", x$sampling_fraction, "\n")
  cat("Records:", x$n_records, "\n\n")

  cat("Sample Uniques:\n")
  cat("  Count:", x$n_sample_uniques,
      sprintf("(%.1f%%)", x$pct_sample_uniques), "\n\n")

  if (!is.null(x$comparison)) {
    cat("Method Comparison:\n")
    for (nm in names(x$comparison)) {
      res <- x$comparison[[nm]]
      cat(sprintf("  %-8s  P(pop unique | sample unique) = %.4f  |  global risk = %.4f  |  est. pop uniques = %.1f\n",
                  nm, max(res$risk_per_record), res$global_risk,
                  res$n_population_uniques_est))
    }
    cat("\n")
  } else {
    cat("Risk Estimate:\n")
    p_val <- max(x$risk_per_record)
    cat("  P(population unique | sample unique):", sprintf("%.4f", p_val), "\n")
    cat("  Global risk (mean across all records):", sprintf("%.4f", x$global_risk), "\n")
    cat("  Estimated population uniques:", sprintf("%.1f", x$n_population_uniques_est), "\n\n")
  }

  cat("Privacy Assessment:")
  if (x$privacy_pass) {
    cat(" PASS\n")
    cat("  Population uniqueness risk is within acceptable bounds (<= 0.1).\n")
  } else {
    cat(" WARNING\n")
    cat("  Elevated population uniqueness risk (> 0.1).\n")
    cat("  Sample uniques may also be unique in the population.\n")
  }

  invisible(x)
}


#' Summary method for population_uniqueness objects
#'
#' @param object an object of class "population_uniqueness"
#' @param ... additional arguments (ignored)
#' @export
summary.population_uniqueness <- function(object, ...) {
  # Frequency distribution summary
  freq_vec <- as.integer(object$freq_table)
  freq_dist <- table(factor(pmin(freq_vec, 10), levels = 1:10))
  names(freq_dist) <- c(as.character(1:9), "10+")

  summ <- list(
    global_risk = object$global_risk,
    n_sample_uniques = object$n_sample_uniques,
    pct_sample_uniques = object$pct_sample_uniques,
    n_population_uniques_est = object$n_population_uniques_est,
    n_records = object$n_records,
    n_keys = length(freq_vec),
    freq_distribution = freq_dist,
    mean_freq = mean(freq_vec),
    median_freq = median(freq_vec),
    method = object$method,
    key_vars = object$key_vars,
    sampling_fraction = object$sampling_fraction,
    privacy_pass = object$privacy_pass,
    comparison = object$comparison
  )

  class(summ) <- "summary.population_uniqueness"
  return(summ)
}


#' Print method for summary.population_uniqueness objects
#'
#' @param x an object of class "summary.population_uniqueness"
#' @param ... additional arguments (ignored)
#' @export
print.summary.population_uniqueness <- function(x, ...) {
  cat("Summary: Population Uniqueness Risk Assessment\n")
  cat("===============================================\n")
  cat("Method:", x$method, "\n")
  cat("Key variables:", paste(x$key_vars, collapse = ", "), "\n")
  cat("Sampling fraction:", x$sampling_fraction, "\n\n")

  cat("Dataset Statistics:\n")
  cat("  Total records:", x$n_records, "\n")
  cat("  Unique QI combinations:", x$n_keys, "\n")
  cat("  Mean frequency:", round(x$mean_freq, 2), "\n")
  cat("  Median frequency:", x$median_freq, "\n\n")

  cat("Risk Assessment:\n")
  cat("  Sample uniques:", x$n_sample_uniques,
      sprintf("(%.1f%%)", x$pct_sample_uniques), "\n")
  cat("  Est. population uniques:", sprintf("%.1f", x$n_population_uniques_est), "\n")
  cat("  Global risk:", sprintf("%.4f", x$global_risk), "\n")
  cat("  Privacy:", ifelse(x$privacy_pass, "PASS", "WARNING"), "\n\n")

  if (!is.null(x$comparison)) {
    cat("Method Comparison:\n")
    comp_df <- data.frame(
      method = names(x$comparison),
      p_unique = vapply(x$comparison, function(r) max(r$risk_per_record), numeric(1)),
      global_risk = vapply(x$comparison, function(r) r$global_risk, numeric(1)),
      est_pop_uniques = vapply(x$comparison, function(r) r$n_population_uniques_est, numeric(1))
    )
    comp_df$p_unique <- sprintf("%.4f", comp_df$p_unique)
    comp_df$global_risk <- sprintf("%.4f", comp_df$global_risk)
    comp_df$est_pop_uniques <- sprintf("%.1f", comp_df$est_pop_uniques)
    print(comp_df, row.names = FALSE)
    cat("\n")
  }

  cat("Frequency Distribution (QI combinations):\n")
  print(x$freq_distribution)

  invisible(x)
}


#' Plot method for population_uniqueness objects
#'
#' @param x an object of class "population_uniqueness"
#' @param y not used
#' @param ... additional arguments passed to plotting functions
#' @param which integer, which plot(s) to show: 1 = per-record risk distribution
#'   histogram, 2 = method comparison (only when method="all")
#' @export
plot.population_uniqueness <- function(x, y = NULL, ..., which = 1) {
  show <- rep(FALSE, 2)
  show[which] <- TRUE

  n_plots <- sum(show)
  if (n_plots > 1) {
    op <- par(mfrow = c(1, 2))
    on.exit(par(op))
  }

  if (show[1]) {
    # Per-record risk distribution histogram
    risk_vals <- x$risk_per_record
    n_at_risk <- sum(risk_vals > 0)
    n_safe <- sum(risk_vals == 0)

    if (n_at_risk > 0) {
      hist(risk_vals[risk_vals > 0],
           main = paste0("Population Uniqueness Risk\n(", x$method, ")"),
           xlab = "P(population unique | sample unique)",
           ylab = "Number of Records",
           col = "coral", border = "white", ...)
      abline(v = 0.1, col = "red", lty = 2, lwd = 2)
      legend("topright",
             legend = c(paste("At risk:", n_at_risk),
                        paste("Safe (f_k > 1):", n_safe),
                        "Threshold = 0.1"),
             col = c("coral", NA, "red"),
             pch = c(15, NA, NA),
             lty = c(NA, NA, 2),
             lwd = c(NA, NA, 2),
             cex = 0.8)
    } else {
      barplot(c("At risk" = 0, "Safe" = n_safe),
              main = paste0("Population Uniqueness Risk\n(", x$method, ")"),
              ylab = "Number of Records",
              col = c("coral", "steelblue"), ...)
      legend("topright", "No sample uniques found", cex = 0.9)
    }
  }

  if (show[2]) {
    if (is.null(x$comparison)) {
      # Show frequency distribution instead
      freq_vec <- as.integer(x$freq_table)
      freq_dist <- table(factor(pmin(freq_vec, 10), levels = 1:10))
      names(freq_dist) <- c(as.character(1:9), "10+")
      cols <- ifelse(as.numeric(names(freq_dist)[1:9]) == 1, "coral", "steelblue")
      cols <- c(cols, "steelblue")
      barplot(freq_dist,
              main = "QI Frequency Distribution",
              xlab = "Frequency of QI Combination",
              ylab = "Number of Combinations",
              col = cols, ...)
      legend("topright",
             legend = c("Uniques (f=1)", "Non-unique"),
             fill = c("coral", "steelblue"), cex = 0.8)
    } else {
      # Method comparison barplot
      comp_risks <- vapply(x$comparison, function(r) max(r$risk_per_record), numeric(1))
      names(comp_risks) <- c("Pitman", "Zayatz", "SNB")
      cols <- c("steelblue", "coral", "darkgreen")

      bp <- barplot(comp_risks, names.arg = names(comp_risks),
                    main = "P(pop unique | sample unique)\nMethod Comparison",
                    ylab = "Probability",
                    col = cols,
                    ylim = c(0, max(0.2, max(comp_risks) * 1.3)), ...)
      abline(h = 0.1, col = "red", lty = 2, lwd = 2)
      text(bp, comp_risks + max(comp_risks) * 0.05,
           labels = round(comp_risks, 4), cex = 0.9)
      legend("topright", "threshold = 0.1", col = "red", lty = 2, lwd = 2,
             cex = 0.8)
    }
  }
}
