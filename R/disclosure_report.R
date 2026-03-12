#' Comprehensive Disclosure Risk Report
#'
#' Computes multiple disclosure risk metrics and produces a comprehensive report
#' for data privacy assessment. Combines attribution-based measures
#' (DCAP, TCAP, WEAP, DiSCO, RAPID), distance-based measures (DCR, NNDR, IMS,
#' dRisk, hitting rate, RF Privacy), privacy models (k-anonymity, l-diversity,
#' t-closeness), and membership inference attacks (singling out, linkability,
#' NNAA, DOMIAS). Works for both synthetic and traditionally anonymized data.
#'
#' @param X data frame of original data
#' @param Y data frame of synthetic/anonymized data
#' @param key_vars character vector of quasi-identifier variable names for
#'   attribution-based and privacy model metrics. If NULL, these metrics are skipped.
#' @param target_var character, name of the sensitive target variable for
#'   attribution-based metrics. If NULL, these metrics are skipped.
#' @param holdout data frame of holdout data for distance-based and membership
#'   inference metrics. If NULL, automatically created from X using holdout_fraction.
#' @param holdout_fraction numeric, fraction of X to use as holdout (default: 0.5)
#' @param distance_vars character vector of variables for distance-based metrics.
#'   If NULL, all common variables are used.
#' @param distance_method character, distance method for DCR/NNDR: "gower" or
#'   "euclidean" (default: "gower")
#' @param compute character vector specifying which metrics to compute.
#'   Default is "all". Options: "all", "attribution", "distance", "privacy",
#'   "membership", or specific metric names (see Details).
#' @param rapid_model character, model type for RAPID: "rf" (default), "lm", "cart"
#' @param na.rm logical, remove records with NA values (default: TRUE)
#' @param seed integer, random seed for holdout sampling (default: NULL)
#' @param verbose logical, print progress messages (default: TRUE)
#' @param ... additional arguments passed to methods
#'
#' @return An object of class "disclosure_report" containing:
#' \itemize{
#'   \item results: list of individual metric results
#'   \item summary: data frame with key metrics and pass/fail status
#'   \item overall_risk: character, overall risk assessment ("LOW", "MEDIUM", "HIGH")
#'   \item n_pass: number of metrics that passed
#'   \item n_fail: number of metrics that failed/warned
#'   \item parameters: list of input parameters used
#' }
#'
#' @details
#' This function provides a one-stop assessment of disclosure risk
#' by computing multiple complementary metrics across four families:
#'
#' \strong{Attribution-Based Metrics} (require key_vars and target_var):
#' \itemize{
#'   \item dcap: Overall correct attribution probability
#'   \item tcap: Per-record attribution risk with categories
#'   \item weap: Within equivalence class attribution (on synthetic only)
#'   \item disco: Count of disclosive synthetic records
#'   \item rapid: ML-based attribute inference risk (requires ranger package)
#' }
#'
#' \strong{Distance-Based Metrics} (use holdout comparison):
#' \itemize{
#'   \item dcr: Distance to closest record ratio
#'   \item nndr: Nearest neighbor distance ratio
#'   \item ims: Identical match share (exact copies)
#'   \item drisk: Disclosure risk for continuous variables
#'   \item hitting_rate: Fraction of records within threshold distance
#'   \item rf_privacy: Random forest proximity-based memorization test (requires ranger)
#' }
#'
#' \strong{Privacy Models} (single-dataset, require key_vars):
#' \itemize{
#'   \item kanonymity: Minimum equivalence class size
#'   \item ldiversity: Sensitive attribute diversity per EC (requires target_var)
#'   \item tcloseness: EMD between EC and overall distribution (requires target_var)
#'   \item individual_risk: Frequency-based per-record risk
#' }
#'
#' \strong{Membership Inference} (require holdout):
#' \itemize{
#'   \item singling_out: Predicate-based uniqueness attack (GDPR criterion)
#'   \item linkability: Record linkage attack (GDPR criterion)
#'   \item nnaa: Nearest-neighbor adversarial accuracy / privacy loss
#'   \item domias: Density-based membership inference (requires ranger package)
#' }
#'
#' The overall risk is determined by the number of failed checks:
#' \itemize{
#'   \item Low: All metrics pass
#'   \item Medium: 1-2 metrics fail
#'   \item High: 3+ metrics fail
#' }
#'
#' @seealso \code{\link{dcap}}, \code{\link{tcap}}, \code{\link{weap}},
#'   \code{\link{disco}}, \code{\link{rapid}}, \code{\link{dcr}}, \code{\link{nndr}},
#'   \code{\link{ims}}, \code{\link{drisk}}, \code{\link{hitting_rate}},
#'   \code{\link{rf_privacy}},
#'   \code{\link{kanonymity}}, \code{\link{ldiversity}}, \code{\link{tcloseness}},
#'   \code{\link{individual_risk}}, \code{\link{singling_out}},
#'   \code{\link{linkability}}, \code{\link{nnaa}}, \code{\link{domias}}
#'
#' @author Matthias Templ
#' @family privacy-models
#' @export
#' @examples
#' # Create example data
#' set.seed(123)
#' n <- 200
#' original <- data.frame(
#'   age_group = sample(c("18-30", "31-45", "46-60", "60+"), n, replace = TRUE),
#'   gender = sample(c("M", "F"), n, replace = TRUE),
#'   region = sample(c("N", "S", "E", "W"), n, replace = TRUE),
#'   income = sample(c("low", "medium", "high"), n, replace = TRUE)
#' )
#'
#' # Good synthetic data
#' synthetic <- data.frame(
#'   age_group = sample(c("18-30", "31-45", "46-60", "60+"), n, replace = TRUE),
#'   gender = sample(c("M", "F"), n, replace = TRUE),
#'   region = sample(c("N", "S", "E", "W"), n, replace = TRUE),
#'   income = sample(c("low", "medium", "high"), n, replace = TRUE)
#' )
#'
#' # Generate comprehensive report
#' report <- disclosure_report(
#'   original, synthetic,
#'   key_vars = c("age_group", "gender", "region"),
#'   target_var = "income",
#'   seed = 42
#' )
#'
#' print(report)
#' summary(report)
#' plot(report)
disclosure_report <- function(X, ...) {
  UseMethod("disclosure_report")
}

#' @rdname disclosure_report
#' @export
disclosure_report.synth_pair <- function(X, ...) {
  disclosure_report.default(
    X = X$original,
    Y = X$synthetic,
    key_vars = X$key_vars,
    target_var = X$target_var,
    holdout = X$holdout,
    ...
  )
}

#' @rdname disclosure_report
#' @export
disclosure_report.default <- function(X, Y,
                              key_vars = NULL,
                              target_var = NULL,
                              holdout = NULL,
                              holdout_fraction = 0.5,
                              distance_vars = NULL,
                              distance_method = c("gower", "euclidean"),
                              compute = "all",
                              rapid_model = "rf",
                              na.rm = TRUE,
                              seed = NULL,
                              verbose = TRUE, ...) {

  # Input validation
  if (!is.data.frame(X)) stop("X must be a data frame.")
  if (!is.data.frame(Y)) stop("Y must be a data frame.")

  distance_method <- match.arg(distance_method)

  # Determine which metrics to compute
  attribution_metrics <- c("dcap", "tcap", "weap", "disco", "rapid")
  distance_metrics <- c("dcr", "nndr", "ims", "drisk", "hitting_rate", "rf_privacy")
  privacy_metrics <- c("kanonymity", "ldiversity", "tcloseness", "individual_risk")
  membership_metrics <- c("singling_out", "linkability", "nnaa", "domias")
  all_metrics <- c(attribution_metrics, distance_metrics,
                   privacy_metrics, membership_metrics)

  if ("all" %in% compute) {
    metrics_to_compute <- all_metrics
  } else {
    metrics_to_compute <- character(0)
    if ("attribution" %in% compute) {
      metrics_to_compute <- c(metrics_to_compute, attribution_metrics)
    }
    if ("distance" %in% compute) {
      metrics_to_compute <- c(metrics_to_compute, distance_metrics)
    }
    if ("privacy" %in% compute) {
      metrics_to_compute <- c(metrics_to_compute, privacy_metrics)
    }
    if ("membership" %in% compute) {
      metrics_to_compute <- c(metrics_to_compute, membership_metrics)
    }
    # Also allow specific metric names
    specific <- setdiff(compute, c("all", "attribution", "distance",
                                   "privacy", "membership"))
    if (length(specific) > 0) {
      metrics_to_compute <- unique(c(metrics_to_compute,
                                     intersect(tolower(specific), all_metrics)))
    }
    if (length(metrics_to_compute) == 0) {
      metrics_to_compute <- intersect(tolower(compute), all_metrics)
    }
  }

  # Check if attribution metrics can be computed
  can_compute_attribution <- !is.null(key_vars) && !is.null(target_var)
  if (!can_compute_attribution) {
    metrics_to_compute <- setdiff(metrics_to_compute, attribution_metrics)
    if (verbose && length(intersect(compute, c("all", "attribution", attribution_metrics))) > 0) {
      message("Note: Attribution metrics skipped (key_vars or target_var not provided)")
    }
  }

  # Check if privacy models can be computed (need at least key_vars)
  can_compute_privacy <- !is.null(key_vars)
  if (!can_compute_privacy) {
    metrics_to_compute <- setdiff(metrics_to_compute, privacy_metrics)
    if (verbose && length(intersect(compute, c("all", "privacy", privacy_metrics))) > 0) {
      message("Note: Privacy model metrics skipped (key_vars not provided)")
    }
  }

  # Privacy models needing target_var (ldiversity, tcloseness)
  needs_target <- c("ldiversity", "tcloseness")
  if (is.null(target_var) && any(needs_target %in% metrics_to_compute)) {
    metrics_to_compute <- setdiff(metrics_to_compute, needs_target)
    if (verbose) {
      message("Note: ldiversity/tcloseness skipped (target_var not provided)")
    }
  }

  # Initialize results
  results <- list()
  summary_rows <- list()

  # Set seed for reproducibility
  if (!is.null(seed)) set.seed(seed)

  # Create holdout if needed for distance or membership metrics
  needs_holdout <- c("dcr", "nndr", "rf_privacy", "singling_out", "linkability", "nnaa", "domias")
  if (any(needs_holdout %in% metrics_to_compute)) {
    if (is.null(holdout)) {
      n_holdout <- max(2, floor(nrow(X) * holdout_fraction))
      holdout_idx <- sample(nrow(X), n_holdout)
      holdout_data <- X[holdout_idx, , drop = FALSE]
      train_data <- X[-holdout_idx, , drop = FALSE]
    } else {
      holdout_data <- holdout
      train_data <- X
    }
  }

  # ============================================
  # ATTRIBUTION-BASED METRICS
  # ============================================

  if ("dcap" %in% metrics_to_compute && can_compute_attribution) {
    if (verbose) message("Computing DCAP...")
    tryCatch({
      results$dcap <- dcap(X, Y, key_vars, target_var, na.rm = na.rm)
      risk_ratio <- results$dcap$dcap / results$dcap$baseline
      pass <- !is.na(risk_ratio) && risk_ratio <= 1.5
      summary_rows$dcap <- data.frame(
        metric = "DCAP",
        value = round(results$dcap$dcap, 4),
        reference = round(results$dcap$baseline, 4),
        ratio = round(risk_ratio, 2),
        status = ifelse(is.na(pass), "N/A",
                        ifelse(pass, "PASS", "WARNING")),
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      if (verbose) message("  DCAP failed: ", e$message)
      summary_rows$dcap <<- data.frame(
        metric = "DCAP", value = NA, reference = NA, ratio = NA,
        status = "ERROR", stringsAsFactors = FALSE
      )
    })
  }

  if ("tcap" %in% metrics_to_compute && can_compute_attribution) {
    if (verbose) message("Computing TCAP...")
    tryCatch({
      results$tcap <- tcap(X, Y, key_vars, target_var, na.rm = na.rm)
      risk_ratio <- results$tcap$tcap_mean / results$tcap$baseline
      pass <- !is.na(risk_ratio) && risk_ratio <= 1.5
      summ <- summary(results$tcap)
      summary_rows$tcap <- data.frame(
        metric = "TCAP",
        value = round(results$tcap$tcap_mean, 4),
        reference = round(results$tcap$baseline, 4),
        ratio = round(risk_ratio, 2),
        status = ifelse(is.na(pass), "N/A",
                        ifelse(pass, "PASS", "WARNING")),
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      if (verbose) message("  TCAP failed: ", e$message)
      summary_rows$tcap <<- data.frame(
        metric = "TCAP", value = NA, reference = NA, ratio = NA,
        status = "ERROR", stringsAsFactors = FALSE
      )
    })
  }

  if ("weap" %in% metrics_to_compute && can_compute_attribution) {
    if (verbose) message("Computing WEAP...")
    tryCatch({
      results$weap <- weap(Y, key_vars, target_var, na.rm = na.rm)
      pass <- results$weap$pct_disclosive <= 5
      summary_rows$weap <- data.frame(
        metric = "WEAP",
        value = round(results$weap$weap_mean, 4),
        reference = paste0(results$weap$n_disclosive, " disclosive"),
        ratio = round(results$weap$pct_disclosive, 1),
        status = ifelse(pass, "PASS", "WARNING"),
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      if (verbose) message("  WEAP failed: ", e$message)
      summary_rows$weap <<- data.frame(
        metric = "WEAP", value = NA, reference = NA, ratio = NA,
        status = "ERROR", stringsAsFactors = FALSE
      )
    })
  }

  if ("disco" %in% metrics_to_compute && can_compute_attribution) {
    if (verbose) message("Computing DiSCO...")
    tryCatch({
      results$disco <- disco(X, Y, key_vars, target_var, na.rm = na.rm)
      pass <- is.na(results$disco$disco_ratio) || results$disco$disco_ratio <= 1.5
      summary_rows$disco <- data.frame(
        metric = "DiSCO",
        value = results$disco$n_disco,
        reference = round(results$disco$baseline_disco, 1),
        ratio = round(results$disco$disco_ratio, 2),
        status = ifelse(is.na(results$disco$disco_ratio), "N/A",
                        ifelse(pass, "PASS", "WARNING")),
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      if (verbose) message("  DiSCO failed: ", e$message)
      summary_rows$disco <<- data.frame(
        metric = "DiSCO", value = NA, reference = NA, ratio = NA,
        status = "ERROR", stringsAsFactors = FALSE
      )
    })
  }

  if ("rapid" %in% metrics_to_compute && can_compute_attribution) {
    if (verbose) message("Computing RAPID...")
    # Check if ranger is available
    if (!requireNamespace("ranger", quietly = TRUE)) {
      if (verbose) message("  RAPID skipped: 'ranger' package not installed")
      summary_rows$rapid <- data.frame(
        metric = "RAPID", value = NA, reference = NA, ratio = NA,
        status = "SKIPPED", stringsAsFactors = FALSE
      )
    } else {
      tryCatch({
        results$rapid <- rapid(X, Y, key_vars, target_var,
                               model_type = rapid_model, verbose = FALSE)
        pass <- results$rapid$rapid < 0.15
        summary_rows$rapid <- data.frame(
          metric = "RAPID",
          value = round(results$rapid$rapid, 4),
          reference = paste0(results$rapid$n_at_risk, " at risk"),
          ratio = round(results$rapid$pct_at_risk, 1),
          status = ifelse(pass, "PASS", "WARNING"),
          stringsAsFactors = FALSE
        )
      }, error = function(e) {
        if (verbose) message("  RAPID failed: ", e$message)
        summary_rows$rapid <<- data.frame(
          metric = "RAPID", value = NA, reference = NA, ratio = NA,
          status = "ERROR", stringsAsFactors = FALSE
        )
      })
    }
  }

  # ============================================
  # DISTANCE-BASED METRICS
  # ============================================

  if ("dcr" %in% metrics_to_compute) {
    if (verbose) message("Computing DCR...")
    tryCatch({
      results$dcr <- dcr(train_data, Y, holdout = holdout_data,
                         vars = distance_vars, method = distance_method,
                         na.rm = na.rm)
      summary_rows$dcr <- data.frame(
        metric = "DCR",
        value = round(results$dcr$dcr_ratio, 4),
        reference = paste0(round(100 * results$dcr$dcr_share, 1), "% closer"),
        ratio = round(results$dcr$dcr_ratio, 2),
        status = ifelse(results$dcr$privacy_pass, "PASS", "WARNING"),
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      if (verbose) message("  DCR failed: ", e$message)
      summary_rows$dcr <<- data.frame(
        metric = "DCR", value = NA, reference = NA, ratio = NA,
        status = "ERROR", stringsAsFactors = FALSE
      )
    })
  }

  if ("nndr" %in% metrics_to_compute) {
    if (verbose) message("Computing NNDR...")
    tryCatch({
      results$nndr <- nndr(train_data, Y, holdout = holdout_data,
                          vars = distance_vars, method = distance_method,
                          na.rm = na.rm)
      summary_rows$nndr <- data.frame(
        metric = "NNDR",
        value = round(results$nndr$nndr_ratio, 4),
        reference = paste0(results$nndr$n_suspicious, " suspicious"),
        ratio = round(results$nndr$nndr_ratio, 2),
        status = ifelse(results$nndr$privacy_pass, "PASS", "WARNING"),
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      if (verbose) message("  NNDR failed: ", e$message)
      summary_rows$nndr <<- data.frame(
        metric = "NNDR", value = NA, reference = NA, ratio = NA,
        status = "ERROR", stringsAsFactors = FALSE
      )
    })
  }

  if ("ims" %in% metrics_to_compute) {
    if (verbose) message("Computing IMS...")
    tryCatch({
      results$ims <- ims(X, Y, vars = distance_vars, na.rm = na.rm)
      summary_rows$ims <- data.frame(
        metric = "IMS",
        value = sprintf("%.2f%%", results$ims$ims_pct),
        reference = paste0(results$ims$n_identical, " copies"),
        ratio = results$ims$ims_pct,
        status = ifelse(results$ims$privacy_pass, "PASS", "WARNING"),
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      if (verbose) message("  IMS failed: ", e$message)
      summary_rows$ims <<- data.frame(
        metric = "IMS", value = NA, reference = NA, ratio = NA,
        status = "ERROR", stringsAsFactors = FALSE
      )
    })
  }

  # ============================================
  # ADDITIONAL DISTANCE-BASED METRICS
  # ============================================

  if ("drisk" %in% metrics_to_compute) {
    if (verbose) message("Computing dRisk...")
    tryCatch({
      results$drisk <- drisk(X, Y, vars = distance_vars, method = "both",
                             na.rm = na.rm)
      drisk_val <- if (!is.null(results$drisk$drisk_interval)) {
        results$drisk$drisk_interval
      } else {
        results$drisk$drisk_rmd
      }
      pass <- drisk_val < 0.05
      summary_rows$drisk <- data.frame(
        metric = "dRisk",
        value = round(drisk_val, 4),
        reference = "< 0.05",
        ratio = round(drisk_val, 4),
        status = ifelse(pass, "PASS", "WARNING"),
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      if (verbose) message("  dRisk failed: ", e$message)
      summary_rows$drisk <<- data.frame(
        metric = "dRisk", value = NA, reference = NA, ratio = NA,
        status = "ERROR", stringsAsFactors = FALSE
      )
    })
  }

  if ("hitting_rate" %in% metrics_to_compute) {
    if (verbose) message("Computing Hitting Rate...")
    tryCatch({
      results$hitting_rate <- hitting_rate(X, Y, vars = distance_vars,
                                           method = distance_method,
                                           na.rm = na.rm)
      pass <- results$hitting_rate$rate < 0.05
      summary_rows$hitting_rate <- data.frame(
        metric = "Hitting Rate",
        value = round(results$hitting_rate$rate, 4),
        reference = "< 0.05",
        ratio = round(results$hitting_rate$rate, 4),
        status = ifelse(pass, "PASS", "WARNING"),
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      if (verbose) message("  Hitting Rate failed: ", e$message)
      summary_rows$hitting_rate <<- data.frame(
        metric = "Hitting Rate", value = NA, reference = NA, ratio = NA,
        status = "ERROR", stringsAsFactors = FALSE
      )
    })
  }

  if ("rf_privacy" %in% metrics_to_compute) {
    if (verbose) message("Computing RF Privacy...")
    if (requireNamespace("ranger", quietly = TRUE)) {
      tryCatch({
        results$rf_privacy <- rf_privacy(train_data, Y,
                                          holdout = holdout_data,
                                          vars = distance_vars,
                                          na.rm = na.rm)
        summary_rows$rf_privacy <- data.frame(
          metric = "RF Privacy",
          value = round(results$rf_privacy$max_prox_share, 4),
          reference = paste0("ratio: ", round(results$rf_privacy$max_prox_ratio, 2)),
          ratio = round(results$rf_privacy$max_prox_ratio, 2),
          status = ifelse(results$rf_privacy$privacy_pass, "PASS", "WARNING"),
          stringsAsFactors = FALSE
        )
      }, error = function(e) {
        if (verbose) message("  RF Privacy failed: ", e$message)
        summary_rows$rf_privacy <<- data.frame(
          metric = "RF Privacy", value = NA, reference = NA, ratio = NA,
          status = "ERROR", stringsAsFactors = FALSE
        )
      })
    } else {
      if (verbose) message("  Skipping RF Privacy (ranger not installed)")
      summary_rows$rf_privacy <- data.frame(
        metric = "RF Privacy", value = NA, reference = NA, ratio = NA,
        status = "SKIPPED", stringsAsFactors = FALSE
      )
    }
  }

  # ============================================
  # PRIVACY MODELS (single-dataset)
  # ============================================

  if ("kanonymity" %in% metrics_to_compute && can_compute_privacy) {
    if (verbose) message("Computing k-Anonymity...")
    tryCatch({
      results$kanonymity <- kanonymity(Y, key_vars = key_vars, na.rm = na.rm)
      k_val <- results$kanonymity$k_level
      pass <- k_val >= 5
      summary_rows$kanonymity <- data.frame(
        metric = "k-Anonymity",
        value = k_val,
        reference = ">= 5",
        ratio = k_val,
        status = ifelse(pass, "PASS", ifelse(k_val >= 2, "WARNING", "FAIL")),
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      if (verbose) message("  k-Anonymity failed: ", e$message)
      summary_rows$kanonymity <<- data.frame(
        metric = "k-Anonymity", value = NA, reference = NA, ratio = NA,
        status = "ERROR", stringsAsFactors = FALSE
      )
    })
  }

  if ("ldiversity" %in% metrics_to_compute && can_compute_privacy) {
    if (verbose) message("Computing l-Diversity...")
    tryCatch({
      results$ldiversity <- ldiversity(Y, key_vars = key_vars,
                                       sensitive_var = target_var,
                                       na.rm = na.rm)
      l_val <- results$ldiversity$distinct_l
      pass <- l_val >= 2
      summary_rows$ldiversity <- data.frame(
        metric = "l-Diversity",
        value = l_val,
        reference = ">= 2",
        ratio = l_val,
        status = ifelse(pass, "PASS", "WARNING"),
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      if (verbose) message("  l-Diversity failed: ", e$message)
      summary_rows$ldiversity <<- data.frame(
        metric = "l-Diversity", value = NA, reference = NA, ratio = NA,
        status = "ERROR", stringsAsFactors = FALSE
      )
    })
  }

  if ("tcloseness" %in% metrics_to_compute && can_compute_privacy) {
    if (verbose) message("Computing t-Closeness...")
    tryCatch({
      results$tcloseness <- tcloseness(Y, key_vars = key_vars,
                                       sensitive_var = target_var,
                                       na.rm = na.rm)
      t_val <- results$tcloseness$t_achieved
      pass <- t_val <= 0.2
      summary_rows$tcloseness <- data.frame(
        metric = "t-Closeness",
        value = round(t_val, 4),
        reference = "<= 0.2",
        ratio = round(t_val, 4),
        status = ifelse(pass, "PASS", "WARNING"),
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      if (verbose) message("  t-Closeness failed: ", e$message)
      summary_rows$tcloseness <<- data.frame(
        metric = "t-Closeness", value = NA, reference = NA, ratio = NA,
        status = "ERROR", stringsAsFactors = FALSE
      )
    })
  }

  if ("individual_risk" %in% metrics_to_compute && can_compute_privacy) {
    if (verbose) message("Computing Individual Risk...")
    tryCatch({
      results$individual_risk <- individual_risk(Y, key_vars = key_vars,
                                                 na.rm = na.rm)
      mean_risk <- mean(results$individual_risk$risk, na.rm = TRUE)
      max_risk <- max(results$individual_risk$risk, na.rm = TRUE)
      n_high <- sum(results$individual_risk$risk > 0.1, na.rm = TRUE)
      pass <- mean_risk < 0.1
      summary_rows$individual_risk <- data.frame(
        metric = "Individual Risk",
        value = round(mean_risk, 4),
        reference = paste0(n_high, " high-risk (max=", round(max_risk, 3), ")"),
        ratio = round(mean_risk, 4),
        status = ifelse(pass, "PASS", "WARNING"),
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      if (verbose) message("  Individual Risk failed: ", e$message)
      summary_rows$individual_risk <<- data.frame(
        metric = "Individual Risk", value = NA, reference = NA, ratio = NA,
        status = "ERROR", stringsAsFactors = FALSE
      )
    })
  }

  # ============================================
  # MEMBERSHIP INFERENCE
  # ============================================

  if ("singling_out" %in% metrics_to_compute) {
    if (verbose) message("Computing Singling Out Risk...")
    tryCatch({
      results$singling_out <- singling_out(train_data, Y,
                                           holdout = holdout_data,
                                           n_attacks = 500,
                                           seed = seed, na.rm = na.rm)
      risk_val <- results$singling_out$risk
      pass <- risk_val < 0.1
      summary_rows$singling_out <- data.frame(
        metric = "Singling Out",
        value = round(risk_val, 4),
        reference = "< 0.1",
        ratio = round(risk_val, 4),
        status = ifelse(pass, "PASS", "WARNING"),
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      if (verbose) message("  Singling Out failed: ", e$message)
      summary_rows$singling_out <<- data.frame(
        metric = "Singling Out", value = NA, reference = NA, ratio = NA,
        status = "ERROR", stringsAsFactors = FALSE
      )
    })
  }

  if ("linkability" %in% metrics_to_compute) {
    if (verbose) message("Computing Linkability Risk...")
    tryCatch({
      results$linkability <- linkability(train_data, Y,
                                         holdout = holdout_data,
                                         n_attacks = 500,
                                         seed = seed, na.rm = na.rm)
      risk_val <- results$linkability$risk
      pass <- risk_val < 0.1
      summary_rows$linkability <- data.frame(
        metric = "Linkability",
        value = round(risk_val, 4),
        reference = "< 0.1",
        ratio = round(risk_val, 4),
        status = ifelse(pass, "PASS", "WARNING"),
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      if (verbose) message("  Linkability failed: ", e$message)
      summary_rows$linkability <<- data.frame(
        metric = "Linkability", value = NA, reference = NA, ratio = NA,
        status = "ERROR", stringsAsFactors = FALSE
      )
    })
  }

  if ("nnaa" %in% metrics_to_compute) {
    if (verbose) message("Computing NNAA...")
    tryCatch({
      results$nnaa <- nnaa(train_data, Y,
                           holdout = holdout_data,
                           vars = distance_vars,
                           method = distance_method,
                           na.rm = na.rm, seed = seed)
      pl_val <- results$nnaa$privacy_loss
      pass <- pl_val < 0.05
      summary_rows$nnaa <- data.frame(
        metric = "NNAA",
        value = round(pl_val, 4),
        reference = "< 0.05",
        ratio = round(pl_val, 4),
        status = ifelse(pass, "PASS", "WARNING"),
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      if (verbose) message("  NNAA failed: ", e$message)
      summary_rows$nnaa <<- data.frame(
        metric = "NNAA", value = NA, reference = NA, ratio = NA,
        status = "ERROR", stringsAsFactors = FALSE
      )
    })
  }

  if ("domias" %in% metrics_to_compute) {
    if (verbose) message("Computing DOMIAS...")
    if (!requireNamespace("ranger", quietly = TRUE)) {
      if (verbose) message("  DOMIAS skipped: 'ranger' package not installed")
      summary_rows$domias <- data.frame(
        metric = "DOMIAS", value = NA, reference = NA, ratio = NA,
        status = "SKIPPED", stringsAsFactors = FALSE
      )
    } else {
      tryCatch({
        results$domias <- domias(train_data, Y,
                                 holdout = holdout_data,
                                 vars = distance_vars,
                                 na.rm = na.rm, seed = seed)
        auc_val <- results$domias$auc
        pass <- auc_val < 0.6
        summary_rows$domias <- data.frame(
          metric = "DOMIAS",
          value = round(auc_val, 4),
          reference = "< 0.6",
          ratio = round(auc_val, 4),
          status = ifelse(pass, "PASS", "WARNING"),
          stringsAsFactors = FALSE
        )
      }, error = function(e) {
        if (verbose) message("  DOMIAS failed: ", e$message)
        summary_rows$domias <<- data.frame(
          metric = "DOMIAS", value = NA, reference = NA, ratio = NA,
          status = "ERROR", stringsAsFactors = FALSE
        )
      })
    }
  }

  # ============================================
  # AGGREGATE RESULTS
  # ============================================

  # Combine summary rows
  summary_df <- do.call(rbind, summary_rows)
  rownames(summary_df) <- NULL

  # Count passes and failures
  n_pass <- sum(summary_df$status == "PASS", na.rm = TRUE)
  n_warn <- sum(summary_df$status == "WARNING", na.rm = TRUE)
  n_error <- sum(summary_df$status == "ERROR", na.rm = TRUE)
  n_na <- sum(summary_df$status == "N/A", na.rm = TRUE)

  # Determine overall risk
  if (n_warn == 0 && n_error == 0) {
    overall_risk <- "LOW"
  } else if (n_warn <= 2 && n_error == 0) {
    overall_risk <- "MEDIUM"
  } else {
    overall_risk <- "HIGH"
  }

  # Build report object
  report <- list(
    results = results,
    summary = summary_df,
    overall_risk = overall_risk,
    n_pass = n_pass,
    n_warn = n_warn,
    n_error = n_error,
    n_na = n_na,
    n_total = nrow(summary_df),
    parameters = list(
      n_original = nrow(X),
      n_synthetic = nrow(Y),
      key_vars = key_vars,
      target_var = target_var,
      distance_method = distance_method,
      holdout_fraction = holdout_fraction,
      metrics_computed = names(results)
    )
  )

  class(report) <- "disclosure_report"
  return(report)
}


#' Print method for disclosure_report objects
#'
#' @param x an object of class "disclosure_report"
#' @param ... additional arguments (ignored)
#' @export
print.disclosure_report <- function(x, ...) {
  cat("\n")
  cat("==============================================================\n")
  cat("           DISCLOSURE RISK REPORT\n")
  cat("==============================================================\n\n")

  cat("Data Summary:\n")
  cat("  Original records:", x$parameters$n_original, "\n")
  cat("  Synthetic records:", x$parameters$n_synthetic, "\n")
  if (!is.null(x$parameters$key_vars)) {
    cat("  Key variables:", paste(x$parameters$key_vars, collapse = ", "), "\n")
    cat("  Target variable:", x$parameters$target_var, "\n")
  }
  cat("\n")

  cat("--------------------------------------------------------------\n")
  cat("                    RESULTS SUMMARY\n")
  cat("--------------------------------------------------------------\n\n")

  # Print summary table
  print(x$summary, row.names = FALSE)

  cat("\n--------------------------------------------------------------\n")
  cat("                   OVERALL ASSESSMENT\n")
  cat("--------------------------------------------------------------\n\n")

  cat("  Metrics passed:", x$n_pass, "/", x$n_total, "\n")
  if (x$n_warn > 0) cat("  Metrics with warnings:", x$n_warn, "\n")
  if (x$n_error > 0) cat("  Metrics with errors:", x$n_error, "\n")
  if (x$n_na > 0) cat("  Metrics not applicable:", x$n_na, "\n")
  cat("\n")

  cat("  OVERALL RISK LEVEL:", x$overall_risk, "\n")

  if (x$overall_risk == "LOW") {
    cat("\n  Interpretation: Synthetic data shows good privacy preservation.\n")
    cat("  No significant disclosure risks detected.\n")
  } else if (x$overall_risk == "MEDIUM") {
    cat("\n  Interpretation: Some potential privacy concerns detected.\n")
    cat("  Review the flagged metrics and consider adjustments.\n")
  } else {
    cat("\n  Interpretation: Significant privacy concerns detected.\n")
    cat("  Synthetic data may leak information about original records.\n")
    cat("  Consider revising the synthesis approach.\n")
  }

  cat("\n==============================================================\n")

  invisible(x)
}


#' Summary method for disclosure_report objects
#'
#' @param object an object of class "disclosure_report"
#' @param ... additional arguments (ignored)
#' @export
summary.disclosure_report <- function(object, ...) {
  # Create detailed summary
  summ <- list(
    overall_risk = object$overall_risk,
    n_pass = object$n_pass,
    n_warn = object$n_warn,
    n_total = object$n_total,
    summary_table = object$summary,
    attribution_results = list(),
    distance_results = list(),
    privacy_results = list(),
    membership_results = list()
  )

  # Extract key values from attribution metrics
  if (!is.null(object$results$dcap)) {
    summ$attribution_results$dcap <- list(
      value = object$results$dcap$dcap,
      baseline = object$results$dcap$baseline,
      risk_ratio = object$results$dcap$dcap / object$results$dcap$baseline,
      n_matched = object$results$dcap$n_matched
    )
  }

  if (!is.null(object$results$tcap)) {
    s <- summary(object$results$tcap)
    summ$attribution_results$tcap <- list(
      mean = s$tcap_mean,
      median = s$tcap_median,
      n_high_risk = s$n_high_risk,
      pct_high_risk = s$pct_high_risk
    )
  }

  if (!is.null(object$results$rapid)) {
    summ$attribution_results$rapid <- list(
      value = object$results$rapid$rapid,
      n_at_risk = object$results$rapid$n_at_risk,
      pct_at_risk = object$results$rapid$pct_at_risk,
      model_type = object$results$rapid$model_type
    )
  }

  # Extract key values from distance metrics
  if (!is.null(object$results$dcr)) {
    summ$distance_results$dcr <- list(
      ratio = object$results$dcr$dcr_ratio,
      share = object$results$dcr$dcr_share,
      pass = object$results$dcr$privacy_pass
    )
  }

  if (!is.null(object$results$nndr)) {
    summ$distance_results$nndr <- list(
      ratio = object$results$nndr$nndr_ratio,
      n_suspicious = object$results$nndr$n_suspicious,
      pass = object$results$nndr$privacy_pass
    )
  }

  if (!is.null(object$results$ims)) {
    summ$distance_results$ims <- list(
      pct = object$results$ims$ims_pct,
      n_identical = object$results$ims$n_identical,
      pass = object$results$ims$privacy_pass
    )
  }

  if (!is.null(object$results$rf_privacy)) {
    summ$distance_results$rf_privacy <- list(
      max_prox_share = object$results$rf_privacy$max_prox_share,
      max_prox_ratio = object$results$rf_privacy$max_prox_ratio,
      pass = object$results$rf_privacy$privacy_pass
    )
  }

  if (!is.null(object$results$drisk)) {
    summ$distance_results$drisk <- list(
      interval = object$results$drisk$drisk_interval,
      rmd = object$results$drisk$drisk_rmd
    )
  }

  if (!is.null(object$results$hitting_rate)) {
    summ$distance_results$hitting_rate <- list(
      rate = object$results$hitting_rate$rate
    )
  }

  # Extract privacy model results
  if (!is.null(object$results$kanonymity)) {
    summ$privacy_results$kanonymity <- list(
      k_level = object$results$kanonymity$k_level
    )
  }

  if (!is.null(object$results$ldiversity)) {
    summ$privacy_results$ldiversity <- list(
      distinct_l = object$results$ldiversity$distinct_l
    )
  }

  if (!is.null(object$results$tcloseness)) {
    summ$privacy_results$tcloseness <- list(
      t_achieved = object$results$tcloseness$t_achieved
    )
  }

  if (!is.null(object$results$individual_risk)) {
    summ$privacy_results$individual_risk <- list(
      mean_risk = mean(object$results$individual_risk$risk, na.rm = TRUE),
      max_risk = max(object$results$individual_risk$risk, na.rm = TRUE)
    )
  }

  # Extract membership inference results
  if (!is.null(object$results$singling_out)) {
    summ$membership_results$singling_out <- list(
      risk = object$results$singling_out$risk
    )
  }

  if (!is.null(object$results$linkability)) {
    summ$membership_results$linkability <- list(
      risk = object$results$linkability$risk
    )
  }

  if (!is.null(object$results$nnaa)) {
    summ$membership_results$nnaa <- list(
      privacy_loss = object$results$nnaa$privacy_loss
    )
  }

  if (!is.null(object$results$domias)) {
    summ$membership_results$domias <- list(
      auc = object$results$domias$auc
    )
  }

  summ$parameters <- object$parameters

  class(summ) <- "summary.disclosure_report"
  return(summ)
}


#' Print method for summary.disclosure_report objects
#'
#' @param x an object of class "summary.disclosure_report"
#' @param ... additional arguments (ignored)
#' @export
print.summary.disclosure_report <- function(x, ...) {
  cat("Disclosure Risk Report Summary\n")
  cat("==============================\n\n")

  cat("Overall Risk:", x$overall_risk, "\n")
  cat("Metrics:", x$n_pass, "passed,", x$n_warn, "warnings,",
      "of", x$n_total, "total\n\n")

  if (length(x$attribution_results) > 0) {
    cat("Attribution-Based Metrics:\n")
    if (!is.null(x$attribution_results$dcap)) {
      cat("  DCAP:", round(x$attribution_results$dcap$value, 4),
          "(baseline:", round(x$attribution_results$dcap$baseline, 4),
          ", ratio:", round(x$attribution_results$dcap$risk_ratio, 2), ")\n")
    }
    if (!is.null(x$attribution_results$tcap)) {
      cat("  TCAP mean:", round(x$attribution_results$tcap$mean, 4),
          ", high-risk records:", x$attribution_results$tcap$n_high_risk,
          sprintf("(%.1f%%)", x$attribution_results$tcap$pct_high_risk), "\n")
    }
    if (!is.null(x$attribution_results$rapid)) {
      cat("  RAPID:", round(x$attribution_results$rapid$value, 4),
          ", at-risk records:", x$attribution_results$rapid$n_at_risk,
          sprintf("(%.1f%%)", x$attribution_results$rapid$pct_at_risk), "\n")
    }
    cat("\n")
  }

  if (length(x$distance_results) > 0) {
    cat("Distance-Based Metrics:\n")
    if (!is.null(x$distance_results$dcr)) {
      cat("  DCR ratio:", round(x$distance_results$dcr$ratio, 4),
          ", share:", sprintf("%.1f%%", 100 * x$distance_results$dcr$share),
          ifelse(x$distance_results$dcr$pass, "(PASS)", "(WARNING)"), "\n")
    }
    if (!is.null(x$distance_results$nndr)) {
      cat("  NNDR ratio:", round(x$distance_results$nndr$ratio, 4),
          ", suspicious:", x$distance_results$nndr$n_suspicious,
          ifelse(x$distance_results$nndr$pass, "(PASS)", "(WARNING)"), "\n")
    }
    if (!is.null(x$distance_results$ims)) {
      cat("  IMS:", sprintf("%.2f%%", x$distance_results$ims$pct),
          ", copies:", x$distance_results$ims$n_identical,
          ifelse(x$distance_results$ims$pass, "(PASS)", "(WARNING)"), "\n")
    }
    if (!is.null(x$distance_results$rf_privacy)) {
      cat("  RF Privacy share:", round(x$distance_results$rf_privacy$max_prox_share, 4),
          ", ratio:", round(x$distance_results$rf_privacy$max_prox_ratio, 2),
          ifelse(x$distance_results$rf_privacy$pass, "(PASS)", "(WARNING)"), "\n")
    }
    if (!is.null(x$distance_results$drisk)) {
      cat("  dRisk interval:", round(x$distance_results$drisk$interval, 4),
          ", RMD:", round(x$distance_results$drisk$rmd, 4), "\n")
    }
    if (!is.null(x$distance_results$hitting_rate)) {
      cat("  Hitting Rate:", round(x$distance_results$hitting_rate$rate, 4), "\n")
    }
    cat("\n")
  }

  if (length(x$privacy_results) > 0) {
    cat("Privacy Models:\n")
    if (!is.null(x$privacy_results$kanonymity)) {
      cat("  k-Anonymity: k =", x$privacy_results$kanonymity$k_level, "\n")
    }
    if (!is.null(x$privacy_results$ldiversity)) {
      cat("  l-Diversity: l =", x$privacy_results$ldiversity$distinct_l, "\n")
    }
    if (!is.null(x$privacy_results$tcloseness)) {
      cat("  t-Closeness: t =", round(x$privacy_results$tcloseness$t_achieved, 4), "\n")
    }
    if (!is.null(x$privacy_results$individual_risk)) {
      cat("  Individual Risk: mean =", round(x$privacy_results$individual_risk$mean_risk, 4),
          ", max =", round(x$privacy_results$individual_risk$max_risk, 4), "\n")
    }
    cat("\n")
  }

  if (length(x$membership_results) > 0) {
    cat("Membership Inference:\n")
    if (!is.null(x$membership_results$singling_out)) {
      cat("  Singling Out risk:", round(x$membership_results$singling_out$risk, 4), "\n")
    }
    if (!is.null(x$membership_results$linkability)) {
      cat("  Linkability risk:", round(x$membership_results$linkability$risk, 4), "\n")
    }
    if (!is.null(x$membership_results$nnaa)) {
      cat("  NNAA privacy loss:", round(x$membership_results$nnaa$privacy_loss, 4), "\n")
    }
    if (!is.null(x$membership_results$domias)) {
      cat("  DOMIAS AUC:", round(x$membership_results$domias$auc, 4), "\n")
    }
  }

  invisible(x)
}


#' Plot method for disclosure_report objects
#'
#' @param x an object of class "disclosure_report"
#' @param y not used
#' @param ... additional arguments passed to plotting functions
#' @param which which plot to show: 1 = summary bar chart, 2 = metric details
#' @export
#' @importFrom graphics barplot par text mtext axis box
plot.disclosure_report <- function(x, y = NULL, ..., which = 1) {
  show <- rep(FALSE, 2)
  show[which] <- TRUE

  if (sum(show) > 1) {
    op <- par(mfrow = c(1, 2))
    on.exit(par(op))
  }

  if (show[1]) {
    # Summary bar chart of pass/warn/error
    counts <- c(x$n_pass, x$n_warn, x$n_error)
    names(counts) <- c("Pass", "Warning", "Error")
    colors <- c("forestgreen", "orange", "firebrick")

    barplot(counts,
            main = paste("Disclosure Risk Assessment\nOverall:", x$overall_risk),
            ylab = "Number of Metrics",
            col = colors,
            ylim = c(0, max(counts) * 1.2), ...)

    # Add count labels
    text(x = seq(0.7, by = 1.2, length.out = 3),
         y = counts + 0.3,
         labels = counts,
         cex = 1.2)
  }

  if (show[2]) {
    # Detailed metric values
    summ <- x$summary
    n_metrics <- nrow(summ)

    # Create color-coded bar chart
    colors <- ifelse(summ$status == "PASS", "forestgreen",
                     ifelse(summ$status == "WARNING", "orange",
                            ifelse(summ$status == "ERROR", "firebrick", "gray")))

    # Use ratio values where available, normalized
    values <- as.numeric(summ$ratio)
    values[is.na(values)] <- 0

    # Cap at reasonable maximum for display
    values <- pmin(values, 3)

    barplot(values,
            names.arg = summ$metric,
            main = "Metric Ratios (capped at 3)",
            ylab = "Ratio Value",
            col = colors,
            las = 2, ...)

    # Reference line at 1.0
    abline(h = 1, col = "blue", lty = 2, lwd = 2)
    abline(h = 1.5, col = "red", lty = 3, lwd = 2)
    legend("topright",
           legend = c("Ideal (1.0)", "Threshold (1.5)"),
           col = c("blue", "red"),
           lty = c(2, 3), lwd = 2, cex = 0.8)
  }
}
