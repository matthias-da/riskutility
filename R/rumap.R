#' Multivariate Risk-Utility Map (RU-Map)
#'
#' Computes multiple disclosure risk and data utility measures for synthetic/anonymized
#' data and provides comprehensive visualization tools for multivariate evaluation.
#' This implements the framework from "Beyond the Trade-off Curve" (Thees, Müller, Templ 2026).
#'
#' @param X A synth_pair object, or a data.frame containing the original dataset.
#' @param synthetic A data.frame or named list of data.frames containing synthetic dataset(s).
#'   If a named list, each element is treated as a different synthetic data generator (SDG).
#' @param risk_measures Character vector of risk measures to compute. Options:
#'   "dcap", "tcap", "disco", "ims", "repu", "dcr", "nndr", "rapid", "rf_privacy".
#'   Default includes all except distance-based measures.
#' @param utility_measures Character vector of utility measures to compute. Options:
#'   "pmse", "wasserstein", "hellinger", "energy", "ci_proximity". Default includes all.
#' @param key_vars Character vector of quasi-identifier variables (for attribution-based risk).
#' @param target_var Character string naming the sensitive target variable (for CAP metrics).
#' @param holdout Optional data.frame for distance-based risk measures. If NULL and
#'   distance-based measures are requested, uses holdout_fraction to split original.
#' @param holdout_fraction Numeric, fraction of original to use as holdout if holdout
#'   is not provided. Default 0.2.
#' @param vars Character vector of variables to use for utility measures. If NULL,
#'   uses all common variables.
#' @param cat_vars Character vector of categorical variables for Hellinger distance.
#'   If NULL, auto-detected.
#' @param num_vars Character vector of numeric variables for energy distance.
#'   If NULL, auto-detected.
#' @param normalize Logical, whether to normalize measures to the 0-1 range. Default TRUE.
#' @param seed Integer, random seed for reproducibility.
#' @param na.rm Logical, whether to remove NA values. Default TRUE.
#' @param ... additional arguments passed to methods
#'
#' @return An object of class "rumap" containing:
#' \itemize{
#'   \item risk: data.frame of raw risk measures per SDG
#'   \item utility: data.frame of raw utility measures per SDG
#'   \item normalized: data.frame of normalized measures (if normalize=TRUE)
#'   \item composites: data.frame with composite risk/utility scores
#'   \item pareto: logical vector indicating Pareto-optimal SDGs
#'   \item pareto_sdgs: names of Pareto-optimal SDGs
#'   \item n_sdgs: number of synthetic datasets evaluated
#'   \item sdg_names: names of SDGs
#'   \item risk_measures, utility_measures: measures computed
#'   \item metadata: list of parameters used
#' }
#'
#' @details
#' The rumap function provides a comprehensive framework for evaluating synthetic data
#' by computing multiple risk and utility measures simultaneously. It supports:
#'
#' \strong{Risk Measures:}
#' \itemize{
#'   \item dcap: Differential Correct Attribution Probability
#'   \item tcap: Targeted CAP (mean per-record risk)
#'   \item disco: Disclosive in Synthetic Correct Original
#'   \item rapid: Risk of Attribute Prediction-Induced Disclosure (ML-based)
#'   \item ims: Identical Match Share
#'   \item repu: Replicated Uniques
#'   \item dcr: Distance to Closest Record ratio
#'   \item nndr: Nearest Neighbor Distance Ratio
#'   \item rf_privacy: Random forest proximity-based memorization test (requires ranger)
#' }
#'
#' \strong{Utility Measures:}
#' \itemize{
#'   \item pmse: Propensity Mean Squared Error (from propscore)
#'   \item wasserstein: Mean Wasserstein distance across numeric variables
#'   \item hellinger: Mean Hellinger distance across categorical variables
#'   \item energy: Energy distance for multivariate numeric data
#'   \item ci_proximity: Confidence interval proximity score
#' }
#'
#' All measures are normalized to the 0-1 range with consistent direction:
#' higher risk values = higher disclosure risk (bad),
#' higher utility values = higher utility (good).
#'
#' Pareto-optimal SDGs are identified as those not dominated by any other SDG
#' (cannot improve risk without worsening utility, and vice versa).
#'
#' @seealso \code{\link{plot.rumap}} for visualization options,
#'   \code{\link{dcap}}, \code{\link{tcap}}, \code{\link{disco}}, \code{\link{ims}},
#'   \code{\link{propscore}}, \code{\link{hellinger}}, \code{\link{energy_distance}}
#'
#' @references
#' Thees, O., Müller, R., & Templ, M. (2026). Beyond the Trade-off Curve:
#' Multivariate and Advanced Risk-Utility Maps for Evaluating Anonymized and
#' Synthetic Data. Journal of Official Statistics.
#'
#' @author Matthias Templ
#' @family utility
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' # Create original data
#' original <- data.frame(
#'   age = sample(20:70, 200, replace = TRUE),
#'   gender = sample(c("M", "F"), 200, replace = TRUE),
#'   region = sample(c("N", "S", "E", "W"), 200, replace = TRUE),
#'   income = rnorm(200, 50000, 15000)
#' )
#'
#' # Create synthetic datasets (simulating different SDGs)
#' synth_good <- data.frame(
#'   age = sample(20:70, 200, replace = TRUE),
#'   gender = sample(c("M", "F"), 200, replace = TRUE),
#'   region = sample(c("N", "S", "E", "W"), 200, replace = TRUE),
#'   income = rnorm(200, 50000, 15000)
#' )
#'
#' synth_poor <- data.frame(
#'   age = sample(30:50, 200, replace = TRUE),
#'   gender = sample(c("M", "F"), 200, replace = TRUE, prob = c(0.8, 0.2)),
#'   region = sample(c("N", "S", "E", "W"), 200, replace = TRUE, prob = c(0.7, 0.1, 0.1, 0.1)),
#'   income = rnorm(200, 60000, 20000)
#' )
#'
#' # Compare multiple SDGs
#' result <- rumap(
#'   original = original,
#'   synthetic = list(good_sdg = synth_good, poor_sdg = synth_poor),
#'   key_vars = c("age", "gender", "region"),
#'   target_var = "income",
#'   risk_measures = c("ims", "dcap"),
#'   utility_measures = c("hellinger", "energy", "ci_proximity"),
#'   seed = 42
#' )
#'
#' print(result)
#' summary(result)
#' }
rumap <- function(X, ...) {
  UseMethod("rumap")
}

#' @rdname rumap
#' @export
rumap.synth_pair <- function(X, ...) {
  rumap.default(
    X = X$original,
    synthetic = X$synthetic,
    key_vars = X$key_vars,
    target_var = X$target_var,
    holdout = X$holdout,
    vars = X$vars,
    cat_vars = X$cat_vars,
    num_vars = X$num_vars,
    ...
  )
}

#' @rdname rumap
#' @export
rumap.default <- function(X,
                  synthetic,
                  risk_measures = c("dcap", "tcap", "disco", "rapid", "ims", "repu"),
                  utility_measures = c("pmse", "wasserstein", "hellinger", "energy", "ci_proximity"),
                  key_vars = NULL,
                  target_var = NULL,
                  holdout = NULL,
                  holdout_fraction = 0.2,
                  vars = NULL,
                  cat_vars = NULL,
                  num_vars = NULL,
                  normalize = TRUE,
                  seed = NULL,
                  na.rm = TRUE, ...) {

  original <- X

  # Set seed for reproducibility
  if (!is.null(seed)) set.seed(seed)

  # Input validation
  if (!is.data.frame(original)) {
    stop("original must be a data.frame")
  }

 # Handle synthetic input - convert single data.frame to named list
  if (is.data.frame(synthetic)) {
    synthetic <- list(synthetic = synthetic)
  }

  if (!is.list(synthetic) || length(synthetic) == 0) {
    stop("synthetic must be a data.frame or a non-empty named list of data.frames")
  }

  # Ensure all elements are data.frames
  for (i in seq_along(synthetic)) {
    if (!is.data.frame(synthetic[[i]])) {
      stop(paste("Element", i, "of synthetic is not a data.frame"))
    }
  }

  # Assign names if missing
  if (is.null(names(synthetic))) {
    names(synthetic) <- paste0("SDG_", seq_along(synthetic))
  }

  sdg_names <- names(synthetic)
  n_sdgs <- length(synthetic)

  # Validate risk measures
  valid_risk <- c("dcap", "tcap", "disco", "rapid", "ims", "repu", "dcr", "nndr", "rf_privacy")
  invalid_risk <- setdiff(risk_measures, valid_risk)
  if (length(invalid_risk) > 0) {
    warning(paste("Unknown risk measures ignored:", paste(invalid_risk, collapse = ", ")))
    risk_measures <- intersect(risk_measures, valid_risk)
  }

  # Validate utility measures
  valid_utility <- c("pmse", "wasserstein", "hellinger", "energy", "ci_proximity",
                     "mmd", "tstr", "copula", "tail", "contingency")
  invalid_utility <- setdiff(utility_measures, valid_utility)
  if (length(invalid_utility) > 0) {
    warning(paste("Unknown utility measures ignored:", paste(invalid_utility, collapse = ", ")))
    utility_measures <- intersect(utility_measures, valid_utility)
  }

  # Check required parameters for attribution-based risk
  attr_risk <- intersect(risk_measures, c("dcap", "tcap", "disco", "rapid"))
  if (length(attr_risk) > 0) {
    if (is.null(key_vars)) {
      warning("key_vars not specified; attribution-based risk measures will be skipped")
      risk_measures <- setdiff(risk_measures, attr_risk)
    }
    if (is.null(target_var) && any(c("dcap", "rapid") %in% risk_measures)) {
      skip_target <- intersect(c("dcap", "rapid"), risk_measures)
      warning(paste(paste(skip_target, collapse = ", "),
                    "require target_var; will be skipped"))
      risk_measures <- setdiff(risk_measures, skip_target)
    }
  }

  # Auto-detect variable types if not specified
  common_vars <- Reduce(intersect, lapply(synthetic, names))
  common_vars <- intersect(names(original), common_vars)

  if (is.null(vars)) {
    vars <- common_vars
  }

  if (is.null(cat_vars)) {
    cat_vars <- names(original)[sapply(original[, vars, drop = FALSE],
                                        function(x) is.factor(x) || is.character(x))]
  }

  if (is.null(num_vars)) {
    num_vars <- names(original)[sapply(original[, vars, drop = FALSE], is.numeric)]
  }

  # Create holdout if needed for distance-based measures
  dist_risk <- intersect(risk_measures, c("dcr", "nndr", "rf_privacy"))
  if (length(dist_risk) > 0 && is.null(holdout)) {
    n_holdout <- floor(nrow(original) * holdout_fraction)
    holdout_idx <- sample(nrow(original), n_holdout)
    holdout <- original[holdout_idx, , drop = FALSE]
    train <- original[-holdout_idx, , drop = FALSE]
  } else {
    train <- original
  }

  # Initialize result storage
  risk_results <- data.frame(sdg = sdg_names, stringsAsFactors = FALSE)
  utility_results <- data.frame(sdg = sdg_names, stringsAsFactors = FALSE)

  # Compute measures for each SDG
  for (i in seq_along(synthetic)) {
    sdg_name <- sdg_names[i]
    Y <- synthetic[[i]]

    # ---- RISK MEASURES ----

    # DCAP
    if ("dcap" %in% risk_measures) {
      tryCatch({
        res <- dcap(original, Y, key_vars = key_vars, target_var = target_var)
        risk_results[i, "dcap"] <- res$dcap
      }, error = function(e) {
        warning(paste("dcap failed for", sdg_name, ":", e$message))
        risk_results[i, "dcap"] <<- NA
      })
    }

    # TCAP
    if ("tcap" %in% risk_measures) {
      tryCatch({
        res <- tcap(original, Y, key_vars = key_vars, target_var = target_var)
        risk_results[i, "tcap"] <- res$tcap_mean
      }, error = function(e) {
        warning(paste("tcap failed for", sdg_name, ":", e$message))
        risk_results[i, "tcap"] <<- NA
      })
    }

    # DiSCO
    if ("disco" %in% risk_measures) {
      tryCatch({
        res <- disco(original, Y, key_vars = key_vars, target_var = target_var)
        risk_results[i, "disco"] <- res$pct_disco / 100  # Convert to proportion
      }, error = function(e) {
        warning(paste("disco failed for", sdg_name, ":", e$message))
        risk_results[i, "disco"] <<- NA
      })
    }

    # RAPID
    if ("rapid" %in% risk_measures) {
      tryCatch({
        res <- rapid(original, Y, key_vars = key_vars, target_var = target_var,
                     model_type = "lm", return_all_records = FALSE,
                     store_model = FALSE, verbose = FALSE)
        risk_results[i, "rapid"] <- res$rapid
      }, error = function(e) {
        warning(paste("rapid failed for", sdg_name, ":", e$message))
        risk_results[i, "rapid"] <<- NA
      })
    }

    # IMS
    if ("ims" %in% risk_measures) {
      tryCatch({
        res <- ims(original, Y, vars = vars, na.rm = na.rm)
        risk_results[i, "ims"] <- res$ims
      }, error = function(e) {
        warning(paste("ims failed for", sdg_name, ":", e$message))
        risk_results[i, "ims"] <<- NA
      })
    }

    # RepU
    if ("repu" %in% risk_measures) {
      tryCatch({
        res <- repu(original, Y, vars = vars, uniques_only = TRUE, na.rm = na.rm)
        risk_results[i, "repu"] <- res$ims
      }, error = function(e) {
        warning(paste("repu failed for", sdg_name, ":", e$message))
        risk_results[i, "repu"] <<- NA
      })
    }

    # DCR
    if ("dcr" %in% risk_measures) {
      tryCatch({
        res <- dcr(train, Y, holdout = holdout, vars = vars)
        risk_results[i, "dcr"] <- 1 - res$dcr_share  # Transform: higher = more risk
      }, error = function(e) {
        warning(paste("dcr failed for", sdg_name, ":", e$message))
        risk_results[i, "dcr"] <<- NA
      })
    }

    # NNDR
    if ("nndr" %in% risk_measures) {
      tryCatch({
        res <- nndr(train, Y, holdout = holdout, vars = vars)
        risk_results[i, "nndr"] <- 1 - res$nndr_ratio  # Transform: higher = more risk
      }, error = function(e) {
        warning(paste("nndr failed for", sdg_name, ":", e$message))
        risk_results[i, "nndr"] <<- NA
      })
    }

    # RF Privacy
    if ("rf_privacy" %in% risk_measures) {
      if (requireNamespace("ranger", quietly = TRUE)) {
        tryCatch({
          res <- rf_privacy(train, Y, holdout = holdout, vars = vars,
                            na.rm = na.rm, null_test = FALSE)
          risk_results[i, "rf_privacy"] <- res$max_prox_share
        }, error = function(e) {
          warning(paste("rf_privacy failed for", sdg_name, ":", e$message))
          risk_results[i, "rf_privacy"] <<- NA
        })
      } else {
        risk_results[i, "rf_privacy"] <- NA
      }
    }

    # ---- UTILITY MEASURES ----

    # pMSE (from propscore)
    if ("pmse" %in% risk_measures || "pmse" %in% utility_measures) {
      tryCatch({
        # Build explicit formula using available variables
        form_vars <- intersect(vars, names(Y))
        if (length(form_vars) > 0) {
          form_str <- paste("group ~", paste(form_vars, collapse = " + "))
          form <- as.formula(form_str)
          res <- propscore(original, Y, form = form, na = "remove")
        } else {
          res <- propscore(original, Y, na = "remove")
        }
        # ps_score: lower is better, transform to utility (higher = better)
        pmse_val <- res$ps_score
        utility_results[i, "pmse"] <- 1 - min(1, pmse_val * 10)  # Scale and invert
      }, error = function(e) {
        warning(paste("pmse failed for", sdg_name, ":", e$message))
        utility_results[i, "pmse"] <<- NA
      })
    }

    # Wasserstein
    if ("wasserstein" %in% utility_measures && length(num_vars) > 0) {
      tryCatch({
        w_vals <- numeric(length(num_vars))
        for (j in seq_along(num_vars)) {
          res <- compare_wasserstein(original, Y, num_var = num_vars[j], var_type = "continuous")
          w_vals[j] <- res$wasserstein[1]
        }
        # Normalize by IQR of original
        iqrs <- sapply(num_vars, function(v) IQR(original[[v]], na.rm = TRUE))
        iqrs[iqrs == 0] <- 1
        w_normalized <- w_vals / iqrs
        # Transform to utility (higher = better)
        utility_results[i, "wasserstein"] <- 1 - min(1, mean(w_normalized, na.rm = TRUE))
      }, error = function(e) {
        warning(paste("wasserstein failed for", sdg_name, ":", e$message))
        utility_results[i, "wasserstein"] <<- NA
      })
    }

    # Hellinger
    if ("hellinger" %in% utility_measures && length(cat_vars) > 0) {
      tryCatch({
        res <- hellinger(original, Y, vars = cat_vars, na.rm = na.rm)
        # Transform to utility (higher = better)
        utility_results[i, "hellinger"] <- res$utility_score
      }, error = function(e) {
        warning(paste("hellinger failed for", sdg_name, ":", e$message))
        utility_results[i, "hellinger"] <<- NA
      })
    }

    # Energy distance
    if ("energy" %in% utility_measures && length(num_vars) > 0) {
      tryCatch({
        res <- energy_distance(original, Y, vars = num_vars, standardize = TRUE,
                               n_sample = 500, seed = seed, na.rm = na.rm)
        utility_results[i, "energy"] <- res$utility_score
      }, error = function(e) {
        warning(paste("energy failed for", sdg_name, ":", e$message))
        utility_results[i, "energy"] <<- NA
      })
    }

    # CI proximity
    if ("ci_proximity" %in% utility_measures && length(num_vars) > 0) {
      tryCatch({
        res <- ci_proximity(original, Y, vars = num_vars, na.rm = na.rm)
        utility_results[i, "ci_proximity"] <- res$proximity_mean
      }, error = function(e) {
        warning(paste("ci_proximity failed for", sdg_name, ":", e$message))
        utility_results[i, "ci_proximity"] <<- NA
      })
    }

    # MMD
    if ("mmd" %in% utility_measures && length(num_vars) > 0) {
      tryCatch({
        res <- mmd(original, Y, vars = num_vars, standardize = TRUE,
                   seed = seed, na.rm = na.rm)
        utility_results[i, "mmd"] <- res$utility_score
      }, error = function(e) {
        warning(paste("mmd failed for", sdg_name, ":", e$message))
        utility_results[i, "mmd"] <<- NA
      })
    }

    # TSTR
    if ("tstr" %in% utility_measures && !is.null(target_var)) {
      tryCatch({
        res <- tstr(original, Y, target_var = target_var,
                    seed = seed, na.rm = na.rm)
        utility_results[i, "tstr"] <- res$utility_score
      }, error = function(e) {
        warning(paste("tstr failed for", sdg_name, ":", e$message))
        utility_results[i, "tstr"] <<- NA
      })
    }

    # Copula fidelity
    if ("copula" %in% utility_measures && length(num_vars) >= 2) {
      tryCatch({
        res <- copula_fidelity(original, Y, vars = num_vars, na.rm = na.rm)
        utility_results[i, "copula"] <- res$utility_score
      }, error = function(e) {
        warning(paste("copula_fidelity failed for", sdg_name, ":", e$message))
        utility_results[i, "copula"] <<- NA
      })
    }

    # Tail fidelity
    if ("tail" %in% utility_measures && length(num_vars) > 0) {
      tryCatch({
        res <- tail_fidelity(original, Y, vars = num_vars, na.rm = na.rm)
        utility_results[i, "tail"] <- res$utility_score
      }, error = function(e) {
        warning(paste("tail_fidelity failed for", sdg_name, ":", e$message))
        utility_results[i, "tail"] <<- NA
      })
    }

    # Contingency fidelity
    if ("contingency" %in% utility_measures && length(cat_vars) >= 2) {
      tryCatch({
        res <- contingency_fidelity(original, Y, vars = cat_vars)
        utility_results[i, "contingency"] <- res$utility_score
      }, error = function(e) {
        warning(paste("contingency_fidelity failed for", sdg_name, ":", e$message))
        utility_results[i, "contingency"] <<- NA
      })
    }
  }

  # Get actual computed measures
  risk_cols <- setdiff(names(risk_results), "sdg")
  utility_cols <- setdiff(names(utility_results), "sdg")

  # Normalize measures to [0,1] if requested
  if (normalize) {
    normalized <- data.frame(sdg = sdg_names, stringsAsFactors = FALSE)

    # Min-max normalize risk measures
    for (col in risk_cols) {
      vals <- risk_results[[col]]
      if (all(is.na(vals))) {
        normalized[[col]] <- NA
      } else {
        min_v <- min(vals, na.rm = TRUE)
        max_v <- max(vals, na.rm = TRUE)
        if (max_v > min_v) {
          normalized[[col]] <- (vals - min_v) / (max_v - min_v)
        } else {
          normalized[[col]] <- 0.5  # All same value
        }
      }
    }

    # Min-max normalize utility measures
    for (col in utility_cols) {
      vals <- utility_results[[col]]
      if (all(is.na(vals))) {
        normalized[[col]] <- NA
      } else {
        min_v <- min(vals, na.rm = TRUE)
        max_v <- max(vals, na.rm = TRUE)
        if (max_v > min_v) {
          normalized[[col]] <- (vals - min_v) / (max_v - min_v)
        } else {
          normalized[[col]] <- 0.5
        }
      }
    }
  } else {
    normalized <- NULL
  }

  # Compute composite scores
  composites <- data.frame(sdg = sdg_names, stringsAsFactors = FALSE)

  # Use normalized values if available, otherwise raw
  if (!is.null(normalized)) {
    risk_data <- normalized[, risk_cols, drop = FALSE]
    utility_data <- normalized[, utility_cols, drop = FALSE]
  } else {
    risk_data <- risk_results[, risk_cols, drop = FALSE]
    utility_data <- utility_results[, utility_cols, drop = FALSE]
  }

  composites$risk_mean <- rowMeans(risk_data, na.rm = TRUE)
  composites$utility_mean <- rowMeans(utility_data, na.rm = TRUE)
  composites$risk_sd <- apply(risk_data, 1, sd, na.rm = TRUE)
  composites$utility_sd <- apply(utility_data, 1, sd, na.rm = TRUE)

  # Identify Pareto-optimal SDGs
  # SDG i dominates SDG j if: risk_i <= risk_j AND utility_i >= utility_j with at least one strict
  is_pareto <- identify_pareto(composites$risk_mean, composites$utility_mean)
  pareto_sdgs <- sdg_names[is_pareto]

  # Build result object
  result <- list(
    risk = risk_results,
    utility = utility_results,
    normalized = normalized,
    composites = composites,
    pareto = is_pareto,
    pareto_sdgs = pareto_sdgs,
    n_sdgs = n_sdgs,
    sdg_names = sdg_names,
    risk_measures = risk_cols,
    utility_measures = utility_cols,
    metadata = list(
      key_vars = key_vars,
      target_var = target_var,
      vars = vars,
      cat_vars = cat_vars,
      num_vars = num_vars,
      holdout_fraction = holdout_fraction,
      normalize = normalize,
      n_original = nrow(original),
      seed = seed
    )
  )

  class(result) <- "rumap"
  return(result)
}


#' Identify Pareto-Optimal Solutions
#'
#' @param risk Numeric vector of risk scores (lower = better)
#' @param utility Numeric vector of utility scores (higher = better)
#' @return Logical vector indicating which observations are Pareto-optimal
#' @keywords internal
identify_pareto <- function(risk, utility) {
  n <- length(risk)
  is_pareto <- rep(TRUE, n)

  for (i in 1:n) {
    if (is.na(risk[i]) || is.na(utility[i])) {
      is_pareto[i] <- FALSE
      next
    }

    for (j in 1:n) {
      if (i == j || is.na(risk[j]) || is.na(utility[j])) next

      # j dominates i if: risk_j <= risk_i AND utility_j >= utility_i
      # with at least one strict inequality
      if (risk[j] <= risk[i] && utility[j] >= utility[i]) {
        if (risk[j] < risk[i] || utility[j] > utility[i]) {
          is_pareto[i] <- FALSE
          break
        }
      }
    }
  }

  return(is_pareto)
}


#' Print method for rumap objects
#'
#' @param x an object of class "rumap"
#' @param ... additional arguments (ignored)
#' @export
print.rumap <- function(x, ...) {
  cat("Multivariate Risk-Utility Map (RU-Map)\n")
  cat("======================================\n\n")

  cat("Configuration:\n")
  cat("  Synthetic datasets (SDGs):", x$n_sdgs, "\n")
  cat("  Risk measures:   ", paste(x$risk_measures, collapse = ", "), "\n")
  cat("  Utility measures:", paste(x$utility_measures, collapse = ", "), "\n")
  cat("  Original data size:", x$metadata$n_original, "\n\n")

  cat("Composite Scores:\n")
  df <- x$composites[, c("sdg", "risk_mean", "utility_mean")]
  df$pareto <- ifelse(x$pareto, "*", "")
  names(df) <- c("SDG", "Risk", "Utility", "Pareto")
  df$Risk <- sprintf("%.4f", df$Risk)
  df$Utility <- sprintf("%.4f", df$Utility)
  print(df, row.names = FALSE)

  cat("\n* = Pareto-optimal\n")
  cat("\nPareto-optimal SDGs:", paste(x$pareto_sdgs, collapse = ", "), "\n")

  invisible(x)
}


#' Summary method for rumap objects
#'
#' @param object an object of class "rumap"
#' @param ... additional arguments (ignored)
#' @export
summary.rumap <- function(object, ...) {
  # Compute internal consistency metrics
  risk_data <- object$normalized[, object$risk_measures, drop = FALSE]
  utility_data <- object$normalized[, object$utility_measures, drop = FALSE]

  # Cronbach's alpha (simplified)
  cronbach_alpha <- function(data) {
    data <- data[complete.cases(data), , drop = FALSE]
    if (nrow(data) < 2 || ncol(data) < 2) return(NA)
    k <- ncol(data)
    item_vars <- apply(data, 2, var, na.rm = TRUE)
    total_var <- var(rowSums(data, na.rm = TRUE), na.rm = TRUE)
    if (total_var == 0) return(NA)
    (k / (k - 1)) * (1 - sum(item_vars) / total_var)
  }

  summ <- list(
    n_sdgs = object$n_sdgs,
    sdg_names = object$sdg_names,
    risk_measures = object$risk_measures,
    utility_measures = object$utility_measures,
    risk_raw = object$risk,
    utility_raw = object$utility,
    normalized = object$normalized,
    composites = object$composites,
    pareto = object$pareto,
    pareto_sdgs = object$pareto_sdgs,
    n_pareto = sum(object$pareto),
    alpha_risk = cronbach_alpha(risk_data),
    alpha_utility = cronbach_alpha(utility_data),
    metadata = object$metadata
  )

  class(summ) <- "summary.rumap"
  return(summ)
}


#' Print method for summary.rumap objects
#'
#' @param x an object of class "summary.rumap"
#' @param ... additional arguments (ignored)
#' @export
print.summary.rumap <- function(x, ...) {
  cat("Summary: Multivariate Risk-Utility Map\n")
  cat("======================================\n\n")

  cat("SDGs evaluated:", x$n_sdgs, "\n")
  cat("Pareto-optimal:", x$n_pareto, "(", paste(x$pareto_sdgs, collapse = ", "), ")\n\n")

  cat("Internal Consistency:\n")
  cat("  Cronbach's alpha (risk):   ", sprintf("%.3f", x$alpha_risk), "\n")
  cat("  Cronbach's alpha (utility):", sprintf("%.3f", x$alpha_utility), "\n\n")

  cat("Risk Measures (raw):\n")
  print(x$risk_raw, row.names = FALSE)

  cat("\nUtility Measures (raw):\n")
  print(x$utility_raw, row.names = FALSE)

  cat("\nNormalized Scores (0-1):\n")
  if (!is.null(x$normalized)) {
    print(x$normalized, row.names = FALSE)
  } else {
    cat("  (normalization disabled)\n")
  }

  invisible(x)
}
