#' RAPID: Risk of Attribute Prediction-Induced Disclosure
#'
#' Assesses inferential disclosure risk for a sensitive variable in synthetic data
#' using the RAPID metric. This function trains a predictive model on synthetic data
#' using quasi-identifiers, then evaluates attribute inference risk by scoring the
#' model on the original data.
#'
#' @param X data frame of original data, or a \code{\link{synth_pair}} object
#' @param Y data frame of synthetic data (not needed if X is a synth_pair)
#' @param key_vars character vector of quasi-identifier variable names
#' @param target_var character, name of the sensitive target variable
#' @param model_type character, model type for inference: "rf" (random forest, default),
#'   "lm" (linear model), "cart" (decision tree), "gbm" (gradient boosting), or
#'   "logit" (logistic regression, binary only)
#' @param num_epsilon numeric threshold for continuous attributes. For percentage-based
#'   metrics (default), specify as percentage (e.g., 5 for 5\%). For absolute error,
#'   specify in units of the sensitive attribute.
#' @param num_epsilon_type character, threshold type: "percentage" (default) or "absolute"
#' @param num_error_metric character, error metric for continuous variables:
#'   "symmetric" (default, recommended), "stabilised_relative", or "absolute"
#' @param num_delta numeric, smoothing constant for percentage-based metrics (default: 0.01)
#' @param cat_tau numeric threshold for categorical risk. Interpretation depends on method:
#'   \itemize{
#'     \item For RCS_conditional: ratio threshold (typically 1.25)
#'     \item For RCS_marginal: normalized gain threshold (typically 0.3)
#'     \item For NCE: risk score threshold (typically 0.5-0.7)
#'   }
#' @param cat_eval_method character, method for categorical evaluation:
#'   "RCS_marginal" (default, recommended), "RCS_conditional", or "NCE"
#' @param na_strategy character, strategy for NA values in sensitive attribute:
#'   "constant" (default), "drop", or "median"
#' @param na_constant_value numeric, value for constant NA imputation (default: 0)
#' @param return_all_records logical, if TRUE (default) return all records with risk status;
#'   if FALSE return only at-risk records
#' @param seed integer, random seed for reproducibility (default: NULL)
#' @param verbose logical, print diagnostic messages (default: FALSE)
#' @param ... additional arguments passed to model fitting functions
#'
#' @return An object of class "rapid" containing:
#' \itemize{
#'   \item rapid: the confidence rate (proportion of records at risk)
#'   \item n_at_risk: number of records at risk
#'   \item pct_at_risk: percentage of records at risk
#'   \item method: evaluation method used (categorical) or "numeric"
#'   \item threshold: threshold used (cat_tau or num_epsilon)
#'   \item model_type: model type used
#'   \item model_metrics: model performance (accuracy for categorical, MAE/RMSE for numeric)
#'   \item records: data frame of at-risk records (or all records if return_all_records=TRUE)
#'   \item key_vars, target_var: input parameters
#'   \item call: the function call
#' }
#'
#' @details
#' RAPID (Risk of Attribute Prediction-Induced Disclosure) measures disclosure risk
#' by training a predictive model on synthetic data and evaluating how well it can
#' predict the sensitive attribute in the original data. High prediction accuracy
#' indicates potential disclosure risk.
#'
#' Unlike CAP-based methods (DCAP, TCAP) which use exact or fuzzy matching, RAPID
#' uses machine learning models to capture more complex relationships between
#' quasi-identifiers and sensitive attributes.
#'
#' \strong{For categorical sensitive variables:}
#' Three evaluation methods are available:
#' \itemize{
#'   \item \code{RCS_marginal} (recommended): Measures if the attribute can be inferred
#'     better than the marginal baseline rate. Uses normalized gain.
#'   \item \code{RCS_conditional}: Measures if an observation is an outlier within its
#'     class using class-conditional baseline.
#'   \item \code{NCE}: Normalized Cross-Entropy, measures information leakage.
#' }
#'
#' \strong{For continuous sensitive variables:}
#' Risk is measured by the proportion of records for which prediction errors fall
#' below a threshold. Three error metrics are available:
#' \itemize{
#'   \item \code{symmetric} (recommended): Symmetric percentage error, treats predicted
#'     and true values equally.
#'   \item \code{stabilised_relative}: Stabilised relative error.
#'   \item \code{absolute}: Raw absolute error (use with num_epsilon_type = "absolute").
#' }
#'
#' @section Threshold Selection Guidelines:
#' \strong{Categorical sensitive variables:}
#' \itemize{
#'   \item RCS_marginal: \code{cat_tau = 0.3} (recommended)
#'   \item RCS_conditional: \code{cat_tau = 1.25}
#'   \item NCE: \code{cat_tau = 0.5-0.7}
#' }
#' \strong{Continuous sensitive variables:}
#' \itemize{
#'   \item symmetric: \code{num_epsilon = 5-10} (percentage)
#'   \item absolute: \code{num_epsilon} depends on domain/scale
#' }
#'
#' @section Comparison with CAP Methods:
#' \itemize{
#'   \item \strong{RAPID} captures complex non-linear relationships via ML models
#'   \item \strong{DCAP/TCAP} use exact/fuzzy matching on key combinations
#'   \item Use RAPID when relationships are complex; use CAP when matching is intuitive
#'   \item Both approaches are complementary and should ideally be used together
#' }
#'
#' @references
#' Templ, M. (2024). Beyond the Trade-off Curve: Multivariate Risk-Utility
#' Visualization for Synthetic Data. \emph{Journal of Official Statistics}.
#'
#' @seealso \code{\link{dcap}}, \code{\link{tcap}} for attribution-based metrics,
#'   \code{\link{disclosure_report}} for comprehensive risk assessment
#'
#' @author Oscar Thees, Matthias Templ
#' @export
#' @importFrom stats predict lm as.formula glm binomial model.matrix complete.cases median
#' @examples
#' # Create example data
#' set.seed(123)
#' n <- 200
#' original <- data.frame(
#'   age = sample(20:70, n, replace = TRUE),
#'   gender = factor(sample(c("M", "F"), n, replace = TRUE)),
#'   region = factor(sample(c("N", "S", "E", "W"), n, replace = TRUE)),
#'   income = round(rnorm(n, 50000, 10000))
#' )
#'
#' # Synthetic version with some noise
#' synthetic <- data.frame(
#'   age = original$age + sample(-3:3, n, replace = TRUE),
#'   gender = factor(sample(c("M", "F"), n, replace = TRUE)),
#'   region = factor(sample(c("N", "S", "E", "W"), n, replace = TRUE)),
#'   income = round(original$income * runif(n, 0.9, 1.1))
#' )
#'
#' # RAPID for continuous sensitive variable
#' \donttest{
#' result <- rapid(
#'   X = original,
#'   Y = synthetic,
#'   key_vars = c("age", "gender", "region"),
#'   target_var = "income",
#'   model_type = "rf",
#'   num_epsilon = 10,
#'   verbose = TRUE
#' )
#' print(result)
#' summary(result)
#' }
#'
#' # RAPID for categorical sensitive variable
#' \donttest{
#' original$health <- factor(sample(c("Good", "Fair", "Poor"), n,
#'                                   replace = TRUE, prob = c(0.6, 0.3, 0.1)))
#' synthetic$health <- factor(sample(c("Good", "Fair", "Poor"), n,
#'                                    replace = TRUE, prob = c(0.55, 0.35, 0.1)))
#'
#' result_cat <- rapid(
#'   X = original,
#'   Y = synthetic,
#'   key_vars = c("age", "gender", "region"),
#'   target_var = "health",
#'   model_type = "rf",
#'   cat_tau = 0.3,
#'   cat_eval_method = "RCS_marginal"
#' )
#' print(result_cat)
#' }
#'
#' # ================================================================
#' # Real data example using SD2011 from synthpop package
#' # ================================================================
#' \donttest{
#' if (requireNamespace("synthpop", quietly = TRUE) &&
#'     requireNamespace("ranger", quietly = TRUE)) {
#'
#'   # Load SD2011 Polish Social Diagnosis survey data
#'   data(SD2011, package = "synthpop")
#'
#'   # Select relevant variables and remove NAs for cleaner example
#'   vars <- c("sex", "age", "edu", "marital", "region", "income", "smoke")
#'   original_sd <- SD2011[, vars]
#'   original_sd <- original_sd[complete.cases(original_sd), ]
#'
#'   # Generate synthetic data using synthpop
#'   synth_obj <- synthpop::syn(original_sd, seed = 123, print.flag = FALSE)
#'   synthetic_sd <- synth_obj$syn
#'
#'   # ------------------------------------------------------------------
#'   # Example 1: Continuous sensitive variable (income)
#'   # ------------------------------------------------------------------
#'   # Assess risk of inferring income from demographic quasi-identifiers
#'
#'   rapid_income <- rapid(
#'     X = original_sd,
#'     Y = synthetic_sd,
#'     key_vars = c("sex", "age", "edu", "marital", "region"),
#'     target_var = "income",
#'     model_type = "rf",
#'     num_epsilon = 10,           # 10% error threshold
#'     num_epsilon_type = "percentage",
#'     num_error_metric = "symmetric",
#'     verbose = TRUE
#'   )
#'
#'   print(rapid_income)
#'   summary(rapid_income)
#'   plot(rapid_income, which = 1:2)
#'
#'   # Interpretation:
#'   # - RAPID score close to 0: low risk (model cannot predict income well)
#'   # - RAPID score > 0.15: elevated risk (sensitive attribute predictable)
#'
#'   # ------------------------------------------------------------------
#'   # Example 2: Categorical sensitive variable (smoking status)
#'   # ------------------------------------------------------------------
#'   # Assess risk of inferring smoking behavior from demographics
#'
#'   rapid_smoke <- rapid(
#'     X = original_sd,
#'     Y = synthetic_sd,
#'     key_vars = c("sex", "age", "edu", "marital", "income"),
#'     target_var = "smoke",
#'     model_type = "rf",
#'     cat_tau = 0.3,              # 30% normalized gain threshold
#'     cat_eval_method = "RCS_marginal",
#'     verbose = TRUE
#'   )
#'
#'   print(rapid_smoke)
#'   summary(rapid_smoke)
#'   plot(rapid_smoke)
#'
#'   # Compare with DCAP for the same data
#'   dcap_smoke <- dcap(
#'     X = original_sd,
#'     Y = synthetic_sd,
#'     key_vars = c("sex", "age", "edu", "marital"),
#'     target_var = "smoke"
#'   )
#'   print(dcap_smoke)
#'
#'   # RAPID may detect risks that DCAP misses when relationships are non-linear
#' }
#' }
rapid <- function(X, ...) {
  UseMethod("rapid")
}

#' @rdname rapid
#' @export
rapid.synth_pair <- function(X, ...) {
  if (is.null(X$key_vars)) {
    stop("synth_pair must have 'key_vars' set for rapid()")
  }
  if (is.null(X$target_var)) {
    stop("synth_pair must have 'target_var' set for rapid()")
  }

  rapid.default(
    X = X$original,
    Y = X$synthetic,
    key_vars = X$key_vars,
    target_var = X$target_var,
    ...
  )
}

#' @rdname rapid
#' @export
rapid.default <- function(X, Y,
                          key_vars,
                          target_var,
                          model_type = c("rf", "lm", "cart", "gbm", "logit"),
                          # Numeric-specific
                          num_epsilon = 10,
                          num_epsilon_type = c("percentage", "absolute"),
                          num_error_metric = c("symmetric", "stabilised_relative", "absolute"),
                          num_delta = 0.01,
                          # Categorical-specific
                          cat_tau = 0.3,
                          cat_eval_method = c("RCS_marginal", "RCS_conditional", "NCE"),
                          # NA handling
                          na_strategy = c("constant", "drop", "median"),
                          na_constant_value = 0,
                          # Output options
                          return_all_records = TRUE,
                          seed = NULL,
                          verbose = FALSE,
                          ...) {

  call <- match.call()

  # Argument matching

  model_type <- match.arg(model_type)
  num_epsilon_type <- match.arg(num_epsilon_type)
  num_error_metric <- match.arg(num_error_metric)
  na_strategy <- match.arg(na_strategy)
  cat_eval_method <- match.arg(cat_eval_method)

  # Input validation
  if (!is.data.frame(X)) stop("X must be a data frame.")
  if (!is.data.frame(Y)) stop("Y must be a data frame.")

  all_vars <- c(key_vars, target_var)
  missing_X <- setdiff(all_vars, names(X))
  missing_Y <- setdiff(all_vars, names(Y))
  if (length(missing_X) > 0) {
    stop(paste("Variables missing in X:", paste(missing_X, collapse = ", ")))
  }
  if (length(missing_Y) > 0) {
    stop(paste("Variables missing in Y:", paste(missing_Y, collapse = ", ")))
  }

  # Check model dependencies
  .check_rapid_deps(model_type)

  # Set seed
  if (!is.null(seed)) set.seed(seed)

  # Determine if target is categorical
  target_vec <- X[[target_var]]
  is_categorical <- is.factor(target_vec) || is.character(target_vec)

  # Model compatibility check
  if (model_type == "lm" && is_categorical) {
    stop("Model 'lm' is not suitable for categorical sensitive variables. Use 'rf', 'cart', 'gbm', or 'logit'.")
  }
  if (model_type == "logit" && !is_categorical) {
    stop("Model 'logit' requires a categorical sensitive variable.")
  }
  if (model_type == "logit" && is_categorical && nlevels(as.factor(target_vec)) != 2) {
    stop("Model 'logit' only supports binary classification (2 levels).")
  }

  if (verbose) message("RAPID: Running disclosure risk assessment...")

  # Subset to needed variables
  X_sub <- X[, c(target_var, key_vars), drop = FALSE]
  Y_sub <- Y[, c(target_var, key_vars), drop = FALSE]

  # Handle NA in sensitive variable
  na_result <- .handle_sensitive_na(X_sub, Y_sub, target_var, na_strategy, na_constant_value)
  X_sub <- na_result$original
  Y_sub <- na_result$synthetic

  # Handle NA in quasi-identifiers
  X_sub <- .handle_qi_na(X_sub, key_vars)
  Y_sub <- .handle_qi_na(Y_sub, key_vars)

  # Ensure factor levels are consistent
  Y_sub <- .ensure_levels(Y_sub, X_sub, target_var, key_vars, verbose)

  # Build formula
  formula <- as.formula(paste(target_var, "~", paste(key_vars, collapse = "+")))

  # Fit model on synthetic data
  if (verbose) message("  Fitting ", model_type, " model on synthetic data...")
  fit <- .fit_rapid_model(model_type, formula, Y_sub, X_sub, target_var, ...)

  # Predict on original data
  if (verbose) message("  Predicting on original data...")
  A <- X_sub[[target_var]]
  B <- .predict_rapid(model_type, fit, X_sub, target_var)

  # Compute model performance metrics
  model_metrics <- .compute_rapid_metrics(A, B)

  # Evaluate disclosure risk
  if (verbose) message("  Evaluating disclosure risk...")
  if (is_categorical) {
    eval_result <- .evaluate_categorical_rapid(
      A = A, B = B,
      original_data = X_sub,
      cat_tau = cat_tau,
      cat_eval_method = cat_eval_method,
      target_var = target_var,
      return_all_records = return_all_records
    )
    threshold_used <- cat_tau
    method_used <- cat_eval_method
  } else {
    eval_result <- .evaluate_numeric_rapid(
      A = A, B = B,
      original_data = X_sub,
      num_epsilon = num_epsilon,
      num_epsilon_type = num_epsilon_type,
      num_error_metric = num_error_metric,
      num_delta = num_delta,
      return_all_records = return_all_records
    )
    threshold_used <- num_epsilon
    method_used <- paste0(num_error_metric, "_", num_epsilon_type)
  }

  # Build result object
  result <- list(
    rapid = eval_result$confidence_rate,
    n_at_risk = eval_result$n_at_risk,
    pct_at_risk = eval_result$percentage,
    n_total = nrow(X_sub),
    method = method_used,
    threshold = threshold_used,
    model_type = model_type,
    model_metrics = model_metrics,
    records = eval_result$rows_risk_df,
    key_vars = key_vars,
    target_var = target_var,
    is_categorical = is_categorical,
    error_range = eval_result$error_range,  # NULL for categorical
    call = call
  )

  class(result) <- "rapid"
  return(result)
}


# =============================================================================
# Internal helper functions
# =============================================================================

#' Check required packages for RAPID models
#' @noRd
.check_rapid_deps <- function(model_type) {
  if (model_type == "rf" && !requireNamespace("ranger", quietly = TRUE)) {
    stop("Package 'ranger' required for model_type='rf'. Install with: install.packages('ranger')")
  }
  if (model_type == "cart" && !requireNamespace("rpart", quietly = TRUE)) {
    stop("Package 'rpart' required for model_type='cart'. Install with: install.packages('rpart')")
  }
  if (model_type == "gbm" && !requireNamespace("xgboost", quietly = TRUE)) {
    stop("Package 'xgboost' required for model_type='gbm'. Install with: install.packages('xgboost')")
  }
}

#' Handle NA in sensitive attribute
#' @noRd
.handle_sensitive_na <- function(original, synthetic, target_var, strategy, constant_value) {
  orig_target <- original[[target_var]]
  synth_target <- synthetic[[target_var]]

  has_na_orig <- anyNA(orig_target)
  has_na_synth <- anyNA(synth_target)

  if (!has_na_orig && !has_na_synth) {
    return(list(original = original, synthetic = synthetic))
  }

  if (is.factor(orig_target)) {
    # For factors: add "missing" level
    if (has_na_orig) {
      original[[target_var]] <- addNA(original[[target_var]])
      levels(original[[target_var]])[is.na(levels(original[[target_var]]))] <- "missing"
    }
    if (has_na_synth) {
      synthetic[[target_var]] <- addNA(synthetic[[target_var]])
      levels(synthetic[[target_var]])[is.na(levels(synthetic[[target_var]]))] <- "missing"
    }
  } else if (is.numeric(orig_target)) {
    # For numeric: apply strategy
    if (strategy == "drop") {
      keep_orig <- !is.na(original[[target_var]])
      keep_synth <- !is.na(synthetic[[target_var]])
      original <- original[keep_orig, , drop = FALSE]
      synthetic <- synthetic[keep_synth, , drop = FALSE]
    } else if (strategy == "constant") {
      original[[target_var]][is.na(original[[target_var]])] <- constant_value
      synthetic[[target_var]][is.na(synthetic[[target_var]])] <- constant_value
    } else if (strategy == "median") {
      orig_med <- median(original[[target_var]], na.rm = TRUE)
      synth_med <- median(synthetic[[target_var]], na.rm = TRUE)
      original[[target_var]][is.na(original[[target_var]])] <- orig_med
      synthetic[[target_var]][is.na(synthetic[[target_var]])] <- synth_med
    }
  }

  list(original = original, synthetic = synthetic)
}

#' Handle NA in quasi-identifiers
#' @noRd
.handle_qi_na <- function(data, qi_vars) {
  for (var in qi_vars) {
    if (!var %in% names(data)) next

    if (is.factor(data[[var]]) && anyNA(data[[var]])) {
      data[[var]] <- addNA(data[[var]])
      levels(data[[var]])[is.na(levels(data[[var]]))] <- "missing"
    } else if (is.numeric(data[[var]]) && anyNA(data[[var]])) {
      stop(sprintf("Variable '%s' is numeric and contains NA. Please handle missing values before using RAPID.", var))
    }
  }
  data
}

#' Ensure factor levels are consistent between datasets
#' @noRd
.ensure_levels <- function(synthetic, original, target_var, key_vars, verbose = FALSE) {
  # Convert character to factor if needed
  if (is.character(synthetic[[target_var]])) {
    synthetic[[target_var]] <- factor(synthetic[[target_var]])
  }
  if (is.character(original[[target_var]])) {
    original[[target_var]] <- factor(original[[target_var]])
  }

  # For categorical target, ensure all original levels exist in synthetic

  if (is.factor(original[[target_var]])) {
    orig_levels <- levels(original[[target_var]])
    synth_levels <- levels(synthetic[[target_var]])

    missing_levels <- setdiff(orig_levels, synth_levels)
    if (length(missing_levels) > 0) {
      # Add missing levels to synthetic
      levels(synthetic[[target_var]]) <- c(synth_levels, missing_levels)
    }

    # Also ensure synthetic levels match original order
    synthetic[[target_var]] <- factor(synthetic[[target_var]], levels = orig_levels)
  }

  # Do the same for categorical key variables
  for (var in key_vars) {
    if (is.factor(original[[var]])) {
      orig_levels <- levels(original[[var]])
      if (is.factor(synthetic[[var]])) {
        synth_levels <- levels(synthetic[[var]])
        all_levels <- union(orig_levels, synth_levels)
        synthetic[[var]] <- factor(synthetic[[var]], levels = all_levels)
      }
    }
  }

  synthetic
}

#' Fit predictive model for RAPID
#' @noRd
.fit_rapid_model <- function(model_type, formula, synthetic, original, target_var, ...) {

  # Auto-convert character to factor
  if (is.character(synthetic[[target_var]])) {
    synthetic[[target_var]] <- factor(synthetic[[target_var]])
  }

  switch(model_type,
    lm = {
      stats::lm(formula, data = synthetic, ...)
    },

    rf = {
      use_prob <- is.factor(synthetic[[target_var]])
      ranger::ranger(
        formula = formula,
        data = synthetic,
        probability = use_prob,
        ...
      )
    },

    cart = {
      method_type <- if (is.factor(synthetic[[target_var]])) "class" else "anova"
      rpart::rpart(
        formula = formula,
        data = synthetic,
        method = method_type,
        ...
      )
    },

    gbm = {
      y_raw <- synthetic[[target_var]]
      X <- model.matrix(formula, data = synthetic)[, -1, drop = FALSE]

      if (!is.factor(y_raw)) {
        # Regression
        dtrain <- xgboost::xgb.DMatrix(data = X, label = y_raw)
        params <- list(
          objective = "reg:squarederror",
          eval_metric = "rmse",
          verbosity = 0
        )
        model <- xgboost::xgb.train(
          params = params,
          data = dtrain,
          nrounds = 100,
          verbose = 0
        )
      } else {
        # Classification
        y_factor <- factor(y_raw)
        y_encoded <- as.integer(y_factor) - 1
        K <- length(levels(y_factor))
        dtrain <- xgboost::xgb.DMatrix(data = X, label = y_encoded)

        if (K == 2) {
          params <- list(
            objective = "binary:logistic",
            eval_metric = "logloss",
            verbosity = 0
          )
        } else {
          params <- list(
            objective = "multi:softprob",
            num_class = K,
            eval_metric = "mlogloss",
            verbosity = 0
          )
        }
        model <- xgboost::xgb.train(
          params = params,
          data = dtrain,
          nrounds = 100,
          verbose = 0
        )
      }
      # Store formula for prediction
      model$.__x_formula__ <- formula
      model$.__levels__ <- if (is.factor(y_raw)) levels(y_raw) else NULL
      model
    },

    logit = {
      if (!is.factor(synthetic[[target_var]])) {
        stop("Logit requires a categorical sensitive variable.")
      }
      glm(formula = formula, data = synthetic, family = binomial())
    },

    stop("Unsupported model type. Use 'lm', 'rf', 'cart', 'gbm', or 'logit'.")
  )
}

#' Predict sensitive attribute from fitted model
#' @noRd
.predict_rapid <- function(model_type, fit, original, target_var) {

  switch(model_type,
    lm = {
      stats::predict(fit, newdata = original)
    },

    rf = {
      predict(fit, data = original, type = "response")$predictions
    },

    cart = {
      is_classification <- is.factor(original[[target_var]])
      if (is_classification) {
        probs <- predict(fit, newdata = original, type = "prob")
        probs_df <- as.data.frame(probs)
        levels_response <- levels(original[[target_var]])
        probs_df[, levels_response, drop = FALSE]
      } else {
        predict(fit, newdata = original)
      }
    },

    gbm = {
      X_test <- model.matrix(fit$.__x_formula__, data = original)[, -1, drop = FALSE]
      dtest <- xgboost::xgb.DMatrix(X_test)
      preds <- predict(fit, newdata = dtest)

      y_true <- original[[target_var]]

      if (!is.factor(y_true)) {
        return(preds)
      }

      levels_response <- fit$.__levels__
      K <- length(levels_response)

      if (K == 2) {
        probs_df <- data.frame(class1 = 1 - preds, class2 = preds)
        colnames(probs_df) <- levels_response
        return(probs_df)
      }

      probs_mat <- matrix(preds, ncol = K, byrow = TRUE)
      probs_df <- as.data.frame(probs_mat)
      colnames(probs_df) <- levels_response
      probs_df
    },

    logit = {
      probs <- predict(fit, newdata = original, type = "response")
      levels_response <- levels(original[[target_var]])
      probs_df <- data.frame(class1 = 1 - probs, class2 = probs)
      colnames(probs_df) <- levels_response
      probs_df
    },

    stop("Unsupported model type.")
  )
}

#' Compute model performance metrics
#' @noRd
.compute_rapid_metrics <- function(A, B) {
  if (is.factor(A)) {
    # Classification
    if (is.matrix(B) || is.data.frame(B)) {
      max_idx <- apply(as.matrix(B), 1, which.max)
      pred_class <- colnames(B)[max_idx]
      pred_class <- factor(pred_class, levels = levels(A))
    } else {
      pred_class <- B
    }
    accuracy <- mean(pred_class == A, na.rm = TRUE)
    return(list(accuracy = accuracy))
  } else {
    # Regression
    mae <- mean(abs(A - B), na.rm = TRUE)
    rmse <- sqrt(mean((A - B)^2, na.rm = TRUE))
    rmae <- mae / mean(abs(A), na.rm = TRUE)
    rrmse <- rmse / sd(A, na.rm = TRUE)
    return(list(mae = mae, rmse = rmse, rmae = rmae, rrmse = rrmse))
  }
}

#' Evaluate numeric attribute inference risk
#' @noRd
.evaluate_numeric_rapid <- function(A, B, original_data,
                                    num_epsilon, num_epsilon_type,
                                    num_error_metric, num_delta,
                                    return_all_records) {

  if (!is.numeric(num_epsilon) || num_epsilon <= 0) {
    stop("`num_epsilon` must be a positive number.")
  }
  if (!is.numeric(num_delta) || num_delta <= 0) {
    stop("`num_delta` must be a positive number.")
  }

  if (num_epsilon_type == "percentage") {
    error_values <- switch(num_error_metric,
      symmetric = {
        err <- 2 * abs(A - B) / (abs(A) + abs(B) + 2 * num_delta)
        err * 100
      },
      stabilised_relative = {
        err <- abs(A - B) / (abs(A) + num_delta)
        err * 100
      },
      absolute = {
        stop("Percentage-based threshold is not meaningful for absolute error metric. Use num_epsilon_type = 'absolute'.")
      }
    )
    metric_name <- paste0(num_error_metric, "_error_pct")
  } else {
    error_values <- abs(A - B)
    metric_name <- "absolute_error"
  }

  # For numeric: at_risk when prediction error is BELOW threshold
  # Low error means the model can accurately predict the sensitive attribute = disclosure risk
  at_risk <- error_values < num_epsilon

  result_df <- data.frame(
    original_data,
    original_value = A,
    predicted_value = B,
    error_metric = error_values,
    threshold = num_epsilon,
    threshold_type = num_epsilon_type,
    at_risk = at_risk,
    stringsAsFactors = FALSE
  )
  colnames(result_df)[colnames(result_df) == "error_metric"] <- metric_name

  # Store full error range before filtering (for plotting)
  error_range <- range(error_values, na.rm = TRUE)

  if (!return_all_records) {
    result_df <- result_df[result_df$at_risk, , drop = FALSE]
  }

  list(
    confidence_rate = sum(at_risk) / length(at_risk),
    n_at_risk = sum(at_risk),
    percentage = 100 * sum(at_risk) / length(at_risk),
    rows_risk_df = result_df,
    error_range = error_range
  )
}

#' Evaluate categorical attribute inference risk
#' @noRd
.evaluate_categorical_rapid <- function(A, B, original_data,
                                        cat_tau, cat_eval_method,
                                        target_var, return_all_records) {

  true_labels <- A
  predicted_probs <- B

  if (!is.factor(true_labels)) {
    stop("true_labels must be a factor for categorical evaluation.")
  }
  if (is.data.frame(predicted_probs)) {
    predicted_probs <- as.matrix(predicted_probs)
  }

  # g_i: predicted probability for the true class
  g_i <- predicted_probs[cbind(
    seq_len(nrow(predicted_probs)),
    match(as.character(true_labels), colnames(predicted_probs))
  )]

  if (cat_eval_method == "RCS_conditional") {
    # Class-conditional baseline
    baseline_by_class <- sapply(unique(as.character(true_labels)), function(k) {
      mean(predicted_probs[which(as.character(true_labels) == k), k], na.rm = TRUE)
    })
    names(baseline_by_class) <- unique(as.character(true_labels))
    b_i <- baseline_by_class[as.character(true_labels)]
    r_i <- g_i / b_i
    at_risk <- r_i > cat_tau
    confidence_rate <- sum(at_risk, na.rm = TRUE) / length(at_risk)

    result_df <- data.frame(
      original_data,
      predicted_class = colnames(predicted_probs)[max.col(predicted_probs)],
      true_prob = g_i,
      baseline = b_i,
      relative_score = r_i,
      threshold = cat_tau,
      at_risk = at_risk,
      stringsAsFactors = FALSE
    )

  } else if (cat_eval_method == "RCS_marginal") {
    # Marginal baseline
    marginal_freq <- prop.table(table(original_data[[target_var]]))
    b_i <- as.numeric(marginal_freq[as.character(true_labels)])
    max_improvement <- 1 - b_i
    normalized_gain <- (g_i - b_i) / max_improvement
    at_risk <- normalized_gain > cat_tau
    confidence_rate <- sum(at_risk, na.rm = TRUE) / length(at_risk)

    result_df <- data.frame(
      original_data,
      predicted_class = colnames(predicted_probs)[max.col(predicted_probs)],
      true_prob = g_i,
      baseline = b_i,
      normalized_gain = normalized_gain,
      threshold = cat_tau,
      at_risk = at_risk,
      stringsAsFactors = FALSE
    )

  } else { # NCE
    ce <- -log(pmax(g_i, 1e-10))
    n_classes <- ncol(predicted_probs)
    max_entropy <- log(n_classes)
    normalized_ce <- ce / max_entropy
    risk_score <- 1 - normalized_ce
    at_risk <- risk_score > cat_tau
    confidence_rate <- sum(at_risk, na.rm = TRUE) / length(at_risk)

    result_df <- data.frame(
      original_data,
      predicted_class = colnames(predicted_probs)[max.col(predicted_probs)],
      true_prob = g_i,
      risk_score = risk_score,
      threshold = cat_tau,
      at_risk = at_risk,
      stringsAsFactors = FALSE
    )
  }

  if (!return_all_records) {
    result_df <- result_df[result_df$at_risk, , drop = FALSE]
  }

  list(
    confidence_rate = confidence_rate,
    n_at_risk = sum(at_risk, na.rm = TRUE),
    percentage = 100 * confidence_rate,
    rows_risk_df = result_df
  )
}


# =============================================================================
# S3 Methods
# =============================================================================

#' Print method for rapid objects
#'
#' @param x an object of class "rapid"
#' @param ... additional arguments (ignored)
#' @export
print.rapid <- function(x, ...) {
  cat("RAPID Disclosure Risk Assessment\n")
  cat("=================================\n")
  cat("Model:", x$model_type, "\n")
  cat("Method:", x$method, "\n")
  cat("Key variables:", paste(x$key_vars, collapse = ", "), "\n")
  cat("Target variable:", x$target_var, "\n\n")

  cat("Risk Assessment:\n")
  cat("  RAPID (confidence rate):", sprintf("%.4f", x$rapid), "\n")
  cat("  Records at risk:", x$n_at_risk, "/", x$n_total,
      sprintf("(%.1f%%)", x$pct_at_risk), "\n")
  cat("  Threshold:", x$threshold, "\n\n")

  cat("Model Performance:\n")
  if (x$is_categorical) {
    cat("  Accuracy:", sprintf("%.4f", x$model_metrics$accuracy), "\n")
  } else {
    cat("  MAE:", sprintf("%.2f", x$model_metrics$mae), "\n")
    cat("  RMSE:", sprintf("%.2f", x$model_metrics$rmse), "\n")
  }

  invisible(x)
}

#' Summary method for rapid objects
#'
#' @param object an object of class "rapid"
#' @param ... additional arguments (ignored)
#' @export
summary.rapid <- function(object, ...) {
  summ <- list(
    rapid = object$rapid,
    n_at_risk = object$n_at_risk,
    pct_at_risk = object$pct_at_risk,
    n_total = object$n_total,
    method = object$method,
    threshold = object$threshold,
    model_type = object$model_type,
    model_metrics = object$model_metrics,
    is_categorical = object$is_categorical,
    key_vars = object$key_vars,
    target_var = object$target_var,
    risk_level = if (object$rapid < 0.05) "LOW" else if (object$rapid < 0.15) "MEDIUM" else "HIGH"
  )
  class(summ) <- "summary.rapid"
  summ
}

#' Print method for summary.rapid objects
#'
#' @param x an object of class "summary.rapid"
#' @param ... additional arguments (ignored)
#' @export
print.summary.rapid <- function(x, ...) {
  cat("Summary: RAPID Disclosure Risk Assessment\n")
  cat("==========================================\n\n")

  cat("Configuration:\n")
  cat("  Model:", x$model_type, "\n")
  cat("  Method:", x$method, "\n")
  cat("  Threshold:", x$threshold, "\n")
  cat("  Key variables:", paste(x$key_vars, collapse = ", "), "\n")
  cat("  Target variable:", x$target_var, "\n\n")

  cat("Risk Assessment:\n")
  cat("  RAPID score:", sprintf("%.4f", x$rapid), "\n")
  cat("  Records at risk:", x$n_at_risk, "/", x$n_total,
      sprintf("(%.2f%%)", x$pct_at_risk), "\n")
  cat("  Risk level:", x$risk_level, "\n\n")

  cat("Model Performance:\n")
  if (x$is_categorical) {
    cat("  Accuracy:", sprintf("%.4f (%.1f%%)", x$model_metrics$accuracy,
                              100 * x$model_metrics$accuracy), "\n")
  } else {
    cat("  MAE:", sprintf("%.4f", x$model_metrics$mae), "\n")
    cat("  RMSE:", sprintf("%.4f", x$model_metrics$rmse), "\n")
    cat("  Relative MAE:", sprintf("%.4f", x$model_metrics$rmae), "\n")
    cat("  Relative RMSE:", sprintf("%.4f", x$model_metrics$rrmse), "\n")
  }

  cat("\nInterpretation Guidelines:\n")
  if (x$risk_level == "LOW") {
    cat("  Low risk: Synthetic data provides good privacy protection.\n")
  } else if (x$risk_level == "MEDIUM") {
    cat("  Medium risk: Some records may be at risk. Review the at-risk records.\n")
  } else {
    cat("  High risk: Significant disclosure risk detected. Consider additional protection.\n")
  }

  invisible(x)
}

#' Plot method for rapid objects
#'
#' Creates visualizations for RAPID disclosure risk assessment results.
#'
#' @param x an object of class "rapid"
#' @param y not used
#' @param ... additional arguments passed to plotting functions
#' @param which which plot to show: 1 = risk distribution, 2 = prediction scatter (numeric only)
#' @param xlim character or numeric vector controlling x-axis limits for risk distribution:
#'   \itemize{
#'     \item \code{"full"} (default): show full data range
#'     \item \code{"zero_to_one"}: fixed range from 0 to 1
#'     \item \code{"threshold_to_one"}: range from threshold to 1 (at-risk region only)
#'     \item numeric vector of length 2: custom limits, e.g., \code{c(-0.5, 1)}
#'   }
#' @param annotate logical, if TRUE (default) add annotation showing percentage at risk
#' @param facet_by character vector of variable names to facet by (creates subgroup panels).
#'   Must be columns present in the records data frame. Uses ggplot2 for faceted plots.
#' @param bins integer, number of histogram bins (default: 30)
#'
#' @details
#' For categorical sensitive variables with RCS_marginal method, the plot shows the
#' distribution of normalized gain values. The normalized gain measures how much better
#' the model predicts compared to the marginal baseline:
#' \itemize{
#'   \item Values > 0: model predicts better than baseline (potential disclosure risk)
#'   \item Values = 0: model predicts at baseline rate (no additional risk)
#'   \item Values < 0: model predicts worse than baseline (no risk)
#' }
#'
#' The red dashed line indicates the threshold; records to the right are "at risk".
#'
#' @return For faceted plots, returns a ggplot2 object (invisibly). For base plots,
#'   returns NULL invisibly.
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' n <- 200
#' original <- data.frame(
#'   age = sample(20:60, n, replace = TRUE),
#'   gender = factor(sample(c("M", "F"), n, replace = TRUE)),
#'   health = factor(sample(c("Good", "Fair", "Poor"), n, replace = TRUE))
#' )
#' synthetic <- original
#' synthetic$health <- factor(sample(c("Good", "Fair", "Poor"), n, replace = TRUE))
#'
#' result <- rapid(original, synthetic,
#'                 key_vars = c("age", "gender"),
#'                 target_var = "health",
#'                 model_type = "rf",
#'                 return_all_records = TRUE)
#'
#' # Basic plot
#' plot(result)
#'
#' # Fixed x-axis from 0 to 1
#' plot(result, xlim = "zero_to_one")
#'
#' # Facet by gender
#' plot(result, facet_by = "gender")
#' }
#'
#' @export
#' @importFrom graphics hist abline par legend points text mtext
#' @importFrom ggplot2 ggplot aes geom_histogram geom_vline annotate facet_wrap
#'   labs theme_minimal theme element_text scale_fill_manual after_stat
plot.rapid <- function(x, y = NULL, ...,
                       which = 1,
                       xlim = "full",
                       annotate = TRUE,
                       facet_by = NULL,
                       bins = 30) {

  show <- rep(FALSE, 2)
  show[which] <- TRUE

  records <- x$records

  # Check if we have all records (needed for proper visualization)
  if (nrow(records) == 0) {
    message("No records to plot. Run rapid() with return_all_records = TRUE for full distribution.")
    return(invisible(NULL))
  }

  # Determine score variable and labels based on method
  if (x$is_categorical) {
    if ("normalized_gain" %in% names(records)) {
      score_var <- "normalized_gain"
      main_title <- "Normalized Gain Distribution"
      x_label <- "Normalized Gain"
    } else if ("relative_score" %in% names(records)) {
      score_var <- "relative_score"
      main_title <- "Relative Score Distribution"
      x_label <- "Relative Score"
    } else if ("risk_score" %in% names(records)) {
      score_var <- "risk_score"
      main_title <- "Risk Score Distribution"
      x_label <- "Risk Score"
    } else {
      score_var <- "true_prob"
      main_title <- "True Probability Distribution"
      x_label <- "Probability"
    }
  } else {
    if ("symmetric_error_pct" %in% names(records)) {
      score_var <- "symmetric_error_pct"
      x_label <- "Symmetric Error (%)"
    } else if ("stabilised_relative_error_pct" %in% names(records)) {
      score_var <- "stabilised_relative_error_pct"
      x_label <- "Stabilised Relative Error (%)"
    } else {
      score_var <- "absolute_error"
      x_label <- "Absolute Error"
    }
    main_title <- "Prediction Error Distribution"
  }

  scores <- records[[score_var]]

  # Calculate x-axis limits
  # Use stored error_range (from full data) if available, otherwise use current data range
  score_range <- range(scores, na.rm = TRUE)
  full_range <- if (!is.null(x$error_range)) x$error_range else score_range

  if (is.character(xlim)) {
    xlim_values <- switch(xlim,
      "full" = {
        # For numeric variables, use the full error range from all records
        if (!x$is_categorical && !is.null(x$error_range)) {
          c(0, x$error_range[2])
        } else {
          full_range
        }
      },
      "zero_to_one" = c(0, 1),
      "threshold_to_one" = c(x$threshold, 1),
      full_range  # default to full
    )
  } else if (is.numeric(xlim) && length(xlim) == 2) {
    xlim_values <- xlim
  } else {
    xlim_values <- full_range
  }

  # Calculate at-risk percentage using stored totals (not nrow which may be filtered)
  pct_at_risk <- x$pct_at_risk
  n_at_risk <- x$n_at_risk
  n_total <- x$n_total

  # Check if we have all records or just at-risk records
  have_all_records <- nrow(records) == n_total

  # Warn if we don't have all records (distribution will be incomplete)
  if (!have_all_records) {
    warning("Only at-risk records are available for plotting. ",
            "'Not at Risk' values cannot be shown.\n",
            "To see the complete distribution, run rapid() with return_all_records = TRUE.",
            call. = FALSE)
  }

  # Use ggplot2 for faceted plots
  if (!is.null(facet_by) && show[1]) {
    # Check facet variables exist
    missing_vars <- setdiff(facet_by, names(records))
    if (length(missing_vars) > 0) {
      stop(paste("Facet variables not found in records:", paste(missing_vars, collapse = ", ")))
    }

    # Create facet variable (combine if multiple)
    if (length(facet_by) == 1) {
      records$facet_var <- records[[facet_by]]
      facet_label <- facet_by
    } else {
      records$facet_var <- interaction(records[, facet_by, drop = FALSE], sep = " / ")
      facet_label <- paste(facet_by, collapse = " / ")
    }

    # Calculate per-facet statistics
    facet_stats <- stats::aggregate(
      at_risk ~ facet_var,
      data = records,
      FUN = function(x) sprintf("%.1f%% at risk", 100 * mean(x))
    )
    names(facet_stats)[2] <- "label"

    # Build ggplot
    p <- ggplot2::ggplot(records, ggplot2::aes(x = .data[[score_var]])) +
      ggplot2::geom_histogram(
        ggplot2::aes(fill = at_risk),
        bins = bins,
        color = "white",
        alpha = 0.8
      ) +
      ggplot2::geom_vline(
        xintercept = x$threshold,
        color = "red",
        linetype = "dashed",
        linewidth = 1
      ) +
      ggplot2::scale_fill_manual(
        values = c("FALSE" = "steelblue", "TRUE" = "coral"),
        labels = c("Not at Risk", "At Risk"),
        name = "Status"
      ) +
      ggplot2::facet_wrap(~ facet_var, scales = "free_y") +
      ggplot2::coord_cartesian(xlim = xlim_values) +
      ggplot2::labs(
        title = main_title,
        subtitle = sprintf("Overall: %.1f%% at risk (threshold = %.2f)", pct_at_risk, x$threshold),
        x = x_label,
        y = "Count"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold"),
        strip.text = ggplot2::element_text(face = "bold")
      )

    print(p)
    return(invisible(p))
  }

  # Base R plots (non-faceted)
  if (sum(show) > 1) {
    op <- par(mfrow = c(1, 2))
    on.exit(par(op))
  }

  if (show[1]) {
    # Create histogram breaks that span the xlim range
    break_min <- min(xlim_values[1], score_range[1])
    break_max <- max(xlim_values[2], score_range[2])
    breaks_seq <- seq(break_min, break_max, length.out = bins + 1)

    # Compute histogram
    h <- hist(scores, breaks = breaks_seq, plot = FALSE)

    # Color bars based on whether they're in the at-risk region
    # For categorical: at_risk when score >= threshold (high confidence = risk)
    # For numeric: at_risk when error < threshold (low error = disclosure risk)
    if (x$is_categorical) {
      bar_colors <- ifelse(h$mids >= x$threshold, "coral", "steelblue")
    } else {
      # For numeric: error < threshold means at-risk (accurate predictions = disclosure risk)
      bar_colors <- ifelse(h$mids < x$threshold, "coral", "steelblue")
    }

    # Plot histogram with explicit xlim (don't pass ... to avoid xlim conflict)
    plot(h, main = main_title,
         xlab = x_label, col = bar_colors, border = "white",
         xlim = xlim_values, ylab = "Frequency")
    abline(v = x$threshold, col = "red", lwd = 2, lty = 2)

    # Add annotation
    if (annotate) {
      # Position annotation in upper right
      usr <- par("usr")
      text_x <- usr[2] - 0.05 * (usr[2] - usr[1])
      text_y <- usr[4] - 0.05 * (usr[4] - usr[3])
      text(text_x, text_y,
           labels = sprintf("%.1f%% at risk\n(n = %d / %d)",
                            pct_at_risk, n_at_risk, n_total),
           adj = c(1, 1), cex = 0.9, font = 2)
    }

    # Legend with clearer labels for numeric
    # Only show "Not at Risk" in legend if those records exist
    if (x$is_categorical) {
      legend("topleft",
             legend = c("At Risk", "Not at Risk", "Threshold"),
             fill = c("coral", "steelblue", NA),
             border = c("white", "white", NA),
             col = c(NA, NA, "red"),
             lty = c(NA, NA, 2),
             lwd = c(NA, NA, 2),
             cex = 0.8)
    } else {
      legend("topright",
             legend = c(sprintf("At Risk (error < %.0f%%)", x$threshold),
                        sprintf("Not at Risk (error >= %.0f%%)", x$threshold),
                        "Threshold"),
             fill = c("coral", "steelblue", NA),
             border = c("white", "white", NA),
             col = c(NA, NA, "red"),
             lty = c(NA, NA, 2),
             lwd = c(NA, NA, 2),
             cex = 0.8)
    }
  }

  if (show[2] && !x$is_categorical) {
    # Scatter plot of predicted vs actual with transparency
    # Plot not-at-risk points first (background), then at-risk points on top
    col_at_risk <- adjustcolor("coral", alpha.f = 0.6)
    col_not_risk <- adjustcolor("steelblue", alpha.f = 0.6)

    # Initialize empty plot
    plot(range(records$original_value), range(records$predicted_value),
         type = "n",
         main = "Predicted vs Actual Values",
         xlab = "Actual Value", ylab = "Predicted Value")

    # Plot not-at-risk points first (in background)
    not_at_risk_idx <- !records$at_risk
    if (any(not_at_risk_idx)) {
      points(records$original_value[not_at_risk_idx],
             records$predicted_value[not_at_risk_idx],
             col = col_not_risk, pch = 19, cex = 0.7)
    }

    # Plot at-risk points on top
    at_risk_idx <- records$at_risk
    if (any(at_risk_idx)) {
      points(records$original_value[at_risk_idx],
             records$predicted_value[at_risk_idx],
             col = col_at_risk, pch = 19, cex = 0.7)
    }

    abline(0, 1, col = "darkgray", lwd = 2, lty = 2)

    if (annotate) {
      usr <- par("usr")
      text_x <- usr[1] + 0.05 * (usr[2] - usr[1])
      text_y <- usr[4] - 0.05 * (usr[4] - usr[3])
      text(text_x, text_y,
           labels = sprintf("%.1f%% at risk\n(n = %d / %d)",
                            pct_at_risk, n_at_risk, n_total),
           adj = c(0, 1), cex = 0.9, font = 2)
    }

    legend("bottomright",
           legend = c(sprintf("At Risk (error < %.0f%%)", x$threshold),
                      sprintf("Not at Risk (error >= %.0f%%)", x$threshold),
                      "Perfect Prediction"),
           col = c(col_at_risk, col_not_risk, "darkgray"),
           pch = c(19, 19, NA), lty = c(NA, NA, 2), lwd = c(NA, NA, 2),
           cex = 0.8)
  }

  invisible(NULL)
}
