#' RAPID: Risk of Attribute Prediction-Induced Disclosure
#'
#' Assesses the inferential disclosure risk of a sensitive variable in synthetic
#' or anonymized data. An attacker model is trained on the released data to
#' predict the sensitive variable from the quasi-identifiers and is then applied
#' to the original records. A record is at risk when the attacker is confident
#' about its true value well beyond what the marginal distribution of the
#' sensitive variable already gives away. RAPID is the share of such records.
#'
#' @param X data frame of original data, or a \code{\link{synth_pair}} object
#' @param Y data frame of synthetic data (not needed if X is a synth_pair)
#' @param key_vars character vector of quasi-identifier variable names
#' @param target_var character, name of the sensitive target variable
#' @param model_type character, the attacker model: "rf" (random forest,
#'   default), "lm" (linear model), "cart" (decision tree), "gbm" (gradient
#'   boosting), or "logit" (logistic regression, two classes only)
#' @param num_eval_method character, rule for a numeric target: "band"
#'   (default) compares the probability the attacker puts into the tolerance band
#'   around the true value with the share of original values in that band;
#'   "point" flags a record when the error of the point prediction is below the
#'   tolerance (the rule of riskutility <= 0.2.0)
#' @param num_error_metric character, error metric that defines the tolerance
#'   band: "symmetric" (default), "stabilised_relative", or "absolute"
#' @param num_epsilon numeric tolerance. For the two relative metrics it is a
#'   proportion (default 0.1, that is 10 percent). For the absolute metric it is
#'   given in the units of the target and has no default. A relative tolerance of
#'   1 or more is read as percent, with a warning.
#' @param num_delta numeric, smoothing constant of the relative metrics (default: 0.01)
#' @param num_tau numeric threshold on the normalized gain of a numeric target
#'   under the band rule (default: 0.3)
#' @param num_dist character, how the attacker's predictive distribution of a
#'   numeric target is evaluated: "auto" (default) uses the exact probability
#'   where the model provides one (normal or lognormal for "lm", the training
#'   values in the leaf for "cart") and a grid of quantiles for "rf" and "gbm";
#'   "grid" uses \code{n_quantiles} quantiles for all models
#' @param num_scale character, working scale of the linear model for a numeric
#'   target: "identity", "log", or "auto" (default), which takes the log scale
#'   when the error metric is relative, all synthetic target values are positive
#'   and the log model has the higher Gaussian likelihood on the original scale.
#'   On the log scale the point prediction is the predictive median.
#' @param n_quantiles integer, number of quantiles of a grid distribution (default: 100)
#' @param cat_eval_method character, rule for a categorical target:
#'   "RCS_marginal" (default, normalized gain over the marginal class share),
#'   "RCS_conditional" (ratio to the class-conditional mean probability), or
#'   "NCE" (normalized cross-entropy)
#' @param cat_tau numeric threshold for a categorical target. The default
#'   depends on the method: 0.3 for RCS_marginal, 1 for RCS_conditional, 0.5 for NCE.
#' @param na_strategy character, strategy for NA values in a numeric sensitive
#'   attribute: "constant" (default), "drop", or "median". Missing values of a
#'   categorical variable become the level "missing".
#' @param na_constant_value numeric, value for constant NA imputation (default: 0)
#' @param n_perm integer, number of permutations of the true values that give the
#'   chance level (default: 10). The attacker is kept fixed, each record's
#'   prediction is scored against a randomly assigned true value and the baseline
#'   of that value. 0 skips the permutations.
#' @param threshold_rule character, \code{"fixed"} (default) flags with the
#'   threshold of the rule; \code{"permutation"} takes the \code{perm_level}
#'   quantile of all permuted scores as threshold. Only for rules with a gain
#'   over a baseline (band, RCS_marginal) and with \code{n_perm >= 1}.
#' @param perm_level numeric in (0, 1), quantile of the permuted scores that is
#'   the threshold under \code{threshold_rule = "permutation"} (default: 0.95)
#' @param store_model logical, if TRUE store the fitted model object in the result
#'   (needed for \code{plot(result, which = 4)} variable importance plot). Default FALSE.
#' @param store_predictions logical, if TRUE (default) the result keeps the
#'   attacker's predictions, so that it can be tested and rescored without
#'   fitting the model again
#' @param seed integer, random seed for reproducibility (default: NULL). The fit
#'   and the permutations are each started from the seed.
#' @param verbose logical, print diagnostic messages (default: FALSE)
#' @param ... additional arguments passed to the model fitting function
#'   (\code{ranger::ranger}, \code{stats::lm}, \code{rpart::rpart},
#'   \code{xgboost::xgb.train} parameters and \code{nrounds}, \code{stats::glm}).
#'   For "cart", \code{prune = FALSE} keeps the tree as grown. The old argument
#'   names \code{original_data}, \code{synthetic_data}, \code{quasi_identifiers}
#'   and \code{sensitive_attribute} are accepted as aliases of the first four
#'   arguments. A misspelled argument of \code{rapid()} is an error.
#' @param num_epsilon_type deprecated. "percentage" reads \code{num_epsilon} as
#'   percent, "absolute" selects the absolute error metric.
#' @param return_all_records deprecated and ignored: the result always holds all
#'   records. Use \code{as.data.frame(x, at_risk_only = TRUE)} for the flagged ones.
#'
#' @return An object of class "rapid" containing:
#' \itemize{
#'   \item rapid: the proportion of records at risk
#'   \item n_at_risk, pct_at_risk, n_total: number and percentage of records at
#'     risk, number of scored records
#'   \item rule: the scoring rule; method and threshold are copies of its
#'     method and threshold
#'   \item score, baseline, at_risk: per record, the score that is compared with
#'     the threshold, the marginal baseline (NULL for NCE) and the flag
#'   \item records: the prepared original data with the columns of the rule
#'     (band rule: predicted_value, band_lo, band_hi, g, b, r; point rule:
#'     predicted_value, error, b; RCS_marginal: predicted_class, g, b, r;
#'     RCS_conditional: predicted_class, g, b, ratio; NCE: predicted_class, g,
#'     nce) and at_risk, for all records. Here g is the attacker's confidence,
#'     b the baseline and r the normalized gain (g - b) / (1 - b). A column whose
#'     name is already used by a variable gets the suffix ".rapid".
#'   \item attack: the attacker's predictions (NULL if store_predictions = FALSE)
#'   \item chance, chance_rates, excess: the chance level (mean share flagged
#'     when the true values are permuted; NA if n_perm = 0), the share of each
#'     permutation, and rapid - chance
#'   \item n_perm, threshold_rule, perm_level: the permutation settings
#'   \item model: the fitted model object (only if store_model = TRUE)
#'   \item model_type, attacker: the model and a line with the settings that
#'     drive its result (number of trees and node size, number of leaves, scale)
#'   \item model_metrics: model performance (accuracy for categorical, MAE/RMSE for numeric)
#'   \item formula, key_vars, target_var, is_categorical, na, call, version
#' }
#'
#' @details
#' RAPID (Risk of Attribute Prediction-Induced Disclosure) measures disclosure risk
#' by training a predictive model on synthetic data and evaluating how well it can
#' predict the sensitive attribute in the original data.
#'
#' Unlike CAP-based methods (DCAP, TCAP) which use exact or fuzzy matching, RAPID
#' uses machine learning models to capture more complex relationships between
#' quasi-identifiers and sensitive attributes.
#'
#' \strong{For categorical sensitive variables:}
#' \itemize{
#'   \item \code{RCS_marginal} (default): the attacker's probability of the true
#'     class, g, is compared with the share of that class in the original data, b.
#'     The record is flagged when the normalized gain (g - b) / (1 - b) exceeds
#'     \code{cat_tau}.
#'   \item \code{RCS_conditional}: ratio of g to its mean among the records of
#'     the same class.
#'   \item \code{NCE}: Normalized Cross-Entropy, measures information leakage.
#' }
#'
#' \strong{For continuous sensitive variables:}
#' The band rule (default) is the continuous counterpart of \code{RCS_marginal}.
#' The tolerance band of a record holds all values whose error against the true
#' value is below \code{num_epsilon}. g is the probability the attacker's
#' predictive distribution puts into the band, b the share of original values in
#' it, and the record is flagged when (g - b) / (1 - b) exceeds \code{num_tau}.
#' Where the band holds all original values (b = 1) the gain is 0. The point
#' rule flags a record when the error of the point prediction is below the
#' tolerance. It also flags records whose value is easy to guess without any
#' model, which is why it is no longer the default.
#'
#' The band is defined by one of three error metrics:
#' \itemize{
#'   \item \code{symmetric} (recommended): symmetric relative error, treats predicted
#'     and true values equally.
#'   \item \code{stabilised_relative}: stabilised relative error.
#'   \item \code{absolute}: absolute error in the units of the target.
#' }
#'
#' \strong{Attacker models:} The result depends on how sharp the attacker's
#' predictive distribution is, so the settings that drive it are part of the
#' result (\code{attacker}) and of the printed output. The tree of "cart" is
#' grown with \code{cp = 0}, \code{minbucket = 5}, \code{xval = 10} and pruned at
#' the smallest cross-validated error; the three values can be overridden and
#' \code{prune = FALSE} keeps the tree as grown. Cross-validation is random, so
#' set \code{seed} for a reproducible tree. To reproduce the results of the RAPID
#' paper and of the RAPID package use \code{num_dist = "grid"},
#' \code{num_scale = "identity"} and, for "cart", \code{prune = FALSE, cp = 0.01,
#' minsplit = 20, minbucket = 7} (the defaults of rpart).
#'
#' @section Changes in 0.3.0:
#' The defaults changed: numeric targets use the band rule with \code{num_tau =
#' 0.3} and a tolerance of 10 percent (before: the point rule), categorical
#' targets use \code{RCS_marginal} with \code{cat_tau = 0.3} (before:
#' \code{RCS_conditional} with \code{cat_tau = 1}). A relative tolerance is a
#' proportion (before: percent). The old rules remain available through
#' \code{num_eval_method = "point"} and \code{cat_eval_method =
#' "RCS_conditional"}.
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
#' @family rapid
#' @export
#' @importFrom stats predict lm as.formula glm binomial model.matrix complete.cases median
#' @examples
#' # Create example data
#' set.seed(123)
#' n <- 200
#' original <- data.frame(
#'   age = sample(20:70, n, replace = TRUE),
#'   gender = factor(sample(c("M", "F"), n, replace = TRUE)),
#'   region = factor(sample(c("N", "S", "E", "W"), n, replace = TRUE))
#' )
#' original$income <- round(20000 + 600 * original$age + rnorm(n, 0, 5000))
#'
#' # Synthetic version with some noise
#' synthetic <- data.frame(
#'   age = original$age + sample(-3:3, n, replace = TRUE),
#'   gender = factor(sample(c("M", "F"), n, replace = TRUE)),
#'   region = factor(sample(c("N", "S", "E", "W"), n, replace = TRUE)),
#'   income = round(original$income * runif(n, 0.9, 1.1))
#' )
#'
#' # RAPID for a continuous sensitive variable (band rule, linear model)
#' result <- rapid(
#'   X = original,
#'   Y = synthetic,
#'   key_vars = c("age", "gender", "region"),
#'   target_var = "income",
#'   model_type = "lm"
#' )
#' result
#' head(as.data.frame(result, at_risk_only = TRUE))
#'
#' \donttest{
#' # Random forest attacker
#' result_rf <- rapid(original, synthetic,
#'                    key_vars = c("age", "gender", "region"),
#'                    target_var = "income", model_type = "rf", seed = 1)
#' summary(result_rf)
#'
#' # The point rule of earlier versions, with a tolerance of 10 percent
#' rapid(original, synthetic,
#'       key_vars = c("age", "gender", "region"), target_var = "income",
#'       model_type = "rf", num_eval_method = "point", num_epsilon = 0.1)
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
#'   model_type = "rf"
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
#'     num_epsilon = 0.1,          # tolerance band of 10 percent
#'     num_error_metric = "symmetric",
#'     verbose = TRUE
#'   )
#'
#'   print(rapid_income)
#'   summary(rapid_income)
#'   plot(rapid_income, which = 1)
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
  if (missing(X)) {
    # Without X the data arrive under their alias names. The default method is
    # called with the expressions of the caller, so that it can store the call.
    cl <- sys.call()
    cl[[1L]] <- rapid.default
    return(eval(cl, parent.frame()))
  }
  UseMethod("rapid")
}

#' @rdname rapid
#' @export
rapid.synth_pair <- function(X, ...) {
  # The default method unpacks the pair (data, key_vars, target_var)
  NextMethod()
}

#' @rdname rapid
#' @export
rapid.default <- function(X, Y,
                          key_vars,
                          target_var,
                          model_type = c("rf", "lm", "cart", "gbm", "logit"),
                          # Numeric-specific
                          num_eval_method = c("band", "point"),
                          num_error_metric = c("symmetric", "stabilised_relative", "absolute"),
                          num_epsilon = NULL,
                          num_delta = 0.01,
                          num_tau = 0.3,
                          num_dist = c("auto", "grid", "exact"),
                          num_scale = c("auto", "identity", "log"),
                          n_quantiles = 100,
                          # Categorical-specific
                          cat_eval_method = c("RCS_marginal", "RCS_conditional", "NCE"),
                          cat_tau = NULL,
                          # NA handling
                          na_strategy = c("constant", "drop", "median"),
                          na_constant_value = 0,
                          # Chance level and threshold
                          n_perm = 10,
                          threshold_rule = c("fixed", "permutation"),
                          perm_level = 0.95,
                          # Output options
                          store_model = FALSE,
                          store_predictions = TRUE,
                          seed = NULL,
                          verbose = FALSE,
                          ...,
                          num_epsilon_type,
                          return_all_records) {

  call <- match.call()
  call[[1L]] <- quote(rapid)
  given <- names(call)[-1L]

  # Choices
  model_type <- .rapid_choice(model_type, .rapid_learner_names, "model_type")
  num_dist <- .rapid_choice(num_dist, c("auto", "grid", "exact"), "num_dist")
  num_scale <- .rapid_choice(num_scale, c("auto", "identity", "log"), "num_scale")
  na_strategy <- .rapid_choice(na_strategy, c("constant", "drop", "median"), "na_strategy")

  # Data arguments (aliases, synth_pair) and arguments for the model
  inputs <- .rapid_inputs(
    if (missing(X)) NULL else X, if (missing(Y)) NULL else Y,
    if (missing(key_vars)) NULL else key_vars,
    if (missing(target_var)) NULL else target_var,
    list(...)
  )
  .rapid_check_dots(inputs$dots, setdiff(names(formals(rapid.default)), "..."))
  if (!missing(return_all_records)) {
    .rapid_warn(c("'return_all_records' is deprecated and ignored: the object always ",
                  "keeps all records. Use as.data.frame(x, at_risk_only = TRUE) for the ",
                  "flagged records."), "deprecated")
  }

  if (!is.null(seed)) set.seed(seed)
  if (verbose) message("RAPID: Running disclosure risk assessment...")

  # Common column types, missing values and factor levels in both data sets
  prepared <- .rapid_prepare_data(inputs$X, inputs$Y, inputs$key_vars, inputs$target_var,
                                  na_strategy, na_constant_value)
  is_categorical <- prepared$is_categorical

  # Scoring rule and permutation settings
  scoring <- .rapid_scoring(
    is_categorical, given,
    num_eval_method = num_eval_method, num_error_metric = num_error_metric,
    num_epsilon = num_epsilon,
    num_epsilon_type = if (missing(num_epsilon_type)) NULL else num_epsilon_type,
    num_delta = num_delta, num_tau = num_tau,
    cat_eval_method = cat_eval_method, cat_tau = cat_tau,
    n_perm = n_perm, threshold_rule = threshold_rule, perm_level = perm_level
  )

  # Attacker
  learner <- .rapid_learner(model_type, inputs$dots)
  .learner_check(learner, is_categorical, num_dist,
                 num_scale = if (is_categorical) "auto" else num_scale)
  if (!learner$unseen_levels) .rapid_check_levels(prepared$unseen, learner$name)

  # Fit on the synthetic data, predict the original records
  if (verbose) message("  Fitting ", learner$name, " model on synthetic data...")
  attack <- .rapid_attack(prepared$X, prepared$Y, inputs$key_vars, inputs$target_var,
                          learner, num_dist, num_scale, scoring$rule$metric, n_quantiles)

  if (verbose) message("  Evaluating disclosure risk...")
  # The permutations start from the seed as well, so that update() without a
  # refit draws the same permutations as a fresh call
  if (!is.null(seed)) set.seed(seed)
  .rapid_assess(
    attack, scoring$rule, data = prepared$X,
    meta = list(key_vars = inputs$key_vars, target_var = inputs$target_var,
                na = list(strategy = na_strategy, constant = na_constant_value,
                          dropped = prepared$dropped),
                call = call),
    store_model = isTRUE(store_model), store_predictions = isTRUE(store_predictions),
    n_perm = scoring$n_perm, threshold_rule = scoring$threshold_rule,
    perm_level = scoring$perm_level
  )
}

# Arguments of rapid() that change only how the attack is scored. update()
# rescores the stored attack when it is given nothing else.
.rapid_scoring_args <- c("num_eval_method", "num_error_metric", "num_epsilon", "num_delta",
                         "num_tau", "cat_eval_method", "cat_tau", "n_perm", "threshold_rule",
                         "perm_level", "seed")

#' Scoring rule and permutation settings from the arguments of rapid()
#'
#' Shared by `rapid.default()` and `update.rapid()`, so that a rescored object
#' gets the rule a fresh call would get (choices, tolerance resolution, default
#' thresholds, deprecation warnings, warnings on arguments of the other type of
#' target).
#'
#' @param is_categorical Type of the target.
#' @param given Names of the arguments given in the call.
#' @return `list(rule, n_perm, threshold_rule, perm_level)`.
#' @noRd
.rapid_scoring <- function(is_categorical, given, num_eval_method, num_error_metric,
                           num_epsilon, num_epsilon_type, num_delta, num_tau,
                           cat_eval_method, cat_tau, n_perm, threshold_rule, perm_level) {
  num_eval_method <- .rapid_choice(num_eval_method, c("band", "point"), "num_eval_method")
  num_error_metric <- .rapid_choice(num_error_metric, .rapid_metrics, "num_error_metric")
  cat_eval_method <- .rapid_choice(cat_eval_method,
                                   c("RCS_marginal", "RCS_conditional", "NCE"),
                                   "cat_eval_method")
  threshold_rule <- .rapid_choice(threshold_rule, c("fixed", "permutation"), "threshold_rule")

  # Arguments that belong to the other type of target
  other <- grep(if (is_categorical) "^(num_|n_quantiles$)" else "^cat_", given, value = TRUE)
  if (length(other) > 0) {
    .rapid_warn(c("Ignored for a ", if (is_categorical) "categorical" else "numeric",
                  " target: ", paste(other, collapse = ", "), "."), "ignored_argument")
  }

  rule <- .rapid_make_rule(
    is_categorical,
    num = list(method = num_eval_method, metric = num_error_metric, epsilon = num_epsilon,
               epsilon_type = num_epsilon_type, delta = num_delta, tau = num_tau),
    cat = list(method = cat_eval_method, tau = cat_tau)
  )
  .rapid_check_perm(n_perm, threshold_rule, perm_level, rule)
  list(rule = rule, n_perm = n_perm, threshold_rule = threshold_rule, perm_level = perm_level)
}

#' Prepared original data of a rapid object
#'
#' The record table is the prepared data followed by the columns of the rule.
#' @noRd
.rapid_data <- function(x) {
  x$records[, seq_len(ncol(x$records) - length(x$record_cols)), drop = FALSE]
}

#' Score the stored attack of a rapid object under a new rule
#'
#' @param x A `rapid` object with stored predictions.
#' @param rule A `rapid_rule`.
#' @param n_perm,threshold_rule,perm_level Permutation settings.
#' @param call The call stored in the result.
#' @noRd
.rapid_rescore <- function(x, rule, n_perm, threshold_rule, perm_level, call) {
  .rapid_require(x, c("version", "predictions"))
  attack <- x$attack
  attack["model"] <- list(x$model)
  .rapid_assess(
    attack, rule, data = .rapid_data(x),
    meta = list(key_vars = x$key_vars, target_var = x$target_var, na = x$na, call = call),
    store_model = !is.null(x$model), store_predictions = TRUE,
    n_perm = n_perm, threshold_rule = threshold_rule, perm_level = perm_level
  )
}


# =============================================================================
# The attack and the rapid object
# =============================================================================

#' Fit the attacker on the synthetic data and predict the original records
#'
#' @param X,Y Prepared original and synthetic data (`.rapid_prepare_data()`).
#' @return A list of class `rapid_attack` with `pred` (probability matrix, or
#'   `list(point, dist)`), `y` (true values), `learner`, `model`, `formula`,
#'   `is_categorical` and `attacker` (the line from `.learner_describe()`).
#' @noRd
.rapid_attack <- function(X, Y, key_vars, target_var, learner, num_dist = "auto",
                          num_scale = "auto", metric = "symmetric", n_quantiles = 100) {
  # The formula lives in the global environment, so that the object does not
  # keep the data of the calling frame alive.
  formula <- stats::as.formula(paste(target_var, "~", paste(key_vars, collapse = " + ")),
                               env = globalenv())
  y <- X[[target_var]]
  is_categorical <- is.factor(y)
  if (is.null(metric)) metric <- "symmetric"

  fit <- .learner_fit(learner, formula, Y, target_var, num_dist, num_scale, metric)
  point <- .learner_predict(learner, fit, X, target_var)
  pred <- if (is_categorical) {
    point
  } else {
    list(point = point,
         dist = .learner_pdist(learner, fit, X, Y, target_var, num_dist, n_quantiles))
  }
  structure(
    list(pred = pred, y = y, learner = learner, model = fit, formula = formula,
         is_categorical = is_categorical, attacker = .learner_describe(learner, fit)),
    class = "rapid_attack"
  )
}

#' Score an attack under a rule and assemble the rapid object
#'
#' The only place where a `rapid` object is put together.
#'
#' @param attack A `rapid_attack`.
#' @param rule A `rapid_rule`.
#' @param data The prepared original data (one row per scored record).
#' @param meta `list(key_vars, target_var, na, call)`.
#' @param n_perm,threshold_rule,perm_level Permutations of the true values for
#'   the chance level, and whether their `perm_level` quantile is the threshold.
#' @noRd
.rapid_assess <- function(attack, rule, data, meta, store_model = FALSE,
                          store_predictions = TRUE, n_perm = 10,
                          threshold_rule = "fixed", perm_level = 0.95) {
  .rapid_check_perm(n_perm, threshold_rule, perm_level, rule)
  # The permutations come first: under a permutation threshold they set the flag
  null <- if (n_perm > 0) .rapid_null(attack, rule, n_perm)
  if (threshold_rule == "permutation") {
    rule <- .rule_set_threshold(
      rule, stats::quantile(null, perm_level, names = FALSE, na.rm = TRUE))
  }
  chance_rates <- if (is.null(null)) numeric(0) else
    apply(null, 2, function(s) mean(.rule_flag(rule, s)))
  chance <- if (n_perm > 0) mean(chance_rates) else NA_real_

  sc <- .rule_score(rule, attack$pred, attack$y, attack$y)
  score <- sc[[rule$score_var]]
  at_risk <- .rule_flag(rule, score)
  n <- length(at_risk)
  rec <- .rapid_records(data, attack, sc, at_risk)

  # The model is kept once and only on request
  model <- attack$model
  attack["model"] <- list(NULL)

  validate_rapid(new_rapid(list(
    rapid = sum(at_risk) / n,
    n_at_risk = sum(at_risk),
    pct_at_risk = 100 * sum(at_risk) / n,
    n_total = n,
    chance = chance,
    chance_rates = chance_rates,
    excess = sum(at_risk) / n - chance,
    n_perm = n_perm,
    threshold_rule = threshold_rule,
    perm_level = perm_level,
    rule = rule,
    method = rule$method,
    threshold = rule$threshold,
    score = score,
    baseline = sc$b,
    at_risk = at_risk,
    records = rec$records,
    record_cols = rec$cols,
    attack = if (store_predictions) attack,
    model = if (store_model) model,
    model_type = attack$learner$name,
    attacker = attack$attacker,
    model_metrics = .compute_rapid_metrics(
      attack$y, if (attack$is_categorical) attack$pred else attack$pred$point),
    formula = attack$formula,
    key_vars = meta$key_vars,
    target_var = meta$target_var,
    is_categorical = attack$is_categorical,
    na = meta$na,
    call = meta$call
  )))
}

#' Check the settings of the permutations
#' @noRd
.rapid_check_perm <- function(n_perm, threshold_rule, perm_level, rule) {
  if (!is.numeric(n_perm) || length(n_perm) != 1 || is.na(n_perm) || n_perm < 0 ||
      n_perm != round(n_perm)) {
    .rapid_abort("'n_perm' must be a single whole number (0 skips the permutations).",
                 "input")
  }
  if (!is.numeric(perm_level) || length(perm_level) != 1 || is.na(perm_level) ||
      perm_level <= 0 || perm_level >= 1) {
    .rapid_abort("'perm_level' must be a single number between 0 and 1.", "input")
  }
  if (threshold_rule == "permutation") {
    if (n_perm < 1) {
      .rapid_abort("A permutation threshold needs 'n_perm' of at least 1.", "input")
    }
    if (!isTRUE(rule$gain)) {
      .rapid_abort(c("A permutation threshold is only available for rules with a gain over ",
                     "a baseline (num_eval_method = \"band\", cat_eval_method = ",
                     "\"RCS_marginal\"); the rule '", rule$method, "' has none."), "input")
    }
  }
  invisible(TRUE)
}

# Record table: the prepared data, the true value (numeric targets), the columns
# of the rule and the flag. A column whose name is taken by a variable gets the
# suffix ".rapid"; `cols` maps the plain names to the names in the table.
.rapid_records <- function(data, attack, sc, at_risk) {
  extra <- c(if (!attack$is_categorical) list(original_value = attack$y),
             as.list(sc), list(at_risk = at_risk))
  plain <- names(extra)
  actual <- ifelse(plain %in% names(data), paste0(plain, ".rapid"), plain)
  names(extra) <- actual
  records <- cbind(data, as.data.frame(extra, stringsAsFactors = FALSE))
  list(records = records, cols = stats::setNames(actual, plain))
}

# A record column by its plain name
.rapid_col <- function(x, name) x$records[[x$record_cols[[name]]]]

.rapid_fields <- c("rapid", "n_at_risk", "pct_at_risk", "n_total", "chance", "chance_rates",
                   "excess", "n_perm", "threshold_rule", "perm_level", "rule", "method",
                   "threshold", "score", "baseline", "at_risk", "records", "record_cols",
                   "attack", "model", "model_type", "attacker", "model_metrics", "formula",
                   "key_vars", "target_var", "is_categorical", "na", "call", "version")

#' Constructor of the rapid class
#' @param fields Named list of the fields of the object.
#' @noRd
new_rapid <- function(fields) {
  fields$version <- utils::packageVersion("riskutility")
  structure(fields, class = "rapid")
}

#' Validator of the rapid class
#' @return The object, if it is consistent; a classed error otherwise.
#' @noRd
validate_rapid <- function(x) {
  missing_fields <- setdiff(.rapid_fields, names(x))
  if (length(missing_fields) > 0) {
    .rapid_abort(c("Not a valid rapid object: missing field(s) ",
                   paste(missing_fields, collapse = ", "), "."), "object")
  }
  if (!inherits(x$rule, "rapid_rule")) {
    .rapid_abort("Not a valid rapid object: 'rule' is not a scoring rule.", "object")
  }
  n <- x$n_total
  sizes <- c(score = length(x$score), at_risk = length(x$at_risk),
             records = if (is.data.frame(x$records)) nrow(x$records) else -1L)
  if (!is.numeric(n) || length(n) != 1 || any(sizes != n)) {
    .rapid_abort(c("Not a valid rapid object: score, at_risk and records must have ",
                   "n_total entries."), "object")
  }
  x
}

#' Check what a method on a rapid object needs
#'
#' @param x A `rapid` object.
#' @param what Any of `"version"` (created by this version of the class),
#'   `"predictions"` (stored attack) and `"model"` (stored model).
#' @noRd
.rapid_require <- function(x, what = "version") {
  if ("version" %in% what && is.null(x$rule)) {
    .rapid_abort(c("This rapid object was created by riskutility <= 0.2.0 and cannot be ",
                   "used with this version; re-run rapid()."), "object")
  }
  if ("predictions" %in% what && is.null(x$attack)) {
    .rapid_abort(c("The object does not hold the attacker's predictions; re-run rapid() ",
                   "with store_predictions = TRUE."), "object")
  }
  if ("model" %in% what && is.null(x$model)) {
    .rapid_abort(c("The object does not hold the fitted model; re-run rapid() with ",
                   "store_model = TRUE."), "object")
  }
  invisible(x)
}

#' Compute model performance metrics
#' @noRd
.compute_rapid_metrics <- function(A, B) {
  if (is.factor(A)) {
    # Classification
    pred_class <- factor(colnames(B)[max.col(as.matrix(B), ties.method = "first")],
                         levels = levels(A))
    return(list(accuracy = mean(pred_class == A, na.rm = TRUE)))
  }
  # Regression
  mae <- mean(abs(A - B), na.rm = TRUE)
  rmse <- sqrt(mean((A - B)^2, na.rm = TRUE))
  rmae <- mae / mean(abs(A), na.rm = TRUE)
  rrmse <- rmse / sd(A, na.rm = TRUE)
  list(mae = mae, rmse = rmse, rmae = rmae, rrmse = rrmse)
}


# =============================================================================
# Internal plot helpers
# =============================================================================

#' Plot QI Variable Importance Barplot
#'
#' @param x A rapid object (must have model stored via store_model=TRUE)
#' @param top_n Maximum number of variables to display
#' @return Invisible NULL
#' @noRd
.plot_qi_importance <- function(x, top_n = 15) {
  .rapid_require(x, c("version", "model"))
  learner <- if (!is.null(x$attack)) x$attack$learner else .rapid_learner(x$model_type)
  imp <- .learner_importance(learner, x$model)

  if (is.null(imp) || length(imp) == 0) {
    message("Variable importance not available for model type: ", x$model_type)
    return(invisible(NULL))
  }

  # Sort and limit
  imp <- sort(imp, decreasing = TRUE)
  if (length(imp) > top_n) {
    imp <- imp[1:top_n]
  }

  # Normalize to 0-1 scale
  imp <- imp / max(imp)

  # Reverse for bottom-to-top display (largest at top)
  imp <- rev(imp)
  n <- length(imp)
  y_pos <- seq_len(n)

  # Adjust margins for long labels
  max_label_len <- max(nchar(names(imp)))
  left_margin <- min(max(8, max_label_len * 0.6), 20)

  old_par <- par(mar = c(5, left_margin, 4, 2))
  on.exit(par(old_par))

  # Lollipop chart: segments from 0 + points
  plot(imp, y_pos, type = "n",
       xlim = c(0, 1.05),
       ylim = c(0.5, n + 0.5),
       yaxt = "n",
       xlab = "Relative Importance",
       ylab = "",
       main = paste0("QI Variable Importance (", x$model_type, ")"))
  segments(0, y_pos, imp, y_pos, col = "gray60", lwd = 1.5)
  points(imp, y_pos, pch = 19, col = "steelblue", cex = 1.3)
  axis(2, at = y_pos, labels = names(imp), las = 1, cex.axis = 0.85)

  invisible(NULL)
}

#' Discretize QI data for attribution GLM
#'
#' Converts numeric and high-cardinality factor variables to quantile bins
#' to prevent quasi-complete separation in logistic regression.
#' @param qi_data data frame of QI columns
#' @param max_levels integer, factors with more levels are rebinned (default 10)
#' @param n_bins integer, number of quantile bins for numeric variables (default 4)
#' @return list with qi_data (discretized), level_map, contrast_list
#' @noRd
.prepare_attribution_data <- function(qi_data, at_risk = NULL, max_levels = 10,
                                      n_bins = 4) {
  # Adaptive binning: reduce bins when at-risk count is low to prevent separation
  # Aim for >= 5 expected at-risk records per bin
  if (!is.null(at_risk)) {
    n_at_risk <- sum(at_risk)
    n_minority <- min(n_at_risk, length(at_risk) - n_at_risk)
    n_bins <- max(2L, min(n_bins, as.integer(floor(n_minority / 5))))
  }

  level_map <- list()
  contrast_list <- list()

  for (qi in names(qi_data)) {
    col <- qi_data[[qi]]

    if (is.numeric(col)) {
      # Discretize numeric into quantile bins
      breaks <- unique(stats::quantile(col, probs = seq(0, 1, length.out = n_bins + 1),
                                       na.rm = TRUE))
      if (length(breaks) <= 2) {
        qi_data[[qi]] <- as.factor(col)
      } else {
        qi_data[[qi]] <- cut(col, breaks = breaks, include.lowest = TRUE,
                             dig.lab = 3)
      }
    } else {
      qi_data[[qi]] <- as.factor(qi_data[[qi]])
    }

    # Rebin factors with too many levels
    if (is.factor(qi_data[[qi]]) && nlevels(qi_data[[qi]]) > max_levels) {
      col_num <- as.numeric(qi_data[[qi]])
      reb <- max(2L, min(n_bins, nlevels(qi_data[[qi]])))
      breaks <- unique(stats::quantile(col_num,
                                       probs = seq(0, 1, length.out = reb + 1),
                                       na.rm = TRUE))
      qi_data[[qi]] <- cut(col_num, breaks = breaks, include.lowest = TRUE,
                           labels = paste0("G", seq_len(length(breaks) - 1)))
    }

    # Set up sum-to-zero contrasts and label mapping
    lvls <- levels(qi_data[[qi]])
    n_levels <- length(lvls)
    if (n_levels > 1) {
      contrast_list[[qi]] <- stats::contr.sum(n_levels)
      for (i in seq_len(n_levels - 1)) {
        coef_name <- paste0(qi, i)
        level_map[[coef_name]] <- paste0(qi, ": ", lvls[i])
      }
      level_map[[paste0(qi, "_ref")]] <- paste0(qi, ": ", lvls[n_levels])
    }
  }

  list(qi_data = qi_data, level_map = level_map, contrast_list = contrast_list)
}

#' Fit attribution GLM and extract coefficients with reference levels
#' @noRd
.fit_attribution_glm <- function(formula, qi_data, at_risk, contrast_list, level_map) {
  model_data <- cbind(at_risk = as.numeric(at_risk), qi_data)

  fit <- tryCatch({
    suppressWarnings(
      stats::glm(formula, data = model_data, family = stats::binomial(),
                 contrasts = contrast_list)
    )
  }, error = function(e) {
    warning("Could not fit attribution model: ", e$message, call. = FALSE)
    return(NULL)
  })

  if (is.null(fit)) return(NULL)

  coefs <- summary(fit)$coefficients
  V <- stats::vcov(fit)
  coef_df <- data.frame(
    term = rownames(coefs),
    estimate = coefs[, "Estimate"],
    se = coefs[, "Std. Error"]
  )

  # Remove intercept
  coef_df <- coef_df[coef_df$term != "(Intercept)", , drop = FALSE]
  if (nrow(coef_df) == 0) return(NULL)

  # Add implied reference-level coefficients (sum-to-zero: beta_ref = -sum(others))
  for (qi in names(contrast_list)) {
    lvls <- levels(qi_data[[qi]])
    n_levels <- length(lvls)
    if (n_levels < 2) next

    est_terms <- paste0(qi, seq_len(n_levels - 1))
    est_terms <- est_terms[est_terms %in% coef_df$term]
    if (length(est_terms) == 0) next

    ref_est <- -sum(coef_df$estimate[coef_df$term %in% est_terms])

    V_sub <- V[est_terms, est_terms, drop = FALSE]
    ref_var <- sum(V_sub)
    ref_se <- if (ref_var > 0) sqrt(ref_var) else NA_real_

    ref_label <- level_map[[paste0(qi, "_ref")]]
    if (is.null(ref_label)) ref_label <- paste0(qi, ": ", lvls[n_levels])

    coef_df <- rbind(coef_df, data.frame(
      term = paste0(qi, "_ref"),
      estimate = ref_est,
      se = ref_se
    ))
  }

  # Improve term labels
  coef_df$label <- vapply(coef_df$term, function(t) {
    if (t %in% names(level_map)) return(level_map[[t]])
    t_clean <- t
    for (qi in names(qi_data)) {
      if (is.factor(qi_data[[qi]])) {
        lvls <- levels(qi_data[[qi]])
        for (i in seq_along(lvls)) {
          t_clean <- gsub(paste0("\\b", qi, i, "\\b"), lvls[i], t_clean)
        }
      }
    }
    gsub(":", " x ", t_clean)
  }, character(1))

  coef_df
}

#' Draw forest-plot style dot chart from coefficient data frame
#' @noRd
.draw_forest_plot <- function(coef_df, main = "QI Attribution (Effect Coding)",
                              xlab = "Log-odds contribution to at-risk") {
  n <- nrow(coef_df)
  y_pos <- seq_len(n)

  colors <- ifelse(coef_df$estimate > 0, "firebrick", "steelblue")

  valid_se <- is.finite(coef_df$se)
  ci_lo <- coef_df$estimate - 1.96 * coef_df$se
  ci_hi <- coef_df$estimate + 1.96 * coef_df$se

  # Base x-axis range on point estimates; include CIs only if they are reasonable
  # (SE < 10 is a practical threshold; larger indicates quasi-separation)
  reasonable <- valid_se & coef_df$se < 10
  pts_range <- range(coef_df$estimate, na.rm = TRUE)
  if (any(reasonable)) {
    x_range <- range(c(pts_range, ci_lo[reasonable], ci_hi[reasonable]))
  } else {
    x_range <- pts_range
  }
  x_pad <- max(diff(x_range) * 0.15, 0.05)
  x_lim <- c(x_range[1] - x_pad, x_range[2] + x_pad)

  max_label_len <- max(nchar(coef_df$label))
  left_margin <- min(max(8, max_label_len * 0.5), 20)

  old_par <- par(mar = c(5, left_margin, 4, 2))
  on.exit(par(old_par))

  plot(coef_df$estimate, y_pos, type = "n",
       xlim = x_lim, ylim = c(0.5, n + 0.5),
       yaxt = "n", xlab = xlab, ylab = "", main = main)
  abline(v = 0, lty = 2, col = "gray40")

  # Draw CIs: only for reasonable SEs, truncated to plot limits
  if (any(reasonable)) {
    ci_lo_draw <- pmax(ci_lo[reasonable], x_lim[1])
    ci_hi_draw <- pmin(ci_hi[reasonable], x_lim[2])
    ci_width <- abs(ci_hi_draw - ci_lo_draw)
    draw <- ci_width > .Machine$double.eps * 100
    if (any(draw)) {
      suppressWarnings(
        arrows(ci_lo_draw[draw], y_pos[reasonable][draw],
               ci_hi_draw[draw], y_pos[reasonable][draw],
               angle = 90, code = 3, length = 0.03, col = "gray50")
      )
    }
  }

  points(coef_df$estimate, y_pos, pch = 19, col = colors, cex = 1.2)
  axis(2, at = y_pos, labels = coef_df$label, las = 1, cex.axis = 0.8)
}

#' Plot QI Attribution Analysis using Effect Coding
#'
#' @param x A rapid object
#' @param formula Optional formula for the attribution model
#' @param top_n Maximum number of terms to display
#' @return Invisible NULL
#' @noRd
.plot_qi_attribution <- function(x, formula = NULL, top_n = 20) {
  records <- x$records

  if (is.null(records) || !("at_risk" %in% names(records))) {
    .rapid_abort("Records with at_risk column required.", "object")
  }

  if (all(records$at_risk) || !any(records$at_risk)) {
    message("Cannot fit attribution model: all records have the same at_risk status.")
    return(invisible(NULL))
  }

  key_vars <- x$key_vars

  # Discretize numeric / high-cardinality variables
  prep <- .prepare_attribution_data(records[, key_vars, drop = FALSE],
                                     at_risk = records$at_risk)

  # Build formula: main effects only
  if (is.null(formula)) {
    formula <- stats::as.formula(paste("at_risk ~", paste(key_vars, collapse = " + ")))
  }

  coef_df <- .fit_attribution_glm(formula, prep$qi_data, records$at_risk,
                                  prep$contrast_list, prep$level_map)
  if (is.null(coef_df) || nrow(coef_df) == 0) {
    message("No coefficients to display.")
    return(invisible(NULL))
  }

  # Limit and sort
  if (nrow(coef_df) > top_n) {
    coef_df <- coef_df[order(abs(coef_df$estimate), decreasing = TRUE), ]
    coef_df <- coef_df[1:top_n, ]
  }
  coef_df <- coef_df[order(coef_df$estimate), ]

  .draw_forest_plot(coef_df, main = "QI Attribution (Effect Coding)")

  invisible(NULL)
}

#' Plot QI Interaction Effects
#'
#' Fits a logistic regression with all pairwise interactions between QI variables
#' and shows only the interaction coefficients.
#' @param x A rapid object
#' @param top_n Maximum number of interaction terms to display
#' @return Invisible NULL
#' @noRd
.plot_qi_interactions <- function(x, top_n = 20) {
  records <- x$records

  if (is.null(records) || !("at_risk" %in% names(records))) {
    .rapid_abort("Records with at_risk column required.", "object")
  }

  if (all(records$at_risk) || !any(records$at_risk)) {
    message("Cannot fit interaction model: all records have the same at_risk status.")
    return(invisible(NULL))
  }

  key_vars <- x$key_vars
  if (length(key_vars) < 2) {
    message("Interaction plot requires at least 2 key variables.")
    return(invisible(NULL))
  }

  # Discretize
  prep <- .prepare_attribution_data(records[, key_vars, drop = FALSE],
                                     at_risk = records$at_risk)

  # Formula with all pairwise interactions
  pairs <- utils::combn(key_vars, 2, simplify = FALSE)
  interaction_terms <- vapply(pairs, function(p) paste(p, collapse = ":"), character(1))
  formula <- stats::as.formula(
    paste("at_risk ~", paste(c(key_vars, interaction_terms), collapse = " + "))
  )

  coef_df <- .fit_attribution_glm(formula, prep$qi_data, records$at_risk,
                                  prep$contrast_list, prep$level_map)
  if (is.null(coef_df) || nrow(coef_df) == 0) {
    message("No coefficients to display.")
    return(invisible(NULL))
  }

  # Keep only interaction terms (those containing " x " in the label, or ":" in the term)
  is_interaction <- grepl(":", coef_df$term)
  coef_df <- coef_df[is_interaction, , drop = FALSE]

  if (nrow(coef_df) == 0) {
    message("No interaction coefficients to display.")
    return(invisible(NULL))
  }

  if (nrow(coef_df) > top_n) {
    coef_df <- coef_df[order(abs(coef_df$estimate), decreasing = TRUE), ]
    coef_df <- coef_df[1:top_n, ]
  }
  coef_df <- coef_df[order(coef_df$estimate), ]

  .draw_forest_plot(coef_df, main = "QI Interaction Effects (Effect Coding)")

  invisible(NULL)
}

#' Plot Threshold Sensitivity Curve
#'
#' @param x A rapid object
#' @param tau_range Thresholds for rules that flag scores above the threshold
#' @param epsilon_range Tolerances for the point rule
#' @return Invisible NULL
#' @noRd
.plot_threshold_sensitivity <- function(x, tau_range = seq(0, 1, by = 0.05),
                                        epsilon_range = seq(0, 1, by = 0.05)) {
  scores <- x$score[is.finite(x$score)]

  if (length(scores) == 0) {
    message("No records available for threshold sensitivity plot.")
    return(invisible(NULL))
  }

  if (x$rule$direction == "below") {
    # The score is an error; a record is at risk below the tolerance
    thresholds <- epsilon_range
    if (identical(epsilon_range, seq(0, 1, by = 0.05)) && max(scores) > 1) {
      # Errors in the units of the target (absolute metric)
      thresholds <- seq(0, max(scores), length.out = 21)
    }
    rapid_values <- vapply(thresholds, function(eps) mean(scores < eps), numeric(1))
    xlab_text <- "Tolerance (epsilon)"
    label <- "epsilon"
  } else {
    thresholds <- tau_range
    rapid_values <- vapply(thresholds, function(tau) mean(scores > tau), numeric(1))
    xlab_text <- "Threshold (tau)"
    label <- "tau"
  }

  plot(thresholds, rapid_values,
       type = "l", lwd = 2, col = "blue",
       xlab = xlab_text,
       ylab = "RAPID (Proportion at Risk)",
       main = "Threshold Sensitivity Curve",
       ylim = c(0, 1))

  points(x$threshold, x$rapid, pch = 19, col = "red", cex = 1.5)
  text(x$threshold, x$rapid,
       paste0(label, " = ", round(x$threshold, 3), "\nRAPID = ", round(x$rapid, 3)),
       pos = 4, col = "red")
  graphics::grid()

  invisible(NULL)
}

#' Plot method for rapid objects
#'
#' Creates visualizations for RAPID disclosure risk assessment results.
#' Six plot types are available via the \code{which} parameter.
#'
#' @param x an object of class "rapid"
#' @param y not used
#' @param ... additional arguments passed to plotting functions
#' @param which integer indicating which plot to show:
#'   \itemize{
#'     \item \code{1}: Risk score distribution histogram (default)
#'     \item \code{2}: Prediction scatter plot (numeric targets only)
#'     \item \code{3}: Threshold sensitivity curve (RAPID vs threshold range)
#'     \item \code{4}: QI variable importance lollipop chart (requires \code{store_model = TRUE})
#'     \item \code{5}: QI attribution main effects (sum-to-zero logistic regression)
#'     \item \code{6}: QI interaction effects (pairwise interactions, requires >= 2 key variables)
#'   }
#' @param xlim character or numeric vector controlling x-axis limits for risk distribution
#'   (plots 1 and 2 only):
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
#' @param tau_range numeric vector of tau thresholds for the sensitivity curve
#'   (plot 3, rules that flag scores above the threshold). Default
#'   \code{seq(0, 1, by = 0.05)}.
#' @param epsilon_range numeric vector of epsilon thresholds for the sensitivity curve
#'   (plot 3, point rule). Default \code{seq(0, 1, by = 0.05)}; for errors in
#'   the units of the target (absolute metric) the range of the errors is used.
#' @param formula optional formula for QI attribution plot (plot 5), e.g.,
#'   \code{at_risk ~ age * gender}. If NULL, uses main effects of all key variables.
#' @param top_n integer, maximum number of variables/terms to show in importance (plot 4),
#'   attribution (plot 5), and interaction (plot 6) plots. Default 15.
#'
#' @details
#' \strong{Plot 1 -- Risk Distribution:}
#' Histogram of risk scores colored by at-risk status. The red dashed line indicates
#' the threshold. Higher scores indicate higher risk, except under the point rule
#' for numeric targets, where lower prediction errors indicate higher risk
#' (accurate predictions imply disclosure).
#'
#' \strong{Plot 2 -- Prediction Scatter (numeric only):}
#' Scatter plot of predicted vs actual values, colored by at-risk status.
#'
#' \strong{Plot 3 -- Threshold Sensitivity Curve:}
#' Shows how RAPID (proportion at risk) changes across a range of threshold values.
#' The current threshold is marked with a red point. Useful for understanding how
#' sensitive the risk assessment is to threshold choice.
#'
#' \strong{Plot 4 -- QI Variable Importance:}
#' Lollipop chart of relative variable importance from the attacker model.
#' Requires the model to be stored via \code{rapid(..., store_model = TRUE)}.
#' Supports random forest, CART, XGBoost, and linear/logistic models.
#'
#' \strong{Plot 5 -- QI Attribution (Main Effects):}
#' Sum-to-zero contrast logistic regression on at-risk status. Shows which
#' quasi-identifier values contribute most to disclosure risk. Numeric variables
#' are discretized into quartile bins. All factor levels are shown including the
#' implied reference level. Positive coefficients (firebrick) increase risk;
#' negative coefficients (steelblue) decrease risk. Whiskers show 95\% CI.
#'
#' \strong{Plot 6 -- QI Interaction Effects:}
#' Same logistic regression framework as plot 5, but fits all pairwise interactions
#' and displays only the interaction coefficients. Requires at least 2 key variables.
#'
#' @return For faceted plots, returns a ggplot2 object (invisibly). For base plots,
#'   returns NULL invisibly.
#'
#' @examples
#' \donttest{
#' set.seed(42)
#' n <- 500
#' original <- data.frame(
#'   age = sample(20:60, n, replace = TRUE),
#'   gender = factor(sample(c("M", "F"), n, replace = TRUE)),
#'   education = factor(sample(c("High", "Medium", "Low"), n, replace = TRUE))
#' )
#' original$health <- factor(ifelse(
#'   original$age > 40 & original$gender == "M",
#'   sample(c("Poor", "Fair", "Good"), n, replace = TRUE, prob = c(0.6, 0.3, 0.1)),
#'   sample(c("Poor", "Fair", "Good"), n, replace = TRUE, prob = c(0.1, 0.3, 0.6))
#' ))
#' synthetic <- original
#' synthetic$health <- original$health
#' idx <- sample(n, n * 0.3)
#' synthetic$health[idx] <- factor(
#'   sample(c("Good", "Fair", "Poor"), length(idx), replace = TRUE)
#' )
#'
#' result <- rapid(original, synthetic,
#'                 key_vars = c("age", "gender", "education"),
#'                 target_var = "health",
#'                 model_type = "rf",
#'                 store_model = TRUE)
#'
#' # Risk distribution histogram (default)
#' plot(result)
#'
#' # Threshold sensitivity curve
#' plot(result, which = 3)
#'
#' # QI variable importance
#' plot(result, which = 4)
#'
#' # QI attribution (main effects)
#' plot(result, which = 5)
#'
#' # QI interaction effects
#' plot(result, which = 6)
#'
#' # Facet by gender
#' plot(result, facet_by = "gender")
#' }
#'
#' @export
#' @importFrom graphics hist abline par legend points text mtext arrows grid segments axis
#' @importFrom stats glm binomial contr.sum as.formula coef vcov quantile
#' @importFrom utils combn modifyList
#' @importFrom ggplot2 ggplot aes geom_histogram geom_vline annotate facet_wrap
#'   labs theme_minimal theme element_text scale_fill_manual after_stat
plot.rapid <- function(x, y = NULL, ...,
                       which = 1,
                       xlim = "full",
                       annotate = TRUE,
                       facet_by = NULL,
                       bins = 30,
                       tau_range = seq(0, 1, by = 0.05),
                       epsilon_range = seq(0, 1, by = 0.05),
                       formula = NULL,
                       top_n = 15) {

  .rapid_require(x, "version")

  # Dispatch to specialized plot helpers for which = 3, 4, 5, 6
  if (3 %in% which) {
    return(.plot_threshold_sensitivity(x, tau_range = tau_range,
                                       epsilon_range = epsilon_range))
  }
  if (4 %in% which) {
    return(.plot_qi_importance(x, top_n = top_n))
  }
  if (5 %in% which) {
    return(.plot_qi_attribution(x, formula = formula, top_n = top_n))
  }
  if (6 %in% which) {
    return(.plot_qi_interactions(x, top_n = top_n))
  }

  show <- rep(FALSE, 2)
  show[which] <- TRUE

  records <- x$records

  # Score column and labels come from the rule
  score_var <- x$record_cols[[x$rule$score_var]]
  x_label <- x$rule$score_label
  main_title <- paste(x_label, "Distribution")
  # TRUE: a record is at risk below the threshold (the score is an error)
  below <- x$rule$direction == "below"

  keep <- is.finite(x$score)
  if (!any(keep)) {
    message("No finite scores to plot.")
    return(invisible(NULL))
  }
  if (!all(keep)) records <- records[keep, , drop = FALSE]
  scores <- x$score[keep]

  # Calculate x-axis limits
  score_range <- range(scores)
  full_range <- if (below) c(0, score_range[2]) else score_range

  if (is.character(xlim)) {
    xlim_values <- switch(xlim,
      "full" = full_range,
      "zero_to_one" = c(0, 1),
      "threshold_to_one" = c(x$threshold, 1),
      full_range  # default to full
    )
  } else if (is.numeric(xlim) && length(xlim) == 2) {
    xlim_values <- xlim
  } else {
    xlim_values <- full_range
  }

  pct_at_risk <- x$pct_at_risk
  n_at_risk <- x$n_at_risk
  n_total <- x$n_total
  risk_labels <- paste0(c("At Risk (", "Not at Risk ("),
                        if (below) c("< ", ">= ") else c("> ", "<= "),
                        format(x$threshold), ")")

  # Use ggplot2 for faceted plots
  if (!is.null(facet_by) && show[1]) {
    # Check facet variables exist
    missing_vars <- setdiff(facet_by, names(records))
    if (length(missing_vars) > 0) {
      .rapid_abort(c("Facet variables not found in records: ",
                     paste(missing_vars, collapse = ", "), "."), "input")
    }

    # Create facet variable (combine if multiple)
    if (length(facet_by) == 1) {
      records$facet_var <- records[[facet_by]]
    } else {
      records$facet_var <- interaction(records[, facet_by, drop = FALSE], sep = " / ")
    }
    records$at_risk <- x$at_risk[keep]

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
        labels = c("FALSE" = "Not at Risk", "TRUE" = "At Risk"),
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
    if (break_max <= break_min) break_max <- break_min + 1
    breaks_seq <- seq(break_min, break_max, length.out = bins + 1)

    # Compute histogram
    h <- hist(scores, breaks = breaks_seq, plot = FALSE)

    # Color bars based on whether they're in the at-risk region
    bar_colors <- if (below) {
      # Low error means an accurate prediction, hence disclosure risk
      ifelse(h$mids < x$threshold, "coral", "steelblue")
    } else {
      ifelse(h$mids >= x$threshold, "coral", "steelblue")
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

    legend(if (below) "topright" else "topleft",
           legend = c(risk_labels, "Threshold"),
           fill = c("coral", "steelblue", NA),
           border = c("white", "white", NA),
           col = c(NA, NA, "red"),
           lty = c(NA, NA, 2),
           lwd = c(NA, NA, 2),
           cex = 0.8)
  }

  if (show[2] && !x$is_categorical) {
    # Scatter plot of predicted vs actual with transparency
    # Plot not-at-risk points first (background), then at-risk points on top
    col_at_risk <- adjustcolor("coral", alpha.f = 0.6)
    col_not_risk <- adjustcolor("steelblue", alpha.f = 0.6)
    actual <- .rapid_col(x, "original_value")
    predicted <- .rapid_col(x, "predicted_value")

    # Initialize empty plot
    plot(range(actual), range(predicted),
         type = "n",
         main = "Predicted vs Actual Values",
         xlab = "Actual Value", ylab = "Predicted Value")

    # Plot not-at-risk points first (in background)
    if (any(!x$at_risk)) {
      points(actual[!x$at_risk], predicted[!x$at_risk],
             col = col_not_risk, pch = 19, cex = 0.7)
    }

    # Plot at-risk points on top
    if (any(x$at_risk)) {
      points(actual[x$at_risk], predicted[x$at_risk],
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
           legend = c("At Risk", "Not at Risk", "Perfect Prediction"),
           col = c(col_at_risk, col_not_risk, "darkgray"),
           pch = c(19, 19, NA), lty = c(NA, NA, 2), lwd = c(NA, NA, 2),
           cex = 0.8)
  }

  invisible(NULL)
}
