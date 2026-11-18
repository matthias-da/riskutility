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
#'   boosting), "logit" (logistic regression, two classes only), or "lookup" (one
#'   fully grown tree that reads off the synthetic records sharing the key
#'   values; read its excess over the chance level, not RAPID itself). A
#'   user-defined attacker built with \code{\link{rapid_learner}()} is accepted
#'   as well.
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
#'   "grid" uses \code{n_quantiles} quantiles for all models; "exact" uses the
#'   exact probability also for "rf" (the training values in the terminal nodes
#'   of the forest) and is not available for "gbm"
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
#' @param holdout data frame of original records that were not used to generate
#'   the synthetic data (default: NULL; the holdout of a \code{synth_pair} is
#'   used when not given). It is prepared like \code{X} and scored with the same
#'   attacker and rule (a permutation threshold as found on \code{X}), with
#'   baselines from the values of \code{X}. RAPID on \code{X} minus RAPID on
#'   the holdout is the risk specific to this release.
#' @param store_model logical, if TRUE store the fitted model object and the
#'   prepared synthetic data in the result (needed for \code{predict()} and for
#'   the variable importance plot \code{plot(result, which = 4)}). Default FALSE.
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
#'   \item predictor: what \code{predict()} needs besides the model: the
#'     attacker, the prepared synthetic data and the distribution settings
#'     (only if store_model = TRUE)
#'   \item holdout: \code{list(rapid, chance, excess, n, records)} of the
#'     holdout records (NULL without holdout)
#'   \item release_specific: rapid - holdout$rapid (NULL without holdout)
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
                          model_type = c("rf", "lm", "cart", "gbm", "logit", "lookup"),
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
                          # Records the attacker has not seen
                          holdout = NULL,
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
  if (!inherits(model_type, "rapid_learner")) {
    model_type <- .rapid_choice(model_type, .rapid_learner_names, "model_type")
  }
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

  # Holdout: an explicit argument wins over the one of a synth_pair
  if (is.null(holdout)) holdout <- inputs$holdout
  if (!is.null(holdout)) {
    holdout <- .rapid_prepare_new(holdout, prepared$Y, inputs$key_vars, inputs$target_var,
                                  list(strategy = na_strategy, constant = na_constant_value,
                                       value = prepared$na_value),
                                  learner, "holdout")$data
  }

  # Fit on the synthetic data, predict the original (and holdout) records
  if (verbose) message("  Fitting ", learner$name, " model on synthetic data...")
  attack <- .rapid_attack(prepared$X, prepared$Y, inputs$key_vars, inputs$target_var,
                          learner, num_dist, num_scale, scoring$rule$metric, n_quantiles,
                          holdout = holdout)

  if (verbose) message("  Evaluating disclosure risk...")
  # The permutations start from the seed as well, so that update() without a
  # refit draws the same permutations as a fresh call
  if (!is.null(seed)) set.seed(seed)
  .rapid_assess(
    attack, scoring$rule, data = prepared$X,
    meta = list(key_vars = inputs$key_vars, target_var = inputs$target_var,
                na = list(strategy = na_strategy, constant = na_constant_value,
                          value = prepared$na_value, dropped = prepared$dropped),
                call = call),
    store_model = isTRUE(store_model), store_predictions = isTRUE(store_predictions),
    n_perm = scoring$n_perm, threshold_rule = scoring$threshold_rule,
    perm_level = scoring$perm_level, seed = seed
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
.rapid_data <- function(x, records = x$records) {
  records[, seq_len(ncol(records) - length(x$record_cols)), drop = FALSE]
}

#' Score the stored attack of a rapid object under a new rule
#'
#' @param x A `rapid` object with stored predictions.
#' @param rule A `rapid_rule`.
#' @param n_perm,threshold_rule,perm_level Permutation settings.
#' @param call The call stored in the result.
#' @param seed Seed set before the permutations of a holdout.
#' @noRd
.rapid_rescore <- function(x, rule, n_perm, threshold_rule, perm_level, call, seed = NULL) {
  .rapid_require(x, c("version", "predictions"))
  attack <- x$attack
  attack["model"] <- list(x$model)
  attack["train"] <- list(x$predictor$train)
  if (!is.null(attack$holdout)) {
    attack$holdout$data <- .rapid_data(x, x$holdout$records)
  }
  .rapid_assess(
    attack, rule, data = .rapid_data(x),
    meta = list(key_vars = x$key_vars, target_var = x$target_var, na = x$na, call = call),
    store_model = !is.null(x$model), store_predictions = TRUE,
    n_perm = n_perm, threshold_rule = threshold_rule, perm_level = perm_level, seed = seed
  )
}


# =============================================================================
# The attack and the rapid object
# =============================================================================

#' Fit the attacker on the synthetic data and predict the original records
#'
#' @param X,Y Prepared original and synthetic data (`.rapid_prepare_data()`).
#' @param holdout Prepared holdout records (`.rapid_prepare_new()`) or `NULL`.
#' @return A list of class `rapid_attack` with `pred` (probability matrix, or
#'   `list(point, dist)`), `y` (true values), `learner`, `model`, `formula`,
#'   `is_categorical`, `attacker` (the line from `.learner_describe()`),
#'   `train` (the synthetic data), `num_dist`, `n_quantiles` and `holdout`
#'   (`list(pred, y, data)` or `NULL`).
#' @noRd
.rapid_attack <- function(X, Y, key_vars, target_var, learner, num_dist = "auto",
                          num_scale = "auto", metric = "symmetric", n_quantiles = 100,
                          holdout = NULL) {
  # The formula lives in the global environment, so that the object does not
  # keep the data of the calling frame alive.
  formula <- stats::as.formula(paste(target_var, "~", paste(key_vars, collapse = " + ")),
                               env = globalenv())
  y <- X[[target_var]]
  is_categorical <- is.factor(y)
  if (is.null(metric)) metric <- "symmetric"

  fit <- .learner_fit(learner, formula, Y, target_var, num_dist, num_scale, metric)
  predict_on <- function(data) {
    .rapid_predict(learner, fit, data, Y, target_var, is_categorical, num_dist, n_quantiles)
  }
  structure(
    list(pred = predict_on(X), y = y, learner = learner, model = fit, formula = formula,
         is_categorical = is_categorical, attacker = .learner_describe(learner, fit),
         note = learner$note, train = Y, num_dist = num_dist, n_quantiles = n_quantiles,
         holdout = if (!is.null(holdout)) {
           list(pred = predict_on(holdout), y = holdout[[target_var]], data = holdout)
         }),
    class = "rapid_attack"
  )
}

#' The attacker's prediction for prepared records
#'
#' @param point_only If TRUE, a numeric target gets only the point prediction.
#' @return A probability matrix (categorical target), `list(point, dist)`
#'   (numeric target), or the point predictions.
#' @noRd
.rapid_predict <- function(learner, fit, data, train, target_var, is_categorical,
                           num_dist, n_quantiles, point_only = FALSE) {
  point <- .learner_predict(learner, fit, data, target_var)
  if (is_categorical || point_only) return(point)
  list(point = point,
       dist = .learner_pdist(learner, fit, data, train, target_var, num_dist, n_quantiles))
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
#' @param seed Seed set again before the permutations of the holdout, so that a
#'   holdout equal to the original data gives the main result.
#' @noRd
.rapid_assess <- function(attack, rule, data, meta, store_model = FALSE,
                          store_predictions = TRUE, n_perm = 10,
                          threshold_rule = "fixed", perm_level = 0.95, seed = NULL) {
  .rapid_check_perm(n_perm, threshold_rule, perm_level, rule)
  # The permutations come first: under a permutation threshold they set the flag
  null <- if (n_perm > 0) .rapid_null(attack, rule, n_perm)
  if (threshold_rule == "permutation") {
    rule <- .rule_set_threshold(
      rule, stats::quantile(null, perm_level, names = FALSE, na.rm = TRUE))
  }
  chance_rates <- .rapid_chance_rates(rule, null)
  chance <- if (n_perm > 0) mean(chance_rates) else NA_real_

  main <- .rapid_score_records(rule, attack$pred, attack$y, attack$y, data,
                               attack$is_categorical)
  sc <- main$sc
  score <- main$score
  at_risk <- main$at_risk
  n <- length(at_risk)
  rec <- main$rec

  # The holdout is scored under the rule as calibrated on the original records
  holdout <- NULL
  if (!is.null(attack$holdout)) {
    if (!is.null(seed)) set.seed(seed)
    holdout <- .rapid_score_holdout(attack, rule, n_perm)
    attack$holdout$data <- NULL
  }

  # The model and the training data are kept once and only on request
  model <- attack$model
  predictor <- if (store_model) {
    list(learner = attack$learner, train = attack$train, num_dist = attack$num_dist,
         n_quantiles = attack$n_quantiles)
  }
  attack["model"] <- list(NULL)
  attack["train"] <- list(NULL)

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
    predictor = predictor,
    holdout = holdout,
    release_specific = if (!is.null(holdout)) sum(at_risk) / n - holdout$rapid,
    seed = seed,
    model_type = attack$learner$name,
    attacker = attack$attacker,
    note = attack$note,
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

#' Share flagged in each permutation
#' @param null Matrix from `.rapid_null()`, or `NULL` (no permutations).
#' @noRd
.rapid_chance_rates <- function(rule, null) {
  if (is.null(null)) return(numeric(0))
  apply(null, 2, function(s) mean(.rule_flag(rule, s)))
}

#' Score records under a rule and build their record table
#'
#' Shared by the original records, the holdout and `predict()`.
#'
#' @param pred,y The attacker's prediction and the true values.
#' @param ref_y The original values the baseline is taken from.
#' @param data The prepared records.
#' @return `list(sc, score, at_risk, rec)`; `rec` from `.rapid_records()`.
#' @noRd
.rapid_score_records <- function(rule, pred, y, ref_y, data, is_categorical) {
  sc <- .rule_score(rule, pred, y, ref_y)
  score <- sc[[rule$score_var]]
  at_risk <- .rule_flag(rule, score)
  rec <- .rapid_records(data, list(is_categorical = is_categorical, y = y), sc, at_risk)
  list(sc = sc, score = score, at_risk = at_risk, rec = rec)
}

#' Score the holdout of an attack
#'
#' The holdout is scored with the rule of the main result, its baseline and its
#' permutations use the original values as reference.
#' @return `list(rapid, chance, excess, n, records)`.
#' @noRd
.rapid_score_holdout <- function(attack, rule, n_perm) {
  ho <- attack$holdout
  null <- if (n_perm > 0) .rapid_null(ho, rule, n_perm, ref_y = attack$y)
  chance <- if (n_perm > 0) mean(.rapid_chance_rates(rule, null)) else NA_real_
  s <- .rapid_score_records(rule, ho$pred, ho$y, attack$y, ho$data, attack$is_categorical)
  n <- length(s$at_risk)
  share <- sum(s$at_risk) / n
  list(rapid = share, chance = chance, excess = share - chance, n = n,
       records = s$rec$records)
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
                   "attack", "model", "predictor", "holdout", "release_specific",
                   "model_type", "attacker", "model_metrics", "formula", "key_vars",
                   "target_var", "is_categorical", "na", "call", "version")

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
#'   `"predictions"` (stored attack), `"model"` (stored model) and `"predictor"` (the
#'   learner and training data `predict()` needs; stored with the model).
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
  if (("model" %in% what && is.null(x$model)) ||
      ("predictor" %in% what && is.null(x$predictor))) {
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
