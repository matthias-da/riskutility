#' RAPID: Risk of Attribute Prediction-Induced Disclosure
#'
#' Estimates the share of original records whose sensitive value an attacker can
#' infer from a released (synthetic or anonymized) data set. The attacker fits a
#' model on the release that predicts the sensitive variable from the
#' quasi-identifiers, and applies it to the quasi-identifiers of the original
#' records. A record is at risk when the attacker puts clearly more probability on
#' its true value than the marginal distribution of the sensitive variable already
#' does. RAPID is the share of such records. Every result also reports the chance
#' level, the share of records flagged when the same predictions are paired with
#' randomly reassigned true values, and the excess of RAPID over it.
#'
#' @param X the original data (a data frame), or a \code{\link{synth_pair}}
#'   that holds the original and the synthetic data, the key variables, the
#'   target and optionally a holdout
#' @param Y data frame of the released (synthetic) data the attacker is fitted
#'   on (not needed if \code{X} is a \code{synth_pair})
#' @param key_vars character vector of the quasi-identifiers, the variables the
#'   attacker knows about the records
#' @param target_var character, name of the sensitive variable. A factor,
#'   character or logical variable is a categorical target, a numeric variable a
#'   continuous one.
#' @param model_type the attacker: \code{"rf"} (random forest from ranger,
#'   default), \code{"lm"} (linear model, continuous targets only),
#'   \code{"cart"} (classification or regression tree from rpart, pruned by
#'   cross-validation), \code{"gbm"} (gradient boosting from xgboost),
#'   \code{"logit"} (logistic regression, two classes only) or \code{"lookup"}
#'   (one fully grown tree that reads off the synthetic records with the same
#'   key values, read by its excess). A user-defined attacker built with
#'   \code{\link{rapid_learner}()} is accepted as well. See the sections
#'   "Predictive distribution per attacker" and "Which attacker answers which
#'   question".
#' @param num_eval_method character, rule for a continuous target:
#'   \code{"band"} (default) or \code{"point"}, the rule of riskutility 0.2.0.
#'   See the section "The measure".
#' @param num_error_metric character, error metric that defines the tolerance
#'   band: \code{"symmetric"} (default), \code{"stabilised_relative"} or
#'   \code{"absolute"}
#' @param num_epsilon numeric, the tolerance. For the two relative metrics it is
#'   a proportion (default 0.1, that is 10 percent). For the absolute metric it
#'   is given in the units of the target and has no default. A relative
#'   tolerance of 1 or more is read as percent, with a warning of class
#'   \code{"rapid_warning_legacy_tolerance"}.
#' @param num_delta numeric, the smoothing constant \eqn{\delta}{delta} of the
#'   relative metrics (default 0.01)
#' @param num_tau numeric, threshold \eqn{\tau}{tau} on the normalized gain of a
#'   continuous target under the band rule (default 0.3). The point rule ignores
#'   it: there the tolerance is the threshold.
#' @param num_dist character, how the predictive distribution of a continuous
#'   target is represented: \code{"auto"} (default), \code{"grid"} or
#'   \code{"exact"}. See the section "Predictive distribution per attacker".
#' @param num_scale character, working scale of the linear model for a
#'   continuous target: \code{"auto"} (default), \code{"identity"} or
#'   \code{"log"}. \code{"auto"} takes the log scale when three conditions hold:
#'   the error metric is relative, all observed synthetic target values are
#'   positive, and the log model has the higher Gaussian log-likelihood on the
#'   original scale (the log model pays the Jacobian, the sum of the logged
#'   target values). Otherwise it takes the identity. Imputed missing values
#'   take no part in this choice, and the log model leaves out imputed values
#'   that are not positive. \code{"log"} with a zero or negative
#'   synthetic value is an error. On the log scale the point prediction is the
#'   predictive median. Other attackers ignore \code{num_scale}, with a warning
#'   when it is given.
#' @param n_quantiles integer, number of quantiles \eqn{K} of a grid
#'   distribution (default 100)
#' @param cat_eval_method character, rule for a categorical target:
#'   \code{"RCS_marginal"} (default), \code{"RCS_conditional"} or \code{"NCE"}.
#'   See the section "The measure".
#' @param cat_tau numeric, threshold for a categorical target. The default
#'   depends on the rule: 0.3 for \code{RCS_marginal}, 1 for
#'   \code{RCS_conditional}, 0.5 for \code{NCE}.
#' @param na_strategy character, treatment of missing values of a continuous
#'   target: \code{"constant"} (default, replaced by \code{na_constant_value}),
#'   \code{"drop"} (the records are removed) or \code{"median"} (replaced by the
#'   median of the original values in \code{X} and of the synthetic values in
#'   \code{Y}). A missing value of a categorical target or a categorical
#'   quasi-identifier becomes the level \code{"missing"}. A numeric
#'   quasi-identifier with missing values is an error. When the target of
#'   \code{X}, \code{Y}, the holdout or the \code{newdata} of \code{predict()}
#'   has missing values, a warning of class
#'   \code{"rapid_warning_missing_target"} gives their number per data set and
#'   their treatment. Under \code{"constant"} and \code{"median"}, and for a
#'   categorical target, the records of \code{X} stay in the denominator of
#'   RAPID. They are scored against the replaced value, so they are rarely
#'   flagged and RAPID tends to be lower than on the observed records. Under
#'   \code{"drop"} RAPID is the share among the records with an observed target.
#'   The linear model chooses its scale under \code{num_scale = "auto"} on the
#'   observed synthetic values.
#' @param na_constant_value numeric, the value of \code{na_strategy =
#'   "constant"} (default 0)
#' @param n_perm integer, number of permutations of the true values that give
#'   the chance level (default 10). 0 skips them, and \code{chance} and
#'   \code{excess} are then \code{NA}.
#' @param threshold_rule character, \code{"fixed"} (default) flags with the
#'   threshold of the rule. \code{"permutation"} replaces it by the
#'   \code{perm_level} quantile of all permuted scores. It needs \code{n_perm >=
#'   1} and a rule with a gain over a baseline (\code{"band"},
#'   \code{"RCS_marginal"}).
#' @param perm_level numeric between 0 and 1, the quantile of the permuted scores
#'   used as threshold under \code{threshold_rule = "permutation"} (default 0.95)
#' @param holdout data frame of original records that were not used to generate
#'   the release (default \code{NULL}, the holdout of a \code{synth_pair} is used
#'   when there is one). They are prepared like \code{X} and scored with the
#'   same attacker and rule (a permutation threshold as found on \code{X}). Their
#'   baselines and permutations take the values of \code{X} as reference. A
#'   missing target value of a holdout record gets the value a missing value of
#'   \code{X} got (\code{na_constant_value}, or the median of \code{X}), or the
#'   record is dropped under \code{na_strategy = "drop"}. A holdout with factor
#'   levels that occur neither in \code{X} nor in \code{Y} is refused with an
#'   error of class \code{"rapid_error_levels"}.
#' @param store_model logical. If \code{TRUE}, the result keeps the fitted model
#'   (field \code{model}) and, in the field \code{predictor}, the attacker, the
#'   prepared synthetic data and the distribution settings. \code{predict()} and
#'   the importance plot \code{plot(x, which = 4)} need both. The prepared
#'   synthetic data make the object about as large as \code{Y}. Default
#'   \code{FALSE}.
#' @param store_predictions logical. If \code{TRUE} (default), the result keeps
#'   the attacker's predictions (field \code{attack}), so that
#'   \code{\link{rapid_test}()}, \code{\link{rapid_threshold_select}()},
#'   \code{update()} and the null outline of \code{plot()} work without a refit.
#' @param seed integer, random seed (default \code{NULL}: the random stream is
#'   left alone). The fit and the permutations are each started from it.
#' @param verbose logical, print progress messages (default \code{FALSE})
#' @param ... arguments passed to the model: to \code{ranger::ranger()} (rf,
#'   lookup), \code{stats::lm()}, \code{rpart::rpart()}, \code{stats::glm()},
#'   or the booster parameters and \code{nrounds} of \code{xgboost::xgb.train()}
#'   (gbm). They override the defaults of the attacker. For cart, \code{prune =
#'   FALSE} keeps the tree as grown. The old argument names
#'   \code{original_data}, \code{synthetic_data}, \code{quasi_identifiers} and
#'   \code{sensitive_attribute} are accepted as aliases of the first four
#'   arguments. An argument whose name is a near miss of an argument of
#'   \code{rapid()} is an error.
#' @param num_epsilon_type deprecated. \code{"percentage"} reads
#'   \code{num_epsilon} as percent, \code{"absolute"} selects the absolute error
#'   metric. A warning of class \code{"rapid_warning_deprecated"} says what to
#'   use instead.
#' @param return_all_records deprecated and ignored: the result always holds all
#'   records. Use \code{as.data.frame(x, at_risk_only = TRUE)} for the flagged
#'   ones.
#'
#' @return An object of class \code{"rapid"}, a list with the fields below. The
#'   methods for it are listed in the section "What you can do with the result".
#' \describe{
#'   \item{\code{rapid}}{RAPID, the share of scored records at risk.}
#'   \item{\code{n_at_risk}, \code{pct_at_risk}, \code{n_total}}{Number and
#'     percentage of records at risk, number of scored records.}
#'   \item{\code{chance}}{The chance level, the mean share flagged over the
#'     permutations (\code{NA} if \code{n_perm = 0}).}
#'   \item{\code{chance_rates}}{The share flagged in each permutation.}
#'   \item{\code{excess}}{\code{rapid - chance}.}
#'   \item{\code{n_perm}, \code{threshold_rule}, \code{perm_level}}{The
#'     permutation settings.}
#'   \item{\code{rule}}{The scoring rule (an internal object with its own
#'     \code{format()} method).}
#'   \item{\code{method}, \code{threshold}}{Name and threshold of the rule. Under
#'     the point rule the threshold is the tolerance, under a permutation
#'     threshold it is the permutation quantile.}
#'   \item{\code{score}}{Per record, the score that is compared with the
#'     threshold: \eqn{r_i} (band, \code{RCS_marginal}), the error (point), the
#'     ratio (\code{RCS_conditional}) or the NCE score.}
#'   \item{\code{baseline}}{Per record, the baseline \eqn{b_i} of the rule
#'     (\code{NULL} for \code{NCE}).}
#'   \item{\code{at_risk}}{Per record, the flag.}
#'   \item{\code{records}}{Data frame with one row per scored record: the
#'     prepared target and quasi-identifiers, for a continuous target the true
#'     value \code{original_value}, the columns of the rule and \code{at_risk}.
#'     The columns of the rule are \code{predicted_value, band_lo, band_hi, g, b,
#'     r} (band), \code{predicted_value, error, b} (point), \code{predicted_class,
#'     g, b, r} (\code{RCS_marginal}), \code{predicted_class, g, b, ratio}
#'     (\code{RCS_conditional}) and \code{predicted_class, g, nce} (\code{NCE}).
#'     A column whose name is taken by a variable gets the suffix
#'     \code{".rapid"}.}
#'   \item{\code{record_cols}}{Named character vector that maps the plain column
#'     names of the rule to their names in \code{records}.}
#'   \item{\code{attack}}{The attacker's predictions on the original (and
#'     holdout) records, \code{NULL} if \code{store_predictions = FALSE}.}
#'   \item{\code{model}}{The fitted model (\code{NULL} unless \code{store_model =
#'     TRUE}).}
#'   \item{\code{predictor}}{What \code{predict()} needs besides the model: the
#'     attacker, the prepared synthetic data and the distribution settings
#'     (\code{NULL} unless \code{store_model = TRUE}).}
#'   \item{\code{holdout}}{\code{list(rapid, chance, excess, n, records)} of the
#'     holdout records (\code{NULL} without holdout).}
#'   \item{\code{release_specific}}{\code{rapid - holdout$rapid}, the part of
#'     RAPID that is specific to the records used for the release (\code{NULL}
#'     without holdout).}
#'   \item{\code{seed}}{The seed of the call, or \code{NULL}.}
#'   \item{\code{model_type}}{Name of the attacker.}
#'   \item{\code{attacker}}{One line with the settings that drive the result,
#'     for example \code{"rf: 500 trees, min.node.size 5"}, \code{"cart: 12
#'     leaves, pruned by cross-validation"} or \code{"lm: log scale (chosen
#'     automatically)"}.}
#'   \item{\code{note}}{A line the print adds for the lookup attacker, otherwise
#'     \code{NULL}.}
#'   \item{\code{model_metrics}}{Accuracy of the predicted class (categorical
#'     target), or \code{mae}, \code{rmse} and their relative forms \code{rmae}
#'     (MAE over the mean absolute value) and \code{rrmse} (RMSE over the
#'     standard deviation) of the point prediction (continuous target).}
#'   \item{\code{formula}, \code{key_vars}, \code{target_var},
#'     \code{is_categorical}}{The attack scenario.}
#'   \item{\code{na}}{The treatment of missing target values:
#'     \code{list(strategy, constant, value, dropped)}, where \code{value} is the
#'     value a missing continuous target got and \code{dropped} the rows of
#'     \code{X} removed under \code{"drop"}.}
#'   \item{\code{call}}{The call, used by \code{update()}.}
#'   \item{\code{version}}{The version of riskutility that created the object.}
#' }
#'
#' @section Arguments by group:
#' \describe{
#'   \item{Data}{\code{X}, \code{Y}, \code{key_vars}, \code{target_var},
#'     \code{holdout}, and the treatment of missing values \code{na_strategy},
#'     \code{na_constant_value}.}
#'   \item{Attacker}{\code{model_type}, the model arguments and \code{prune} in
#'     \code{...}, \code{num_dist}, \code{num_scale}, \code{n_quantiles}.}
#'   \item{Continuous target}{\code{num_eval_method}, \code{num_error_metric},
#'     \code{num_epsilon}, \code{num_delta}, \code{num_tau}, and the deprecated
#'     \code{num_epsilon_type}.}
#'   \item{Categorical target}{\code{cat_eval_method}, \code{cat_tau}.}
#'   \item{Chance level}{\code{n_perm}, \code{threshold_rule},
#'     \code{perm_level}, \code{seed} (which also seeds the fit).}
#'   \item{Storage}{\code{store_model}, \code{store_predictions},
#'     \code{verbose}, and the deprecated \code{return_all_records}.}
#' }
#' Arguments of the other type of target are ignored with a warning of class
#' \code{"rapid_warning_ignored_argument"}.
#'
#' @section The measure:
#' Record \eqn{i} of the original data has the true value \eqn{y_i} and the
#' quasi-identifiers \eqn{x_i}. The attacker's model, fitted on the release,
#' gives a predictive distribution of the target for \eqn{x_i}. Three numbers are
#' computed per record. The confidence \eqn{g_i} is the probability the attacker
#' puts on the true value. The baseline \eqn{b_i} is the probability an attacker
#' who ignores the quasi-identifiers and guesses from the original values would
#' put on it. The normalized gain is
#' \deqn{r_i = \frac{g_i - b_i}{1 - b_i},}{r_i = (g_i - b_i) / (1 - b_i),}
#' with \eqn{r_i = 0} where \eqn{b_i = 1}. A gain of 0 means the attacker does no
#' better than the marginal distribution, a gain of 1 means it is certain of the
#' true value, and \eqn{r_i = 0.3} means it covers 30 percent of the distance
#' between the baseline and certainty. The record is at risk when \eqn{r_i >
#' \tau}{r_i > tau}, and RAPID is the share of records at risk.
#'
#' \strong{Categorical target} (\code{cat_eval_method = "RCS_marginal"},
#' \eqn{\tau}{tau} = \code{cat_tau}). \eqn{g_i} is the predicted probability of
#' the true class and \eqn{b_i} the share of that class among the \eqn{n}
#' original values, \eqn{b_i = n^{-1} \sum_j I(y_j = y_i)}{b_i = (1/n) sum_j
#' I(y_j = y_i)}. Two older rules are available. \code{"RCS_conditional"}
#' replaces \eqn{b_i} by the mean of \eqn{g_j} over the records \eqn{j} of the
#' same class and flags when the ratio \eqn{g_i / b_i} exceeds \code{cat_tau}
#' (default 1). \code{"NCE"} flags when \eqn{1 + \log(g_i) / \log(C)}{1 +
#' log(g_i) / log(C)} exceeds \code{cat_tau} (default 0.5), where \eqn{C} is
#' the number of classes and \eqn{g_i} is bounded below by \eqn{10^{-10}}{1e-10}.
#'
#' \strong{Continuous target, band rule} (\code{num_eval_method = "band"},
#' \eqn{\tau}{tau} = \code{num_tau}). The true class is replaced by the tolerance
#' band \eqn{B_i = \{v : e(y_i, v) < \varepsilon\}}{B_i = {v : e(y_i, v) <
#' epsilon}}, the open interval of all values within the tolerance
#' \eqn{\varepsilon}{epsilon} (\code{num_epsilon}) of the true value. The error
#' \eqn{e} is set by \code{num_error_metric}, with the smoothing constant
#' \eqn{\delta}{delta} (\code{num_delta}):
#' \describe{
#'   \item{\code{"symmetric"}}{\eqn{e(y, v) = 2|y - v| / (|y| + |v| +
#'     2\delta)}{e(y, v) = 2|y - v| / (|y| + |v| + 2 delta)}}
#'   \item{\code{"stabilised_relative"}}{\eqn{e(y, v) = |y - v| / (|y| +
#'     \delta)}{e(y, v) = |y - v| / (|y| + delta)}}
#'   \item{\code{"absolute"}}{\eqn{e(y, v) = |y - v|}}
#' }
#' \eqn{g_i} is the probability the attacker's predictive
#' distribution puts into \eqn{B_i}, and \eqn{b_i = n^{-1} \sum_j I(y_j \in
#' B_i)}{b_i = (1/n) sum_j I(y_j in B_i)} is the share of original values inside
#' the band. With the symmetric error and a tolerance of 10 percent, a true
#' income of 50000 has the band from about 45200 to 55300. The band depends only
#' on the true value and the tolerance, so \eqn{b_i} is the same for every
#' attacker. Near the centre of the distribution \eqn{b_i} is large, and a close
#' guess there earns little gain.
#'
#' \strong{Continuous target, point rule} (\code{num_eval_method = "point"}). A
#' record is at risk when the error of the point prediction is below the
#' tolerance, \eqn{e(y_i, \hat y_i) < \varepsilon}{e(y_i, yhat_i) < epsilon}.
#' The rule has no baseline. It flags records near the centre of the
#' distribution even when the release carries no information about the target
#' (see the examples), which is why it is no longer the default. The column
#' \code{b} of its record table holds the band baseline for comparison.
#'
#' @section Predictive distribution per attacker:
#' A categorical target needs class probabilities, which every attacker
#' provides. A continuous target under the band rule needs the probability the
#' attacker puts into an interval. \code{num_dist} chooses how it is computed.
#' Under \code{"grid"} every attacker gives \eqn{K} = \code{n_quantiles}
#' quantiles \eqn{q_{i1}, \ldots, q_{iK}}{q_i1, ..., q_iK} at the levels
#' \eqn{(k - 0.5)/K}{(k - 0.5)/K}, and \eqn{g_i} is the share of them strictly
#' inside the band. \code{"exact"} computes the probability from the model's own
#' distribution. \code{"auto"} (default) takes the exact form where the attacker
#' has one and the grid otherwise.
#' \describe{
#'   \item{lm}{Normal distribution with mean \eqn{\hat y_i}{yhat_i} and standard
#'     deviation \eqn{\sqrt{se_i^2 + \hat\sigma^2}}{sqrt(se_i^2 + sigma^2)} (the
#'     standard error of the fit and the residual standard deviation) on the
#'     working scale, lognormal on the log scale. \code{"auto"} and
#'     \code{"exact"} use the distribution function, \code{"grid"} its
#'     quantiles.}
#'   \item{cart}{The synthetic target values in the record's leaf.
#'     \code{"auto"} and \code{"exact"} count the share of them inside the band,
#'     \code{"grid"} uses their quantiles (type 1, so that point masses such as
#'     zeros stay).}
#'   \item{rf}{\code{"auto"} and \code{"grid"} use the quantiles of a quantile
#'     regression forest (Meinshausen, 2006) from ranger, which draws one value
#'     per leaf and tree, so \eqn{g_i} carries Monte Carlo noise.
#'     \code{"exact"} uses the weights of Meinshausen (2006): per tree, each
#'     synthetic value in the record's terminal node has weight one over the node
#'     size, all values in the node count (in-bag or not), and the mass is the
#'     mean over the trees.}
#'   \item{gbm}{A second model with the same settings and the quantile loss
#'     estimates the \eqn{K} quantiles, crossing quantiles are sorted. Only the
#'     grid is available (\code{"auto"} gives it).}
#'   \item{lookup}{The synthetic values in the record's leaf of one fully grown
#'     tree, exact (\code{"auto"} or \code{"exact"}, no grid).}
#'   \item{user-defined}{The quantiles of \code{predict_quantiles} of
#'     \code{\link{rapid_learner}()} (\code{"auto"} or \code{"grid"}).}
#' }
#' The point rule uses only the point prediction, which for lm on the log scale
#' is the predictive median.
#'
#' @section Chance level:
#' Some records pass any fixed threshold by chance. An attacker with
#' concentrated predictions, such as a forest with small leaves, lands near some
#' true values even when the quasi-identifiers say nothing about the target. To
#' measure this share, the attacker and its predictions are kept fixed, the true
#' values are permuted \code{n_perm} times, and each record's prediction is
#' scored against the value it receives, with the baseline of that value. The
#' chance level \code{chance} is the share flagged, averaged over the
#' permutations. Each permutation scores all \eqn{n} records, so ten
#' permutations give a stable mean. The print states that RAPID is not
#' distinguishable from the chance level when the excess is at most twice the
#' standard deviation of the shares of the permutations.
#'
#' The excess \code{rapid - chance} is the part of RAPID that comes from the
#' link between the quasi-identifiers and the target. An excess of 0.2 means
#' that 20 percent of the records are flagged in addition to those that the
#' same predictions flag when they are matched with true values at random,
#' that is, the release lets the attacker infer the sensitive value of one
#' record in five beyond chance. \code{\link{rapid_test}()} gives a p-value for the excess.
#'
#' @section Choosing tolerance and threshold:
#' The tolerance \eqn{\varepsilon}{epsilon} says when a guess counts as a
#' disclosure. For the relative metrics it is a proportion (\code{num_epsilon =
#' 0.1} is 10 percent), which suits targets such as income whose spread grows
#' with their level. A tolerance of 10 percent is narrow for a target that
#' varies over orders of magnitude and wide for one with a small relative
#' spread, so the same tolerance does not make RAPID comparable across targets.
#' The absolute metric takes the tolerance in the units of the target (for
#' example 2 years of age) and has no default.
#'
#' The threshold \eqn{\tau}{tau} is the share of the possible gain an attacker
#' must reach. The paper uses 0.3 as reference for both target types. Two
#' data-driven choices exist. \code{threshold_rule = "permutation"} takes as
#' threshold the \code{perm_level} quantile of all permuted scores, so that
#' about 5 percent of the records or fewer (\code{perm_level = 0.95}) pass it
#' when the predictions carry no information. \code{\link{rapid_threshold_select}()}
#' scans a grid of thresholds and returns the largest one at which RAPID is
#' above the \code{1 - alpha} quantile of its permutation distribution.
#'
#' @section Which attacker answers which question:
#' RAPID depends on the attacker, and the attackers answer two different
#' questions. Generalising attackers, the linear model on its fitting scale and
#' the tree pruned by cross-validation, learn how the target depends on the
#' quasi-identifiers. Their RAPID measures inference from the quasi-identifiers,
#' which any release that keeps this dependence allows, for records in the
#' original data and for others alike. They do not flag the records of an exact
#' copy when the target does not depend on the quasi-identifiers.
#'
#' The lookup attacker (\code{model_type = "lookup"}) is one tree of ranger
#' grown on all synthetic records with node size 1, trying all
#' quasi-identifiers at every split. It reads off the synthetic records that
#' share a record's key values and so measures whether the release reproduces
#' the records' values. Its predictions are concentrated, so its chance level is
#' high, and only its excess can be read. The print adds a line that says so.
#'
#' The random forest lies between the two, and where depends on ranger's node
#' size \code{min.node.size} (default 5 for a continuous target, 10 for class
#' probabilities), which \code{rapid()} leaves unchanged and reports in the
#' field \code{attacker}. The table gives RAPID / chance level for the band rule
#' (\eqn{\tau}{tau} = 0.3, tolerance 10 percent) as computed with Oscar Thees'
#' implementation of the band method in a review of 18 November 2026: n = 1000,
#' lognormal target, mean of 3 replications. The first column is an exact copy
#' of the original as the release with a target independent of the
#' quasi-identifiers, so every flag beyond chance is a reproduced record. The
#' second is an independent synthetic sample with a strong dependence of the
#' target on the quasi-identifiers (\eqn{\rho}{rho} = 0.99 in the review). The
#' review reports node size 3 only with the sentence that it flags 100 percent
#' of the copied records, as node size 1 does.
#' \tabular{lll}{
#'   \code{min.node.size} \tab exact copy, no dependence \tab independent sample,
#'     \eqn{\rho}{rho} = 0.99 \cr
#'   1 \tab 1.000 / 0.074 \tab 0.517 / 0.076 \cr
#'   3 \tab 100 percent, chance level not reported \tab not reported \cr
#'   5 (ranger default) \tab 0.736 / 0.051 \tab 0.497 / 0.073 \cr
#'   10 \tab 0.071 / 0.004 \tab 0.443 / 0.065 \cr
#'   20 \tab 0.003 / 0.000 \tab 0.320 / 0.049 \cr
#' }
#' Small nodes detect copied records and large nodes do not, while the
#' inference from the quasi-identifiers falls slowly with the node size. To
#' tell the two apart on real data, pass records that were not used to
#' generate the release as \code{holdout}. RAPID on them measures inference
#' from the quasi-identifiers alone, and \code{release_specific}, RAPID minus
#' holdout RAPID, the part that is specific to the records behind the release.
#' \code{\link{rapid_suite}()} runs several attackers at once, by default rf,
#' cart, gbm and lookup.
#'
#' The tree of \code{"cart"} is grown with \code{cp = 0}, \code{minbucket = 5}
#' and \code{xval = 10} and pruned at the complexity with the smallest
#' cross-validated error. Each of the three values can be set through
#' \code{...}, a \code{control} object replaces all three, and \code{prune =
#' FALSE} keeps the tree as grown. Pruning needs \code{xval >= 1}. The
#' cross-validation draws random numbers, so the tree depends on \code{seed}.
#' The argument \code{prune} is ignored by the other attackers, with a warning.
#' The linear model chooses its scale by the rule under \code{num_scale}.
#'
#' @section What you can do with the result:
#' \describe{
#'   \item{\code{print(x)}, \code{format(x)}}{The attacker line, the rule,
#'     RAPID, chance level and excess, the holdout result.}
#'   \item{\code{summary(x)}}{RAPID with its Wilson interval, the spread of
#'     the flags over rare and common values, model performance.}
#'   \item{\code{as.data.frame(x)}}{The record table, with \code{at_risk_only =
#'     TRUE} the flagged records.}
#'   \item{\code{nobs(x)}}{The number of scored records.}
#'   \item{\code{confint(x)}}{Wilson, Wald, Clopper-Pearson or bootstrap
#'     interval for RAPID.}
#'   \item{\code{update(x, ...)}}{The object with changed arguments, rescored
#'     without a refit where only the rule changes.}
#'   \item{\code{predict(x, newdata)}}{The attacker applied to new records
#'     (needs \code{store_model = TRUE}).}
#'   \item{\code{plot(x, which)}}{Six plot types, see \code{\link{plot.rapid}}.}
#'   \item{\code{rapid_test(x)}}{Permutation test and Wald test of RAPID against
#'     the chance level, see \code{\link{rapid_test}}.}
#'   \item{\code{rapid_threshold_select(x)}}{Threshold chosen by permutation,
#'     see \code{\link{rapid_threshold_select}}.}
#'   \item{\code{rapid_suite(x, ...)}}{Several fitted objects combined into a
#'     model suite, see \code{\link{rapid_suite}}.}
#' }
#' The first eight are documented together in \code{\link{rapid-methods}}.
#'
#' @section Changes in 0.3.0:
#' Calls of riskutility 0.2.0 and the call that gives the same rule now:
#' \describe{
#'   \item{Continuous target, default}{\code{rapid(X, Y, kv, tv)} is now
#'     \code{rapid(X, Y, kv, tv, num_eval_method = "point", num_epsilon =
#'     0.1)}.}
#'   \item{Tolerance in percent}{\code{num_epsilon = 10} is now
#'     \code{num_epsilon = 0.1}.}
#'   \item{Absolute tolerance}{\code{num_epsilon_type = "absolute", num_epsilon
#'     = 2} is now \code{num_error_metric = "absolute", num_epsilon = 2}.}
#'   \item{Categorical target, default}{\code{rapid(X, Y, kv, tv)} is now
#'     \code{rapid(X, Y, kv, tv, cat_eval_method = "RCS_conditional")}.}
#'   \item{Linear model}{\code{model_type = "lm"} is now \code{model_type =
#'     "lm", num_scale = "identity"}.}
#'   \item{Tree}{\code{model_type = "cart"} is now \code{model_type = "cart",
#'     prune = FALSE, cp = 0.01, minsplit = 20, minbucket = 7}.}
#'   \item{Flagged records only}{\code{return_all_records = FALSE} is now
#'     \code{as.data.frame(x, at_risk_only = TRUE)}.}
#'   \item{Old argument names}{\code{rapid(original_data = X, synthetic_data =
#'     Y, quasi_identifiers = kv, sensitive_attribute = tv)} still works and is
#'     the same as \code{rapid(X, Y, kv, tv)}.}
#' }
#' The risk levels LOW, MEDIUM and HIGH of \code{summary()} are gone. Read
#' the excess over the chance level and \code{\link{rapid_test}()} instead.
#'
#' @section Reproducing the paper and the RAPID package:
#' The RAPID paper (Templ, Thees and Müller, 2026) and the RAPID package of
#' Oscar Thees use the band rule with the defaults above and take three
#' settings that differ from the defaults here:
#' \describe{
#'   \item{\code{num_dist = "grid"}}{The predictive distribution of every
#'     attacker as 100 quantiles.}
#'   \item{\code{num_scale = "identity"}}{The linear model on the scale of the
#'     target.}
#'   \item{\code{prune = FALSE, cp = 0.01, minsplit = 20, minbucket = 7}}{The
#'     tree with the defaults of rpart, not pruned. \code{minsplit} must be
#'     given: with \code{minbucket = 7} alone, rpart sets \code{minsplit = 21}.}
#' }
#' With these settings the lm and cart results agree with those of the RAPID
#' package. The rf and gbm results agree up to the random draws of the two
#' implementations.
#'
#' @references
#' Templ, M., Thees, O. and Müller, R. (2026). RAPID: Risk of Attribute
#' Prediction-Induced Disclosure. Manuscript under review at EPJ Data Science.
#'
#' Meinshausen, N. (2006). Quantile Regression Forests. \emph{Journal of Machine
#' Learning Research}, 7, 983--999.
#'
#' @seealso \code{\link{rapid-methods}}, \code{\link{plot.rapid}},
#'   \code{\link{rapid_test}}, \code{\link{rapid_threshold_select}},
#'   \code{\link{rapid_suite}}, \code{\link{rapid_learner}},
#'   \code{\link{rapid_qi_contribution}}, \code{\link{rapid_synthesizer_cv}},
#'   \code{\link{dcap}} and \code{\link{tcap}} for attribution by matching,
#'   \code{\link{disclosure_report}}
#'
#' @author Oscar Thees, Matthias Templ
#' @family rapid
#' @export
#' @importFrom stats predict lm as.formula glm binomial model.matrix complete.cases median
#' @examples
#' set.seed(1)
#' n <- 300
#' X <- data.frame(
#'   age = sample(20:70, n, replace = TRUE),
#'   sex = factor(sample(c("F", "M"), n, replace = TRUE)),
#'   region = factor(sample(c("N", "S", "E", "W"), n, replace = TRUE))
#' )
#' X$income <- round(exp(9 + 0.03 * X$age + 0.3 * (X$sex == "M") + rnorm(n, 0, 0.1)))
#' X$smoker <- factor(ifelse(runif(n) < ifelse(X$age < 40, 0.6, 0.15), "yes", "no"))
#' # A release that keeps the dependence on the quasi-identifiers
#' Y <- X
#' Y$income <- round(Y$income * exp(rnorm(n, 0, 0.05)))
#' kv <- c("age", "sex", "region")
#'
#' # Continuous target: band rule, linear model on the scale it chooses
#' r <- rapid(X, Y, key_vars = kv, target_var = "income", model_type = "lm", seed = 1)
#' r
#' head(as.data.frame(r, at_risk_only = TRUE))
#'
#' # A release that says nothing about income: the band rule flags nothing, the
#' # point rule of 0.2.0 still flags records, as many as it does by chance
#' Y0 <- Y
#' Y0$income <- sample(Y0$income)
#' rapid(X, Y0, key_vars = kv, target_var = "income", model_type = "lm", seed = 1)
#' rapid(X, Y0, key_vars = kv, target_var = "income", model_type = "lm",
#'       num_eval_method = "point", seed = 1)
#'
#' # Categorical target with a pruned tree
#' if (requireNamespace("rpart", quietly = TRUE)) {
#'   rapid(X, Y, key_vars = kv, target_var = "smoker", model_type = "cart", seed = 1)
#' }
#'
#' # Records 201 to 300 were not used for the release
#' rh <- rapid(X[1:200, ], Y[1:200, ], key_vars = kv, target_var = "income",
#'             model_type = "lm", holdout = X[201:300, ], seed = 1)
#' rh
#' rh$release_specific
#'
#' \donttest{
#' # Forest and lookup attacker; read the lookup attacker by its excess
#' rapid(X, Y, key_vars = kv, target_var = "income", model_type = "rf", seed = 1)
#' rapid(X, Y, key_vars = kv, target_var = "income", model_type = "lookup", seed = 1)
#'
#' # Settings of the RAPID paper and the RAPID package
#' rapid(X, Y, key_vars = kv, target_var = "income", model_type = "lm",
#'       num_dist = "grid", num_scale = "identity", seed = 1)
#' if (requireNamespace("rpart", quietly = TRUE)) {
#'   rapid(X, Y, key_vars = kv, target_var = "income", model_type = "cart",
#'         num_dist = "grid", prune = FALSE, cp = 0.01, minsplit = 20, minbucket = 7,
#'         seed = 1)
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
