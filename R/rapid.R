#' RAPID: Risk of Attribute Prediction-Induced Disclosure
#'
#' @description
#' Estimates the share of original records whose sensitive value an attacker can
#' infer from a released (synthetic or anonymized) data set. The attacker fits a
#' model on the release that predicts the sensitive variable from the
#' quasi-identifiers and applies it to the original records. A record is at risk
#' when the attacker puts clearly more probability on its true value than the
#' marginal distribution of the sensitive variable does. RAPID is the share of
#' such records.
#'
#' A first call needs four arguments: the original data \code{X}, the release
#' \code{Y}, the quasi-identifiers \code{key_vars} and the sensitive variable
#' \code{target_var}. The other arguments have defaults. The first ones to
#' decide on are the attacker (\code{model_type}), for a numeric target the
#' tolerance (\code{num_epsilon}), and \code{seed}.
#'
#' Read the result by three numbers. \code{rapid} is the share of records at
#' risk. \code{chance} is the share at risk when the same predictions are paired
#' with randomly reassigned true values. \code{excess}, their difference, is the
#' part that comes from the link between the quasi-identifiers and the sensitive
#' variable. See the section "Reading the result" and
#' \code{vignette("rapid", package = "riskutility")}.
#'
#' @param X the original data, a data frame with one row per record, or a
#'   \code{\link{synth_pair}} that holds the original and the synthetic data,
#'   the quasi-identifiers, the target and optionally a holdout
#' @param Y data frame of the released (synthetic or anonymized) data the
#'   attacker is fitted on. Not needed if \code{X} is a \code{synth_pair}. It
#'   need not have the same number of rows as \code{X}.
#' @param key_vars character vector with the names of the quasi-identifiers,
#'   the variables the attacker knows about the records. They must be in
#'   \code{X} and \code{Y} and of the same kind in both, numeric or categorical
#'   (factor, character or logical). Other types, such as dates, are an error.
#'   The names must be syntactic, because they enter a model formula.
#'   \code{target_var} must not be among them.
#' @param target_var character, name of the sensitive variable, in \code{X} and
#'   \code{Y} and of the same kind in both. A factor, character or logical
#'   variable is a categorical target, and each distinct value is a class (an
#'   ordered factor is treated as unordered). A numeric variable is a numeric
#'   target and is scored with a tolerance band. Convert class codes that are
#'   stored as numbers (0 and 1, or 1 to 5) to a factor first. Otherwise they
#'   are scored as quantities.
#' @param model_type the attacker, one of the names below or an object from
#'   \code{\link{rapid_learner}()}. RAPID depends on this choice, see the section
#'   "Attackers".
#'   \describe{
#'     \item{\code{"rf"}}{Random forest from ranger with its defaults (500
#'       trees). The default.}
#'     \item{\code{"lm"}}{Linear model. Numeric target only.}
#'     \item{\code{"cart"}}{Classification or regression tree from rpart, pruned
#'       by cross-validation. Needs the package rpart.}
#'     \item{\code{"gbm"}}{Gradient boosting from xgboost with 100 rounds and
#'       the default booster parameters. Needs the package xgboost.}
#'     \item{\code{"logit"}}{Logistic regression. Categorical target with
#'       exactly two classes.}
#'     \item{\code{"lookup"}}{One fully grown tree that reads off the synthetic
#'       records with the same key values. Read it by its excess.}
#'   }
#'   A missing package is an error of class \code{"rapid_error_dependency"}.
#'   \code{"lm"}, \code{"logit"} and user-defined attackers cannot predict a
#'   record with a level of a quasi-identifier that \code{Y} does not contain.
#'   They stop with an error of class \code{"rapid_error_levels"}. The
#'   tree-based attackers accept such records.
#' @param num_eval_method character, rule for a numeric target: \code{"band"}
#'   (default) or \code{"point"}, the rule of riskutility 0.2.0. Numeric target
#'   only. See the section "The measure".
#' @param num_error_metric character, the error \eqn{e(y, v)} between the true
#'   value \eqn{y} and a guess \eqn{v} that defines the tolerance band. Numeric
#'   target only. \code{"symmetric"} (default) is the difference relative to the
#'   mean size of the two values, \code{"stabilised_relative"} the difference
#'   relative to the true value, \code{"absolute"} the difference in the units
#'   of the target. Take a relative metric for a target whose spread grows with
#'   its level, such as income, and the absolute metric for a target with a
#'   natural unit of closeness, such as age in years. The formulas are in the
#'   section "The measure".
#' @param num_epsilon numeric, the tolerance \eqn{\varepsilon}{epsilon}: how
#'   close a guess must be to the true value to count as a disclosure. Numeric
#'   target only. For the two relative metrics it is a proportion, above 0 and
#'   below 1 (default 0.1, that is 10 percent). For the absolute metric it is in
#'   the units of the target, has no default and must be given. A relative
#'   tolerance of 1 or more is read as percent, with a warning of class
#'   \code{"rapid_warning_legacy_tolerance"}. A larger tolerance widens the
#'   band, which raises the confidence and the baseline. Under the point rule
#'   the tolerance is also the threshold.
#' @param num_delta numeric, the constant \eqn{\delta}{delta} that the two
#'   relative metrics add to the denominator, in the units of the target
#'   (default 0.01, at least 0). Numeric target only. It keeps the error finite
#'   for true values at or near zero. Around a true value of 0 the band reaches
#'   about \eqn{\varepsilon \delta}{epsilon times delta} to each side. For
#'   values far above \code{num_delta} it has no visible effect. Choose it small
#'   against the typical size of the target, for example 0.001 for a share
#'   between 0 and 1. The absolute metric does not use it.
#' @param num_tau numeric, threshold \eqn{\tau}{tau} on the normalized gain of a
#'   numeric target. Band rule only. The default 0.3 flags a record when the
#'   attacker covers more than 30 percent of the distance between the baseline
#'   and certainty. Allowed values are 0 or more and less than 1, and any other
#'   value is an error of class \code{"rapid_error_input"}. The point rule
#'   ignores \code{num_tau}, with a warning when it is given, because there the
#'   tolerance is the threshold.
#' @param num_dist character, how the predictive distribution of a numeric
#'   target is represented under the band rule: \code{"auto"} (default, exact
#'   where the attacker has an exact form, a grid of \code{n_quantiles}
#'   quantiles otherwise), \code{"grid"} or \code{"exact"}. A value the attacker
#'   does not offer is an error of class \code{"rapid_error_model"}: gbm has no
#'   exact form, lookup no grid. See the section "Attackers".
#' @param num_scale character, working scale of the linear model:
#'   \code{"auto"} (default), \code{"identity"} or \code{"log"}. Only for
#'   \code{model_type = "lm"} and a numeric target. Other attackers ignore it,
#'   with a warning when it is given. \code{"auto"} takes the log scale when the
#'   error metric is relative, all observed synthetic target values are positive
#'   and the log model has the higher likelihood. \code{"log"} with a zero or
#'   negative synthetic value is an error. On the log scale the point prediction
#'   is the predictive median. See the section "Attackers".
#' @param n_quantiles integer, number of quantiles \eqn{K} of a grid
#'   distribution (default 100, a grid needs at least 2). Numeric target, band
#'   rule only. It is used by \code{"rf"} and \code{"gbm"}, by user-defined
#'   attackers, and by \code{"lm"} and \code{"cart"} under \code{num_dist =
#'   "grid"}. It has no effect where the distribution is exact. A larger value
#'   gives a finer confidence \eqn{g_i} (steps of \eqn{1/K}) and a larger
#'   object.
#' @param cat_eval_method character, rule for a categorical target:
#'   \code{"RCS_marginal"} (default, the gain over the marginal share of the
#'   class), \code{"RCS_conditional"} (the rule of riskutility 0.2.0, the ratio
#'   to the mean confidence within the class) or \code{"NCE"} (normalized
#'   cross-entropy). Categorical target only. See the section "The measure".
#' @param cat_tau numeric, threshold for a categorical target. Categorical
#'   target only. The default depends on the rule: 0.3 for \code{RCS_marginal},
#'   1 for \code{RCS_conditional}, 0.5 for \code{NCE}. Under \code{RCS_marginal}
#'   it is a threshold on a gain, with the meaning and the range of
#'   \code{num_tau}. A value outside the range is an error of class
#'   \code{"rapid_error_input"}. For \code{cat_tau = 1} the message names the
#'   rule of 0.2.0, which had this default.
#' @param na_strategy character, treatment of missing values of a numeric
#'   target: \code{"drop"} (default) removes the records, \code{"constant"}
#'   replaces the values by \code{na_constant_value} (the default of 0.2.0),
#'   \code{"median"} replaces them by the median of their data set. It has no
#'   effect on a categorical target or on the quasi-identifiers. See the section
#'   "Missing values".
#' @param na_constant_value numeric, the value of \code{na_strategy =
#'   "constant"} (default 0). Numeric target only.
#' @param n_perm integer, number of permutations of the true values that give
#'   the chance level (default 10, a whole number of 0 or more). 0 skips them,
#'   and \code{chance} and \code{excess} are then \code{NA}. The print can say
#'   that RAPID is not distinguishable from the chance level only with 2 or more
#'   permutations.
#' @param threshold_rule character, \code{"fixed"} (default) flags with the
#'   threshold of the rule. \code{"permutation"} replaces it by the
#'   \code{perm_level} quantile of all permuted scores. Band rule and
#'   \code{RCS_marginal} only, the rules with a gain over a baseline. It needs
#'   \code{n_perm >= 1}.
#' @param perm_level numeric between 0 and 1, the quantile of the permuted scores
#'   used as threshold (default 0.95). Only under \code{threshold_rule =
#'   "permutation"}. Under a fixed threshold it is ignored, with a warning when
#'   it is given.
#' @param holdout data frame of original records that were not used to generate
#'   the release, with the columns \code{key_vars} and \code{target_var}
#'   (default \code{NULL}, the holdout of a \code{synth_pair} is used when there
#'   is one). They are scored with the same attacker and rule as \code{X}, under
#'   a permutation threshold with the threshold found on \code{X}. Their
#'   baselines and permutations take the values of \code{X} as reference. The
#'   result is in the field \code{holdout}, and \code{release_specific} is RAPID
#'   minus the holdout RAPID. Missing values are treated as described in the
#'   section "Missing values". A holdout with factor levels that occur neither
#'   in \code{X} nor in \code{Y} is refused with an error of class
#'   \code{"rapid_error_levels"}. This includes a missing value of a categorical
#'   target when neither \code{X} nor \code{Y} has one.
#' @param store_model logical. If \code{TRUE}, the result keeps the fitted model
#'   (field \code{model}) and, in the field \code{predictor}, the attacker, the
#'   prepared synthetic data and the distribution settings. \code{predict()} and
#'   the importance plot \code{plot(x, which = 4)} need both. The prepared
#'   synthetic data make the object at least as large as \code{Y}, and a
#'   random forest is much larger than its training data: with \code{"rf"} the
#'   stored model can take many times the size of \code{Y}. Default
#'   \code{FALSE}.
#' @param store_predictions logical. If \code{TRUE} (default), the result keeps
#'   the attacker's predictions (field \code{attack}), so that
#'   \code{\link{rapid_test}()}, \code{\link{rapid_threshold_select}()},
#'   \code{update()} and the permutation null of \code{plot()} work without a
#'   refit. The predictions of a numeric target include the predictive
#'   distribution. For a grid distribution (the default of \code{"rf"} and
#'   \code{"gbm"}) they take \code{8 * n_quantiles} bytes per record, 0.8 kB
#'   with the default, which can be many times the size of the data.
#'   \code{FALSE} keeps the object small. \code{update()} then refits,
#'   \code{rapid_test()} and \code{rapid_threshold_select()} stop with an error
#'   of class \code{"rapid_error_object"}, and the plots leave out the
#'   permutation null.
#' @param seed integer, random seed, or \code{NULL} (default). With a seed,
#'   \code{set.seed(seed)} is called before the fit and again before the
#'   permutations. The result is then reproducible, and \code{update()} without
#'   a refit draws the same permutations. The random stream of the session is
#'   not restored afterwards. It continues from this seed. With \code{NULL} the
#'   stream is used as it stands. In a simulation loop, leave \code{seed = NULL}
#'   and call \code{set.seed()} once before the loop. Otherwise the draws after
#'   every call repeat. See the section "Reproducibility".
#' @param verbose logical, print progress messages (default \code{FALSE})
#' @param ... arguments for the model of the attacker. They override its
#'   defaults and go to \code{ranger::ranger()} (rf, lookup),
#'   \code{stats::lm()}, \code{rpart::rpart()} (cart), \code{stats::glm()}
#'   (logit), or the \code{fit} function of a user-defined attacker. For gbm
#'   they are booster parameters of \code{xgboost::xgb.train()} such as
#'   \code{max_depth} and \code{eta}, and \code{nrounds}. Examples:
#'   \code{num.trees = 1000}, \code{min.node.size = 1}, \code{cp = 0.01}. For
#'   cart, \code{prune = FALSE} keeps the tree as grown. The old argument names
#'   \code{original_data}, \code{synthetic_data}, \code{quasi_identifiers} and
#'   \code{sensitive_attribute} are accepted as aliases of the first four
#'   arguments. A name within a typing error of an argument of \code{rapid()}
#'   (for example \code{num_taus} or \code{nperm}) is an error of class
#'   \code{"rapid_error_input"}. Every other unknown name is passed to the
#'   model. \code{tau}, \code{epsilon}, \code{tolerance} and
#'   \code{n_permutations} are not arguments of \code{rapid()}: the threshold is
#'   \code{num_tau} or \code{cat_tau}, the tolerance \code{num_epsilon}, the
#'   number of permutations \code{n_perm}. Given by these names they reach the
#'   model, and the result uses the defaults. ranger and \code{lm()} then give
#'   a warning of their own (ranger: "Unused arguments: tau"), rpart stops with
#'   an error, and xgboost takes them without a message. Check the line "Rule:"
#'   of the print.
#' @param num_epsilon_type deprecated. Numeric target only.
#'   \code{"percentage"} reads \code{num_epsilon} as percent, \code{"absolute"}
#'   selects the absolute error metric. A warning of class
#'   \code{"rapid_warning_deprecated"} says what to use instead.
#' @param return_all_records deprecated and ignored: the result always holds all
#'   records. Use \code{as.data.frame(x, at_risk_only = TRUE)} for the records
#'   at risk.
#'
#' @return An object of class \code{"rapid"}, a list. Start with \code{rapid},
#'   \code{chance} and \code{excess}, then look at \code{records}. The methods
#'   for the object are listed in the section "Methods and companion functions".
#'
#' \strong{The result}
#' \describe{
#'   \item{\code{rapid}}{RAPID, the share of scored records at risk.}
#'   \item{\code{n_at_risk}, \code{pct_at_risk}, \code{n_total}}{Number and
#'     percentage of records at risk, number of scored records.}
#'   \item{\code{chance}}{The chance level, the mean share at risk over the
#'     permutations (\code{NA} if \code{n_perm = 0}).}
#'   \item{\code{chance_rates}}{The share at risk in each permutation.}
#'   \item{\code{excess}}{\code{rapid - chance}.}
#'   \item{\code{holdout}}{\code{list(rapid, chance, excess, n, records)} of the
#'     holdout records (\code{NULL} without holdout).}
#'   \item{\code{release_specific}}{\code{rapid - holdout$rapid}, the part of
#'     RAPID that is specific to the records used for the release (\code{NULL}
#'     without holdout). No interval is computed. With the flags treated as
#'     independent, as \code{confint()} does, its standard error is
#'     \code{sqrt(x$rapid * (1 - x$rapid) / x$n_total + x$holdout$rapid * (1 -
#'     x$holdout$rapid) / x$holdout$n)}. In the holdout example below (200 and
#'     100 records) it is 0.053. \code{rapid_test()}, \code{confint()} and
#'     \code{summary()} work on the records of \code{X} and do not use the
#'     holdout.}
#'   \item{\code{attacker}}{One line with the settings that drive the result,
#'     for example \code{"rf: 500 trees, min.node.size 5"}, \code{"cart: 12
#'     leaves, pruned by cross-validation"} or \code{"lm: log scale (chosen
#'     automatically)"}.}
#'   \item{\code{note}}{A line the print adds for the lookup attacker, otherwise
#'     \code{NULL}.}
#'   \item{\code{model_metrics}}{Accuracy of the predicted class (categorical
#'     target), or \code{mae}, \code{rmse} and their relative forms \code{rmae}
#'     (MAE over the mean absolute value) and \code{rrmse} (RMSE over the
#'     standard deviation) of the point prediction (numeric target).}
#' }
#' \strong{Per record}
#' \describe{
#'   \item{\code{records}}{Data frame with one row per scored record. The row
#'     names are those of \code{X}, also after records were removed, so they
#'     link a row back to \code{X}. The columns are the prepared target and
#'     quasi-identifiers, for a numeric target the true value again as
#'     \code{original_value}, the columns of the rule and \code{at_risk}.
#'     \code{original_value} repeats the target column under a fixed name, so
#'     that code does not depend on the name of the target. The columns of the
#'     rule are \code{predicted_value, band_lo, band_hi, g, b, r} (band),
#'     \code{predicted_value, error, b} (point), \code{predicted_class, g, b, r}
#'     (\code{RCS_marginal}), \code{predicted_class, g, b, ratio}
#'     (\code{RCS_conditional}) and \code{predicted_class, g, nce} (\code{NCE}).
#'     Here \code{g} is the confidence, \code{b} the baseline, \code{r} the
#'     gain, and \code{band_lo}, \code{band_hi} the limits of the tolerance
#'     band. \code{predicted_value} is the point prediction of the attacker,
#'     which the band rule does not use. It is the fitted value for lm (the
#'     predictive median on the log scale), the mean of the synthetic target
#'     values in the record's leaf for cart and lookup, the mean over the trees
#'     for rf, the prediction of the boosted model for gbm, and the value of
#'     \code{predict_point} for a user-defined attacker. A column whose name is
#'     taken by a variable gets the suffix \code{".rapid"}.}
#'   \item{\code{score}}{Per record, the score that is compared with the
#'     threshold: \eqn{r_i} (band, \code{RCS_marginal}), the error (point), the
#'     ratio (\code{RCS_conditional}) or the NCE score.}
#'   \item{\code{baseline}}{Per record, the baseline \eqn{b_i} of the rule
#'     (\code{NULL} for \code{NCE}).}
#'   \item{\code{at_risk}}{Per record, the flag.}
#'   \item{\code{record_cols}}{Named character vector that maps the plain column
#'     names of the rule to their names in \code{records}.}
#' }
#' \strong{Settings of the call}
#' \describe{
#'   \item{\code{method}, \code{threshold}}{Name and threshold of the rule. Under
#'     the point rule the threshold is the tolerance, under a permutation
#'     threshold it is the permutation quantile.}
#'   \item{\code{rule}}{The scoring rule (an internal object with its own
#'     \code{format()} method).}
#'   \item{\code{n_perm}, \code{threshold_rule}, \code{perm_level}}{The
#'     permutation settings.}
#'   \item{\code{model_type}}{Name of the attacker.}
#'   \item{\code{formula}, \code{key_vars}, \code{target_var},
#'     \code{is_categorical}}{The attack scenario.}
#'   \item{\code{na}}{The treatment of missing target values:
#'     \code{list(strategy, constant, value, dropped)}. \code{value} is the
#'     value a missing numeric target got (\code{NULL} under \code{"drop"} and
#'     for a categorical target), \code{dropped} the rows of \code{X} removed
#'     under \code{"drop"}.}
#'   \item{\code{seed}}{The seed of the call, or \code{NULL}.}
#'   \item{\code{call}}{The call, used by \code{update()}.}
#'   \item{\code{version}}{The version of riskutility that created the object.}
#' }
#' \strong{Kept for the methods}
#' \describe{
#'   \item{\code{attack}}{The attacker's predictions on the original (and
#'     holdout) records, \code{NULL} if \code{store_predictions = FALSE}.}
#'   \item{\code{model}}{The fitted model (\code{NULL} unless \code{store_model =
#'     TRUE}).}
#'   \item{\code{predictor}}{What \code{predict()} needs besides the model: the
#'     attacker, the prepared synthetic data and the distribution settings
#'     (\code{NULL} unless \code{store_model = TRUE}).}
#' }
#'
#' @section Reading the result:
#' The first example below (numeric target, default attacker, \code{seed = 1})
#' prints three lines of numbers.
#' \describe{
#'   \item{\code{RAPID = 0.120: 36 of 300 records at risk}}{For 36 original
#'     records the gain of the attacker over the marginal distribution is above
#'     the threshold of 0.3.}
#'   \item{\code{Chance level = 0.012 (10 permutations)}}{About 1 percent of the
#'     records pass the same threshold when the true values are permuted, that
#'     is, when the predictions carry no information about the record.}
#'   \item{\code{Excess = 0.108}}{The release lets this attacker infer the
#'     sensitive value of about 11 percent of the records beyond chance.}
#' }
#' These points help to read the numbers.
#' \itemize{
#'   \item Look at the excess first. It is the share of records the attacker
#'     infers beyond random pairing. This matters most when the chance level is
#'     high, as for a categorical target with few classes and for the lookup
#'     attacker. \code{\link{rapid_test}()} says whether the excess differs
#'     from 0. It does not say whether the excess is large.
#'   \item RAPID depends on the attacker. On the same data the linear model
#'     gives 0.760 where the forest gives 0.120. Report the attacker with the
#'     number, and compare attackers with \code{\link{rapid_suite}()}.
#'   \item Ask where the excess comes from. With a holdout, the excess on the
#'     holdout records is inference from the quasi-identifiers. It applies to
#'     records in the original data and to others alike, because the release
#'     keeps the dependence of the target on the quasi-identifiers.
#'     \code{release_specific} is what the records behind the release carry in
#'     addition.
#'   \item Without a holdout, compare a generalizing attacker (\code{"lm"},
#'     \code{"cart"}) with the lookup attacker. An excess under the lookup
#'     attacker that the linear model does not show points to reproduced
#'     records.
#'   \item \code{confint()} gives an interval for RAPID. The print says when
#'     RAPID is not distinguishable from the chance level.
#'   \item \code{as.data.frame(x, at_risk_only = TRUE)} lists the records at
#'     risk with their confidence \code{g}, baseline \code{b} and gain \code{r}.
#' }
#' There is no general cut-off. Choose the tolerance and the threshold for a
#' sensitive variable before looking at the result. A small excess under the
#' attackers tried does not show that the release is safe.
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
#' \tau}{r_i > tau}, and RAPID is the share of records at risk. Such a record is
#' also called flagged.
#'
#' \strong{Categorical target} (\code{cat_eval_method = "RCS_marginal"},
#' \eqn{\tau}{tau} = \code{cat_tau}). \eqn{g_i} is the predicted probability of
#' the true class and \eqn{b_i} the share of that class among the \eqn{n}
#' original values, \eqn{b_i = n^{-1} \sum_j I(y_j = y_i)}{b_i = (1/n) sum_j
#' I(y_j = y_i)}. Two older rules are available. \code{"RCS_conditional"}
#' replaces \eqn{b_i} by the mean of \eqn{g_j} over the records \eqn{j} of the
#' same class and flags when the ratio \eqn{g_i / b_i} exceeds \code{cat_tau}
#' (default 1). \code{"NCE"}, the normalized cross-entropy, flags when \eqn{1 +
#' \log(g_i) / \log(C)}{1 + log(g_i) / log(C)} exceeds \code{cat_tau} (default
#' 0.5), where \eqn{C} is the number of classes and \eqn{g_i} is bounded below
#' by \eqn{10^{-10}}{1e-10}.
#'
#' \strong{Numeric target, band rule} (\code{num_eval_method = "band"},
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
#' attacker. Near the center of the distribution \eqn{b_i} is large, and a close
#' guess there earns little gain.
#'
#' \strong{Numeric target, point rule} (\code{num_eval_method = "point"}). A
#' record is at risk when the error of the point prediction is below the
#' tolerance, \eqn{e(y_i, \hat y_i) < \varepsilon}{e(y_i, yhat_i) < epsilon}.
#' The rule has no baseline. It flags records near the center of the
#' distribution even when the release carries no information about the target
#' (see the examples), which is why it is no longer the default. The column
#' \code{b} of its record table holds the band baseline for comparison.
#'
#' Arguments of the other type of target are ignored with a warning of class
#' \code{"rapid_warning_ignored_argument"}, and so are \code{num_tau} under the
#' point rule and \code{perm_level} under a fixed threshold.
#'
#' \strong{Terms.} The \emph{prepared} data are the columns \code{key_vars} and
#' \code{target_var}, with character and logical variables turned into factors,
#' the same factor levels in \code{X} and \code{Y}, and missing values treated
#' as described in the section "Missing values". The \emph{rule} (field
#' \code{method}) says how a record is scored and when it is at risk. A record
#' \emph{at risk} is the same as a flagged record. The \emph{attacker} is the
#' prediction model (\code{model_type}). The \emph{quasi-identifiers} are the
#' key variables (\code{key_vars}).
#'
#' @section Chance level:
#' Some records pass any fixed threshold by chance. An attacker with
#' concentrated predictions, such as a forest with small leaves, lands near some
#' true values even when the quasi-identifiers say nothing about the target. To
#' measure this share, the attacker and its predictions are kept fixed, the true
#' values are permuted \code{n_perm} times, and each record's prediction is
#' scored against the value it receives, with the baseline of that value. The
#' chance level \code{chance} is the share at risk, averaged over the
#' permutations. For a categorical target with few classes the chance level is
#' high by construction. A record is at risk in every permutation that assigns
#' it a class whose predicted probability passes the threshold for that class,
#' and with two classes this happens often. RAPID alone then says little. Read
#' the excess.
#'
#' Each permutation scores all \eqn{n} records, so ten permutations give a
#' stable mean. The print states that RAPID is not distinguishable from the
#' chance level when the excess is at most twice the standard deviation of the
#' shares of the permutations.
#'
#' The excess \code{rapid - chance} is the part of RAPID that comes from the
#' link between the quasi-identifiers and the target. An excess of 0.2 means
#' that 20 percent of the records are at risk in addition to those that the
#' same predictions flag when they are matched with true values at random,
#' that is, the release lets the attacker infer the sensitive value of one
#' record in five beyond chance. \code{\link{rapid_test}()} gives a p-value for
#' the excess. See also the vignette section "Chance level and excess".
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
#' when the predictions carry no information.
#' \code{\link{rapid_threshold_select}()} scans a grid of thresholds and returns
#' the strictest one at which RAPID is above the \code{1 - alpha} quantile of
#' its permutation distribution. Its help page and the vignette section
#' "Choosing the threshold" say how to read the selected threshold.
#'
#' @section Attackers:
#' RAPID depends on the attacker. For a report, run
#' \code{\link{rapid_suite}()} and give the largest excess with its attacker
#' (fields \code{excess_max} and \code{worst_model_excess}). Its default
#' attackers are rf, cart, gbm and lookup. For a numeric target add the linear
#' model with \code{model_types = c("lm", "rf", "cart", "gbm", "lookup")}.
#'
#' The attackers answer two different questions. Generalizing attackers, the
#' linear model and the tree pruned by cross-validation, learn how the target
#' depends on the quasi-identifiers. Their RAPID measures inference from the
#' quasi-identifiers, which any release that keeps this dependence allows, for
#' records in the original data and for others alike. They do not flag the
#' records of an exact copy when the target does not depend on the
#' quasi-identifiers.
#'
#' The lookup attacker (\code{model_type = "lookup"}) is one tree of ranger
#' grown on all synthetic records with node size 1, trying all
#' quasi-identifiers at every split. It reads off the synthetic records that
#' share a record's key values and so measures whether the release reproduces
#' the records' values. Its predictions are concentrated, so its chance level is
#' high, and only its excess can be read. The print adds a line that says so.
#'
#' The random forest, the default, takes both types of target and numeric and
#' categorical quasi-identifiers. It lies between the two kinds of attackers,
#' and its single value does not say whether it comes from inference or from
#' reproduced records. Where it lies depends on settings of ranger that
#' \code{rapid()} leaves at their defaults: the node size \code{min.node.size}
#' (5 for a numeric target, 10 for class probabilities, reported in the field
#' \code{attacker}) and \code{mtry}. When all quasi-identifiers are tried at
#' every split, small nodes detect copied records and large nodes do not. On an
#' exact copy of the example data with a permuted target, the forest with
#' \code{mtry = 3} gave a RAPID of 0.840 at node size 1 and of 0.007 at node
#' size 20. With the default \code{mtry} it flagged no record at either node
#' size. The forest also depends on the seed. Over the seeds 1 to 10 its RAPID
#' on the example data lay between 0.083 and 0.160, and the interval of
#' \code{confint()} does not include this variation.
#'
#' To tell inference and reproduction apart on real data, pass records that
#' were not used to generate the release as \code{holdout}. RAPID on them
#' measures inference from the quasi-identifiers alone, and
#' \code{release_specific}, RAPID minus holdout RAPID, the part that is specific
#' to the records behind the release. The vignette sections "Inference and
#' reproduction: which attacker answers which question", "Model suite and worst
#' case" and "Holdout and release-specific risk" show these comparisons.
#'
#' The tree of \code{"cart"} is grown with \code{cp = 0}, \code{minbucket = 5}
#' and \code{xval = 10} and pruned at the complexity with the smallest
#' cross-validated error. Each of the three values can be set through
#' \code{...}, a \code{control} object replaces all three, and \code{prune =
#' FALSE} keeps the tree as grown. Pruning needs \code{xval >= 1}. The
#' cross-validation draws random numbers, so the size of the tree depends on
#' \code{seed}. The field \code{attacker} gives the number of leaves. On the
#' example data the tree had 25 to 31 leaves over the seeds 1 to 10, with a
#' RAPID between 0.740 and 0.757. The argument \code{prune} is ignored by the
#' other attackers, with a warning.
#'
#' The linear model works on the scale set by \code{num_scale}. Under
#' \code{"auto"} it takes the log scale when three conditions hold: the error
#' metric is relative, all observed synthetic target values are positive, and
#' the log model has the higher Gaussian log-likelihood on the original scale
#' (the log model pays the Jacobian, the sum of the logged target values).
#' Otherwise it takes the identity. Imputed missing values take no part in this
#' choice, and the log model leaves out imputed values that are not positive.
#' The field \code{attacker} names the scale.
#'
#' \strong{Predictive distribution.} A categorical target needs class
#' probabilities, which every attacker provides. A numeric target under the band
#' rule needs the probability the attacker puts into an interval.
#' \code{num_dist} chooses how it is computed. Under \code{"grid"} every
#' attacker gives \eqn{K} = \code{n_quantiles} quantiles \eqn{q_{i1}, \ldots,
#' q_{iK}}{q_i1, ..., q_iK} at the levels \eqn{(k - 0.5)/K}{(k - 0.5)/K}, and
#' \eqn{g_i} is the share of them strictly inside the band. \code{"exact"}
#' computes the probability from the model's own distribution. \code{"auto"}
#' (default) takes the exact form where the attacker has one and the grid
#' otherwise.
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
#'     \code{"exact"} uses the weights of Meinshausen (2006) and removes this
#'     noise, at the cost of a larger object.}
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
#' @section Missing values:
#' \describe{
#'   \item{Numeric target}{\code{na_strategy} decides. Under \code{"drop"}
#'     (default) the records of \code{X}, \code{Y} and the holdout with a missing
#'     target are removed. RAPID is the share among the remaining records of
#'     \code{X}, \code{nobs()} gives their number, \code{x$na$dropped} the
#'     removed rows of \code{X}, and the record table keeps the row names of
#'     \code{X}. Under \code{"constant"} the values are replaced by
#'     \code{na_constant_value}. Under \code{"median"} they are replaced by the
#'     median of the observed values, of \code{X} in \code{X} and of \code{Y} in
#'     \code{Y}. A holdout record gets the value a record of \code{X} got. Under
#'     both, the records stay in the denominator of RAPID. They are scored
#'     against the replaced value, so they are rarely at risk and RAPID tends to
#'     be lower than on the observed records. The linear model chooses its scale
#'     under \code{num_scale = "auto"} on the observed synthetic values.}
#'   \item{Categorical target}{A missing value becomes the class
#'     \code{"missing"} in \code{X} and \code{Y}, whatever \code{na_strategy}
#'     says, and the records stay in the denominator. A holdout record or a
#'     record of the \code{newdata} of \code{predict()} with a missing target
#'     is refused with an error of class \code{"rapid_error_levels"} when
#'     neither \code{X} nor \code{Y} has a missing target value. Remove such
#'     records first.}
#'   \item{Categorical quasi-identifier}{A missing value becomes the level
#'     \code{"missing"}. The attackers \code{"lm"}, \code{"logit"} and
#'     user-defined ones stop with an error of class
#'     \code{"rapid_error_levels"} when \code{X} has this level and \code{Y}
#'     does not. The tree-based attackers accept it.}
#'   \item{Numeric quasi-identifier}{A missing value is an error of class
#'     \code{"rapid_error_input"}. Impute or remove these records before the
#'     call.}
#' }
#' When the target of \code{X}, \code{Y}, the holdout or the \code{newdata} of
#' \code{predict()} has missing values, a warning of class
#' \code{"rapid_warning_missing_target"} gives their number per data set and
#' their treatment.
#'
#' @section Reproducibility:
#' A result is reproducible when every random step has a seed.
#' \describe{
#'   \item{\code{rapid()}}{The fit of the forest, the cross-validation of the
#'     pruned tree, and the permutations of every attacker are random. All
#'     start from \code{seed}.}
#'   \item{\code{update()}, \code{plot()}}{The permutations are random. They
#'     start from the seed stored in the object. Without one they are new
#'     draws.}
#'   \item{\code{rapid_test()}, \code{rapid_threshold_select()} on an
#'     object}{The permutations are random. They take the \code{seed} of
#'     these functions. The seed of the object is not used.}
#'   \item{\code{confint(method = "bootstrap")}}{The resampling is random. It
#'     takes the \code{seed} of \code{confint()}.}
#'   \item{\code{rapid_suite()}, \code{rapid_qi_contribution()}}{Every fit and
#'     its permutations are random. Each fit starts from the same \code{seed}.}
#'   \item{\code{rapid_synthesizer_cv()}}{The folds, the synthesizer and the
#'     fits are random. Its \code{seed} covers all three.}
#' }
#' \code{set.seed()} before a call without \code{seed} works as well. The linear
#' model has no random fit, only its permutations are random. A forest with a
#' seed gave the same result with 1 and with 4 threads (\code{num.threads}).
#'
#' @section Methods and companion functions:
#' \describe{
#'   \item{\code{print(x)}, \code{format(x)}}{The attacker line, the rule,
#'     RAPID, chance level and excess, the holdout result.}
#'   \item{\code{summary(x)}}{RAPID with its Wilson interval, the spread of
#'     the flags over rare and common values, model performance.}
#'   \item{\code{as.data.frame(x)}}{The record table, with \code{at_risk_only =
#'     TRUE} the records at risk.}
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
#' The first seven are documented together in \code{\link{rapid-methods}}, the
#' plots in \code{\link{plot.rapid}}. \code{\link{rapid_qi_contribution}()} and
#' \code{\link{rapid_synthesizer_cv}()} start from the data and fit their own
#' models. They do not take a \code{rapid} object.
#'
#' @section Changes in 0.3.0:
#' Calls of riskutility 0.2.0 and the call that gives the same rule now:
#' \describe{
#'   \item{Numeric target, default}{\code{rapid(X, Y, kv, tv)} is now
#'     \code{rapid(X, Y, kv, tv, num_eval_method = "point", num_epsilon =
#'     0.1)}.}
#'   \item{Tolerance in percent}{\code{num_epsilon = 10} is now
#'     \code{num_epsilon = 0.1}.}
#'   \item{Absolute tolerance}{\code{num_epsilon_type = "absolute", num_epsilon
#'     = 2} is now \code{num_error_metric = "absolute", num_epsilon = 2}.}
#'   \item{Categorical target, default}{\code{rapid(X, Y, kv, tv)} is now
#'     \code{rapid(X, Y, kv, tv, cat_eval_method = "RCS_conditional")}.}
#'   \item{Categorical target, threshold 1}{\code{rapid(X, Y, kv, tv, cat_tau =
#'     1)} is now an error, because the default rule \code{RCS_marginal} takes a
#'     threshold below 1. The 0.2.0 rule is \code{rapid(X, Y, kv, tv, cat_tau = 1,
#'     cat_eval_method = "RCS_conditional")}.}
#'   \item{Linear model}{\code{model_type = "lm"} is now \code{model_type =
#'     "lm", num_scale = "identity"}.}
#'   \item{Tree}{\code{model_type = "cart"} is now \code{model_type = "cart",
#'     prune = FALSE, cp = 0.01, minsplit = 20, minbucket = 7}.}
#'   \item{Missing values of a numeric target}{\code{rapid(X, Y, kv, tv)}
#'     is now \code{rapid(X, Y, kv, tv, na_strategy = "constant")}.}
#'   \item{Records at risk only}{\code{return_all_records = FALSE} is now
#'     \code{as.data.frame(x, at_risk_only = TRUE)}.}
#'   \item{Old argument names}{\code{rapid(original_data = X, synthetic_data =
#'     Y, quasi_identifiers = kv, sensitive_attribute = tv)} still works and is
#'     the same as \code{rapid(X, Y, kv, tv)}.}
#' }
#' For the random forest no attacker setting changed. With the same random
#' stream (\code{set.seed()} before the call in 0.2.0, the same
#' \code{set.seed()} or \code{seed} in 0.3.0) and \code{num_eval_method =
#' "point", num_epsilon = 0.1}, the call reproduces the value of 0.2.0.
#'
#' The risk levels LOW, MEDIUM and HIGH of \code{summary()} are gone. Read
#' the excess over the chance level and \code{\link{rapid_test}()} instead.
#' Several fields and record columns of 0.2.0 are renamed or removed, and
#' reading them gives \code{NULL} without an error. The confidence is the
#' column \code{g} (0.2.0: \code{true_prob}), the baseline \code{b}
#' (\code{baseline}), the score of \code{RCS_conditional} \code{ratio}
#' (\code{relative_score}) and the error of the point rule \code{error}, a
#' proportion (\code{symmetric_error_pct}, in percent). The full list is in
#' the file NEWS.md and in the vignette section "Migration from 0.2.0".
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
#' @seealso \code{vignette("rapid", package = "riskutility")} for a worked
#'   example. \code{\link{rapid-methods}} and \code{\link{plot.rapid}} for the
#'   methods, \code{\link{rapid_test}}, \code{\link{rapid_threshold_select}},
#'   \code{\link{rapid_suite}}, \code{\link{rapid_learner}},
#'   \code{\link{rapid_qi_contribution}} and \code{\link{rapid_synthesizer_cv}}
#'   for the companion functions, \code{\link{dcap}} and \code{\link{tcap}} for
#'   attribution by matching, \code{\link{disclosure_report}}.
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
#' # Numeric target with the defaults: random forest, band rule, tolerance
#' # 10 percent, threshold 0.3
#' r <- rapid(X, Y, key_vars = kv, target_var = "income", seed = 1)
#' r
#' # RAPID is the share of records at risk. The chance level is the share at
#' # risk when the true values are permuted. Read the excess, their difference.
#' c(rapid = r$rapid, chance = r$chance, excess = r$excess)
#'
#' # The records at risk: tolerance band, confidence g, baseline b, gain r
#' head(as.data.frame(r, at_risk_only = TRUE), 3)
#'
#' # Interval for RAPID, and a stricter threshold without a new fit
#' confint(r)
#' update(r, num_tau = 0.5)
#'
#' # The result depends on the attacker. The linear model learns how income
#' # depends on age and sex and flags more records than the forest.
#' r_lm <- rapid(X, Y, key_vars = kv, target_var = "income", model_type = "lm",
#'               seed = 1)
#' r_lm
#'
#' # A release that says nothing about income: the band rule flags nothing, the
#' # point rule of 0.2.0 still flags records, as many as it does by chance
#' Y0 <- Y
#' Y0$income <- sample(Y0$income)
#' rapid(X, Y0, key_vars = kv, target_var = "income", model_type = "lm", seed = 1)
#' rapid(X, Y0, key_vars = kv, target_var = "income", model_type = "lm",
#'       num_eval_method = "point", seed = 1)
#'
#' # A tolerance in the units of the target: age within 5 years
#' rapid(X, Y, key_vars = c("sex", "income"), target_var = "age",
#'       num_error_metric = "absolute", num_epsilon = 5, seed = 1)
#'
#' # Categorical target. With two classes many records pass the threshold by
#' # chance, so the chance level is high and the excess is the number to read.
#' rapid(X, Y, key_vars = kv, target_var = "smoker", seed = 1)
#'
#' # Records 201 to 300 were not used for the release. RAPID on them is about as
#' # high as on the records behind the release. The attacker infers income from
#' # the quasi-identifiers and does not need the records themselves.
#' rh <- rapid(X[1:200, ], Y[1:200, ], key_vars = kv, target_var = "income",
#'             model_type = "lm", holdout = X[201:300, ], seed = 1)
#' rh
#' rh$release_specific
#'
#' # The same call as the first one, with the data in a synth_pair
#' sp <- synth_pair(X, Y, key_vars = kv, target_var = "income")
#' rapid(sp, seed = 1)
#'
#' \donttest{
#' # Lookup attacker: read it by its excess
#' rapid(X, Y, key_vars = kv, target_var = "income", model_type = "lookup", seed = 1)
#'
#' # Categorical target with a pruned tree
#' if (requireNamespace("rpart", quietly = TRUE)) {
#'   rapid(X, Y, key_vars = kv, target_var = "smoker", model_type = "cart", seed = 1)
#' }
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
                          na_strategy = c("drop", "constant", "median"),
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
  na_strategy <- .rapid_choice(na_strategy, c("drop", "constant", "median"), "na_strategy")

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
# rescores the stored attack when it is given nothing else. `seed` is not one of
# them: it seeds the fit as well.
.rapid_scoring_args <- c("num_eval_method", "num_error_metric", "num_epsilon", "num_delta",
                         "num_tau", "cat_eval_method", "cat_tau", "n_perm", "threshold_rule",
                         "perm_level")

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

  # Arguments that the rule or the threshold rule does not use
  tau_name <- if (is_categorical) "cat_tau" else "num_tau"
  used <- .rapid_rule_spec(rule$method)$threshold_arg
  if (tau_name %in% given && !identical(used, tau_name)) {
    .rapid_warn(c("Argument '", tau_name, "' is ignored: the rule \"", rule$method,
                  "\" has no gain threshold, its threshold is '", used, "'."),
                "ignored_argument")
  }
  if ("perm_level" %in% given && threshold_rule == "fixed") {
    .rapid_warn(c("Argument 'perm_level' is ignored: it sets the threshold under ",
                  "threshold_rule = \"permutation\" only, and the threshold here is fixed."),
                "ignored_argument")
  }
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
  # The permutation settings were checked by .rapid_scoring()
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
#'
#' The one place that turns permuted scores into shares, used by `rapid()`,
#' `rapid_test()`, `rapid_threshold_select()` and the plots.
#' @param null Matrix from `.rapid_null()`, or `NULL` (no permutations).
#' @noRd
.rapid_chance_rates <- function(rule, null) {
  if (is.null(null)) return(numeric(0))
  flags <- .rule_flag(rule, null)
  dim(flags) <- dim(null)
  colMeans(flags)
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
  if (!is.numeric(n_perm) || length(n_perm) != 1 || !is.finite(n_perm) || n_perm < 0 ||
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
                   "attack", "model", "predictor", "holdout", "release_specific", "seed",
                   "model_type", "attacker", "note", "model_metrics", "formula", "key_vars",
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
