#' Record Linkage Risk After Perturbation
#'
#' Measures targeted re-identification risk by linking each original record to
#' the most similar record(s) in a perturbed dataset using quasi-identifiers.
#' Supports deterministic (distance-based), probabilistic (Fellegi-Sunter),
#' PRAM (transition-matrix), and predictive (propensity-score-based) linkage
#' methods.
#'
#' @section Attacker scenario:
#' Targeted record linkage with membership knowledge and exact quasi-identifier knowledge:
#' \itemize{
#'   \item \strong{Membership knowledge:} the attacker knows the target individual has a
#'   record in the released dataset.
#'   \item \strong{Exact QI knowledge:} the attacker knows the target's quasi-identifiers
#'   exactly for the variables specified in \code{key}.
#'   \item \strong{Linking rule:} the attacker applies a strategy (see \code{strategy})
#'   based on distances, probabilistic scores, or transition probabilities.
#' }
#'
#' Ground truth for scoring is taken either from row alignment (\code{truth="row"})
#' or a shared identifier column (\code{truth="id"}). Note that \code{truth} affects only the
#' \emph{evaluation} of success; it is not an attacker capability.
#'
#' @section Distance measure:
#' For the deterministic method, record linkage is based on a weighted Gower-type
#' distance over the quasi-identifiers. For an original record \eqn{i} and a candidate
#' record \eqn{j} in the perturbed data, the distance is
#' \deqn{
#' d(i,j) =
#' \frac{\sum_{v \in \mathcal{V}_{ij}} w_v \, \delta_v(i,j)}
#'      {\sum_{v \in \mathcal{V}_{ij}} w_v},
#' \qquad 0 \le d(i,j) \le 1,
#' }
#' where \eqn{w_v \ge 0} are variable weights and
#' \eqn{\mathcal{V}_{ij}} denotes the subset of variables contributing
#' to the comparison (depending on missing-value handling).
#'
#' Variable-specific dissimilarities are defined as:
#' \itemize{
#'   \item Numeric / ordinal variables:
#'   \deqn{
#'   \delta_v(i,j) =
#'   \frac{|x_{iv} - x^{anon}_{jv}|}
#'        {\max(x_v) - \min(x_v)}
#'   }
#'   (absolute difference scaled to the unit interval).
#'
#'   \item Nominal variables:
#'   \deqn{
#'   \delta_v(i,j) =
#'   \begin{cases}
#'   0, & \text{if } x_{iv} = x^{anon}_{jv}, \\
#'   1, & \text{otherwise}.
#'   \end{cases}
#'   }
#' }
#'
#' The resulting distance is bounded in \eqn{[0,1]}, where 0 indicates
#' complete agreement and 1 maximal disagreement on the contributing
#' quasi-identifiers. If \code{na_anon = "ignore"}, variables with missing
#' values in either record are excluded from the numerator and denominator
#' for that pairwise comparison.
#'
#' @section Probabilistic method:
#' The Fellegi-Sunter probabilistic method estimates match (m) and
#' non-match (u) probabilities for each variable, then computes
#' log-likelihood ratios to score candidate pairs. The m-probabilities
#' represent the probability of agreement among true matched pairs, while
#' u-probabilities represent agreement among random non-matched pairs.
#'
#' User-supplied \code{m_probs} and \code{u_probs} override the supervised
#' estimation. The \code{fs_threshold} controls the minimum log-likelihood
#' ratio for inclusion in the candidate set.
#'
#' @section PRAM method:
#' For data protected via Post-Randomization Method (PRAM), risk is computed
#' directly from transition matrices. For each original record, the probability
#' of producing each anonymized record is computed as the product of per-variable
#' transition probabilities. This avoids distance computation entirely.
#'
#' @section Predictive method:
#' When \code{method = "predictive"}, a propensity model is fit on the stacked
#' original + anonymized data to predict record origin (original = 1,
#' synthetic = 0) from the quasi-identifiers. Records are then linked in the
#' 1D propensity-score space. The propensity score
#' \eqn{e(x) = P(\text{original} \mid \text{QIs} = x)} projects the multivariate
#' QI space into a single dimension that maximally separates the two datasets.
#'
#' When \code{pred_se = TRUE} (default) and \code{pred_model = "logit"}, the
#' bandwidth for kernel weighting is derived from the Fisher information
#' matrix of the logistic regression, giving each record a bandwidth equal to
#' its prediction standard error. This follows the principle from
#' Gaffert et al. (2015): records where the model is uncertain get wider
#' kernels, while confidently classified records get narrow kernels (higher
#' individual risk).
#'
#' @section RF method:
#' When \code{method = "rf"}, a random forest (ranger) is trained on the
#' stacked original + anonymized records to learn a supervised similarity.
#' Terminal-node co-occurrence across trees yields a proximity matrix in \eqn{[0, 1]}.
#' Each original record is linked to the anonymized record with highest
#' proximity.  Per-block RF is used by default; set \code{rf_global = TRUE}
#' to train a single global forest and restrict proximity lookups within blocks.
#' Blocks with fewer than 2 records per class fall back to a global RF
#' second pass.  Bijective matching is supported via the existing
#' \code{matching = "bijective"} mechanism with \code{maximize = TRUE}.
#'
#' Note: for the RF method, \code{risk} is the maximum proximity (similarity
#' in \eqn{[0, 1]}) to any candidate, not a re-identification probability.
#' The \code{d_min} and \code{d_true} fields store proximity values (higher
#' = more similar), which is inverted from the distance convention used by
#' other methods (lower = closer).
#'
#' @section RBRL method:
#' When \code{method = "rbrl"}, rank-based record linkage is used following
#' Muralidhar & Domingo-Ferrer (2016). Each numeric/ordinal variable is
#' replaced by its normalized rank (within each dataset independently),
#' yielding values in \eqn{[0, 1]}. Linkage distance is the weighted mean absolute
#' rank difference across variables. Nominal variables use exact matching
#' (0 if equal, 1 if not), consistent with Gower distance.
#'
#' RBRL is robust to monotone perturbations (additive noise, rounding,
#' top/bottom coding) because rank order is preserved under such
#' transformations. This makes it particularly useful for assessing risk
#' when the anonymization method preserves variable ordering.
#' All \code{strategy}, \code{risk_weighting}, \code{matching},
#' and \code{block} parameters work with RBRL.
#'
#' @section Mahalanobis method:
#' When \code{method = "mahalanobis"}, distances are computed using the
#' Mahalanobis metric, which accounts for correlations between variables.
#' Records that appear 'close' in Gower distance may be 'far' in
#' Mahalanobis distance if their difference goes against the data's
#' correlation structure, and vice versa.
#'
#' By default (\code{robust = TRUE}), the covariance matrix is estimated
#' using the Minimum Covariance Determinant (MCD) via
#' \code{robustbase::covMcd()}, with fallback to \code{MASS::cov.rob()}
#' then classical \code{cov()} if MCD fails. The covariance is always
#' estimated from the original data regardless of \code{direction},
#' modelling an attacker who knows the population covariance structure
#' (Templ & Meindl, 2008).
#'
#' For mixed data, numeric/ordinal variables use Mahalanobis distance
#' (normalized by a chi-squared threshold) while nominal variables use
#' Gower-style exact matching. The two components are combined as a
#' weighted average, where the weight is the proportion of each variable
#' type. For purely nominal data, use \code{method = "deterministic"}
#' instead.
#'
#' Ordinal factors are converted to integer codes (1, 2, 3, ...),
#' assuming approximately equal spacing between levels.
#'
#' @section OT matching:
#' When \code{matching = "ot"}, an entropy-regularized optimal transport
#' plan replaces the hard assignment of bijective matching.
#' The Sinkhorn-Knopp algorithm computes a transport plan \eqn{P}
#' that minimizes
#' \deqn{
#' \sum_{ij} P_{ij} C_{ij} + \varepsilon \sum_{ij} P_{ij} \log P_{ij}
#' }
#' subject to marginal constraints, where \eqn{C} is the cost matrix
#' from the chosen \code{method} and \eqn{\varepsilon} is the
#' regularization strength (\code{ot_epsilon}).
#'
#' Risk for record \eqn{i} is defined as
#' \eqn{P_{i,\text{true}} \times n_{\text{query}}}, where
#' \eqn{P_{i,\text{true}}} is the transport weight from query record
#' \eqn{i} to its true match. Under uniform random matching this gives
#' \eqn{1/n_{\text{search}}} (the same baseline as independent matching).
#'
#' OT matching produces \emph{continuous} risk scores in \eqn{[0, 1]},
#' interpolating between independent matching (smooth) and bijective
#' matching (hard assignment). As \eqn{\varepsilon \to 0}, the transport
#' plan converges to the Hungarian assignment; as
#' \eqn{\varepsilon \to \infty}, it approaches a uniform plan.
#' If \code{ot_epsilon = NULL} (default), epsilon is auto-calibrated
#' from the cost matrix as \code{median(C[C > 0]) * 0.05}.
#'
#' @section Softmax risk weighting:
#' When \code{risk_weighting = "softmax"}, closer candidates receive higher
#' attribution probability via:
#' \deqn{w_j = \frac{\exp(-\kappa \cdot d_j)}{\sum_k \exp(-\kappa \cdot d_k)}}
#' The temperature \code{kappa} is auto-calibrated from the distance range
#' if not supplied. This replaces the uniform \eqn{1/|candidate\_set|} risk.
#'
#' @section Kernel risk weighting:
#' When \code{risk_weighting = "kernel"}, donor probabilities are computed
#' via kernel-weighted distances. Each candidate \eqn{j} receives weight
#' proportional to a kernel function evaluated at its scaled distance:
#' \deqn{w_j = K(d_j / h) \Big/ \sum_k K(d_k / h)}
#' where \eqn{h} is a bandwidth and \eqn{K} is one of:
#' \itemize{
#'   \item \strong{gaussian} (default): \eqn{K(u) = \exp(-u^2/2)}
#'   \item \strong{epanechnikov}: \eqn{K(u) = \max(0, 1 - u^2) \cdot 3/4}
#'   \item \strong{tricube}: \eqn{K(u) = \max(0, (1 - |u|^3)^3) \cdot 70/81}
#' }
#' With compact kernels (Epanechnikov, tricube), candidates beyond the
#' bandwidth receive exactly zero weight and are effectively excluded.
#' The bandwidth is auto-selected via Silverman's rule if not supplied.
#'
#' @section Direction:
#' By default (\code{direction = "forward"}) the function loops over original
#' records and finds matches in the anonymized data. This quantifies how easily
#' each individual in the population can be re-identified.
#'
#' With \code{direction = "reverse"} the loop runs over anonymized records,
#' searching for matches in the original data. This gives a risk profile of
#' the data to be released: each row in \code{per_record} corresponds to a
#' released record and its probability of being correctly linked back to the
#' original. Deterministic, probabilistic, and predictive methods are symmetric
#' (distance/agreement does not depend on direction), so only the evaluation
#' perspective changes. The PRAM method is asymmetric: the transition matrix
#' \eqn{P(\text{output} \mid \text{input})} is looked up with swapped roles.
#'
#' @param X data.frame or \code{\link{synth_pair}} object. Original microdata.
#' @param x_anon data.frame. Perturbed/anonymized microdata.
#' @param key character. Names of quasi-identifier variables used for linkage.
#' @param method character. Linkage method: \code{"deterministic"} (default),
#'   \code{"probabilistic"} (Fellegi-Sunter), \code{"pram"} (transition matrix),
#'   \code{"predictive"} (propensity-score-based), \code{"rf"}
#'   (random forest proximity-based; requires \pkg{ranger}),
#'   \code{"rbrl"} (rank-based record linkage), or
#'   \code{"mahalanobis"} (Mahalanobis distance with robust covariance).
#' @param direction character. Direction of the linkage attack:
#'   \code{"forward"} (default) loops over original records and searches in the
#'   anonymized data, answering "how safe is each original individual?";
#'   \code{"reverse"} loops over anonymized records and searches in the
#'   original data, answering "how disclosive is each released record?".
#' @param risk_weighting character. How to weight candidates: \code{"uniform"}
#'   (default, 1/|set|), \code{"softmax"} (exponential distance-weighting),
#'   or \code{"kernel"} (kernel-weighted donor probabilities).
#' @param kappa numeric or NULL. Temperature parameter for softmax weighting.
#'   If NULL (default), auto-calibrated as \code{2 / range(distances)}.
#' @param bandwidth numeric or NULL. Bandwidth for kernel weighting.
#'   If NULL (default), auto-selected via Silverman's rule on the candidate
#'   distances.
#' @param kernel character. Kernel function for kernel weighting:
#'   \code{"gaussian"} (default), \code{"epanechnikov"}, or \code{"tricube"}.
#' @param truth character. How to define the true match for scoring:
#'   one of \code{"row"} (default) or \code{"id"}.
#'   \code{"row"} requires equal row counts and assumes row \eqn{i}
#'   in \code{X} corresponds to row \eqn{i} in \code{x_anon}.
#'   When datasets have unequal sizes, use \code{truth = "id"} with
#'   a shared identifier column.
#' @param id character or NULL. If \code{truth="id"}, name of an identifier column
#'   present in both \code{X} and \code{x_anon} that uniquely defines the true match.
#' @param type named character vector or NULL. Optional per-key type override:
#'   values in \code{c("numeric","ordinal","nominal")}. If NULL, inferred from \code{X}.
#' @param weights named numeric vector or NULL. Optional nonnegative weights for keys.
#'   If NULL, equal weights are used.
#' @param na_anon character. How to handle NA in the perturbed data during comparison:
#'   \code{"ignore"} (default), \code{"match"}, or \code{"mismatch"}.
#'   \itemize{
#'     \item \code{"ignore"} removes the variable from the pairwise distance denominator if
#'       either side is NA.
#'     \item \code{"match"} treats NA in \code{x_anon} as a perfect match (distance contribution 0).
#'     \item \code{"mismatch"} treats NA in \code{x_anon} as a full mismatch (distance contribution 1).
#'   }
#' @param strategy character. Adversary strategy variant:
#'   \code{"nearest"}, \code{"threshold"}, \code{"topk"}, \code{"topk_threshold"},
#'   or \code{"nearest_threshold"}.
#' @param k integer or NULL. Used when \code{strategy} is \code{"topk"} or \code{"topk_threshold"}.
#'   Ties at the cut-off distance may yield more than \code{k} candidates.
#' @param threshold numeric, list, or NULL. If numeric, absolute distance cutoff in \eqn{[0,1]}.
#'   If list, supports \code{list(type="quantile", p=...)} where \code{p} is in (0,1).
#'   The quantile is computed per record within its candidate set (and within-block if blocking is used).
#' @param block character or NULL. Optional subset of \code{key} used for exact blocking.
#'   Distances are computed only within blocks defined by these variables.
#' @param m_probs named numeric vector or NULL. User-supplied m-probabilities
#'   for the probabilistic method (probability of agreement given true match).
#' @param u_probs named numeric vector or NULL. User-supplied u-probabilities
#'   for the probabilistic method (probability of agreement given non-match).
#' @param fs_threshold numeric or NULL. Log-likelihood ratio threshold for the
#'   probabilistic method. Records with LR below this are excluded. If NULL,
#'   defaults to 0 (include all records with positive evidence of match).
#' @param pram_matrix named list of matrices. Transition matrices for the PRAM method.
#'   Each element is a square matrix where entry (i,j) is P(output=j | input=i).
#'   Names must correspond to variables in \code{key}.
#' @param pred_model character. Model for the predictive method:
#'   \code{"logit"} (default, logistic regression) or \code{"rf"} (random forest).
#' @param pred_se logical. For the predictive method with \code{pred_model="logit"},
#'   use the Fisher-information prediction SE as kernel bandwidth (default TRUE).
#'   Ignored for \code{"rf"} or when \code{bandwidth} is user-supplied.
#' @param return_matches logical. If TRUE, returns candidate indices per record (may be memory-heavy).
#' @param matching character. Matching mode: \code{"independent"} (default) scores
#'   each record independently (classical DBRL), \code{"bijective"} enforces
#'   one-to-one assignment via the Hungarian algorithm (GDBRL), and \code{"ot"}
#'   uses entropy-regularized optimal transport (Sinkhorn) for soft global
#'   assignment producing continuous risk. Bijective matching requires the
#'   \pkg{clue} package. See Herranz, Nin, Rodriguez & Tassa (2016).
#' @param risk_threshold numeric. Threshold for classifying records as "high risk"
#'   and for the \code{privacy_pass} flag (default 0.1). Records with risk above
#'   this threshold are counted in \code{n_high_risk} and \code{pct_high_risk}.
#'   The \code{privacy_pass} flag is TRUE when \code{mean_risk <= risk_threshold}.
#' @param n_trees integer. Number of trees for \code{method = "rf"} (default 500).
#'   Must be at least 10.
#' @param rf_global logical. For \code{method = "rf"} with blocking: if
#'   \code{TRUE}, train a single global forest on all data and restrict
#'   proximity computation within blocks. If \code{FALSE} (default), train
#'   a separate forest per block (with automatic global fallback for blocks
#'   that are too small).
#' @param robust logical. For \code{method = "mahalanobis"}: if \code{TRUE}
#'   (default), use the Minimum Covariance Determinant (MCD) estimator
#'   via \code{robustbase::covMcd()} for robust covariance estimation,
#'   with automatic fallback to \code{MASS::cov.rob()} or classical
#'   \code{cov()} if MCD fails. If \code{FALSE}, use classical covariance.
#'   Ignored for other methods.
#' @param ot_epsilon numeric or NULL. Regularization strength for OT matching.
#'   Smaller values produce sharper (more bijective-like) assignments;
#'   larger values produce smoother (more uniform) risk. If \code{NULL}
#'   (default), auto-calibrated as \code{median(C[C > 0]) * 0.05}.
#'   Ignored unless \code{matching = "ot"}.
#' @param ot_max_iter integer. Maximum Sinkhorn iterations for OT matching
#'   (default 100). Ignored unless \code{matching = "ot"}.
#' @param ... additional arguments passed to methods.
#' @author Roman Mueller and Matthias Templ
#'
#' @return An object of class \code{"recordLinkageRisk"}: a list with components
#'   \describe{
#'     \item{per_record}{data.frame with \code{n_query} rows and columns:
#'       \code{risk} (numeric, re-identification risk),
#'       \code{cand_n} (integer, candidate set size),
#'       \code{true_in_set} (logical),
#'       \code{d_true} (numeric, distance/score to true match),
#'       \code{d_min} (numeric, minimum distance),
#'       \code{d_rank} (integer, rank of true match among candidates),
#'       \code{risk_band} (factor with levels \code{"very_low"}, \code{"low"},
#'         \code{"moderate"}, \code{"high"}, \code{"very_high"},
#'         \code{"unique_match"}).
#'       When \code{matching = "bijective"}, an additional column
#'       \code{bijective_assigned} gives the search-side row index
#'       assigned by the Hungarian algorithm. In bijective mode,
#'       \code{risk} and \code{bijective_assigned} reflect the global
#'       one-to-one assignment, while \code{d_true}, \code{d_min}, and
#'       \code{d_rank} retain their per-record independent-scoring values
#'       from the main loop (useful for diagnostics).
#'       When \code{matching = "ot"}, \code{risk} is the continuous
#'       transport-weighted risk from the Sinkhorn plan and
#'       \code{d_rank} is the rank of the true match by transport weight.}
#'     \item{overall}{list with aggregate statistics including \code{risk_gini}
#'       (Gini coefficient of risk concentration).}
#'     \item{var_importance}{named numeric vector of per-variable importance.}
#'     \item{direction}{character, \code{"forward"} or \code{"reverse"}.}
#'     \item{matches}{(optional) list of integer vectors with candidate indices.}
#'   }
#'   When \code{direction = "forward"}, \code{per_record} has one row per
#'   original record; when \code{"reverse"}, one row per anonymized record.
#'
#'   Use \code{\link{top_at_risk}}, \code{\link{risk_by_group}},
#'   \code{\link{merge_per_record}}, and \code{\link{inspect_record}} for
#'   post-hoc per-record analysis.
#'
#' @examples
#' set.seed(1)
#' n <- 100
#' x <- data.frame(
#'   age    = sample(18:80, n, TRUE),
#'   sex    = factor(sample(c("f","m"), n, TRUE)),
#'   region = factor(sample(paste0("R",1:5), n, TRUE)),
#'   stringsAsFactors = FALSE
#' )
#'
#' # Simple perturbation: swap age within region
#' x_anon <- x
#' for (r in levels(x$region)) {
#'   idx <- which(x$region == r)
#'   x_anon$age[idx] <- sample(x$age[idx])
#' }
#'
#' # Deterministic nearest-neighbor (default)
#' res1 <- recordLinkage(x, x_anon, key = c("age","sex","region"))
#' print(res1)
#' summary(res1)
#'
#' # Reverse direction: risk per released record
#' res1r <- recordLinkage(x, x_anon, key = c("age","sex","region"),
#'                        direction = "reverse")
#' print(res1r)
#'
#' \donttest{
#' # Softmax distance-weighted risk
#' res2 <- recordLinkage(x, x_anon, key = c("age","sex","region"),
#'                       risk_weighting = "softmax")
#' print(res2)
#'
#' # Kernel-weighted risk (Gaussian kernel, auto bandwidth)
#' res2b <- recordLinkage(x, x_anon, key = c("age","sex","region"),
#'                        risk_weighting = "kernel")
#' print(res2b)
#'
#' # Kernel weighting with Epanechnikov kernel
#' res2c <- recordLinkage(x, x_anon, key = c("age","sex","region"),
#'                        risk_weighting = "kernel",
#'                        kernel = "epanechnikov")
#' print(res2c)
#'
#' # Probabilistic (Fellegi-Sunter) method
#' res3 <- recordLinkage(x, x_anon, key = c("age","sex","region"),
#'                       method = "probabilistic")
#' print(res3)
#'
#' # With synth_pair
#' pair <- synth_pair(x, x_anon, key_vars = c("age","sex","region"))
#' res4 <- recordLinkage(pair)
#' print(res4)
#'
#' # Predictive method (propensity-score-based)
#' res5 <- recordLinkage(x, x_anon, key = c("age","sex","region"),
#'                       method = "predictive")
#' print(res5)
#' plot(res5, which = 4)  # propensity distributions
#'
#' # Plot risk distribution
#' plot(res1)
#' plot(res1, which = 2)
#' }
#'
#' @references
#' Gower, J. C. (1971). A general coefficient of similarity and some of its
#' properties. Biometrics, 27(4), 857-871.
#'
#' Fellegi, I. P., & Sunter, A. B. (1969). A theory for record linkage.
#' Journal of the American Statistical Association, 64(328), 1183-1210.
#'
#' Domingo-Ferrer, J., & Torra, V. (2002). Validating distance-based record
#' linkage with probabilistic record linkage. Lecture Notes in Computer
#' Science, 2504, 207-215. Springer.
#'
#' Gaffert, P., Meinfelder, F., & Bosch, V. (2015). Towards an MI-proper
#' Predictive Mean Matching. Discussion Paper, University of Bamberg.
#'
#' Pagliuca, D., & Seri, G. (1999). Some results of individual ranking
#' method on the system of enterprise accounts annual survey (Esprit SDC Project,
#' Deliverable MI-3/D2). Unpublished technical report.
#'
#' Herranz, J., Nin, J., Rodriguez, P., & Tassa, T. (2016). Revisiting
#' distance-based record linkage for privacy-preserving release of statistical
#' datasets. Data & Knowledge Engineering, 100, 78-93.
#'
#' @seealso \code{\link{individual_risk}}, \code{\link{dcr}},
#'   \code{\link{nndr}}
#' @family privacy-models
#' @importFrom stats quantile mad
#' @importFrom graphics hist abline legend par plot points
#' @export
recordLinkage <- function(X, ...) UseMethod("recordLinkage")

#' @rdname recordLinkage
#' @export
recordLinkage.synth_pair <- function(X, ...) {
    if (is.null(X$key_vars))
        stop("synth_pair must have 'key_vars' for recordLinkage()", call. = FALSE)
    recordLinkage.default(
        X = X$original,
        x_anon = X$synthetic,
        key = X$key_vars,
        ...
    )
}

#' @rdname recordLinkage
#' @export
recordLinkage.default <- function(X,
                                  x_anon,
                                  key,
                                  method = c("deterministic", "probabilistic",
                                             "pram", "predictive", "rf",
                                             "rbrl", "mahalanobis"),
                                  direction = c("forward", "reverse"),
                                  risk_weighting = c("uniform", "softmax",
                                                     "kernel"),
                                  kappa = NULL,
                                  bandwidth = NULL,
                                  kernel = c("gaussian", "epanechnikov",
                                             "tricube"),
                                  truth = c("row", "id"),
                                  id = NULL,
                                  type = NULL,
                                  weights = NULL,
                                  na_anon = c("ignore", "match", "mismatch"),
                                  strategy = c("nearest", "threshold", "topk",
                                               "topk_threshold",
                                               "nearest_threshold"),
                                  k = NULL,
                                  threshold = NULL,
                                  block = NULL,
                                  m_probs = NULL,
                                  u_probs = NULL,
                                  fs_threshold = NULL,
                                  pram_matrix = NULL,
                                  pred_model = c("logit", "rf"),
                                  pred_se = TRUE,
                                  return_matches = FALSE,
                                  matching = c("independent", "bijective", "ot"),
                                  risk_threshold = 0.1,
                                  n_trees = 500L,
                                  rf_global = FALSE,
                                  robust = TRUE,
                                  ot_epsilon = NULL,
                                  ot_max_iter = 100L,
                                  ...) {

    method <- match.arg(method)
    direction <- match.arg(direction)
    risk_weighting <- match.arg(risk_weighting)
    kernel <- match.arg(kernel)
    pred_model <- match.arg(pred_model)
    truth <- match.arg(truth)
    na_anon <- match.arg(na_anon)
    strategy <- match.arg(strategy)
    matching <- match.arg(matching)

    if (matching == "bijective") {
        if (!requireNamespace("clue", quietly = TRUE))
            stop("Package 'clue' is required for bijective matching. ",
                 "Install it with install.packages('clue').", call. = FALSE)
        if (strategy != "nearest")
            message("Note: bijective matching uses the full candidate set; ",
                    "'strategy' applies only to independent-scoring diagnostics.")
        if (risk_weighting != "uniform")
            message("Note: bijective matching produces binary risk (0/1); ",
                    "'risk_weighting' applies only to independent-scoring diagnostics.")
    }

    if (matching == "ot") {
        if (risk_weighting != "uniform")
            message("Note: matching = 'ot' computes risk from the transport plan; ",
                    "'risk_weighting' applies only to independent-scoring diagnostics.")
    }

    if (method == "rf") {
        if (!requireNamespace("ranger", quietly = TRUE)) {
            stop("Package 'ranger' required for recordLinkage(method = 'rf'). ",
                 "Install with install.packages('ranger')", call. = FALSE)
        }
        if (!missing(weights)) {
            message("method = 'rf': weights ignored. ",
                    "RF uses data-driven variable importance instead.")
        }
        if (!missing(strategy) && strategy != "nearest") {
            message("method = 'rf': strategy '", strategy,
                    "' uses proximity-based scores [0,1], ",
                    "not distance-based thresholds.")
        }
    }

    stopifnot(is.data.frame(X), is.data.frame(x_anon))
    stopifnot(is.character(key), length(key) >= 1L)

    if (!all(key %in% names(X))) stop("Some 'key' not in 'X'.")
    if (!all(key %in% names(x_anon))) stop("Some 'key' not in 'x_anon'.")

    # PRAM requires transition matrix
    if (method == "pram") {
        if (is.null(pram_matrix) || !is.list(pram_matrix))
            stop("method='pram' requires 'pram_matrix' (named list of transition matrices).",
                 call. = FALSE)
        missing_vars <- setdiff(key, names(pram_matrix))
        if (length(missing_vars) > 0)
            stop("pram_matrix missing entries for: ",
                 paste(missing_vars, collapse = ", "), call. = FALSE)
    }

    n <- nrow(X)

    # ground truth mapping (forward, for FS estimation) ----
    if (truth == "row") {
        if (nrow(x_anon) != n)
            stop("If truth='row', x_anon must have same number of rows as X.")
        fs_true_idx <- seq_len(n)
    } else {
        if (is.null(id) || !is.character(id) || length(id) != 1L)
            stop("If truth='id', provide a single 'id' column name.")
        if (!(id %in% names(X)) || !(id %in% names(x_anon)))
            stop("'id' must exist in both X and x_anon.")
        fs_true_idx <- match(X[[id]], x_anon[[id]])
        if (anyNA(fs_true_idx)) stop("Some ids in X do not exist in x_anon.")
    }

    # direction-dependent query/search setup ----
    if (direction == "forward") {
        n_query <- n
        query_data <- X
        search_data <- x_anon
        true_idx <- fs_true_idx
    } else {
        n_query <- nrow(x_anon)
        query_data <- x_anon
        search_data <- X
        if (truth == "row") {
            true_idx <- seq_len(n_query)
        } else {
            true_idx <- match(x_anon[[id]], X[[id]])
            if (anyNA(true_idx))
                stop("Some ids in x_anon do not exist in X.")
        }
    }

    # infer / validate types ----
    if (is.null(type)) {
        type <- vapply(key, function(v) {
            xv <- X[[v]]
            if (is.numeric(xv)) "numeric"
            else if (is.ordered(xv)) "ordinal"
            else "nominal"
        }, character(1))
        names(type) <- key
    } else {
        if (is.null(names(type)) || !all(key %in% names(type)))
            stop("'type' must be a named character vector containing all keys.")
        type <- type[key]
        if (!all(type %in% c("numeric", "ordinal", "nominal")))
            stop("Invalid 'type' values.")
    }

    # weights ----
    if (is.null(weights)) {
        weights <- rep(1, length(key))
        names(weights) <- key
    } else {
        if (is.null(names(weights)) || !all(key %in% names(weights)))
            stop("'weights' must be named and contain all keys.")
        weights <- weights[key]
        if (any(!is.finite(weights)) || any(weights < 0))
            stop("'weights' must be nonnegative and finite.")
    }
    wsum <- sum(weights)
    if (wsum <= 0) stop("At least one weight must be > 0.")

    # numeric scaling ranges (for deterministic/probabilistic) ----
    rng <- lapply(key, function(v) {
        if (type[[v]] %in% c("numeric", "ordinal")) {
            allv <- c(X[[v]], x_anon[[v]])
            # ordinal factors: use integer codes for range
            if (is.factor(allv)) allv <- as.integer(allv)
            r <- range(allv, na.rm = TRUE)
            if (!is.finite(r[1]) || !is.finite(r[2]) || r[1] == r[2])
                r <- c(0, 1)
            r
        } else NULL
    })
    names(rng) <- key

    # Fellegi-Sunter: estimate m/u if probabilistic ----
    fs_params <- NULL
    if (method == "probabilistic") {
        if (is.null(m_probs) || is.null(u_probs)) {
            fs_est <- .fs_estimate(X, x_anon, key, type, fs_true_idx, rng)
            if (is.null(m_probs)) m_probs <- fs_est$m
            if (is.null(u_probs)) u_probs <- fs_est$u
        }
        if (is.null(fs_threshold)) fs_threshold <- 0
        fs_params <- list(m_probs = m_probs, u_probs = u_probs,
                          threshold_used = fs_threshold)
    }

    # Predictive: fit propensity model ----
    prop_fit <- NULL
    if (method == "predictive") {
        prop_fit <- .fit_propensity(X, x_anon, key,
                                    pred_model = pred_model, ...)
    }

    # precompute per-variable tolerance for probabilistic method ----
    fs_tol <- NULL
    if (method == "probabilistic") {
        fs_tol <- vapply(key, function(v) {
            if (type[[v]] %in% c("numeric", "ordinal")) {
                r <- rng[[v]]
                span <- r[2] - r[1]
                if (span <= 0) span <- 1
                tol_v <- stats::mad(X[[v]], na.rm = TRUE)
                if (!is.finite(tol_v) || tol_v <= 0)
                    tol_v <- 0.1 * span
                tol_v
            } else NA_real_
        }, numeric(1))
    }

    # RBRL: rank-transform numeric/ordinal variables ----
    rbrl_query <- NULL
    rbrl_search <- NULL
    rbrl_var_dist_acc <- NULL
    if (method == "rbrl") {
        rbrl_query <- query_data[, key, drop = FALSE]
        rbrl_search <- search_data[, key, drop = FALSE]
        for (v in key) {
            if (type[v] %in% c("numeric", "ordinal")) {
                qv <- rbrl_query[[v]]
                sv <- rbrl_search[[v]]
                if (is.factor(qv)) qv <- as.integer(qv)
                if (is.factor(sv)) sv <- as.integer(sv)
                rq <- rank(qv, ties.method = "average", na.last = "keep")
                rs <- rank(sv, ties.method = "average", na.last = "keep")
                nq <- sum(!is.na(rq))
                ns <- sum(!is.na(rs))
                rbrl_query[[v]] <- if (nq > 0) rq / nq else rq
                rbrl_search[[v]] <- if (ns > 0) rs / ns else rs
            }
            # Nominal: leave as-is (exact match in distance computation)
        }
        rbrl_var_dist_acc <- setNames(numeric(length(key)), key)
    }

    # Mahalanobis: prepare covariance from original data ----
    # Always use X (original) regardless of direction, because the attacker
    # model assumes knowledge of the population covariance (Templ & Meindl 2008).
    mahal_prep <- NULL
    if (method == "mahalanobis") {
        mahal_prep <- .mahal_prepare(X, key, type, robust = robust)
    }

    # blocking ----
    if (!is.null(block)) {
        if (!is.character(block) || length(block) < 1L)
            stop("'block' must be a character vector.")
        if (!all(block %in% key))
            stop("'block' must be a subset of 'key'.")
        blk_search <- .make_block_id(search_data, block)
        split_search <- split(seq_len(nrow(search_data)), blk_search)
        blk_query <- .make_block_id(query_data, block)
    } else {
        split_search <- list(all = seq_len(nrow(search_data)))
        blk_query <- rep("all", n_query)
    }

    # storage ----
    risk <- numeric(n_query)
    cand_n <- integer(n_query)
    true_in_set <- logical(n_query)
    d_true <- rep(NA_real_, n_query)
    d_min <- rep(NA_real_, n_query)
    d_rank <- rep(NA_integer_, n_query)
    matches <- if (isTRUE(return_matches)) vector("list", n_query) else NULL
    # Bijective matching: cache full score vectors per record
    score_cache <- if (matching %in% c("bijective", "ot")) vector("list", n_query) else NULL
    # Per-variable distance accumulator (deterministic only)
    var_dist_acc <- if (method == "deterministic") {
        setNames(numeric(length(key)), key)
    } else {
        NULL
    }

    # main loop ----
    if (method != "rf") {
    for (i in seq_len(n_query)) {
        cand <- split_search[[blk_query[i]]]
        if (is.null(cand) || length(cand) == 0L) {
            risk[i] <- 0
            cand_n[i] <- 0L
            true_in_set[i] <- FALSE
            next
        }

        if (method == "pram") {
            # PRAM: use transition probabilities directly
            pram_probs <- .pram_risk(
                x_row = query_data[i, key, drop = FALSE],
                anon_block = search_data[cand, key, drop = FALSE],
                key = key, pram_matrix = pram_matrix,
                reverse = (direction == "reverse")
            )

            tpos <- true_idx[i]

            # Candidate set: all records with nonzero probability
            nonzero <- which(pram_probs > 0)
            if (length(nonzero) == 0L) {
                risk[i] <- 0
                cand_n[i] <- 0L
                true_in_set[i] <- FALSE
                if (isTRUE(return_matches)) matches[[i]] <- integer(0)
                next
            }

            guess <- cand[nonzero]
            cand_n[i] <- length(guess)
            true_in_set[i] <- (tpos %in% guess)

            # d_true for PRAM: raw transition probability of true match
            j_in <- match(tpos, cand)
            if (!is.na(j_in)) {
                d_true[i] <- pram_probs[j_in]
                # Rank by descending probability (rank 1 = highest prob)
                d_rank[i] <- as.integer(rank(-pram_probs, ties.method = "min")[j_in])
            }

            if (!true_in_set[i]) {
                risk[i] <- 0
            } else {
                # Normalize probabilities over candidate set
                norm_probs <- pram_probs[nonzero] / sum(pram_probs[nonzero])
                true_local <- match(tpos, guess)
                risk[i] <- norm_probs[true_local]
            }

            if (isTRUE(return_matches)) matches[[i]] <- guess
            if (!is.null(score_cache))
                score_cache[[i]] <- list(cand = cand, scores = pram_probs,
                                          maximize = TRUE)
            next
        }

        if (method == "predictive") {
            # Distance in propensity-score space
            if (direction == "forward") {
                p_i <- prop_fit$p_original[i]
                p_cand <- prop_fit$p_synthetic[cand]
            } else {
                p_i <- prop_fit$p_synthetic[i]
                p_cand <- prop_fit$p_original[cand]
            }
            di <- abs(p_i - p_cand)

            d_min[i] <- min(di)
            tpos <- true_idx[i]
            j_in <- match(tpos, cand)
            if (!is.na(j_in)) {
                d_true[i] <- di[j_in]
                # Rank by ascending distance (rank 1 = closest)
                d_rank[i] <- as.integer(rank(di, ties.method = "min")[j_in])
            }

            # Use strategy to select guess set
            guess <- .choose_guess_set(
                d = di, cand = cand, strategy = strategy,
                k = k, threshold = threshold
            )

            cand_n[i] <- length(guess)
            if (cand_n[i] == 0L) {
                risk[i] <- 0
                true_in_set[i] <- FALSE
                if (isTRUE(return_matches)) matches[[i]] <- integer(0)
                next
            }

            true_in_set[i] <- (true_idx[i] %in% guess)

            if (!true_in_set[i]) {
                risk[i] <- 0
            } else {
                guess_local_idx <- match(guess, cand)
                d_guess <- di[guess_local_idx]

                # Determine bandwidth: Fisher SE if available
                bw <- bandwidth
                if (is.null(bw) && isTRUE(pred_se) &&
                    pred_model == "logit") {
                    if (direction == "forward" &&
                        !is.null(prop_fit$se_original)) {
                        bw <- prop_fit$se_original[i]
                    } else if (direction == "reverse" &&
                               !is.null(prop_fit$se_synthetic)) {
                        bw <- prop_fit$se_synthetic[i]
                    }
                }

                if (risk_weighting == "kernel") {
                    w <- .kernel_risk(d_guess, bandwidth = bw, kernel = kernel)
                } else if (risk_weighting == "softmax") {
                    w <- .softmax_risk(d_guess, kappa)
                } else {
                    w <- rep(1 / length(d_guess), length(d_guess))
                }
                true_local <- match(true_idx[i], guess)
                risk[i] <- w[true_local]
            }

            if (isTRUE(return_matches)) matches[[i]] <- guess
            if (!is.null(score_cache))
                score_cache[[i]] <- list(cand = cand, scores = di,
                                          maximize = FALSE)
            next
        }

        if (method == "probabilistic") {
            # Fellegi-Sunter: compute log-likelihood ratios
            lr <- .fs_log_lr(
                x_row = query_data[i, key, drop = FALSE],
                anon_block = search_data[cand, key, drop = FALSE],
                key = key, type = type, rng = rng,
                m_probs = m_probs, u_probs = u_probs,
                tol = fs_tol
            )

            d_min[i] <- NA_real_
            tpos <- true_idx[i]
            j_in <- match(tpos, cand)
            if (!is.na(j_in)) {
                d_true[i] <- lr[j_in]
                # Rank by descending LR (rank 1 = best match)
                d_rank[i] <- as.integer(rank(-lr, ties.method = "min")[j_in])
            }

            # Candidate set: records above threshold
            above <- which(lr >= fs_threshold)
            if (length(above) == 0L) {
                risk[i] <- 0
                cand_n[i] <- 0L
                true_in_set[i] <- FALSE
                if (isTRUE(return_matches)) matches[[i]] <- integer(0)
                next
            }

            guess <- cand[above]
            cand_n[i] <- length(guess)
            true_in_set[i] <- (tpos %in% guess)

            if (!true_in_set[i]) {
                risk[i] <- 0
            } else if (risk_weighting == "softmax") {
                # Softmax over LRs (higher LR = higher risk)
                lr_above <- lr[above]
                lr_shifted <- lr_above - max(lr_above)
                w <- exp(lr_shifted)
                w <- w / sum(w)
                true_local <- match(tpos, guess)
                risk[i] <- w[true_local]
            } else if (risk_weighting == "kernel") {
                # Kernel weighting over LRs
                # Transform LRs to distances: negate so higher LR = smaller d
                lr_above <- lr[above]
                pseudo_d <- max(lr_above) - lr_above
                w <- .kernel_risk(pseudo_d, bandwidth, kernel)
                true_local <- match(tpos, guess)
                risk[i] <- w[true_local]
            } else {
                risk[i] <- 1 / cand_n[i]
            }

            if (isTRUE(return_matches)) matches[[i]] <- guess
            if (!is.null(score_cache))
                score_cache[[i]] <- list(cand = cand, scores = lr,
                                          maximize = TRUE)
            next
        }

        if (method == "rbrl") {
            # RBRL: weighted mean absolute rank difference
            di <- numeric(length(cand))
            for (v in key) {
                wv <- weights[v]
                if (type[v] == "nominal") {
                    di <- di + wv * (rbrl_query[[v]][i] != rbrl_search[[v]][cand])
                } else {
                    di <- di + wv * abs(rbrl_query[[v]][i] -
                                        rbrl_search[[v]][cand])
                }
            }
            di <- di / wsum

            d_min[i] <- min(di)
            tpos <- true_idx[i]
            j_in <- match(tpos, cand)
            if (!is.na(j_in)) {
                d_true[i] <- di[j_in]
                d_rank[i] <- as.integer(rank(di, ties.method = "min")[j_in])
            }

            # Per-variable rank distance accumulation
            for (v in key) {
                if (type[v] == "nominal") {
                    rbrl_var_dist_acc[v] <- rbrl_var_dist_acc[v] +
                        min(as.numeric(rbrl_query[[v]][i] !=
                                       rbrl_search[[v]][cand]))
                } else {
                    rbrl_var_dist_acc[v] <- rbrl_var_dist_acc[v] +
                        min(abs(rbrl_query[[v]][i] - rbrl_search[[v]][cand]))
                }
            }

            guess <- .choose_guess_set(
                d = di, cand = cand, strategy = strategy,
                k = k, threshold = threshold
            )

            cand_n[i] <- length(guess)
            if (cand_n[i] == 0L) {
                risk[i] <- 0
                true_in_set[i] <- FALSE
                if (isTRUE(return_matches)) matches[[i]] <- integer(0)
                next
            }

            true_in_set[i] <- (true_idx[i] %in% guess)

            if (!true_in_set[i]) {
                risk[i] <- 0
            } else if (risk_weighting == "softmax") {
                guess_local_idx <- match(guess, cand)
                w <- .softmax_risk(di[guess_local_idx], kappa)
                true_local <- match(true_idx[i], guess)
                risk[i] <- w[true_local]
            } else if (risk_weighting == "kernel") {
                guess_local_idx <- match(guess, cand)
                w <- .kernel_risk(di[guess_local_idx], bandwidth, kernel)
                true_local <- match(true_idx[i], guess)
                risk[i] <- w[true_local]
            } else {
                risk[i] <- 1 / cand_n[i]
            }

            if (isTRUE(return_matches)) matches[[i]] <- guess
            if (!is.null(score_cache))
                score_cache[[i]] <- list(cand = cand, scores = di,
                                          maximize = FALSE)
            next
        }

        if (method == "mahalanobis") {
            # Mahalanobis distance
            di <- .mahal_dist(
                x_row = query_data[i, key, drop = FALSE],
                candidates = search_data[cand, key, drop = FALSE],
                prep = mahal_prep, type = type
            )

            d_min[i] <- min(di)

            tpos <- true_idx[i]
            j_in <- match(tpos, cand)
            if (!is.na(j_in)) {
                d_true[i] <- di[j_in]
                d_rank[i] <- as.integer(rank(di, ties.method = "min")[j_in])
            }

            guess <- .choose_guess_set(
                d = di, cand = cand, strategy = strategy,
                k = k, threshold = threshold
            )

            cand_n[i] <- length(guess)
            if (cand_n[i] == 0L) {
                risk[i] <- 0
                true_in_set[i] <- FALSE
                if (isTRUE(return_matches)) matches[[i]] <- integer(0)
                next
            }

            true_in_set[i] <- (true_idx[i] %in% guess)

            if (!true_in_set[i]) {
                risk[i] <- 0
            } else if (risk_weighting == "softmax") {
                guess_local_idx <- match(guess, cand)
                w <- .softmax_risk(di[guess_local_idx], kappa)
                true_local <- match(true_idx[i], guess)
                risk[i] <- w[true_local]
            } else if (risk_weighting == "kernel") {
                guess_local_idx <- match(guess, cand)
                w <- .kernel_risk(di[guess_local_idx], bandwidth, kernel)
                true_local <- match(true_idx[i], guess)
                risk[i] <- w[true_local]
            } else {
                risk[i] <- 1 / cand_n[i]
            }

            if (isTRUE(return_matches)) matches[[i]] <- guess
            if (!is.null(score_cache))
                score_cache[[i]] <- list(cand = cand, scores = di,
                                          maximize = FALSE)
            next
        }

        # --- Deterministic method ---
        di <- .dist_to_candidates(
            x_row = query_data[i, key, drop = FALSE],
            anon_block = search_data[cand, key, drop = FALSE],
            key = key, type = type, weights = weights, wsum = wsum,
            rng = rng, na_anon = na_anon
        )

        d_min[i] <- min(di)

        # distance to true (only if true is inside candidate set)
        tpos <- true_idx[i]
        j_in <- match(tpos, cand)
        if (!is.na(j_in)) {
            d_true[i] <- di[j_in]
            # Rank by ascending distance (rank 1 = closest)
            d_rank[i] <- as.integer(rank(di, ties.method = "min")[j_in])
        }

        # Accumulate per-variable distances for variable importance
        di_var <- .dist_per_variable(
            x_row = query_data[i, key, drop = FALSE],
            anon_block = search_data[cand, key, drop = FALSE],
            key = key, type = type, weights = weights, rng = rng,
            na_anon = na_anon
        )
        var_dist_acc <- var_dist_acc + di_var

        # attacker guess set
        guess <- .choose_guess_set(
            d = di, cand = cand, strategy = strategy,
            k = k, threshold = threshold
        )

        cand_n[i] <- length(guess)
        if (cand_n[i] == 0L) {
            risk[i] <- 0
            true_in_set[i] <- FALSE
            if (isTRUE(return_matches)) matches[[i]] <- integer(0)
            next
        }

        true_in_set[i] <- (true_idx[i] %in% guess)

        if (!true_in_set[i]) {
            risk[i] <- 0
        } else if (risk_weighting == "softmax") {
            # Softmax distance-weighted risk
            guess_local_idx <- match(guess, cand)
            w <- .softmax_risk(di[guess_local_idx], kappa)
            true_local <- match(true_idx[i], guess)
            risk[i] <- w[true_local]
        } else if (risk_weighting == "kernel") {
            # Kernel-weighted donor risk
            guess_local_idx <- match(guess, cand)
            w <- .kernel_risk(di[guess_local_idx], bandwidth, kernel)
            true_local <- match(true_idx[i], guess)
            risk[i] <- w[true_local]
        } else {
            risk[i] <- 1 / cand_n[i]
        }

        if (isTRUE(return_matches)) matches[[i]] <- guess
        if (!is.null(score_cache))
            score_cache[[i]] <- list(cand = cand, scores = di,
                                      maximize = FALSE)
    }
    } else if (method == "rf") {

    # RF method: proximity-based record linkage ----
    if (n_query + nrow(search_data) > 10000 && is.null(block))
        message("RF proximity with >10,000 records and no blocking variable ",
                "may be slow. Consider setting 'block'.")
    rf_var_importance <- NULL

    .rf_process_block <- function(q_idx, s_idx, cross_prox) {
        for (r in seq_along(q_idx)) {
            qi <- q_idx[r]
            prox_vec <- cross_prox[r, ]
            cand_n[qi] <<- length(s_idx)
            best_col <- which.max(prox_vec)
            risk[qi] <<- prox_vec[best_col]

            # Truth evaluation
            tpos <- true_idx[qi]
            if (!is.na(tpos) && tpos > 0L) {
                t_col <- match(tpos, s_idx)
                if (!is.na(t_col)) {
                    true_in_set[qi] <<- TRUE
                    d_true[qi] <<- prox_vec[t_col]
                    d_min[qi] <<- max(prox_vec)
                    d_rank[qi] <<- sum(prox_vec >= prox_vec[t_col])
                }
            }

            if (!is.null(score_cache)) {
                score_cache[[qi]] <<- list(cand = s_idx, scores = prox_vec,
                                           maximize = TRUE)
            }

            if (isTRUE(return_matches)) {
                matches[[qi]] <<- s_idx[best_col]
            }
        }
    }

    if (is.null(block) || length(split_search) == 1) {
        # No blocking: single RF
        block_res <- .rf_linkage_block(query_data, search_data, key,
                                       n_trees = n_trees, ...)
        s_idx <- seq_len(nrow(search_data))
        .rf_process_block(seq_len(n_query), s_idx, block_res$cross_prox)
        rf_var_importance <- block_res$var_importance

    } else {
        # Blocked matching
        all_blocks <- names(split_search)
        rf_var_importance <- NULL
        small_blocks <- 0L
        total_blocks <- length(all_blocks)
        fallback_blocks <- character(0)
        importance_list <- list()
        block_sizes <- integer(0)

        if (rf_global) {
            # Global RF: train once, restrict proximity within blocks
            global_rf <- .rf_proximity(query_data[, key, drop = FALSE],
                                       search_data[, key, drop = FALSE],
                                       n_trees = n_trees, importance = TRUE,
                                       ...)
            rf_var_importance <- global_rf$importance
            n_q_all <- nrow(query_data)

            for (blk in all_blocks) {
                s_idx <- split_search[[blk]]
                q_idx <- which(blk_query == blk)
                if (length(q_idx) == 0 || length(s_idx) == 0) next

                tn_q <- q_idx
                tn_s <- n_q_all + s_idx
                cross_prox <- .proximity_from_nodes(
                    global_rf$terminal_nodes, tn_s, tn_q
                )
                .rf_process_block(q_idx, s_idx, cross_prox)
            }

        } else {
            # Per-block RF
            for (blk in all_blocks) {
                s_idx <- split_search[[blk]]
                q_idx <- which(blk_query == blk)
                if (length(q_idx) == 0 || length(s_idx) == 0) next

                min_per_class <- min(length(q_idx), length(s_idx))

                if (min_per_class < 2) {
                    small_blocks <- small_blocks + 1L
                    fallback_blocks <- c(fallback_blocks, blk)
                    next
                }

                block_res <- .rf_linkage_block(
                    query_data[q_idx, , drop = FALSE],
                    search_data[s_idx, , drop = FALSE],
                    key, n_trees = n_trees, ...
                )
                .rf_process_block(q_idx, s_idx, block_res$cross_prox)
                importance_list[[blk]] <- block_res$var_importance
                block_sizes <- c(block_sizes, length(q_idx) + length(s_idx))
            }

            # Second pass: global RF fallback for blocks with < 2 per class
            if (length(fallback_blocks) > 0) {
                global_rf_fb <- .rf_proximity(
                    query_data[, key, drop = FALSE],
                    search_data[, key, drop = FALSE],
                    n_trees = n_trees, importance = TRUE, ...
                )
                n_q_all <- nrow(query_data)

                for (blk in fallback_blocks) {
                    s_idx <- split_search[[blk]]
                    q_idx <- which(blk_query == blk)
                    if (length(q_idx) == 0 || length(s_idx) == 0) next

                    cross_prox <- .proximity_from_nodes(
                        global_rf_fb$terminal_nodes,
                        n_q_all + s_idx, q_idx
                    )
                    .rf_process_block(q_idx, s_idx, cross_prox)
                }
            }

            # Aggregate variable importance (weighted by block size)
            if (length(importance_list) > 0) {
                imp_mat <- do.call(rbind, importance_list)
                rf_var_importance <- colSums(imp_mat * block_sizes) /
                    sum(block_sizes)
            }

            if (small_blocks > 0) {
                message(small_blocks, " of ", total_blocks,
                        " blocks have < 2 records per class. ",
                        "Consider rf_global = TRUE or coarser blocking.")
            }
        }
    }
    }

    # bijective matching override ----
    bijective_assigned <- NULL
    if (matching == "bijective") {
        bij <- .solve_bijective(score_cache, true_idx, n_query,
                                 split_search, blk_query)
        bijective_assigned <- bij$assigned_search
        risk <- bij$risk
        cand_n <- ifelse(bijective_assigned > 0L, 1L, 0L)
        true_in_set <- (bijective_assigned == true_idx)
        if (isTRUE(return_matches)) {
            matches <- lapply(bijective_assigned, function(a) {
                if (a > 0L) a else integer(0)
            })
        }
    }

    if (matching == "ot") {
        ot_res <- .solve_ot(score_cache, true_idx, n_query,
                             split_search, blk_query,
                             epsilon = ot_epsilon, max_iter = ot_max_iter)
        risk <- ot_res$risk
        d_rank <- ot_res$d_rank
        transport_plans <- ot_res$transport_plans
    }

    # risk_band ----
    risk_band <- cut(risk,
        breaks = c(-Inf, 0.05, 0.1, 0.2, 0.5, 1 - .Machine$double.eps, Inf),
        labels = c("very_low", "low", "moderate", "high", "very_high",
                    "unique_match"),
        right = TRUE)

    # var_importance ----
    if (method == "deterministic") {
        var_importance <- var_dist_acc / n_query
    } else if (method == "probabilistic" && !is.null(fs_params)) {
        var_importance <- log(fs_params$m_probs / fs_params$u_probs)
        names(var_importance) <- key
    } else if (method == "predictive" && !is.null(prop_fit)) {
        if (pred_model == "logit" && !is.null(prop_fit$model)) {
            # Use absolute model coefficients (excluding intercept)
            coefs <- stats::coef(prop_fit$model)[-1]
            var_importance <- .map_coefs_to_vars(names(coefs), key,
                                                  abs(coefs))
        } else if (pred_model == "rf" && !is.null(prop_fit$model) &&
                   !is.null(prop_fit$model$variable.importance)) {
            vi_raw <- prop_fit$model$variable.importance
            var_importance <- .map_coefs_to_vars(names(vi_raw), key,
                                                  vi_raw)
        } else {
            var_importance <- setNames(rep(NA_real_, length(key)), key)
        }
    } else if (method == "pram") {
        # Perturbation strength: mean off-diagonal probability
        var_importance <- setNames(numeric(length(key)), key)
        for (v in key) {
            tm <- pram_matrix[[v]]
            var_importance[v] <- 1 - mean(diag(tm))
        }
    } else if (method == "rf") {
        var_importance <- if (!is.null(rf_var_importance)) {
            rf_var_importance
        } else {
            setNames(rep(NA_real_, length(key)), key)
        }
    } else if (method == "rbrl") {
        # Mean per-variable rank distance to best match
        var_importance <- rbrl_var_dist_acc / n_query
    } else if (method == "mahalanobis") {
        # Diagonal of precision matrix, normalized to proportions (scale-free).
        # diag(Sigma^{-1})_j = 1/Var(x_j|x_{-j}): larger = more discriminating.
        vi <- setNames(rep(NA_real_, length(key)), key)
        if (!is.null(mahal_prep) && length(mahal_prep$numeric_keys) > 0) {
            diag_inv <- diag(mahal_prep$cov_inv)
            diag_inv <- diag_inv / sum(diag_inv)
            for (j in seq_along(mahal_prep$numeric_keys)) {
                vi[mahal_prep$numeric_keys[j]] <- diag_inv[j]
            }
        }
        var_importance <- vi
    } else {
        var_importance <- setNames(rep(NA_real_, length(key)), key)
    }

    overall <- list(
        mean_risk = mean(risk),
        max_risk = max(risk),
        n_high_risk = sum(risk > risk_threshold),
        pct_high_risk = mean(risk > risk_threshold),
        risk_quantiles = as.numeric(
            stats::quantile(risk, probs = c(0, .25, .5, .75, 1))
        ),
        pct_risk_gt0 = mean(risk > 0),
        pct_true_in_set = mean(true_in_set),
        mean_candidate_size = mean(cand_n),
        candidate_size_quantiles = as.numeric(
            stats::quantile(cand_n, probs = c(0, .25, .5, .75, 1))
        ),
        risk_gini = .gini(risk)
    )
    names(overall$risk_quantiles) <- c("min", "q25", "median", "q75", "max")
    names(overall$candidate_size_quantiles) <- c("min", "q25", "median",
                                                  "q75", "max")

    per_rec <- data.frame(
        risk = risk,
        cand_n = cand_n,
        true_in_set = true_in_set,
        d_true = d_true,
        d_min = d_min,
        d_rank = d_rank,
        risk_band = risk_band
    )
    if (!is.null(bijective_assigned))
        per_rec$bijective_assigned <- bijective_assigned

    out <- list(
        per_record = per_rec,
        overall = overall,
        privacy_pass = overall$mean_risk <= risk_threshold,
        n_original = nrow(X),
        n_synthetic = nrow(x_anon),
        n_query = n_query,
        direction = direction,
        key_vars = key,
        method = method,
        settings = list(
            key = key,
            truth = truth,
            id = id,
            type = type,
            weights = weights,
            na_anon = na_anon,
            strategy = strategy,
            k = k,
            threshold = threshold,
            block = block,
            direction = direction,
            risk_weighting = risk_weighting,
            kappa = kappa,
            bandwidth = bandwidth,
            kernel = kernel,
            pred_model = pred_model,
            pred_se = pred_se,
            matching = matching,
            risk_threshold = risk_threshold,
            n_trees = if (method == "rf") n_trees else NULL,
            rf_global = if (method == "rf") rf_global else NULL,
            robust = if (method == "mahalanobis") robust else NULL,
            ot_epsilon = if (matching == "ot") ot_epsilon else NULL,
            ot_max_iter = if (matching == "ot") ot_max_iter else NULL
        )
    )
    if (isTRUE(return_matches)) out$matches <- matches
    if (matching == "ot" && exists("transport_plans"))
        out$transport_plans <- transport_plans
    if (!is.null(fs_params)) out$fs_params <- fs_params
    if (method == "predictive" && !is.null(prop_fit)) {
        out$propensity_info <- list(
            pred_model = pred_model,
            pred_se = pred_se,
            mean_propensity_original = mean(prop_fit$p_original),
            mean_propensity_synthetic = mean(prop_fit$p_synthetic),
            propensity_original = prop_fit$p_original,
            propensity_synthetic = prop_fit$p_synthetic
        )
    }
    if (method == "pram") {
        out$pram_info <- list(
            variables_used = key,
            mean_transition_prob = mean(vapply(key, function(v) {
                mean(diag(pram_matrix[[v]]))
            }, numeric(1)))
        )
    }
    out$var_importance <- var_importance

    class(out) <- "recordLinkageRisk"
    out
}


# ---------- per-record risk helpers ----------

#' Return Highest-Risk Records
#'
#' Returns the top-n riskiest records from a \code{recordLinkageRisk} object
#' with their per-record diagnostics.
#'
#' @param x object of class \code{"recordLinkageRisk"}.
#' @param n integer, number of records to return (default 10).
#' @param data optional data frame (original or synthetic depending on
#'   direction) whose key-variable columns are appended.
#' @param ... ignored.
#' @return A data frame with per-record diagnostics for the top-n riskiest
#'   records, including a \code{record_id} column.
#' @family privacy-models
#' @author Matthias Templ, Oscar Thees
#' @export
top_at_risk <- function(x, ...) UseMethod("top_at_risk")

#' @rdname top_at_risk
#' @export
top_at_risk.recordLinkageRisk <- function(x, n = 10, data = NULL, ...) {
    pr <- x$per_record
    n <- min(n, nrow(pr))
    idx <- order(pr$risk, decreasing = TRUE)[seq_len(n)]
    out <- pr[idx, , drop = FALSE]
    out$record_id <- idx
    if (!is.null(data)) {
        out <- cbind(out, data[idx, x$key_vars, drop = FALSE])
    }
    rownames(out) <- NULL
    out
}


#' Aggregate Risk by Group
#'
#' Computes per-group risk statistics from a \code{recordLinkageRisk} object.
#'
#' @param x object of class \code{"recordLinkageRisk"}.
#' @param group a vector of group labels (same length as number of records),
#'   or a single column name if \code{data} is provided.
#' @param data optional data frame from which to extract the grouping column.
#' @param ... ignored.
#' @return A data frame with columns \code{mean_risk}, \code{max_risk},
#'   \code{n}, \code{n_high}, \code{pct_high}, and the grouping variable,
#'   sorted by \code{mean_risk} descending.
#' @family privacy-models
#' @author Matthias Templ, Oscar Thees
#' @export
risk_by_group <- function(x, ...) UseMethod("risk_by_group")

#' @rdname risk_by_group
#' @export
risk_by_group.recordLinkageRisk <- function(x, group, data = NULL, ...) {
    if (is.character(group) && length(group) == 1L && !is.null(data)) {
        g <- data[[group]]
        gname <- group
    } else {
        g <- group
        gname <- "group"
    }
    stopifnot(length(g) == nrow(x$per_record))
    thresh <- x$settings$risk_threshold
    risk <- x$per_record$risk
    res <- tapply(risk, g, function(r) {
        c(mean_risk = mean(r), max_risk = max(r), n = length(r),
          n_high = sum(r > thresh), pct_high = 100 * mean(r > thresh))
    })
    out <- as.data.frame(do.call(rbind, res))
    out[[gname]] <- rownames(out)
    rownames(out) <- NULL
    out[order(out$mean_risk, decreasing = TRUE), ]
}


#' Merge Per-Record Risks Back to Data
#'
#' Joins per-record risk diagnostics back to the original (or synthetic)
#' data frame.
#'
#' @param x object of class \code{"recordLinkageRisk"}.
#' @param data data frame with the same number of rows as \code{x$per_record}.
#' @param ... ignored.
#' @return A data frame combining \code{data} and \code{x$per_record}.
#' @family privacy-models
#' @author Matthias Templ, Oscar Thees
#' @export
merge_per_record <- function(x, ...) UseMethod("merge_per_record")

#' @rdname merge_per_record
#' @export
merge_per_record.recordLinkageRisk <- function(x, data, ...) {
    stopifnot(nrow(data) == nrow(x$per_record))
    cbind(data, x$per_record)
}


#' Inspect a Single Record's Linkage Detail
#'
#' Returns detailed linkage diagnostics for a single record, including
#' candidate IDs and (optionally) the data for the query record and its
#' candidates. Requires \code{recordLinkage(..., return_matches = TRUE)}.
#'
#' @param x object of class \code{"recordLinkageRisk"}.
#' @param i integer, the record index to inspect.
#' @param data_orig optional data frame of the original records (query side).
#' @param data_anon optional data frame of the anonymized records (search side).
#' @param ... ignored.
#' @return An object of class \code{"inspect_record"} containing:
#'   \describe{
#'     \item{record_id}{the record index}
#'     \item{risk}{re-identification risk}
#'     \item{d_rank}{rank of the true match among candidates}
#'     \item{risk_band}{categorical risk band}
#'     \item{d_true, d_min}{distances for the true match and closest candidate}
#'     \item{n_candidates}{number of candidates}
#'     \item{true_in_set}{whether the true match is in the candidate set}
#'     \item{candidate_ids}{vector of candidate indices}
#'     \item{query_record}{optional: data for the query record}
#'     \item{candidate_records}{optional: data for the candidate records}
#'   }
#' @family privacy-models
#' @author Matthias Templ, Oscar Thees
#' @export
inspect_record <- function(x, ...) UseMethod("inspect_record")

#' @rdname inspect_record
#' @export
inspect_record.recordLinkageRisk <- function(x, i, data_orig = NULL,
                                              data_anon = NULL, ...) {
    if (is.null(x$matches))
        stop("'inspect_record()' requires recordLinkage(..., return_matches = TRUE)")
    if (!is.numeric(i) || length(i) != 1L || i < 1L || i > nrow(x$per_record))
        stop("'i' must be an integer between 1 and ", nrow(x$per_record))
    pr <- x$per_record[i, ]
    candidates <- x$matches[[i]]
    out <- list(
        record_id = i,
        risk = pr$risk,
        d_rank = pr$d_rank,
        risk_band = as.character(pr$risk_band),
        d_true = pr$d_true,
        d_min = pr$d_min,
        n_candidates = pr$cand_n,
        true_in_set = pr$true_in_set,
        candidate_ids = candidates
    )
    if (!is.null(pr$bijective_assigned))
        out$bijective_assigned <- pr$bijective_assigned
    if (!is.null(data_orig))
        out$query_record <- data_orig[i, , drop = FALSE]
    if (!is.null(data_anon) && length(candidates) > 0)
        out$candidate_records <- data_anon[candidates, , drop = FALSE]
    class(out) <- "inspect_record"
    out
}

#' @rdname inspect_record
#' @param x object of class \code{"inspect_record"}.
#' @export
print.inspect_record <- function(x, ...) {
    cat("Record Linkage Inspection: Record", x$record_id, "\n")
    cat(strrep("=", 45), "\n\n")
    cat(sprintf("  Risk:        %.4f (%s)\n", x$risk, x$risk_band))
    cat(sprintf("  True rank:   %s\n",
                if (is.na(x$d_rank)) "N/A" else as.character(x$d_rank)))
    cat(sprintf("  d_true:      %s\n",
                if (is.na(x$d_true)) "N/A" else sprintf("%.4f", x$d_true)))
    cat(sprintf("  d_min:       %s\n",
                if (is.na(x$d_min)) "N/A" else sprintf("%.4f", x$d_min)))
    cat(sprintf("  Candidates:  %d\n", x$n_candidates))
    cat(sprintf("  True in set: %s\n", x$true_in_set))
    if (!is.null(x$query_record)) {
        cat("\nQuery record:\n")
        print(x$query_record)
    }
    if (!is.null(x$candidate_records)) {
        cat("\nCandidate records:\n")
        print(x$candidate_records)
    }
    invisible(x)
}


# ---------- internal helpers ----------

#' Solve bijective (one-to-one) assignment via the Hungarian algorithm
#'
#' Builds a cost matrix per blocking group and solves the LSAP so that each
#' query record is assigned to at most one search record. Risk is binary:
#' 1 if assigned to the true match, 0 otherwise.
#'
#' @param score_cache list of length n_query.
#'   Each element is \code{list(cand, scores, maximize)}.
#' @param true_idx integer vector of true search-side indices.
#' @param n_query integer, number of query records.
#' @param split_search named list of search-side indices per block.
#' @param blk_query character vector of block labels per query record.
#' @return list with \code{assigned_search} (integer) and \code{risk} (numeric).
#' @keywords internal
.solve_bijective <- function(score_cache, true_idx, n_query,
                              split_search, blk_query) {
    assigned <- integer(n_query)
    risk_out <- numeric(n_query)

    blocks <- unique(blk_query)
    for (blk in blocks) {
        q_idx <- which(blk_query == blk)
        s_idx <- split_search[[blk]]
        if (is.null(s_idx) || length(s_idx) == 0L || length(q_idx) == 0L)
            next

        nq <- length(q_idx)
        ns <- length(s_idx)

        # Determine direction from first non-NULL cache entry in block
        maximize <- FALSE
        for (qi in q_idx) {
            if (!is.null(score_cache[[qi]])) {
                maximize <- score_cache[[qi]]$maximize
                break
            }
        }

        # Build cost matrix: rows = query records, cols = search records
        BIG <- 1e12
        dim_n <- max(nq, ns)
        cost <- matrix(BIG, nrow = dim_n, ncol = dim_n)

        for (r in seq_along(q_idx)) {
            qi <- q_idx[r]
            sc <- score_cache[[qi]]
            if (is.null(sc)) next

            # Map cached candidate indices to column positions
            col_pos <- match(sc$cand, s_idx)
            valid <- !is.na(col_pos)
            if (!any(valid)) next

            if (maximize) {
                # Transform to minimization: cost = max_score - score
                max_s <- max(sc$scores[valid])
                cost[r, col_pos[valid]] <- max_s - sc$scores[valid]
            } else {
                cost[r, col_pos[valid]] <- sc$scores[valid]
            }
        }

        # Solve LSAP (linear sum assignment problem)
        sol <- clue::solve_LSAP(cost)

        # Map solution back
        for (r in seq_along(q_idx)) {
            qi <- q_idx[r]
            assigned_col <- sol[r]
            if (assigned_col <= ns && cost[r, assigned_col] < BIG) {
                assigned[qi] <- s_idx[assigned_col]
                risk_out[qi] <- if (assigned[qi] == true_idx[qi]) 1 else 0
            }
            # else: assigned to dummy column -> remains 0
        }
    }

    list(assigned_search = assigned, risk = risk_out)
}


#' @keywords internal
.gini <- function(x) {
    x <- sort(x)
    n <- length(x)
    if (n == 0L || sum(x) == 0) return(0)
    2 * sum(x * seq_len(n)) / (n * sum(x)) - (n + 1) / n
}

#' @keywords internal
# Map model coefficient names to key variables (handles prefix collisions)
.map_coefs_to_vars <- function(coef_names, key, values) {
    out <- setNames(numeric(length(key)), key)
    # Sort key vars longest-first to avoid prefix collisions
    # (e.g. "age_group" before "age")
    key_sorted <- key[order(nchar(key), decreasing = TRUE)]
    claimed <- logical(length(coef_names))
    for (v in key_sorted) {
        idx <- which(!claimed & startsWith(coef_names, v))
        if (length(idx) > 0L) {
            out[v] <- mean(values[idx])
            claimed[idx] <- TRUE
        }
    }
    out
}

#' @keywords internal
.make_block_id <- function(df, vars) {
    tmp <- df[vars]
    for (v in vars) {
        tmp[[v]] <- ifelse(is.na(tmp[[v]]), "<NA>", as.character(tmp[[v]]))
    }
    do.call(paste, c(tmp, sep = "\r"))
}

#' @keywords internal
.dist_to_candidates <- function(x_row, anon_block,
                                key, type, weights,
                                wsum, rng, na_anon) {
    m <- nrow(anon_block)
    acc <- numeric(m)
    denom <- rep(wsum, m)

    for (v in key) {
        w <- weights[[v]]
        xv <- x_row[[v]]
        av <- anon_block[[v]]

        # NA default: "match" -> 0, "mismatch" -> 1, "ignore" -> 0
        dv <- rep(if (na_anon == "mismatch") 1 else 0, m)
        ok <- !is.na(av) & !is.na(xv)

        if (type[[v]] %in% c("numeric", "ordinal")) {
            r <- rng[[v]]
            span <- (r[2] - r[1])
            if (span <= 0) span <- 1
            # ordinal factors: use integer codes for arithmetic
            xv_num <- if (is.factor(xv)) as.integer(xv) else xv
            av_num <- if (is.factor(av)) as.integer(av) else av
            dv[ok] <- abs(av_num[ok] - xv_num) / span
        } else { # nominal
            # character cast avoids "level sets of factors are different"
            dv[ok] <- ifelse(as.character(av[ok]) == as.character(xv), 0, 1)
        }

        if (na_anon == "ignore") {
            miss <- is.na(av) | is.na(xv)
            denom[miss] <- denom[miss] - w
        }

        acc <- acc + w * dv
    }

    denom[denom <= 0] <- 1
    acc / denom
}

#' @keywords internal
.dist_per_variable <- function(x_row, anon_block, key, type, weights,
                               rng, na_anon) {
    # Returns named numeric: mean weighted distance per variable across candidates
    m <- nrow(anon_block)
    out <- setNames(numeric(length(key)), key)

    for (v in key) {
        w <- weights[[v]]
        xv <- x_row[[v]]
        av <- anon_block[[v]]

        dv <- rep(if (na_anon == "mismatch") 1 else 0, m)
        ok <- !is.na(av) & !is.na(xv)

        if (type[[v]] %in% c("numeric", "ordinal")) {
            r <- rng[[v]]
            span <- (r[2] - r[1])
            if (span <= 0) span <- 1
            xv_num <- if (is.factor(xv)) as.integer(xv) else xv
            av_num <- if (is.factor(av)) as.integer(av) else av
            dv[ok] <- abs(av_num[ok] - xv_num) / span
        } else {
            dv[ok] <- ifelse(as.character(av[ok]) == as.character(xv), 0, 1)
        }

        # When na_anon="ignore", average only over non-NA pairs
        if (na_anon == "ignore") {
            ok_count <- sum(ok)
            out[v] <- if (ok_count > 0L) w * mean(dv[ok]) else 0
        } else {
            out[v] <- w * mean(dv)
        }
    }
    out
}

#' @keywords internal
.choose_guess_set <- function(d, cand, strategy, k = NULL,
                              threshold = NULL) {
    # Resolve tau if supplied
    tau <- NULL
    if (!is.null(threshold)) {
        if (is.numeric(threshold) && length(threshold) == 1L) {
            tau <- threshold
        } else if (is.list(threshold) && is.character(threshold$type)) {
            if (threshold$type == "quantile") {
                p <- threshold$p
                if (!is.numeric(p) || length(p) != 1L || p <= 0 || p >= 1)
                    stop("Invalid quantile p.")
                tau <- as.numeric(stats::quantile(d, probs = p,
                                                  names = FALSE, type = 7))
            } else {
                stop("Only threshold$type = 'quantile' is supported.")
            }
        } else {
            stop("Invalid 'threshold' format.")
        }
    }

    if (strategy %in% c("threshold", "topk_threshold",
                         "nearest_threshold") && is.null(tau))
        stop("This strategy requires 'threshold' (numeric or list(type='quantile', ...)).")
    if (strategy %in% c("topk", "topk_threshold")) {
        if (is.null(k)) stop("This strategy requires 'k'.")
        k <- as.integer(k)
        if (!is.finite(k) || k < 1L) stop("'k' must be >= 1.")
    }

    .topk_with_ties <- function(dv, idx, k) {
        ord <- order(dv[idx])
        if (length(ord) <= k) return(idx)
        kth <- dv[idx][ord[k]]
        idx[dv[idx] <= kth]
    }

    if (strategy == "nearest") {
        dmin <- min(d)
        idx <- which(d == dmin)
        return(cand[idx])
    }

    if (strategy == "nearest_threshold") {
        dmin <- min(d)
        if (dmin > tau) return(integer(0))
        idx <- which(d == dmin)
        return(cand[idx])
    }

    if (strategy == "threshold") {
        idx <- which(d <= tau)
        if (length(idx) == 0L) return(integer(0))
        return(cand[idx])
    }

    if (strategy == "topk") {
        idx <- seq_along(d)
        idx <- .topk_with_ties(d, idx, k)
        return(cand[idx])
    }

    if (strategy == "topk_threshold") {
        idx <- which(d <= tau)
        if (length(idx) == 0L) return(integer(0))
        idx <- .topk_with_ties(d, idx, k)
        return(cand[idx])
    }

    stop("Unknown strategy.")
}

#' @keywords internal
.softmax_risk <- function(d, kappa = NULL) {
    if (length(d) == 0L) return(numeric(0))
    if (length(d) == 1L) return(1)
    if (is.null(kappa)) {
        d_range <- diff(range(d))
        if (d_range < .Machine$double.eps) return(rep(1 / length(d), length(d)))
        kappa <- 2 / d_range
    }
    neg_kd <- -kappa * d
    neg_kd <- neg_kd - max(neg_kd)  # numerical stability
    w <- exp(neg_kd)
    w / sum(w)
}

#' @keywords internal
.kernel_risk <- function(d, bandwidth = NULL,
                             kernel = "gaussian") {
    n <- length(d)
    if (n == 0L) return(numeric(0))
    if (n == 1L) return(1)

    # Auto-bandwidth via Silverman's rule
    if (is.null(bandwidth)) {
        s <- stats::sd(d)
        iqr <- stats::IQR(d)
        if (!is.finite(s) || s <= 0) s <- 1
        if (!is.finite(iqr) || iqr <= 0) iqr <- s * 1.34
        bandwidth <- 0.9 * min(s, iqr / 1.34) * n^(-1/5)
        if (bandwidth < .Machine$double.eps)
            bandwidth <- max(d) / 2
        if (bandwidth < .Machine$double.eps)
            return(rep(1 / n, n))
    }

    u <- d / bandwidth

    w <- switch(kernel,
        gaussian = exp(-0.5 * u^2),
        epanechnikov = pmax(0, 1 - u^2) * 0.75,
        tricube = pmax(0, (1 - abs(u)^3)^3) * (70 / 81),
        stop("Unknown kernel: ", kernel)
    )

    wsum <- sum(w)
    if (wsum < .Machine$double.eps) return(rep(1 / n, n))
    w / wsum
}

#' @keywords internal
.fit_propensity <- function(X, x_anon, key, pred_model = "logit", ...) {
    X_key <- X[, key, drop = FALSE]
    A_key <- x_anon[, key, drop = FALSE]
    n_orig <- nrow(X_key)
    n_anon <- nrow(A_key)

    stacked <- rbind(X_key, A_key)
    stacked$.label <- c(rep(1L, n_orig), rep(0L, n_anon))

    # Convert characters to factors for GLM
    for (v in key) {
        if (is.character(stacked[[v]])) stacked[[v]] <- factor(stacked[[v]])
    }
    # Drop unused factor levels
    for (v in key) {
        if (is.factor(stacked[[v]])) stacked[[v]] <- droplevels(stacked[[v]])
    }

    formula <- stats::as.formula(paste(".label ~", paste(key, collapse = " + ")))

    if (pred_model == "logit") {
        fit <- stats::glm(formula, data = stacked, family = stats::binomial(),
                          control = list(maxit = 50))
        p_all <- stats::predict(fit, type = "response")

        # Fisher-information-based SE: sqrt(diag(X %*% vcov %*% t(X)))
        X_mat <- stats::model.matrix(fit)
        V <- stats::vcov(fit)
        # SE on linear predictor (logit) scale
        se_link <- sqrt(pmax(0, rowSums((X_mat %*% V) * X_mat)))
        # Transform to response scale via delta method: se_p = p*(1-p)*se_link
        se_all <- p_all * (1 - p_all) * se_link
        se_all <- pmax(se_all, .Machine$double.eps)

        list(
            p_original = unname(p_all[seq_len(n_orig)]),
            p_synthetic = unname(p_all[n_orig + seq_len(n_anon)]),
            se_original = unname(se_all[seq_len(n_orig)]),
            se_synthetic = unname(se_all[n_orig + seq_len(n_anon)]),
            model = fit
        )

    } else if (pred_model == "rf") {
        if (!requireNamespace("ranger", quietly = TRUE))
            stop("Package 'ranger' required for pred_model='rf'. ",
                 "Install with: install.packages('ranger')", call. = FALSE)

        stacked$.label <- factor(stacked$.label, levels = c("0", "1"))
        user_args <- list(...)
        default_args <- list(formula = formula, data = stacked,
                             probability = TRUE, num.trees = 500,
                             importance = "impurity")
        args <- modifyList(default_args, user_args)
        fit <- do.call(ranger::ranger, args)
        p_all <- stats::predict(fit, data = stacked)$predictions[, "1"]

        list(
            p_original = p_all[seq_len(n_orig)],
            p_synthetic = p_all[n_orig + seq_len(n_anon)],
            se_original = NULL,
            se_synthetic = NULL,
            model = fit
        )
    } else {
        stop("pred_model must be 'logit' or 'rf'.", call. = FALSE)
    }
}

#' @keywords internal
.fs_estimate <- function(X, x_anon, key, type, true_idx, rng,
                   n_sample = 500L) {
    n <- nrow(X)
    nv <- length(key)

    # degenerate case: cannot estimate from a single record
    if (n <= 1L) {
        return(list(
            m = setNames(rep(0.99, nv), key),
            u = setNames(rep(0.01, nv), key)
        ))
    }

    # m-probabilities: agreement rates among true matched pairs
    m <- numeric(nv)
    names(m) <- key
    for (vi in seq_along(key)) {
        v <- key[vi]
        if (type[[v]] %in% c("numeric", "ordinal")) {
            r <- rng[[v]]
            span <- r[2] - r[1]
            if (span <= 0) span <- 1
            tol_v <- stats::mad(X[[v]], na.rm = TRUE)
            if (!is.finite(tol_v) || tol_v <= 0)
                tol_v <- 0.1 * span
            agree <- vapply(seq_len(n), function(i) {
                xv <- X[[v]][i]
                av <- x_anon[[v]][true_idx[i]]
                if (is.na(xv) || is.na(av)) return(NA)
                abs(xv - av) <= tol_v
            }, logical(1))
            m[vi] <- mean(agree, na.rm = TRUE)
        } else {
            agree <- vapply(seq_len(n), function(i) {
                xv <- X[[v]][i]
                av <- x_anon[[v]][true_idx[i]]
                if (is.na(xv) || is.na(av)) return(NA)
                as.character(xv) == as.character(av)
            }, logical(1))
            m[vi] <- mean(agree, na.rm = TRUE)
        }
        if (is.na(m[vi]) || m[vi] <= 0) m[vi] <- 0.01
        if (m[vi] >= 1) m[vi] <- 0.99
    }

    # u-probabilities: agreement rates among random non-matched pairs
    u <- numeric(nv)
    names(u) <- key
    n_pairs <- min(n_sample, n * (n - 1L))
    sample_i <- sample.int(n, n_pairs, replace = TRUE)
    sample_j <- sample.int(nrow(x_anon), n_pairs, replace = TRUE)
    # Avoid true matches
    is_true <- sample_j == true_idx[sample_i]
    if (any(is_true)) {
        sample_j[is_true] <- ((sample_j[is_true]) %% nrow(x_anon)) + 1L
    }

    for (vi in seq_along(key)) {
        v <- key[vi]
        if (type[[v]] %in% c("numeric", "ordinal")) {
            r <- rng[[v]]
            span <- r[2] - r[1]
            if (span <= 0) span <- 1
            tol_v <- stats::mad(X[[v]], na.rm = TRUE)
            if (!is.finite(tol_v) || tol_v <= 0)
                tol_v <- 0.1 * span
            agree <- vapply(seq_len(n_pairs), function(idx) {
                xv <- X[[v]][sample_i[idx]]
                av <- x_anon[[v]][sample_j[idx]]
                if (is.na(xv) || is.na(av)) return(NA)
                abs(xv - av) <= tol_v
            }, logical(1))
            u[vi] <- mean(agree, na.rm = TRUE)
        } else {
            agree <- vapply(seq_len(n_pairs), function(idx) {
                xv <- X[[v]][sample_i[idx]]
                av <- x_anon[[v]][sample_j[idx]]
                if (is.na(xv) || is.na(av)) return(NA)
                as.character(xv) == as.character(av)
            }, logical(1))
            u[vi] <- mean(agree, na.rm = TRUE)
        }
        if (is.na(u[vi]) || u[vi] <= 0) u[vi] <- 0.01
        if (u[vi] >= 1) u[vi] <- 0.99
    }

    list(m = m, u = u)
}

#' @keywords internal
.fs_log_lr <- function(x_row, anon_block, key, type, rng,
                       m_probs, u_probs, tol = NULL) {
    m <- nrow(anon_block)
    lr <- numeric(m)

    for (vi in seq_along(key)) {
        v <- key[vi]
        xv <- x_row[[v]]
        av <- anon_block[[v]]
        mv <- m_probs[v]
        uv <- u_probs[v]

        if (type[[v]] %in% c("numeric", "ordinal")) {
            if (!is.null(tol) && !is.na(tol[v])) {
                tol_v <- tol[v]
            } else {
                r <- rng[[v]]
                span <- r[2] - r[1]
                if (span <= 0) span <- 1
                tol_v <- stats::mad(c(xv, av), na.rm = TRUE)
                if (!is.finite(tol_v) || tol_v <= 0)
                    tol_v <- 0.1 * span
            }
            gamma <- ifelse(!is.na(av) & !is.na(xv),
                            abs(av - xv) <= tol_v, FALSE)
        } else {
            gamma <- ifelse(!is.na(av) & !is.na(xv),
                            as.character(av) == as.character(xv), FALSE)
        }

        lr <- lr + ifelse(gamma,
                          log(mv / uv),
                          log((1 - mv) / (1 - uv)))
    }

    lr
}

#' @keywords internal
.pram_risk <- function(x_row, anon_block, key, pram_matrix,
                       reverse = FALSE) {
    m <- nrow(anon_block)
    log_probs <- numeric(m)

    for (v in key) {
        xv <- as.character(x_row[[v]])
        av <- as.character(anon_block[[v]])
        tm <- pram_matrix[[v]]
        row_names <- rownames(tm)
        col_names <- colnames(tm)

        if (is.null(row_names) || is.null(col_names)) {
            # Assume levels are 1:nrow if no names
            row_names <- as.character(seq_len(nrow(tm)))
            col_names <- as.character(seq_len(ncol(tm)))
        }

        for (j in seq_len(m)) {
            if (reverse) {
                # reverse: query=anon, search=original
                # P(anon_value | original_value) = tm[original, anon]
                ri <- match(av[j], row_names)
                ci <- match(xv, col_names)
            } else {
                # forward: query=original, search=anon
                # P(anon_value | original_value) = tm[original, anon]
                ri <- match(xv, row_names)
                ci <- match(av[j], col_names)
            }
            if (is.na(ri) || is.na(ci)) {
                log_probs[j] <- -Inf
            } else {
                p <- tm[ri, ci]
                log_probs[j] <- log_probs[j] + if (p > 0) log(p) else -Inf
            }
        }
    }

    exp(log_probs)
}


# ---------- S3 methods ----------

#' Print method for recordLinkageRisk
#'
#' @param x object of class \code{"recordLinkageRisk"}.
#' @param ... ignored.
#' @export
print.recordLinkageRisk <- function(x, ...) {
    if (!is.list(x) || is.null(x$overall) || is.null(x$settings)) {
        cat("<invalid recordLinkageRisk object>\n")
        return(invisible(x))
    }

    s <- x$settings
    o <- x$overall
    meth <- x$method

    dir <- if (is.null(x$direction)) "forward" else x$direction

    cat("Record Linkage Risk\n")
    cat("===================\n\n")

    cat("Method:      ", meth, "\n", sep = "")
    cat("Direction:   ", dir, "\n", sep = "")
    mtch <- if (!is.null(s$matching)) s$matching else "independent"
    if (mtch == "bijective")
        cat("Matching:    bijective (Hungarian algorithm)\n")
    wt_label <- s$risk_weighting
    if (s$risk_weighting == "kernel")
        wt_label <- paste0(wt_label, " (", s$kernel, " kernel)")
    cat("Weighting:   ", wt_label, "\n", sep = "")
    cat("Keys:        ", paste(s$key, collapse = ", "), "\n", sep = "")
    rec_label <- if (dir == "reverse") " (risk per synthetic record)" else ""
    cat("Records:     ", x$n_original, " original, ",
        x$n_synthetic, " synthetic", rec_label, "\n", sep = "")
    cat("Truth:       ", s$truth,
        if (!is.null(s$id)) paste0(" (id: ", s$id, ")") else "", "\n", sep = "")

    if (meth %in% c("deterministic", "predictive")) {
        cat("Strategy:    ", s$strategy, "\n", sep = "")
        if (!is.null(s$block))
            cat("Blocking:    ", paste(s$block, collapse = ", "), "\n", sep = "")
        cat("NA handling: ", s$na_anon, "\n", sep = "")
    } else if (!is.null(s$block)) {
        cat("Blocking:    ", paste(s$block, collapse = ", "), "\n", sep = "")
    }

    if (meth == "mahalanobis") {
        cat("Robust MCD:  ", if (isTRUE(s$robust)) "yes" else "no", "\n", sep = "")
    }
    if (meth == "predictive") {
        cat("Pred. model: ", s$pred_model, "\n", sep = "")
        cat("Fisher SE:   ", if (isTRUE(s$pred_se)) "yes" else "no", "\n", sep = "")
        if (!is.null(x$propensity_info)) {
            cat(sprintf("Mean propensity: %.3f (orig) / %.3f (synth)\n",
                        x$propensity_info$mean_propensity_original,
                        x$propensity_info$mean_propensity_synthetic))
        }
    }

    cat("\nRisk Summary\n")
    cat(sprintf("  Mean risk:       %6.4f\n", o$mean_risk))
    cat(sprintf("  Max risk:        %6.4f\n", o$max_risk))
    thresh <- if (!is.null(x$settings$risk_threshold)) x$settings$risk_threshold else 0.1
    cat(sprintf("  High risk (>%.2g): %d (%.1f%%)\n",
                thresh, o$n_high_risk, 100 * o$pct_high_risk))

    cat("\nPrivacy Assessment: ")
    if (x$privacy_pass) {
        cat("PASS\n")
    } else {
        cat("WARNING\n")
        cat(sprintf("  Mean re-identification risk exceeds %.2g.\n", thresh))
    }

    invisible(x)
}


#' Summary method for recordLinkageRisk
#'
#' @param object object of class \code{"recordLinkageRisk"}.
#' @param ... ignored.
#' @export
summary.recordLinkageRisk <- function(object, ...) {
    risk <- object$per_record$risk

    # Risk quantiles
    risk_quantiles <- stats::quantile(risk,
        probs = c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1))

    # Risk bands (from per_record if available, else recompute)
    if (!is.null(object$per_record$risk_band)) {
        tb <- table(object$per_record$risk_band)
        risk_bands <- c(
            unique_match = unname(tb["unique_match"]),
            very_high    = unname(tb["very_high"]),
            high         = unname(tb["high"]),
            moderate     = unname(tb["moderate"]),
            low          = unname(tb["low"]),
            very_low     = unname(tb["very_low"])
        )
        risk_bands[is.na(risk_bands)] <- 0L
    } else {
        risk_bands <- c(
            unique_match = sum(risk == 1),
            very_high    = sum(risk > 0.5 & risk < 1),
            high         = sum(risk > 0.2 & risk <= 0.5),
            moderate     = sum(risk > 0.1 & risk <= 0.2),
            low          = sum(risk > 0.05 & risk <= 0.1),
            very_low     = sum(risk <= 0.05)
        )
    }

    # Distance-to-true-match statistics
    d_true <- object$per_record$d_true
    d_true_valid <- d_true[!is.na(d_true)]
    d_true_stats <- if (length(d_true_valid) > 0) {
        c(mean = mean(d_true_valid),
          median = stats::median(d_true_valid),
          sd = stats::sd(d_true_valid),
          min = min(d_true_valid),
          max = max(d_true_valid))
    } else {
        NULL
    }

    summ <- list(
        method = object$method,
        direction = if (is.null(object$direction)) "forward" else object$direction,
        matching = if (!is.null(object$settings$matching)) object$settings$matching
                   else "independent",
        risk_weighting = object$settings$risk_weighting,
        kernel = object$settings$kernel,
        bandwidth = object$settings$bandwidth,
        key_vars = object$key_vars,
        n_original = object$n_original,
        n_synthetic = object$n_synthetic,
        mean_risk = object$overall$mean_risk,
        max_risk = object$overall$max_risk,
        risk_quantiles = risk_quantiles,
        risk_bands = risk_bands,
        d_true_stats = d_true_stats,
        pct_true_in_set = object$overall$pct_true_in_set,
        mean_candidate_size = object$overall$mean_candidate_size,
        risk_gini = object$overall$risk_gini,
        var_importance = object$var_importance,
        privacy_pass = object$privacy_pass,
        fs_params = object$fs_params,
        propensity_info = object$propensity_info
    )

    class(summ) <- "summary.recordLinkageRisk"
    summ
}


#' Print method for summary.recordLinkageRisk
#'
#' @param x object of class \code{"summary.recordLinkageRisk"}.
#' @param ... ignored.
#' @export
print.summary.recordLinkageRisk <- function(x, ...) {
    dir <- if (is.null(x$direction)) "forward" else x$direction

    cat("Summary: Record Linkage Risk\n")
    cat("============================\n\n")

    wt_label <- x$risk_weighting
    if (x$risk_weighting == "kernel")
        wt_label <- paste0(wt_label, " (", x$kernel, " kernel)")
    mtch <- if (!is.null(x$matching)) x$matching else "independent"
    cat("Method:", x$method, "| Direction:", dir,
        "| Weighting:", wt_label, "\n")
    if (mtch == "bijective")
        cat("Matching: bijective (Hungarian algorithm)\n")
    cat("Key variables:", paste(x$key_vars, collapse = ", "), "\n")
    cat("Records:", x$n_original, "original,", x$n_synthetic, "synthetic\n\n")

    cat("Risk Distribution:\n")
    rq <- x$risk_quantiles
    nms <- c("Min", "P5", "Q1", "Median", "Q3", "P95", "Max")
    for (i in seq_along(rq)) {
        cat(sprintf("  %-8s %7.4f\n", paste0(nms[i], ":"), rq[i]))
    }

    cat("\nRisk Bands:\n")
    bn <- c("Unique match (=1)", "Very high (>0.5)",
            "High (0.2-0.5)", "Moderate (0.1-0.2)",
            "Low (0.05-0.1)", "Very low (<0.05)")
    n_denom <- if (dir == "reverse") x$n_synthetic else x$n_original
    for (i in seq_along(x$risk_bands)) {
        pct <- 100 * x$risk_bands[i] / n_denom
        cat(sprintf("  %-22s %5d (%5.1f%%)\n",
                    bn[i], x$risk_bands[i], pct))
    }

    if (!is.null(x$d_true_stats)) {
        d_label <- if (x$method == "probabilistic") "Log-LR for True Match:"
                   else "Distance to True Match:"
        cat("\n", d_label, "\n", sep = "")
        cat(sprintf("  Mean:    %7.4f\n", x$d_true_stats["mean"]))
        cat(sprintf("  Median:  %7.4f\n", x$d_true_stats["median"]))
        cat(sprintf("  SD:      %7.4f\n", x$d_true_stats["sd"]))
    }

    if (!is.null(x$fs_params)) {
        cat("\nFellegi-Sunter Parameters:\n")
        cat("  Variable    m_prob   u_prob\n")
        for (v in names(x$fs_params$m_probs)) {
            cat(sprintf("  %-12s %5.3f    %5.3f\n",
                        v, x$fs_params$m_probs[v], x$fs_params$u_probs[v]))
        }
    }

    if (!is.null(x$propensity_info)) {
        cat("\nPropensity Model:\n")
        cat("  Model:", x$propensity_info$pred_model, "\n")
        cat(sprintf("  Mean propensity: %.4f (orig) / %.4f (synth)\n",
                    x$propensity_info$mean_propensity_original,
                    x$propensity_info$mean_propensity_synthetic))
    }

    if (!is.null(x$var_importance) && !all(is.na(x$var_importance))) {
        vi_label <- switch(x$method,
            deterministic = "Variable Importance (mean weighted distance):",
            probabilistic = "Variable Importance (log-LR on agreement):",
            predictive    = "Variable Importance (model coefficients):",
            pram          = "Variable Importance (perturbation strength):",
            mahalanobis   = "Variable Importance (precision matrix proportion):",
            "Variable Importance:")
        cat("\n", vi_label, "\n", sep = "")
        for (v in names(x$var_importance)) {
            cat(sprintf("  %-12s %7.4f\n", v, x$var_importance[v]))
        }
    }

    cat(sprintf("\nMean candidate size: %.1f\n", x$mean_candidate_size))
    cat(sprintf("True match in candidate set: %.1f%%\n",
                100 * x$pct_true_in_set))
    if (!is.null(x$risk_gini)) {
        cat(sprintf("Risk concentration (Gini): %.4f\n", x$risk_gini))
    }
    cat("Privacy:", ifelse(x$privacy_pass, "PASS", "WARNING"), "\n")

    invisible(x)
}


#' Plot method for recordLinkageRisk
#'
#' @param x object of class \code{"recordLinkageRisk"}.
#' @param y not used.
#' @param ... additional arguments passed to plotting functions.
#' @param which integer, which plot(s) to show:
#'   1 = Risk distribution histogram,
#'   2 = Distance/score vs risk scatterplot,
#'   3 = Per-variable importance (all methods),
#'   4 = Propensity score distributions (predictive method only),
#'   5 = Risk band barplot,
#'   6 = True-match rank distribution,
#'   7 = Risk by group barplot (requires \code{group}),
#'   8 = Lorenz curve of risk concentration
#' @param group optional grouping vector or column name (with \code{data})
#'   for \code{which = 7}.
#' @param data optional data frame for extracting \code{group} column.
#' @importFrom graphics hist abline legend par plot plot.new points text barplot lines polygon segments axis box
#' @export
plot.recordLinkageRisk <- function(x, y = NULL, ..., which = 1,
                                    group = NULL, data = NULL) {
    n_types <- 8L
    show <- rep(FALSE, n_types)
    show[which] <- TRUE

    dir <- if (is.null(x$direction)) "forward" else x$direction
    dir_suffix <- if (dir == "reverse") " (reverse)" else ""
    thresh <- if (!is.null(x$settings$risk_threshold))
        x$settings$risk_threshold else 0.1

    n_plots <- sum(show)
    if (n_plots > 1) {
        ncol <- min(n_plots, 2)
        nrow <- ceiling(n_plots / ncol)
        op <- par(mfrow = c(nrow, ncol))
        on.exit(par(op))
    }

    if (show[1]) {
        # Risk distribution histogram
        risk <- x$per_record$risk
        hist(risk,
             main = paste0("Record Linkage Risk Distribution", dir_suffix),
             xlab = "Re-identification Risk",
             ylab = "Number of Records",
             col = "steelblue", border = "white", ...)
        abline(v = thresh, col = "red", lty = 2, lwd = 2)
        abline(v = mean(risk), col = "darkblue", lty = 3, lwd = 1.5)
        legend("topright",
               legend = c(paste0("threshold = ", thresh),
                           paste0("mean = ", round(mean(risk), 3))),
               col = c("red", "darkblue"),
               lty = c(2, 3), lwd = c(2, 1.5), cex = 0.8)
    }

    if (show[2]) {
        # Distance vs risk scatter
        d <- x$per_record$d_true
        risk <- x$per_record$risk
        valid <- !is.na(d)
        if (sum(valid) > 0) {
            xlab2 <- if (x$method == "pram") "Transition Probability"
                     else if (x$method == "probabilistic") "Log-LR (True Match)"
                     else "Distance to True Match"
            plot(d[valid], risk[valid],
                 main = paste0(xlab2, " vs Risk", dir_suffix),
                 xlab = xlab2,
                 ylab = "Re-identification Risk",
                 pch = 16, col = adjustcolor("steelblue", 0.5),
                 ...)
            abline(h = thresh, col = "red", lty = 2, lwd = 2)
        } else {
            plot.new()
            text(0.5, 0.5, "No distance data available",
                 cex = 1.2)
        }
    }

    if (show[3]) {
        # Per-variable importance / discriminative power
        if (!is.null(x$fs_params)) {
            # Probabilistic: grouped agree/disagree barplot
            m <- x$fs_params$m_probs
            u <- x$fs_params$u_probs
            log_lr_agree <- log(m / u)
            log_lr_disagree <- log((1 - m) / (1 - u))

            vars <- names(m)
            mat <- rbind(agree = log_lr_agree, disagree = log_lr_disagree)
            colnames(mat) <- vars

            barplot(mat, beside = TRUE,
                    main = paste0("Per-Variable Log-Likelihood Ratios",
                                  dir_suffix),
                    ylab = "log(LR)",
                    col = c("steelblue", "coral"),
                    las = 2, ...)
            abline(h = 0, col = "gray40", lty = 1)
            legend("topright",
                   legend = c("Agree", "Disagree"),
                   fill = c("steelblue", "coral"),
                   cex = 0.8)
        } else if (!is.null(x$var_importance) &&
                   !all(is.na(x$var_importance))) {
            vi <- x$var_importance
            vi_label <- switch(x$method,
                deterministic = "Mean Weighted Distance",
                predictive    = "Absolute Model Coefficient",
                pram          = "Perturbation Strength (1 - diag)",
                rf            = "RF Impurity Importance",
                mahalanobis   = "Precision Matrix Proportion",
                "Importance")
            barplot(vi[order(vi, decreasing = TRUE)],
                    main = paste0("Variable Importance", dir_suffix),
                    ylab = vi_label,
                    col = "steelblue", las = 2, ...)
        } else {
            plot.new()
            text(0.5, 0.5,
                 "No variable importance available",
                 cex = 1.2)
        }
    }

    if (show[4]) {
        if (!is.null(x$propensity_info)) {
            p_orig <- x$propensity_info$propensity_original
            p_synth <- x$propensity_info$propensity_synthetic
            d_orig <- stats::density(p_orig, from = 0, to = 1)
            d_synth <- stats::density(p_synth, from = 0, to = 1)
            ylim <- range(c(d_orig$y, d_synth$y))
            plot(d_orig, main = paste0("Propensity Score Distributions", dir_suffix),
                 xlab = "Propensity Score P(original | QIs)",
                 ylab = "Density", col = "steelblue", lwd = 2,
                 ylim = ylim)
            lines(d_synth, col = "coral", lwd = 2)
            legend("topright",
                   legend = c("Original", "Synthetic"),
                   col = c("steelblue", "coral"),
                   lwd = 2, cex = 0.8)
        } else {
            plot.new()
            text(0.5, 0.5,
                 "Plot 4 requires method='predictive'",
                 cex = 1.2)
        }
    }

    if (show[5]) {
        # Risk band barplot
        if (!is.null(x$per_record$risk_band)) {
            tb <- table(x$per_record$risk_band)
            # Reverse so unique_match is on top
            tb_rev <- rev(tb)
            cols <- rev(c("#2166ac", "#67a9cf", "#d1e5f0",
                          "#fddbc7", "#ef8a62", "#b2182b"))
            bp <- barplot(tb_rev,
                    main = paste0("Risk Band Distribution", dir_suffix),
                    ylab = "Number of Records",
                    col = cols, las = 2,
                    names.arg = gsub("_", " ", names(tb_rev)),
                    ...)
            # Add count labels on bars
            nonzero <- tb_rev > 0
            if (any(nonzero)) {
                text(bp[nonzero], tb_rev[nonzero],
                     labels = tb_rev[nonzero],
                     pos = 3, cex = 0.8)
            }
        } else {
            plot.new()
            text(0.5, 0.5, "No risk_band data available", cex = 1.2)
        }
    }

    if (show[6]) {
        # True-match rank distribution
        dr <- x$per_record$d_rank
        valid <- !is.na(dr)
        if (sum(valid) > 0) {
            dr_valid <- dr[valid]
            max_rank <- max(dr_valid)
            brks <- seq(0.5, max_rank + 0.5, by = 1)
            hist(dr_valid, breaks = brks,
                 main = paste0("True-Match Rank Distribution", dir_suffix),
                 xlab = "Rank of True Match Among Candidates",
                 ylab = "Number of Records",
                 col = "steelblue", border = "white", ...)
            abline(v = 1, col = "red", lty = 2, lwd = 2)
            pct_rank1 <- 100 * mean(dr_valid == 1L)
            legend("topright",
                   legend = paste0("Rank 1: ", round(pct_rank1, 1), "%"),
                   col = "red", lty = 2, lwd = 2, cex = 0.8)
        } else {
            plot.new()
            text(0.5, 0.5, "No rank data available\n(no true matches found)",
                 cex = 1.2)
        }
    }

    if (show[7]) {
        # Risk by group barplot
        if (!is.null(group)) {
            rg <- risk_by_group(x, group = group, data = data)
            # Find the group column name (last column that isn't a stat)
            gcol <- setdiff(names(rg),
                            c("mean_risk", "max_risk", "n", "n_high",
                              "pct_high"))
            labels <- rg[[gcol[1]]]
            bp <- barplot(rg$mean_risk,
                    names.arg = labels,
                    main = paste0("Mean Risk by Group", dir_suffix),
                    ylab = "Mean Re-identification Risk",
                    col = "steelblue", las = 2, ...)
            abline(h = thresh, col = "red", lty = 2, lwd = 2)
            # Add error bars showing max_risk
            segments(bp, rg$mean_risk, bp, rg$max_risk,
                     col = "gray30", lwd = 1.5)
            points(bp, rg$max_risk, pch = 4, col = "gray30", cex = 0.8)
            legend("topright",
                   legend = c(paste0("threshold = ", thresh), "max risk"),
                   col = c("red", "gray30"),
                   lty = c(2, NA), pch = c(NA, 4),
                   lwd = c(2, NA), cex = 0.8)
        } else {
            plot.new()
            text(0.5, 0.5,
                 "Plot 7 requires 'group' argument",
                 cex = 1.2)
        }
    }

    if (show[8]) {
        # Lorenz curve of risk concentration
        risk <- sort(x$per_record$risk)
        n <- length(risk)
        cum_risk <- cumsum(risk) / sum(risk)
        cum_pop <- seq_len(n) / n
        gini_val <- x$overall$risk_gini

        plot(c(0, cum_pop), c(0, cum_risk),
             type = "l", lwd = 2, col = "steelblue",
             main = paste0("Lorenz Curve of Risk", dir_suffix),
             xlab = "Cumulative Share of Records",
             ylab = "Cumulative Share of Risk",
             xlim = c(0, 1), ylim = c(0, 1), ...)
        # Diagonal (perfect equality)
        lines(c(0, 1), c(0, 1), lty = 2, col = "gray50", lwd = 1.5)
        # Shade area between diagonal and Lorenz curve
        polygon(c(0, cum_pop, 1, 0),
                c(0, cum_risk, 1, 0),
                col = adjustcolor("steelblue", 0.15), border = NA)
        if (!is.null(gini_val)) {
            legend("topleft",
                   legend = c("Lorenz curve",
                              "Perfect equality",
                              paste0("Gini = ", round(gini_val, 3))),
                   col = c("steelblue", "gray50", NA),
                   lty = c(1, 2, NA), lwd = c(2, 1.5, NA),
                   cex = 0.8)
        }
    }
}
