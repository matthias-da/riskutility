#' Record Linkage Risk
#'
#' Measures targeted re-identification risk by linking each original record to
#' the most similar record(s) in a perturbed dataset using quasi-identifiers.
#' Supports distance-based (Gower distance) and probabilistic
#' (Fellegi-Sunter) linkage methods.
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
#' For the distance-based method, record linkage is based on a weighted Gower-type
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
#' non-match (u) probabilities for each variable using exact agreement,
#' then computes log-likelihood ratios to score candidate pairs.
#' Per-record risk is the Skinner (2008) posterior identification probability:
#' \deqn{P(\text{match} \mid \gamma) = \frac{\Lambda \cdot p}{\Lambda \cdot p + (1 - p)}}
#' where \eqn{\Lambda = m(\gamma)/u(\gamma)} is the likelihood ratio and
#' \eqn{p = 1/|B|} is the closed-world prior (one true match among the
#' \eqn{|B|} candidates in blocking group \eqn{B}; without blocking
#' \eqn{|B| = n}). This is directly interpretable as the probability that
#' the attacker's link is correct, given the observed agreement pattern.
#'
#' User-supplied \code{m_probs} and \code{u_probs} (via \code{rl_control()})
#' override the supervised estimation. NA pairs always contribute 0 to the
#' log-likelihood ratio (variable dropped). Direction, \code{na_anon}, and
#' \code{risk_weighting} do not apply to this method.
#'
#' \strong{Matching mode:}
#' Use \code{matching = "independent"} (default) for classical per-record risk
#' (DBRL). Bijective matching is not available for the probabilistic method;
#' use \code{method = "distance-based"} if one-to-one assignment is required.
#'
#' @section Softmax risk weighting:
#' When \code{risk_weighting = "softmax"} (set via \code{rl_control()}), closer
#' candidates receive higher attribution probability via:
#' \deqn{w_j = \frac{\exp(-\kappa \cdot d_j)}{\sum_k \exp(-\kappa \cdot d_k)}}
#' The temperature \code{kappa} is auto-calibrated from the distance range
#' if not supplied. This replaces the uniform \eqn{1/|candidate\_set|} risk.
#'
#' @section Harmonic rank weighting:
#' When \code{risk_weighting = "harmonic"} (set via \code{rl_control()}),
#' candidates are weighted by the reciprocal of their distance rank:
#' \deqn{w_j = \frac{1/r_j}{\sum_k 1/r_k}}
#' where \eqn{r_j} is the rank of candidate \eqn{j} by ascending distance
#' (rank 1 = closest). Ties receive the minimum rank. Unlike softmax, harmonic
#' weighting requires no tuning parameter and depends only on the ordering of
#' candidates, not their absolute distances.
#'
#' @section Direction:
#' Only applies to \code{method = "distance-based"}.
#' With \code{direction = "forward"} (default) the function loops over original
#' records and finds nearest matches in the anonymized data. This corresponds
#' to a prosecutor or nosy-neighbour attacker who holds external knowledge about
#' a target individual and searches the released file for that individual.
#'
#' With \code{direction = "reverse"} the loop runs over anonymized records and
#' finds nearest matches in the original data. This is the direction used by the
#' classical DBRL literature (Domingo-Ferrer and Torra, 2003; Herranz et al.,
#' 2016) and by \pkg{sdcMicro}'s \code{dRiskRMD}, and corresponds to a
#' journalist-type attacker who holds the released file and a reference database
#' and tries to identify each released record in that reference.
#' The two directions can produce different risk estimates because the
#' nearest-neighbour relation is not symmetric.
#'
#' @param X data.frame or \code{\link{synth_pair}} object. Original microdata.
#' @param x_anon data.frame. Perturbed/anonymized microdata.
#' @param key character. Names of quasi-identifier variables used for linkage.
#' @param method character. Linkage method: \code{"distance-based"} (default,
#'   weighted Gower distance) or \code{"probabilistic"} (Fellegi-Sunter
#'   log-likelihood ratios).
#' @param direction character. Direction of the linkage attack:
#'   \code{"forward"} (default) loops over original records and searches in the
#'   anonymized data (prosecutor / nosy-neighbour scenario);
#'   \code{"reverse"} loops over anonymized records and searches in the
#'   original data (journalist scenario).
#'   Only applies to \code{method = "distance-based"}; probabilistic linkage
#'   always runs forward.
#' @param truth character. How to define the true match for scoring:
#'   one of \code{"row"} (default) or \code{"id"}.
#'   \code{"row"} requires equal row counts and assumes row \eqn{i}
#'   in \code{X} corresponds to row \eqn{i} in \code{x_anon}.
#'   When datasets have unequal sizes, use \code{truth = "id"} with
#'   a shared identifier column.
#' @param id character or NULL. If \code{truth="id"}, name of an identifier column
#'   present in both \code{X} and \code{x_anon} that uniquely defines the true match.
#' @param block character or NULL. Optional subset of \code{key} used for exact blocking.
#'   Distances are computed only within blocks defined by these variables.
#' @param return_matches logical. If TRUE, returns candidate indices per record (may be memory-heavy).
#' @param matching character. Matching mode: \code{"independent"} (default) scores
#'   each record independently (classical DBRL), \code{"bijective"} enforces
#'   one-to-one assignment via the Hungarian algorithm (GDBRL).
#'   Bijective matching requires the \pkg{clue} package and is only supported
#'   for \code{method = "distance-based"}.
#'   See Herranz, Nin, Rodriguez & Tassa (2016).
#' @param risk_threshold numeric. Threshold for classifying records as "high risk"
#'   and for the \code{privacy_pass} flag (default 0.1). Records with risk above
#'   this threshold are counted in \code{n_high_risk} and \code{pct_high_risk}.
#'   The \code{privacy_pass} flag is TRUE when \code{mean_risk <= risk_threshold}.
#' @param control an \code{rl_control} object from \code{\link{rl_control}()},
#'   carrying method-specific tuning parameters. Use
#'   \code{rl_control(weights, type, risk_weighting, kappa, strategy, k)} for
#'   distance-based tuning and \code{rl_control(m_probs, u_probs)} for
#'   probabilistic tuning. Irrelevant parameters for the chosen method produce
#'   a warning and are ignored.
#' @param ... additional arguments passed to methods.
#' @author Matthias Templ and Roman Müller
#'
#' @return An object of class \code{"recordLinkageRisk"}: a list with components
#'   \describe{
#'     \item{per_record}{data.frame with \code{n_query} rows. Columns common
#'       to both methods: \code{risk} (numeric, re-identification risk),
#'       \code{cand_n} (integer; for distance-based, size of the attacker's
#'         guess set; for probabilistic, size of the blocking group used in
#'         the Skinner prior),
#'       \code{true_in_set} (logical),
#'       \code{risk_band} (ordered factor with levels \code{"very_low"},
#'         \code{"low"}, \code{"moderate"}, \code{"high"}, \code{"very_high"},
#'         \code{"unique_match"}).
#'       Distance-based only: \code{d_true} (Gower distance to true match),
#'       \code{d_min} (minimum distance to any candidate),
#'       \code{d_rank} (rank of true match by ascending distance).
#'       With \code{matching = "bijective"}, an additional
#'       \code{bijective_assigned} column gives the Hungarian-algorithm
#'       assignment; \code{d_true}, \code{d_min}, \code{d_rank} retain
#'       independent-scoring diagnostics.
#'       Probabilistic only: \code{lr_true} (log-likelihood ratio of the true
#'       match), \code{lr_rank} (rank of true match by descending log-LR).}
#'     \item{overall}{list with aggregate statistics including \code{risk_gini}
#'       (Gini coefficient of risk concentration).}
#'     \item{var_importance}{named numeric vector of per-variable mean weighted
#'       distance (distance-based only; absent for probabilistic).}
#'     \item{direction}{character, \code{"forward"} or \code{"reverse"}
#'       (distance-based only; always \code{"forward"} for probabilistic).}
#'     \item{matches}{(optional) list of integer vectors with candidate indices.}
#'   }
#'   For distance-based: when \code{direction = "forward"}, \code{per_record}
#'   has one row per original record; when \code{"reverse"}, one row per
#'   anonymized record. Probabilistic always runs forward.
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
#'   region = factor(sample(paste0("R",1:5), n, TRUE))
#' )
#'
#' # Simple perturbation: swap age within region
#' x_anon <- x
#' for (r in levels(x$region)) {
#'   idx <- which(x$region == r)
#'   x_anon$age[idx] <- sample(x$age[idx])
#' }
#'
#' # Distance-based nearest-neighbor (default)
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
#' # Softmax distance-weighted risk (via rl_control)
#' res2 <- recordLinkage(x, x_anon, key = c("age","sex","region"),
#'                       control = rl_control(risk_weighting = "softmax"))
#' print(res2)
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
#' Herranz, J., Nin, J., Rodriguez, P., & Tassa, T. (2015). Revisiting
#' distance-based record linkage for privacy-preserving release of statistical
#' datasets. Data & Knowledge Engineering, 100, 78-93.
#'
#' @seealso \code{\link{individual_risk}}, \code{\link{dcr}},
#'   \code{\link{nndr}}
#' @family privacy-models
#' @importFrom stats quantile
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

#' Control Parameters for Record Linkage
#'
#' Creates a control object that carries method-specific tuning parameters for
#' \code{\link{recordLinkage}()}. Distance-based and probabilistic parameters
#' are grouped together so that \code{recordLinkage()} can warn when parameters
#' irrelevant to the chosen method are supplied.
#'
#' @section Distance-based parameters:
#' \describe{
#'   \item{\code{weights}}{Named nonnegative numeric vector of variable weights.
#'     If \code{NULL}, equal weights are used.}
#'   \item{\code{type}}{Named character vector overriding per-variable types:
#'     values in \code{c("numeric","ordinal","nominal")}. If \code{NULL},
#'     types are inferred from the data.}
#'   \item{\code{risk_weighting}}{How to weight candidates within the guess set:
#'     \code{"uniform"} (default), \code{"softmax"} (exponential
#'     distance-weighting), or \code{"harmonic"} (rank-based reciprocal
#'     weighting, no tuning parameter).}
#'   \item{\code{kappa}}{Temperature for softmax weighting. If \code{NULL}
#'     (default), auto-calibrated as \code{2 / range(distances)}.}
#'   \item{\code{strategy}}{Adversary candidate-set strategy: \code{"nearest"}
#'     (default) or \code{"topk"}.}
#'   \item{\code{k}}{Candidate set size for \code{strategy = "topk"}.}
#' }
#'
#' @section Probabilistic parameters:
#' \describe{
#'   \item{\code{m_probs}}{Named numeric vector of m-probabilities (probability
#'     of agreement given true match). If \code{NULL}, estimated from data.}
#'   \item{\code{u_probs}}{Named numeric vector of u-probabilities (probability
#'     of agreement given non-match). If \code{NULL}, estimated from data by
#'     enumerating all distinct non-matching pairs when \eqn{n(n-1) \le 500},
#'     or sampling 500 pairs without replacement otherwise.}
#'   \item{\code{breaks}}{Named list of break points for discretizing numeric
#'     variables before probabilistic scoring. Each element is a sorted numeric
#'     vector of cut points (excluding \code{-Inf}/\code{Inf}); the variable is
#'     replaced by the integer bin index before computing agreement. Only used
#'     when \code{method = "probabilistic"}. Example:
#'     \code{breaks = list(age = c(30, 60), income = c(2000, 5000))}.}
#' }
#'
#' @section Distance-based-only parameters:
#' \describe{
#'   \item{\code{na_anon}}{How a missing quasi-identifier value in the anonymized
#'     data is handled during distance computation: \code{"ignore"} (default,
#'     variable dropped from that pairwise comparison), \code{"match"} (NA
#'     counts as agreement, distance contribution 0), or \code{"mismatch"} (NA
#'     counts as full disagreement, distance contribution 1). Ignored for
#'     \code{method = "probabilistic"}, which always drops NA pairs.}
#' }
#'
#' @param weights named numeric vector or NULL.
#' @param type named character vector or NULL.
#' @param risk_weighting character, one of \code{"uniform"} (default),
#'   \code{"softmax"}, \code{"harmonic"}.
#' @param kappa numeric or NULL.
#' @param strategy character, one of \code{"nearest"} (default), \code{"topk"}.
#' @param k integer or NULL.
#' @param m_probs named numeric vector or NULL.
#' @param u_probs named numeric vector or NULL.
#' @param breaks named list of numeric break vectors or NULL.
#' @param na_anon character, one of \code{"ignore"} (default), \code{"match"},
#'   \code{"mismatch"}.
#'
#' @return An object of class \code{"rl_control"}.
#' @seealso \code{\link{recordLinkage}}
#' @export
rl_control <- function(weights = NULL,
                       type = NULL,
                       risk_weighting = c("uniform", "softmax", "harmonic"),
                       kappa = NULL,
                       strategy = c("nearest", "topk"),
                       k = NULL,
                       m_probs = NULL,
                       u_probs = NULL,
                       breaks = NULL,
                       na_anon = c("ignore", "match", "mismatch")) {
    risk_weighting <- match.arg(risk_weighting)
    strategy       <- match.arg(strategy)
    na_anon        <- match.arg(na_anon)
    structure(
        list(weights = weights, type = type,
             risk_weighting = risk_weighting, kappa = kappa,
             strategy = strategy, k = k,
             m_probs = m_probs, u_probs = u_probs,
             breaks = breaks, na_anon = na_anon),
        class = "rl_control"
    )
}


#' @rdname recordLinkage
#' @export
recordLinkage.default <- function(X,
                                  x_anon,
                                  key,
                                  method = c("distance-based", "probabilistic"),
                                  direction = c("forward", "reverse"),
                                  truth = c("row", "id"),
                                  id = NULL,
                                  block = NULL,
                                  return_matches = FALSE,
                                  matching = c("independent", "bijective"),
                                  risk_threshold = 0.1,
                                  control = rl_control(),
                                  ...) {

    stopifnot(!missing(x_anon))
    stopifnot(!missing(key))

    method    <- match.arg(method)
    direction <- match.arg(direction)
    truth     <- match.arg(truth)
    matching  <- match.arg(matching)

    if (!inherits(control, "rl_control"))
        stop("'control' must be an object created by rl_control().", call. = FALSE)

    # Extract method-specific parameters from control ----
    weights        <- control$weights
    type           <- control$type
    risk_weighting <- control$risk_weighting
    kappa          <- control$kappa
    strategy       <- control$strategy
    k              <- control$k
    m_probs        <- control$m_probs
    u_probs        <- control$u_probs
    breaks         <- control$breaks
    na_anon        <- control$na_anon

    # Validate method/matching combinations ----
    if (method == "probabilistic") {
        if (matching == "bijective")
            stop("'matching = \"bijective\"' is not supported for ",
                 "method = \"probabilistic\": it would discard the Skinner ",
                 "posterior. Use matching = \"independent\".", call. = FALSE)
        if (direction != "forward") {
            warning("'direction' is only relevant for method = \"distance-based\". ",
                    "Probabilistic linkage always runs forward (original as query).",
                    call. = FALSE)
            direction <- "forward"
        }
        if (na_anon != "ignore")
            warning("'na_anon' in rl_control() is only used for ",
                    "method = \"distance-based\". ",
                    "The probabilistic method always ignores NA pairs.",
                    call. = FALSE)
        if (risk_weighting != "uniform")
            warning("'risk_weighting' is ignored for method = \"probabilistic\", ",
                    "which always reports the Skinner (2008) posterior.",
                    call. = FALSE)
        if (!is.null(weights) || !is.null(type) || !is.null(kappa) ||
            strategy != "nearest" || !is.null(k))
            warning("Distance-based control parameters (weights, type, kappa, ",
                    "strategy, k) are ignored for method = \"probabilistic\".",
                    call. = FALSE)
        if (!is.null(breaks)) {
            if (!is.list(breaks) || is.null(names(breaks)))
                stop("'breaks' must be a named list.", call. = FALSE)
            bad <- names(breaks)[!names(breaks) %in% key]
            if (length(bad))
                stop("'breaks' contains names not in 'key': ",
                     paste(bad, collapse = ", "), call. = FALSE)
            for (nm in names(breaks)) {
                b <- breaks[[nm]]
                if (!is.numeric(b) || length(b) < 1L || anyNA(b))
                    stop("breaks[['", nm, "']] must be a non-empty numeric vector ",
                         "without NAs.", call. = FALSE)
                breaks[[nm]] <- sort(b)
            }
        }
    }
    if (method == "distance-based" && (!is.null(m_probs) || !is.null(u_probs)))
        warning("'m_probs' and 'u_probs' are ignored for ",
                "method = \"distance-based\".", call. = FALSE)
    if (method == "distance-based" && !is.null(breaks))
        warning("'breaks' is only used for method = \"probabilistic\" and is ",
                "ignored for method = \"distance-based\".", call. = FALSE)

    if (matching == "bijective") {
        if (!requireNamespace("clue", quietly = TRUE))
            stop("Package 'clue' is required for bijective matching. ",
                 "Install it with install.packages('clue').", call. = FALSE)
        if (strategy != "nearest")
            message("Note: bijective matching uses the full candidate set; ",
                    "'strategy' applies only to independent-scoring diagnostics.")
        if (risk_weighting != "uniform")
            message("Note: bijective matching produces binary risk (0/1); ",
                    "'risk_weighting' has no effect.")
    }

    stopifnot(is.data.frame(X), is.data.frame(x_anon))
    stopifnot(is.character(key), length(key) >= 1L)

    if (!all(key %in% names(X))) stop("Some 'key' not in 'X'.")
    if (!all(key %in% names(x_anon))) stop("Some 'key' not in 'x_anon'.")

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

    # distance-based setup: types, weights, numeric scaling ranges ----
    type <- weights <- wsum <- rng <- NULL
    if (method == "distance-based") {
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

        rng <- lapply(key, function(v) {
            if (type[[v]] %in% c("numeric", "ordinal")) {
                allv <- c(X[[v]], x_anon[[v]])
                if (is.factor(allv)) allv <- as.integer(allv)
                r <- range(allv, na.rm = TRUE)
                if (!is.finite(r[1]) || !is.finite(r[2]) || r[1] == r[2])
                    r <- c(0, 1)
                r
            } else NULL
        })
        names(rng) <- key
    }

    # Fellegi-Sunter: estimate m/u if probabilistic ----
    fs_params <- NULL
    if (method == "probabilistic") {
        if (is.null(m_probs) || is.null(u_probs)) {
            fs_est <- .fs_estimate(X, x_anon, key, fs_true_idx, breaks = breaks)
            if (is.null(m_probs)) m_probs <- fs_est$m
            if (is.null(u_probs)) u_probs <- fs_est$u
        }
        # Validate and clamp m/u to the open interval (0,1) so the
        # log-likelihood ratios stay finite. Estimated values are already
        # clamped; this guards user-supplied 0/1 (which would give +/-Inf/NaN).
        .chk_prob <- function(p, nm) {
            if (is.null(names(p)) || !all(key %in% names(p)))
                stop("'", nm, "' must be a named vector containing all keys.",
                     call. = FALSE)
            p <- p[key]
            if (any(!is.finite(p)))
                stop("'", nm, "' must be finite.", call. = FALSE)
            bad <- p <= 0 | p >= 1
            if (any(bad)) {
                warning("'", nm, "' values outside (0,1) clamped to ",
                        "[0.01, 0.99]: ", paste(key[bad], collapse = ", "),
                        call. = FALSE)
                p <- pmin(pmax(p, 0.01), 0.99)
            }
            p
        }
        m_probs <- .chk_prob(m_probs, "m_probs")
        u_probs <- .chk_prob(u_probs, "u_probs")
        fs_params <- list(m_probs = m_probs, u_probs = u_probs)
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
    matches <- if (isTRUE(return_matches)) vector("list", n_query) else NULL
    if (method == "distance-based") {
        d_true <- rep(NA_real_, n_query)
        d_min  <- rep(NA_real_, n_query)
        d_rank <- rep(NA_integer_, n_query)
        score_cache <- if (matching == "bijective") vector("list", n_query) else NULL
        var_dist_acc <- setNames(numeric(length(key)), key)
    } else {
        lr_true <- rep(NA_real_, n_query)
        lr_rank <- rep(NA_integer_, n_query)
    }

    # main loop ----
    for (i in seq_len(n_query)) {
        cand <- split_search[[blk_query[i]]]
        if (is.null(cand) || length(cand) == 0L) {
            risk[i] <- 0
            cand_n[i] <- 0L
            true_in_set[i] <- FALSE
            next
        }

        if (method == "probabilistic") {
            # Fellegi-Sunter: compute log-likelihood ratios (exact agreement, NAs ignored)
            lr <- .fs_log_lr(
                x_row = query_data[i, key, drop = FALSE],
                anon_block = search_data[cand, key, drop = FALSE],
                key = key, m_probs = m_probs, u_probs = u_probs,
                breaks = breaks
            )

            tpos <- true_idx[i]
            j_in <- match(tpos, cand)
            if (!is.na(j_in)) {
                lr_true[i] <- lr[j_in]
                # Rank by descending LR (rank 1 = best match)
                lr_rank[i] <- as.integer(rank(-lr, ties.method = "min")[j_in])
                # FIXME: replace Skinner (2008) pairwise posterior with the
                # normalized multinomial posterior LR_true / sum(LR_all).
                # Skinner computes P(match | gamma_true) using only the true
                # match's LR and a flat prior 1/|cand|; it ignores how well the
                # true match's LR dominates competitors, so it systematically
                # underestimates risk when the LR discriminates (e.g. two male
                # candidates with similar age give ~25% instead of ~50%).
                # The normalized posterior risk = exp(lr[j_in]) / sum(exp(lr))
                # correctly concentrates mass on the top-ranked candidates and
                # reduces to 1/N when all LRs are equal.
                p_prior  <- 1 / length(cand)
                lr_ratio <- exp(lr[j_in])
                risk[i]  <- lr_ratio * p_prior / (lr_ratio * p_prior + (1 - p_prior))
            }
            # else: risk[i] stays 0 (true match not in search block)

            cand_n[i] <- length(cand)
            true_in_set[i] <- !is.na(j_in)

            if (isTRUE(return_matches)) matches[[i]] <- cand
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
            k = k
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
        } else if (risk_weighting == "harmonic") {
            # Harmonic rank weighting: w_j = (1/r_j) / sum_k(1/r_k)
            guess_local_idx <- match(guess, cand)
            w <- .harmonic_risk(di[guess_local_idx])
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

    # risk_band ----
    risk_band <- cut(risk,
        breaks = c(-Inf, 0.05, 0.1, 0.2, 0.5, 1 - .Machine$double.eps, Inf),
        labels = c("very_low", "low", "moderate", "high", "very_high",
                    "unique_match"),
        right = TRUE)

    # var_importance: distance-based only ----
    var_importance <- if (method == "distance-based") var_dist_acc / n_query else NULL

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

    if (method == "distance-based") {
        per_rec <- data.frame(
            risk = risk, cand_n = cand_n, true_in_set = true_in_set,
            d_true = d_true, d_min = d_min, d_rank = d_rank,
            risk_band = risk_band
        )
        if (!is.null(bijective_assigned))
            per_rec$bijective_assigned <- bijective_assigned
    } else {
        per_rec <- data.frame(
            risk = risk, cand_n = cand_n, true_in_set = true_in_set,
            lr_true = lr_true, lr_rank = lr_rank,
            risk_band = risk_band
        )
    }

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
            block = block,
            direction = direction,
            matching = matching,
            risk_threshold = risk_threshold,
            control = control
        )
    )
    if (isTRUE(return_matches)) out$matches <- matches
    if (!is.null(fs_params))      out$fs_params <- fs_params
    if (!is.null(var_importance)) out$var_importance <- var_importance

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
        risk_band = as.character(pr$risk_band),
        n_candidates = pr$cand_n,
        true_in_set = pr$true_in_set,
        candidate_ids = candidates
    )
    if (x$method == "distance-based") {
        out$d_true <- pr$d_true
        out$d_min  <- pr$d_min
        out$d_rank <- pr$d_rank
        if (!is.null(pr$bijective_assigned))
            out$bijective_assigned <- pr$bijective_assigned
    } else {
        out$lr_true <- pr$lr_true
        out$lr_rank <- pr$lr_rank
    }
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
    if (!is.null(x$d_rank)) {
        cat(sprintf("  True rank:   %s\n",
                    if (is.na(x$d_rank)) "N/A" else as.character(x$d_rank)))
        cat(sprintf("  d_true:      %s\n",
                    if (is.na(x$d_true)) "N/A" else sprintf("%.4f", x$d_true)))
        cat(sprintf("  d_min:       %s\n",
                    if (is.na(x$d_min)) "N/A" else sprintf("%.4f", x$d_min)))
    } else {
        cat(sprintf("  True rank:   %s\n",
                    if (is.null(x$lr_rank) || is.na(x$lr_rank)) "N/A"
                    else as.character(x$lr_rank)))
        cat(sprintf("  lr_true:     %s\n",
                    if (is.null(x$lr_true) || is.na(x$lr_true)) "N/A"
                    else sprintf("%.4f", x$lr_true)))
    }
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
.choose_guess_set <- function(d, cand, strategy, k = NULL) {
    if (strategy == "topk") {
        if (is.null(k)) stop("strategy 'topk' requires 'k'.", call. = FALSE)
        k <- as.integer(k)
        if (!is.finite(k) || k < 1L) stop("'k' must be >= 1.", call. = FALSE)
        ord <- order(d)
        if (length(ord) <= k) return(cand)
        kth <- d[ord[k]]
        return(cand[d <= kth])
    }
    dmin <- min(d)
    cand[d == dmin]
}

#' @keywords internal
.harmonic_risk <- function(d) {
    n <- length(d)
    if (n == 0L) return(numeric(0))
    if (n == 1L) return(1)
    r <- rank(d, ties.method = "min")
    w <- 1 / r
    w / sum(w)
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
.discretize <- function(x, brks) {
    as.integer(cut(x, breaks = c(-Inf, brks, Inf), right = TRUE,
                   include.lowest = TRUE))
}

#' @keywords internal
.fs_estimate <- function(X, x_anon, key, true_idx, n_sample = 500L,
                         breaks = NULL) {
    n <- nrow(X)
    nv <- length(key)

    if (n <= 1L) {
        return(list(
            m = setNames(rep(0.99, nv), key),
            u = setNames(rep(0.01, nv), key)
        ))
    }

    # Pre-discretize numeric variables that have user-supplied breaks
    X_est     <- X
    x_anon_est <- x_anon
    for (v in intersect(names(breaks), key)) {
        X_est[[v]]      <- .discretize(X[[v]],      breaks[[v]])
        x_anon_est[[v]] <- .discretize(x_anon[[v]], breaks[[v]])
    }

    # m-probabilities: exact agreement rates among true matched pairs
    m <- numeric(nv)
    names(m) <- key
    for (vi in seq_along(key)) {
        v <- key[vi]
        agree <- vapply(seq_len(n), function(i) {
            xv <- X_est[[v]][i]
            av <- x_anon_est[[v]][true_idx[i]]
            if (is.na(xv) || is.na(av)) return(NA)
            as.character(xv) == as.character(av)
        }, logical(1))
        m[vi] <- mean(agree, na.rm = TRUE)
        if (is.na(m[vi]) || m[vi] <= 0) m[vi] <- 0.01
        if (m[vi] >= 1) m[vi] <- 0.99
    }

    # u-probabilities: exact agreement rates among non-matched pairs.
    # Full enumeration when all distinct pairs fit within n_sample; otherwise
    # sample without replacement (no duplicate pairs).
    u <- numeric(nv)
    names(u) <- key
    n_anon   <- nrow(x_anon)
    n_total  <- n * (n_anon - 1L)   # distinct non-matching ordered pairs
    if (n_total <= n_sample) {
        # Enumerate all (i, j) pairs where j != true_idx[i]
        all_j  <- seq_len(n_anon)
        pairs  <- do.call(rbind, lapply(seq_len(n), function(i) {
            js <- all_j[all_j != true_idx[i]]
            cbind(i = i, j = js)
        }))
        sample_i <- pairs[, "i"]
        sample_j <- pairs[, "j"]
    } else {
        idx      <- sample.int(n_total, n_sample, replace = FALSE)
        sample_i <- ((idx - 1L) %/% (n_anon - 1L)) + 1L
        offset   <- ((idx - 1L) %% (n_anon - 1L)) + 1L
        sample_j <- ifelse(offset < true_idx[sample_i], offset,
                           offset + 1L)
    }
    n_pairs <- length(sample_i)

    for (vi in seq_along(key)) {
        v <- key[vi]
        agree <- vapply(seq_len(n_pairs), function(idx) {
            xv <- X_est[[v]][sample_i[idx]]
            av <- x_anon_est[[v]][sample_j[idx]]
            if (is.na(xv) || is.na(av)) return(NA)
            as.character(xv) == as.character(av)
        }, logical(1))
        u[vi] <- mean(agree, na.rm = TRUE)
        if (is.na(u[vi]) || u[vi] <= 0) u[vi] <- 0.01
        if (u[vi] >= 1) u[vi] <- 0.99
    }

    list(m = m, u = u)
}

#' @keywords internal
.fs_log_lr <- function(x_row, anon_block, key, m_probs, u_probs,
                       breaks = NULL) {
    m <- nrow(anon_block)
    lr <- numeric(m)

    for (vi in seq_along(key)) {
        v <- key[vi]
        xv <- x_row[[v]]
        av <- anon_block[[v]]

        if (!is.null(breaks[[v]])) {
            xv <- .discretize(xv, breaks[[v]])
            av <- .discretize(av, breaks[[v]])
        }

        mv <- m_probs[v]
        uv <- u_probs[v]
        na <- is.na(av) | is.na(xv)

        # Exact agreement for all variable types; NA pairs drop out (contribute 0)
        gamma <- ifelse(!na, as.character(av) == as.character(xv), FALSE)
        contrib <- ifelse(gamma, log(mv / uv), log((1 - mv) / (1 - uv)))
        contrib[na] <- 0
        lr <- lr + contrib
    }

    lr
}


# ---------- S3 methods ----------

#' Print method for recordLinkageRisk
#'
#' @param x object of class \code{"recordLinkageRisk"}.
#' @param ... ignored.
#' @return The input object, invisibly.
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
    cat("Weighting:   ", s$control$risk_weighting, "\n", sep = "")
    cat("Keys:        ", paste(s$key, collapse = ", "), "\n", sep = "")
    rec_label <- if (dir == "reverse") " (risk per synthetic record)" else ""
    cat("Records:     ", x$n_original, " original, ",
        x$n_synthetic, " synthetic", rec_label, "\n", sep = "")
    cat("Truth:       ", s$truth,
        if (!is.null(s$id)) paste0(" (id: ", s$id, ")") else "", "\n", sep = "")

    if (meth == "distance-based") {
        cat("Strategy:    ", s$control$strategy, "\n", sep = "")
        if (!is.null(s$block))
            cat("Blocking:    ", paste(s$block, collapse = ", "), "\n", sep = "")
        cat("NA handling: ", s$control$na_anon, "\n", sep = "")
    } else if (!is.null(s$block)) {
        cat("Blocking:    ", paste(s$block, collapse = ", "), "\n", sep = "")
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
#' @return A list of summary statistics for the corresponding object.
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

    # True-match score statistics (d_true for distance-based, lr_true for probabilistic)
    score_col <- if (object$method == "distance-based") "d_true" else "lr_true"
    score_vals <- object$per_record[[score_col]]
    score_valid <- score_vals[!is.na(score_vals)]
    score_stats <- if (length(score_valid) > 0) {
        c(mean = mean(score_valid),
          median = stats::median(score_valid),
          sd = stats::sd(score_valid),
          min = min(score_valid),
          max = max(score_valid))
    } else {
        NULL
    }

    summ <- list(
        method = object$method,
        direction = if (is.null(object$direction)) "forward" else object$direction,
        matching = if (!is.null(object$settings$matching)) object$settings$matching
                   else "independent",
        key_vars = object$key_vars,
        n_original = object$n_original,
        n_synthetic = object$n_synthetic,
        mean_risk = object$overall$mean_risk,
        max_risk = object$overall$max_risk,
        risk_quantiles = risk_quantiles,
        risk_bands = risk_bands,
        score_stats = score_stats,
        pct_true_in_set = object$overall$pct_true_in_set,
        mean_candidate_size = object$overall$mean_candidate_size,
        risk_gini = object$overall$risk_gini,
        var_importance = object$var_importance,
        privacy_pass = object$privacy_pass,
        fs_params = object$fs_params
    )

    class(summ) <- "summary.recordLinkageRisk"
    summ
}


#' Print method for summary.recordLinkageRisk
#'
#' @param x object of class \code{"summary.recordLinkageRisk"}.
#' @param ... ignored.
#' @return The input object, invisibly.
#' @export
print.summary.recordLinkageRisk <- function(x, ...) {
    dir <- if (is.null(x$direction)) "forward" else x$direction

    cat("Summary: Record Linkage Risk\n")
    cat("============================\n\n")

    mtch <- if (!is.null(x$matching)) x$matching else "independent"
    if (x$method == "distance-based") {
        cat("Method:", x$method, "| Direction:", dir, "\n")
        if (mtch == "bijective")
            cat("Matching: bijective (Hungarian algorithm)\n")
    } else {
        cat("Method:", x$method, "\n")
    }
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

    if (!is.null(x$score_stats)) {
        s_label <- if (x$method == "probabilistic") "Log-LR for True Match:"
                   else "Distance to True Match:"
        cat("\n", s_label, "\n", sep = "")
        cat(sprintf("  Mean:    %7.4f\n", x$score_stats["mean"]))
        cat(sprintf("  Median:  %7.4f\n", x$score_stats["median"]))
        cat(sprintf("  SD:      %7.4f\n", x$score_stats["sd"]))
    }

    if (!is.null(x$fs_params)) {
        cat("\nFellegi-Sunter Parameters:\n")
        cat("  Variable    m_prob   u_prob\n")
        for (v in names(x$fs_params$m_probs)) {
            cat(sprintf("  %-12s %5.3f    %5.3f\n",
                        v, x$fs_params$m_probs[v], x$fs_params$u_probs[v]))
        }
    }

    if (!is.null(x$var_importance)) {
        cat("\nVariable Importance (mean weighted distance):\n")
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
#'   4 = (unused),
#'   5 = Risk band barplot,
#'   6 = True-match rank distribution,
#'   7 = Risk by group barplot (requires \code{group}),
#'   8 = Lorenz curve of risk concentration
#' @param group optional grouping vector or column name (with \code{data})
#'   for \code{which = 7}.
#' @param data optional data frame for extracting \code{group} column.
#' @importFrom graphics hist abline legend par plot plot.new points text barplot lines polygon segments axis box
#' @return No return value; called for the side effect of producing a plot.
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
        if (x$method == "distance-based") {
            sc <- x$per_record$d_true
            xlab2 <- "Distance to True Match"
        } else {
            sc <- x$per_record$lr_true
            xlab2 <- "Log-LR (True Match)"
        }
        risk <- x$per_record$risk
        valid <- !is.na(sc)
        if (sum(valid) > 0) {
            plot(sc[valid], risk[valid],
                 main = paste0(xlab2, " vs Risk", dir_suffix),
                 xlab = xlab2,
                 ylab = "Re-identification Risk",
                 pch = 16, col = adjustcolor("steelblue", 0.5),
                 ...)
            abline(h = thresh, col = "red", lty = 2, lwd = 2)
        } else {
            plot.new()
            text(0.5, 0.5, "No score data available", cex = 1.2)
        }
    }

    if (show[3]) {
        if (x$method == "probabilistic" && !is.null(x$fs_params)) {
            m <- x$fs_params$m_probs
            u <- x$fs_params$u_probs
            mat <- rbind(agree    = log(m / u),
                         disagree = log((1 - m) / (1 - u)))
            colnames(mat) <- names(m)
            barplot(mat, beside = TRUE,
                    main = "Per-Variable Log-Likelihood Ratios",
                    ylab = "log(LR)",
                    col = c("steelblue", "coral"), las = 2, ...)
            abline(h = 0, col = "gray40", lty = 1)
            legend("topright", legend = c("Agree", "Disagree"),
                   fill = c("steelblue", "coral"), cex = 0.8)
        } else if (!is.null(x$var_importance)) {
            vi <- x$var_importance
            barplot(vi[order(vi, decreasing = TRUE)],
                    main = paste0("Variable Importance", dir_suffix),
                    ylab = "Mean Weighted Distance",
                    col = "steelblue", las = 2, ...)
        } else {
            plot.new()
            text(0.5, 0.5, "No variable importance available", cex = 1.2)
        }
    }

    if (show[4]) {
        plot.new()
        text(0.5, 0.5, "Plot 4 not available", cex = 1.2)
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
        rank_col <- if (x$method == "distance-based") "d_rank" else "lr_rank"
        dr <- x$per_record[[rank_col]]
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
