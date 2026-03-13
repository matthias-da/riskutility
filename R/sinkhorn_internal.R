#' Sinkhorn Optimal Transport Solver
#'
#' Computes the entropy-regularized optimal transport plan between
#' two discrete distributions using the Sinkhorn-Knopp algorithm.
#'
#' Given a cost matrix \eqn{C} (n_orig x n_anon) and regularization
#' \eqn{\varepsilon > 0}, the transport plan \eqn{P} minimizes
#' \deqn{
#' \sum_{ij} P_{ij} C_{ij} + \varepsilon \sum_{ij} P_{ij} \log P_{ij}
#' }
#' subject to marginal constraints. With \eqn{\varepsilon \to 0},
#' the solution approaches the Hungarian (LSAP) assignment.
#'
#' Supports rectangular matrices (unequal dataset sizes).
#' Marginals default to uniform: \eqn{r_i = 1/n_{orig}},
#' \eqn{c_j = 1/n_{anon}}.
#'
#' @param C numeric matrix. Cost matrix (n_orig x n_anon). Non-negative.
#' @param epsilon numeric or NULL. Regularization strength.
#'   If NULL, auto-calibrated as \code{median(C[C > 0]) * 0.05}.
#'   Smaller = sharper assignment (closer to Hungarian);
#'   larger = smoother (closer to uniform).
#' @param max_iter integer. Maximum Sinkhorn iterations (default 100).
#' @param tol numeric. Convergence tolerance on marginal error (default 1e-8).
#' @return numeric matrix. Transport plan P of same dimensions as C.
#'   Entry \eqn{P_{ij}} is the transport weight from original i to anonymized j.
#'   Row sums = 1/n_orig, column sums = 1/n_anon.
#' @keywords internal
.sinkhorn <- function(C, epsilon = NULL, max_iter = 100L, tol = 1e-8) {
  nr <- nrow(C)
  nc <- ncol(C)

  # Auto-calibrate epsilon
  if (is.null(epsilon)) {
    pos <- C[C > 0]
    epsilon <- if (length(pos) > 0) median(pos) * 0.05 else 0.1
  }
  if (epsilon <= 0) epsilon <- 0.01

  # Uniform marginals
  r <- rep(1 / nr, nr)
  cc <- rep(1 / nc, nc)

  # Gibbs kernel: K = exp(-C / epsilon)
  # Stabilize by subtracting row-wise min to avoid underflow
  C_shifted <- C - apply(C, 1, min)
  K <- exp(-C_shifted / epsilon)

  # Initialize scaling vectors
  u <- rep(1, nr)
  v <- rep(1, nc)

  for (iter in seq_len(max_iter)) {
    u <- r / (K %*% v)
    v <- cc / (t(K) %*% u)

    # Replace non-finite values (numerical instability)
    u[!is.finite(u)] <- 1e-10
    v[!is.finite(v)] <- 1e-10

    # Check convergence: marginal error
    row_err <- max(abs(u * (K %*% v) - r))
    if (row_err < tol) break
  }

  # Transport plan
  P <- diag(as.numeric(u)) %*% K %*% diag(as.numeric(v))

  # Ensure non-negative (numerical cleanup)
  P[P < 0] <- 0

  P
}


#' Solve Optimal Transport Matching from Score Cache
#'
#' Post-loop handler parallel to \code{.solve_bijective()}.
#' Builds cost matrices from \code{score_cache}, runs Sinkhorn per block,
#' and extracts continuous risk scores from the transport plan.
#'
#' Risk for record i is defined as:
#' \deqn{risk_i = P_{i, \text{true}} \cdot n_{\text{query}}}
#' where \eqn{P_{i, \text{true}}} is the transport weight from query
#' record i to its true match in the search data. Under uniform random
#' matching, this gives risk = 1/n_search (same baseline as independent
#' matching).
#'
#' @param score_cache list of length n_query. Each element is
#'   \code{list(cand, scores, maximize)}.
#' @param true_idx integer vector of true search-side indices.
#' @param n_query integer, number of query records.
#' @param split_search named list of search-side indices per block.
#' @param blk_query character vector of block labels per query record.
#' @param epsilon numeric or NULL. Regularization parameter for Sinkhorn.
#' @param max_iter integer. Max Sinkhorn iterations (default 100).
#' @return list with \code{risk} (numeric), \code{d_rank} (integer),
#'   \code{transport_plans} (list of matrices per block).
#' @keywords internal
.solve_ot <- function(score_cache, true_idx, n_query,
                       split_search, blk_query,
                       epsilon = NULL, max_iter = 100L) {
  risk_out <- numeric(n_query)
  d_rank_out <- rep(NA_integer_, n_query)
  transport_plans <- list()

  blocks <- unique(blk_query)
  for (blk in blocks) {
    q_idx <- which(blk_query == blk)
    s_idx <- split_search[[blk]]
    if (is.null(s_idx) || length(s_idx) == 0L || length(q_idx) == 0L)
      next

    nq <- length(q_idx)
    ns <- length(s_idx)

    # Determine direction from first non-NULL cache entry
    maximize <- FALSE
    for (qi in q_idx) {
      if (!is.null(score_cache[[qi]])) {
        maximize <- score_cache[[qi]]$maximize
        break
      }
    }

    # Build cost matrix: rows = query records, cols = search records
    # Use worst-case score as fill for missing entries
    all_scores <- unlist(lapply(q_idx, function(qi) {
      if (!is.null(score_cache[[qi]])) score_cache[[qi]]$scores
    }))
    fill_val <- if (length(all_scores) > 0) {
      if (maximize) min(all_scores) else max(all_scores)
    } else 1

    cost <- matrix(fill_val, nrow = nq, ncol = ns)

    for (r in seq_along(q_idx)) {
      qi <- q_idx[r]
      sc <- score_cache[[qi]]
      if (is.null(sc)) next

      col_pos <- match(sc$cand, s_idx)
      valid <- !is.na(col_pos)
      if (!any(valid)) next

      if (maximize) {
        # Transform to minimization: cost = max - score
        cost[r, col_pos[valid]] <- max(sc$scores) - sc$scores[valid]
      } else {
        cost[r, col_pos[valid]] <- sc$scores[valid]
      }
    }

    # Run Sinkhorn
    P <- .sinkhorn(cost, epsilon = epsilon, max_iter = max_iter)
    transport_plans[[blk]] <- P

    # Extract risk: P[i, true_col] * n_query_in_block
    for (r in seq_along(q_idx)) {
      qi <- q_idx[r]
      tpos <- true_idx[qi]
      if (!is.na(tpos) && tpos > 0L) {
        t_col <- match(tpos, s_idx)
        if (!is.na(t_col)) {
          # Normalize: under uniform matching, P[i,j] = 1/(nq*ns)
          # risk = P[i,true] * nq  -> gives 1/ns for uniform
          risk_out[qi] <- P[r, t_col] * nq
          # Rank: how many search records have higher transport weight?
          d_rank_out[qi] <- as.integer(sum(P[r, ] >= P[r, t_col]))
        }
      }
    }
  }

  # Cap risk at 1
  risk_out <- pmin(risk_out, 1)

  list(risk = risk_out, d_rank = d_rank_out,
       transport_plans = transport_plans)
}
