# tests/testthat/test-sinkhorn-internal.R

# --- .sinkhorn() tests ---

test_that(".sinkhorn returns doubly-stochastic-like transport plan", {
  set.seed(1)
  C <- matrix(abs(rnorm(25)), 5, 5)
  P <- .sinkhorn(C, epsilon = 0.1)

  expect_true(is.matrix(P))
  expect_equal(dim(P), c(5, 5))
  expect_true(all(P >= 0))
  # Row sums should be approximately equal (balanced OT)
  expect_true(all(abs(rowSums(P) - 1/5) < 1e-4))
  # Column sums should be approximately equal
  expect_true(all(abs(colSums(P) - 1/5) < 1e-4))
})

test_that(".sinkhorn concentrates mass on low-cost entries", {
  # Diagonal cost = 0, off-diagonal = 1 -> mass should be on diagonal
  C <- matrix(1, 4, 4)
  diag(C) <- 0
  P <- .sinkhorn(C, epsilon = 0.01)

  # Diagonal entries should be much larger than off-diagonal
  expect_true(all(diag(P) > 0.2 * (1/4)))
  # With small epsilon, nearly all mass on diagonal
  expect_true(sum(diag(P)) > 0.9 * sum(P))
})

test_that(".sinkhorn with large epsilon approaches uniform", {
  C <- matrix(c(0, 1, 1, 0), 2, 2)
  P_sharp <- .sinkhorn(C, epsilon = 0.001)
  P_smooth <- .sinkhorn(C, epsilon = 100)

  # Large epsilon: all entries near 1/n = 0.25
  expect_true(all(abs(P_smooth - 0.25) < 0.05))
  # Small epsilon: diagonal entries dominate
  expect_true(P_sharp[1,1] > P_smooth[1,1])
})

test_that(".sinkhorn handles rectangular matrices (n_orig != n_anon)", {
  set.seed(1)
  # 3 original, 5 anonymized
  C <- matrix(abs(rnorm(15)), 3, 5)
  P <- .sinkhorn(C, epsilon = 0.1)

  expect_equal(dim(P), c(3, 5))
  expect_true(all(P >= 0))
  # Row sums = 1/n_orig, col sums = 1/n_anon
  expect_true(all(abs(rowSums(P) - 1/3) < 1e-4))
  expect_true(all(abs(colSums(P) - 1/5) < 1e-4))
})

test_that(".sinkhorn handles identical cost rows gracefully", {
  C <- matrix(rep(c(0, 1, 1, 0.5), each = 4), 4, 4)
  P <- .sinkhorn(C, epsilon = 0.1)

  expect_equal(dim(P), c(4, 4))
  expect_true(all(is.finite(P)))
  expect_true(all(P >= 0))
})

test_that(".sinkhorn converges with discrete (categorical) costs", {
  C <- matrix(c(0, 1, 1, 0,
                1, 0, 1, 1,
                1, 1, 0, 1,
                0, 1, 1, 0), 4, 4, byrow = TRUE)
  P <- .sinkhorn(C, epsilon = 0.05)

  expect_true(all(is.finite(P)))
  expect_true(all(P >= 0))
  # Mass concentrated on diagonal (exact matches)
  expect_true(P[1,1] > P[1,2])
  expect_true(P[2,2] > P[2,1])
})

test_that(".sinkhorn auto-calibrates epsilon when NULL", {
  set.seed(1)
  C <- matrix(abs(rnorm(25)), 5, 5)
  P <- .sinkhorn(C, epsilon = NULL)

  expect_true(is.matrix(P))
  expect_true(all(is.finite(P)))
  expect_true(all(P >= 0))
})

test_that(".sinkhorn respects max_iter", {
  set.seed(1)
  C <- matrix(abs(rnorm(25)), 5, 5)
  # Very few iterations should still return a valid matrix
  P <- .sinkhorn(C, epsilon = 0.1, max_iter = 3L)

  expect_true(is.matrix(P))
  expect_true(all(is.finite(P)))
  expect_true(all(P >= 0))
})

test_that(".sinkhorn matches Hungarian for small epsilon", {
  skip_if_not_installed("clue")
  set.seed(42)
  C <- matrix(abs(rnorm(16)), 4, 4)

  # Hungarian solution
  sol <- clue::solve_LSAP(C)
  hungarian_pairs <- cbind(1:4, as.integer(sol))

  # Sinkhorn with very small epsilon -> hard assignment
  P <- .sinkhorn(C, epsilon = 0.001, max_iter = 500L)

  # The argmax of each row should match Hungarian
  sinkhorn_pairs <- cbind(1:4, apply(P, 1, which.max))
  expect_equal(sinkhorn_pairs, hungarian_pairs)
})

# --- .solve_ot() tests ---

test_that(".solve_ot returns continuous risk from score_cache", {
  score_cache <- list(
    list(cand = 1:4, scores = c(0.0, 0.5, 0.8, 1.0), maximize = FALSE),
    list(cand = 1:4, scores = c(0.6, 0.0, 0.7, 0.9), maximize = FALSE),
    list(cand = 1:4, scores = c(0.9, 0.8, 0.0, 0.3), maximize = FALSE),
    list(cand = 1:4, scores = c(1.0, 0.7, 0.4, 0.0), maximize = FALSE)
  )
  true_idx <- 1:4
  split_search <- list(all = 1:4)
  blk_query <- rep("all", 4)

  res <- .solve_ot(score_cache, true_idx, 4L,
                    split_search, blk_query, epsilon = 0.05)

  expect_true(is.list(res))
  expect_true(is.numeric(res$risk))
  expect_equal(length(res$risk), 4)
  # Risk should be continuous, not just 0/1
  expect_true(any(res$risk > 0 & res$risk < 1))
  # Diagonal has zero cost -> true matches should get high risk
  expect_true(all(res$risk > 0.1))
  # Transport plan should be returned
  expect_true(is.list(res$transport_plans))
})

test_that(".solve_ot with maximize = TRUE (RF proximity)", {
  score_cache <- list(
    list(cand = 1:3, scores = c(0.9, 0.1, 0.2), maximize = TRUE),
    list(cand = 1:3, scores = c(0.2, 0.8, 0.1), maximize = TRUE),
    list(cand = 1:3, scores = c(0.1, 0.3, 0.7), maximize = TRUE)
  )
  true_idx <- 1:3
  split_search <- list(all = 1:3)
  blk_query <- rep("all", 3)

  res <- .solve_ot(score_cache, true_idx, 3L,
                    split_search, blk_query, epsilon = 0.05)

  expect_equal(length(res$risk), 3)
  # High proximity to true -> high risk
  expect_true(all(res$risk > 0.1))
})

test_that(".solve_ot handles blocked data", {
  score_cache <- list(
    list(cand = 1:2, scores = c(0.0, 0.5), maximize = FALSE),
    list(cand = 1:2, scores = c(0.6, 0.0), maximize = FALSE),
    list(cand = 3:4, scores = c(0.0, 0.3), maximize = FALSE),
    list(cand = 3:4, scores = c(0.4, 0.0), maximize = FALSE)
  )
  true_idx <- c(1, 2, 3, 4)
  split_search <- list(A = 1:2, B = 3:4)
  blk_query <- c("A", "A", "B", "B")

  res <- .solve_ot(score_cache, true_idx, 4L,
                    split_search, blk_query, epsilon = 0.05)

  expect_equal(length(res$risk), 4)
  expect_true(all(res$risk > 0))
})

test_that(".solve_ot with unequal sizes", {
  score_cache <- list(
    list(cand = 1:5, scores = c(0.0, 0.5, 0.8, 1.0, 0.3), maximize = FALSE),
    list(cand = 1:5, scores = c(0.6, 0.0, 0.7, 0.9, 0.4), maximize = FALSE),
    list(cand = 1:5, scores = c(0.9, 0.8, 0.0, 0.3, 0.7), maximize = FALSE)
  )
  true_idx <- c(1, 2, 3)
  split_search <- list(all = 1:5)
  blk_query <- rep("all", 3)

  res <- .solve_ot(score_cache, true_idx, 3L,
                    split_search, blk_query, epsilon = 0.1)

  expect_equal(length(res$risk), 3)
  expect_true(all(is.finite(res$risk)))
  expect_true(all(res$risk >= 0))
})

test_that(".solve_ot risk is higher when OT assigns to true match", {
  # Record 1: true=col1, cost 0.0 -> OT assigns here, HIGH risk
  # Record 2: true=col2, but col3 is cheapest (0.0 vs 0.9)
  #   OT prefers 2->3, so P[2,2] is low -> LOW risk
  # Record 3: true=col3, but col2 is cheapest (0.0 vs 0.5)
  #   OT prefers 3->2, so P[3,3] is low -> LOW risk
  score_cache <- list(
    list(cand = 1:3, scores = c(0.0, 0.5, 0.9), maximize = FALSE),
    list(cand = 1:3, scores = c(0.1, 0.9, 0.0), maximize = FALSE),
    list(cand = 1:3, scores = c(0.9, 0.0, 0.5), maximize = FALSE)
  )
  true_idx <- 1:3
  split_search <- list(all = 1:3)
  blk_query <- rep("all", 3)

  res <- .solve_ot(score_cache, true_idx, 3L,
                    split_search, blk_query, epsilon = 0.01)

  # Record 1 correctly linked -> high risk; records 2,3 swapped -> low risk
  expect_true(res$risk[1] > res$risk[2])
  expect_true(res$risk[1] > res$risk[3])
})

# --- Additional edge case tests ---

test_that(".sinkhorn handles 1x1 cost matrix", {
  C <- matrix(0.5, 1, 1)
  P <- .sinkhorn(C, epsilon = 0.1)
  expect_equal(dim(P), c(1, 1))
  expect_equal(P[1, 1], 1, tolerance = 1e-6)
})

test_that(".sinkhorn handles all-zero cost matrix", {
  C <- matrix(0, 4, 4)
  P <- .sinkhorn(C, epsilon = 0.1)
  expect_equal(dim(P), c(4, 4))
  expect_true(all(is.finite(P)))
  # All entries should be uniform: 1/16
  expect_true(all(abs(P - 1/16) < 1e-4))
})

test_that(".sinkhorn handles constant cost matrix", {
  C <- matrix(0.5, 4, 4)
  P <- .sinkhorn(C, epsilon = 0.1)
  expect_true(all(is.finite(P)))
  # All entries should be uniform: 1/16
  expect_true(all(abs(P - 1/16) < 1e-4))
})

test_that(".sinkhorn with epsilon = 0 falls back to 0.01", {
  C <- matrix(c(0, 1, 1, 0), 2, 2)
  expect_silent(P <- .sinkhorn(C, epsilon = 0))
  expect_true(all(is.finite(P)))
  expect_true(all(P >= 0))
})

test_that(".sinkhorn with negative epsilon falls back to 0.01", {
  C <- matrix(c(0, 1, 1, 0), 2, 2)
  expect_silent(P <- .sinkhorn(C, epsilon = -1))
  expect_true(all(is.finite(P)))
})

test_that(".sinkhorn warns on very small epsilon causing non-finite values", {
  # Matrix where row-shift still leaves cost differences that cause underflow
  C <- matrix(c(1, 2, 3, 1.5, 2.5, 3.5, 2, 3, 4), 3, 3)
  expect_warning(.sinkhorn(C, epsilon = 1e-300), "non-finite")
})

test_that(".sinkhorn handles 1x5 rectangular matrix", {
  C <- matrix(c(0.1, 0.5, 0.9, 0.3, 0.7), 1, 5)
  P <- .sinkhorn(C, epsilon = 0.1)
  expect_equal(dim(P), c(1, 5))
  expect_true(all(is.finite(P)))
  expect_equal(sum(P), 1, tolerance = 1e-4)  # single row sums to 1/1 = 1
})

test_that(".solve_ot handles nq > ns (more query than search)", {
  score_cache <- list(
    list(cand = 1:2, scores = c(0.0, 0.5), maximize = FALSE),
    list(cand = 1:2, scores = c(0.6, 0.0), maximize = FALSE),
    list(cand = 1:2, scores = c(0.3, 0.4), maximize = FALSE)
  )
  true_idx <- c(1, 2, NA)
  split_search <- list(all = 1:2)
  blk_query <- rep("all", 3)

  res <- .solve_ot(score_cache, true_idx, 3L,
                    split_search, blk_query, epsilon = 0.1)
  expect_equal(length(res$risk), 3)
  expect_true(all(is.finite(res$risk[1:2])))
  # Record 3 has NA true_idx -> risk stays 0
  expect_equal(res$risk[3], 0)
})

test_that(".solve_ot handles NULL score_cache entries", {
  score_cache <- list(
    list(cand = 1:3, scores = c(0.0, 0.5, 0.8), maximize = FALSE),
    NULL,  # missing entry
    list(cand = 1:3, scores = c(0.9, 0.8, 0.0), maximize = FALSE)
  )
  true_idx <- 1:3
  split_search <- list(all = 1:3)
  blk_query <- rep("all", 3)

  res <- .solve_ot(score_cache, true_idx, 3L,
                    split_search, blk_query, epsilon = 0.1)
  expect_equal(length(res$risk), 3)
  expect_true(all(is.finite(res$risk)))
})

test_that(".solve_ot analytical verification: uniform plan gives risk 1/ns", {
  # With very large epsilon, plan approaches uniform
  # For 3x3: P[i,j] ≈ 1/9, risk = P[i,true] * 3 ≈ 1/3
  score_cache <- list(
    list(cand = 1:3, scores = c(0.1, 0.2, 0.3), maximize = FALSE),
    list(cand = 1:3, scores = c(0.2, 0.1, 0.3), maximize = FALSE),
    list(cand = 1:3, scores = c(0.3, 0.2, 0.1), maximize = FALSE)
  )
  true_idx <- 1:3
  split_search <- list(all = 1:3)
  blk_query <- rep("all", 3)

  res <- .solve_ot(score_cache, true_idx, 3L,
                    split_search, blk_query, epsilon = 1000)
  # With huge epsilon, all risks should approach 1/3
  expect_true(all(abs(res$risk - 1/3) < 0.05))
})
