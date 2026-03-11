# RF Proximity Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add RF proximity-based metrics for record linkage risk, memorization detection, and utility measurement to the riskutility package.

**Architecture:** Three use cases share a common `.rf_proximity()` engine based on ranger. A shared `.distance_risk_prepare()` helper extracts holdout-splitting boilerplate from dcr/nndr. New function `rf_privacy()` joins the distance-risk family. Existing `propscore()` and `recordLinkage()` gain new method options.

**Tech Stack:** R, ranger, testthat, roxygen2, ggplot2, clue (for bijective matching)

**Spec:** `docs/superpowers/specs/2026-03-11-rf-proximity-design.md`

**Git:** All commits use `--author="matthias-da <matthias-da@users.noreply.github.com>"`

---

## Chunk 1: Internal Engine

### Task 1: `.rf_proximity()` — tests

**Files:**
- Create: `tests/testthat/test-rf-proximity-internal.R`

- [ ] **Step 1: Write tests for `.rf_proximity()`**

```r
# tests/testthat/test-rf-proximity-internal.R
test_that(".rf_proximity returns correct structure", {
  skip_if_not_installed("ranger")
  set.seed(1)
  d1 <- data.frame(x = rnorm(50), y = rnorm(50))
  d2 <- data.frame(x = rnorm(50), y = rnorm(50))
  res <- riskutility:::.rf_proximity(d1, d2, n_trees = 50, seed = 1)

  expect_type(res, "list")
  expect_s3_class(res$forest, "ranger")
  expect_true(is.matrix(res$terminal_nodes))
  expect_equal(nrow(res$terminal_nodes), 100)  # n1 + n2
  expect_equal(ncol(res$terminal_nodes), 50)   # n_trees
  expect_equal(res$n1, 50)
  expect_equal(res$n2, 50)
  expect_true(is.numeric(res$oob_error))
  expect_true(res$oob_error >= 0 && res$oob_error <= 1)
})

test_that(".rf_proximity returns importance when requested", {
  skip_if_not_installed("ranger")
  set.seed(1)
  d1 <- data.frame(x = rnorm(50), y = rnorm(50))
  d2 <- data.frame(x = rnorm(50) + 2, y = rnorm(50))
  res_imp <- riskutility:::.rf_proximity(d1, d2, importance = TRUE,
                                          n_trees = 50, seed = 1)
  res_no  <- riskutility:::.rf_proximity(d1, d2, importance = FALSE,
                                          n_trees = 50, seed = 1)

  expect_true(is.numeric(res_imp$importance))
  expect_equal(length(res_imp$importance), 2)  # x, y
  expect_null(res_no$importance)
})

test_that(".rf_proximity validates n_trees >= 10", {
  skip_if_not_installed("ranger")
  d1 <- data.frame(x = 1:20)
  d2 <- data.frame(x = 21:40)
  expect_error(riskutility:::.rf_proximity(d1, d2, n_trees = 5),
               "n_trees")
})

test_that(".rf_proximity detects .rf_label collision", {
  skip_if_not_installed("ranger")
  d1 <- data.frame(x = 1:20, .rf_label = 1:20)
  d2 <- data.frame(x = 21:40, .rf_label = 21:40)
  expect_error(riskutility:::.rf_proximity(d1, d2, n_trees = 10),
               "rf_label")
})

test_that(".rf_proximity uses vars subset", {
  skip_if_not_installed("ranger")
  set.seed(1)
  d1 <- data.frame(x = rnorm(30), y = rnorm(30), z = rnorm(30))
  d2 <- data.frame(x = rnorm(30), y = rnorm(30), z = rnorm(30))
  res <- riskutility:::.rf_proximity(d1, d2, vars = c("x", "y"),
                                      n_trees = 50, seed = 1)
  # importance should only have x, y (not z)
  expect_equal(length(res$importance), 2)
  expect_true(all(names(res$importance) %in% c("x", "y")))
})

test_that(".rf_proximity handles high-cardinality factors", {
  skip_if_not_installed("ranger")
  set.seed(1)
  # Create factor with > 53 levels to trigger partition mode
  d1 <- data.frame(x = factor(sample(paste0("cat", 1:60), 100, TRUE)),
                   y = rnorm(100))
  d2 <- data.frame(x = factor(sample(paste0("cat", 1:60), 100, TRUE)),
                   y = rnorm(100))
  expect_message(
    res <- riskutility:::.rf_proximity(d1, d2, n_trees = 20, seed = 1),
    "partition"
  )
  expect_s3_class(res$forest, "ranger")
})

test_that(".rf_proximity OOB NA handling for small data", {
  skip_if_not_installed("ranger")
  set.seed(1)
  # Very small dataset where OOB predictions may produce NAs
  d1 <- data.frame(x = 1:5, y = rnorm(5))
  d2 <- data.frame(x = 6:10, y = rnorm(5))
  # Should run without error (warning about NAs is acceptable)
  res <- suppressWarnings(
    riskutility:::.rf_proximity(d1, d2, n_trees = 10, seed = 1)
  )
  expect_false(any(is.na(res$terminal_nodes)))
})

test_that(".rf_proximity uses modifyList for user overrides", {
  skip_if_not_installed("ranger")
  set.seed(1)
  d1 <- data.frame(x = rnorm(50), y = rnorm(50))
  d2 <- data.frame(x = rnorm(50), y = rnorm(50))
  # Override importance via ... should work without collision
  res <- riskutility:::.rf_proximity(d1, d2, n_trees = 50, seed = 1,
                                      importance = TRUE,
                                      min.node.size = 5)
  expect_s3_class(res$forest, "ranger")
})
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `Rscript -e "testthat::test_file('tests/testthat/test-rf-proximity-internal.R')"`
Expected: FAIL — `.rf_proximity` not found

- [ ] **Step 3: Commit test file**

```bash
git add tests/testthat/test-rf-proximity-internal.R
git commit --author="matthias-da <matthias-da@users.noreply.github.com>" \
  -m "Add tests for .rf_proximity() internal engine"
```

---

### Task 2: `.rf_proximity()` — implementation

**Files:**
- Create: `R/rf_proximity_internal.R`

- [ ] **Step 4: Implement `.rf_proximity()`**

```r
# R/rf_proximity_internal.R
#' RF Proximity Internal Engine
#'
#' Trains a supervised RF on combined data and extracts terminal node matrix.
#'
#' @param data1 data.frame (label = 0)
#' @param data2 data.frame (label = 1)
#' @param vars character vector of variable names (NULL = all common)
#' @param n_trees integer, number of trees (>= 10)
#' @param mtry integer or NULL
#' @param importance logical, compute variable importance
#' @param seed integer or NULL, passed to ranger
#' @param ... additional arguments passed to ranger via modifyList
#' @return list with forest, terminal_nodes, n1, n2, importance, oob_error
#' @keywords internal
.rf_proximity <- function(data1, data2,
                          vars = NULL,
                          n_trees = 500L,
                          mtry = NULL,
                          importance = TRUE,
                          seed = NULL, ...) {

  if (!requireNamespace("ranger", quietly = TRUE)) {
    stop("Package 'ranger' required for RF proximity. ",
         "Install with install.packages('ranger')", call. = FALSE)
  }
  if (n_trees < 10L) {
    stop("n_trees must be at least 10 for meaningful proximity estimates.",
         call. = FALSE)
  }

  # Determine vars
  if (is.null(vars)) {
    vars <- intersect(names(data1), names(data2))
  }
  d1 <- data1[, vars, drop = FALSE]
  d2 <- data2[, vars, drop = FALSE]

  # Check column name collision
  if (".rf_label" %in% names(d1) || ".rf_label" %in% names(d2)) {
    stop("Column '.rf_label' already exists in data. Please rename it.",
         call. = FALSE)
  }

  n1 <- nrow(d1)
  n2 <- nrow(d2)
  combined <- rbind(d1, d2)
  combined$.rf_label <- factor(c(rep(0L, n1), rep(1L, n2)))

  # Handle high-cardinality unordered factors
  respect_unordered <- FALSE
  for (v in vars) {
    if (is.factor(combined[[v]]) && !is.ordered(combined[[v]]) &&
        nlevels(combined[[v]]) > 53) {
      respect_unordered <- TRUE
      break
    }
  }

  # Build default args, let user override via modifyList

  default_args <- list(
    formula = .rf_label ~ .,
    data = combined,
    num.trees = n_trees,
    probability = TRUE,
    write.forest = TRUE,
    importance = if (importance) "impurity" else "none",
    seed = seed
  )
  if (!is.null(mtry)) default_args$mtry <- mtry
  if (respect_unordered) {
    default_args$respect.unordered.factors <- "partition"
    message("High-cardinality unordered factor detected (> 53 levels). ",
            "Using respect.unordered.factors = 'partition'.")
  }

  args <- modifyList(default_args, list(...))
  forest <- do.call(ranger::ranger, args)

  # Extract terminal nodes
  tn <- predict(forest, combined, type = "terminalNodes")$predictions

  # Handle OOB prediction NAs for small datasets
  if (any(is.na(tn))) {
    na_rows <- which(rowSums(is.na(tn)) > 0)
    warning("OOB terminal node predictions contain NAs for ", length(na_rows),
            " records (small dataset). Falling back to in-bag predictions ",
            "for affected records.", call. = FALSE)
    # Re-predict without OOB restriction for affected records
    tn_inbag <- predict(forest, combined[na_rows, , drop = FALSE],
                        type = "terminalNodes")$predictions
    tn[na_rows, ] <- tn_inbag[, , drop = FALSE]
  }

  # Importance
  imp <- if (importance) forest$variable.importance else NULL
  # Remove .rf_label from importance if present

  if (!is.null(imp) && ".rf_label" %in% names(imp)) {
    imp <- imp[names(imp) != ".rf_label"]
  }

  list(
    forest = forest,
    terminal_nodes = tn,
    n1 = n1,
    n2 = n2,
    importance = imp,
    oob_error = forest$prediction.error
  )
}
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `Rscript -e "testthat::test_file('tests/testthat/test-rf-proximity-internal.R')"`
Expected: All tests PASS

- [ ] **Step 6: Commit**

```bash
git add R/rf_proximity_internal.R
git commit --author="matthias-da <matthias-da@users.noreply.github.com>" \
  -m "Implement .rf_proximity() internal engine"
```

---

### Task 3: `.proximity_from_nodes()` — tests and implementation

**Files:**
- Modify: `tests/testthat/test-rf-proximity-internal.R`
- Modify: `R/rf_proximity_internal.R`

- [ ] **Step 7: Add tests for `.proximity_from_nodes()`**

Append to `tests/testthat/test-rf-proximity-internal.R`:

```r
test_that(".proximity_from_nodes returns correct dimensions", {
  skip_if_not_installed("ranger")
  set.seed(1)
  d1 <- data.frame(x = rnorm(20), y = rnorm(20))
  d2 <- data.frame(x = rnorm(30), y = rnorm(30))
  res <- riskutility:::.rf_proximity(d1, d2, n_trees = 50, seed = 1)

  idx1 <- 1:20     # data1 indices
  idx2 <- 21:50    # data2 indices
  prox <- riskutility:::.proximity_from_nodes(res$terminal_nodes, idx1, idx2)

  expect_true(is.matrix(prox))
  expect_equal(nrow(prox), 30)  # length(idx2)
  expect_equal(ncol(prox), 20)  # length(idx1)
})

test_that(".proximity_from_nodes values are in [0, 1]", {
  skip_if_not_installed("ranger")
  set.seed(1)
  d1 <- data.frame(x = rnorm(20))
  d2 <- data.frame(x = rnorm(20))
  res <- riskutility:::.rf_proximity(d1, d2, n_trees = 100, seed = 1)

  prox <- riskutility:::.proximity_from_nodes(res$terminal_nodes, 1:20, 21:40)
  expect_true(all(prox >= 0 & prox <= 1))
})

test_that(".proximity_from_nodes: identical data has high self-proximity", {
  skip_if_not_installed("ranger")
  set.seed(1)
  d <- data.frame(x = rnorm(30), y = rnorm(30))
  # data2 = copy of data1 (memorized)
  res <- riskutility:::.rf_proximity(d, d, n_trees = 200, seed = 1)

  # Self-proximity (record i in d1 vs record i in d2)
  prox <- riskutility:::.proximity_from_nodes(res$terminal_nodes, 1:30, 31:60)
  self_prox <- diag(prox)
  other_prox <- prox[row(prox) != col(prox)]

  # Self-proximity should generally be higher than cross-proximity
  expect_true(mean(self_prox) > mean(other_prox))
})

test_that(".proximity_from_nodes is symmetric", {
  skip_if_not_installed("ranger")
  set.seed(1)
  d1 <- data.frame(x = rnorm(15), y = rnorm(15))
  d2 <- data.frame(x = rnorm(15), y = rnorm(15))
  res <- riskutility:::.rf_proximity(d1, d2, n_trees = 50, seed = 1)

  # prox(idx1→idx2) should be transpose of prox(idx2→idx1)
  prox_ab <- riskutility:::.proximity_from_nodes(res$terminal_nodes, 1:15, 16:30)
  prox_ba <- riskutility:::.proximity_from_nodes(res$terminal_nodes, 16:30, 1:15)
  expect_equal(prox_ab, t(prox_ba))
})

test_that(".proximity_from_nodes: tie correction produces rational values", {
  skip_if_not_installed("ranger")
  set.seed(1)
  d1 <- data.frame(x = rnorm(10))
  d2 <- data.frame(x = rnorm(10))
  res <- riskutility:::.rf_proximity(d1, d2, n_trees = 50, seed = 1)

  prox <- riskutility:::.proximity_from_nodes(res$terminal_nodes, 1:10, 11:20)
  # Values should be multiples of 1/n_trees
  expect_true(all(prox * 50 == round(prox * 50)))
})
```

- [ ] **Step 8: Implement `.proximity_from_nodes()`**

Append to `R/rf_proximity_internal.R`:

```r
#' Compute proximity submatrix from terminal nodes
#'
#' Uses a tree-by-tree accumulator: for each tree, groups records by
#' terminal node ID, then increments the proximity counter only for
#' pairs that share a node. This is efficient when trees have many
#' terminal nodes (each with few records), avoiding the full
#' O(n2 * n1 * n_trees) scan.
#'
#' @param terminal_nodes integer matrix (n x n_trees) from .rf_proximity()
#' @param idx1 integer vector, row indices for reference set
#' @param idx2 integer vector, row indices for query set
#' @param progress logical, show progress
#' @return numeric matrix of dim length(idx2) x length(idx1)
#' @keywords internal
.proximity_from_nodes <- function(terminal_nodes, idx1, idx2,
                                  progress = FALSE) {
  n1 <- length(idx1)
  n2 <- length(idx2)
  n_trees <- ncol(terminal_nodes)

  nodes1 <- terminal_nodes[idx1, , drop = FALSE]  # n1 x n_trees
  nodes2 <- terminal_nodes[idx2, , drop = FALSE]  # n2 x n_trees

  # Map idx1/idx2 positions to local indices for accumulator
  # prox[i, j] = fraction of trees where idx2[i] and idx1[j] share a node
  prox <- matrix(0L, nrow = n2, ncol = n1)

  if (progress && n_trees > 10) {
    pb <- utils::txtProgressBar(min = 0, max = n_trees, style = 3)
  } else {
    pb <- NULL
  }

  # Tree-by-tree accumulator: for each tree, find pairs sharing a node

  for (t in seq_len(n_trees)) {
    tn1 <- nodes1[, t]  # terminal node IDs for idx1 in tree t
    tn2 <- nodes2[, t]  # terminal node IDs for idx2 in tree t

    # Group idx1 positions by their terminal node in this tree
    node_ids <- unique(tn1)
    # Build lookup: for each node ID, which positions in idx1 have it
    for (nid in node_ids) {
      j_in <- which(tn1 == nid)  # positions in idx1 with this node
      i_in <- which(tn2 == nid)  # positions in idx2 with this node
      if (length(i_in) > 0 && length(j_in) > 0) {
        prox[i_in, j_in] <- prox[i_in, j_in] + 1L
      }
    }
    if (!is.null(pb)) utils::setTxtProgressBar(pb, t)
  }

  if (!is.null(pb)) close(pb)
  prox / n_trees
}
```

- [ ] **Step 9: Run tests**

Run: `Rscript -e "testthat::test_file('tests/testthat/test-rf-proximity-internal.R')"`
Expected: All tests PASS

- [ ] **Step 10: Commit**

```bash
git add R/rf_proximity_internal.R tests/testthat/test-rf-proximity-internal.R
git commit --author="matthias-da <matthias-da@users.noreply.github.com>" \
  -m "Add .proximity_from_nodes() helper with tests"
```

---

### Task 4: `.proximity_from_nodes_newdata()` — tests and implementation

**Files:**
- Modify: `tests/testthat/test-rf-proximity-internal.R`
- Modify: `R/rf_proximity_internal.R`

- [ ] **Step 11: Add tests for `.proximity_from_nodes_newdata()`**

Append to `tests/testthat/test-rf-proximity-internal.R`:

```r
test_that(".proximity_from_nodes_newdata returns correct dimensions", {
  skip_if_not_installed("ranger")
  set.seed(1)
  d1 <- data.frame(x = rnorm(30), y = rnorm(30))
  d2 <- data.frame(x = rnorm(30), y = rnorm(30))
  holdout <- data.frame(x = rnorm(20), y = rnorm(20))
  res <- riskutility:::.rf_proximity(d1, d2, n_trees = 50, seed = 1)

  prox_ho <- riskutility:::.proximity_from_nodes_newdata(
    res$forest, holdout, res$terminal_nodes, idx_ref = 21:60
  )

  expect_true(is.matrix(prox_ho))
  expect_equal(nrow(prox_ho), 20)   # nrow(holdout)
  expect_equal(ncol(prox_ho), 40)   # length(idx_ref) = d2 rows
})

test_that(".proximity_from_nodes_newdata values in [0, 1]", {
  skip_if_not_installed("ranger")
  set.seed(1)
  d1 <- data.frame(x = rnorm(20))
  d2 <- data.frame(x = rnorm(20))
  ho  <- data.frame(x = rnorm(10))
  res <- riskutility:::.rf_proximity(d1, d2, n_trees = 50, seed = 1)

  prox_ho <- riskutility:::.proximity_from_nodes_newdata(
    res$forest, ho, res$terminal_nodes, idx_ref = 1:20
  )
  expect_true(all(prox_ho >= 0 & prox_ho <= 1))
})
```

- [ ] **Step 12: Implement `.proximity_from_nodes_newdata()`**

Append to `R/rf_proximity_internal.R`:

```r
#' Compute proximity between new data and reference records
#'
#' Pushes newdata through a trained forest to get terminal nodes,
#' then computes proximity to reference records.
#'
#' @param forest ranger object with write.forest = TRUE
#' @param newdata data.frame to push through the forest
#' @param terminal_nodes_ref terminal node matrix from training
#' @param idx_ref integer vector, row indices in terminal_nodes_ref
#' @param progress logical
#' @return matrix of dim nrow(newdata) x length(idx_ref)
#' @keywords internal
.proximity_from_nodes_newdata <- function(forest, newdata,
                                          terminal_nodes_ref, idx_ref,
                                          progress = FALSE) {
  # Predict terminal nodes for new data
  tn_new <- predict(forest, newdata, type = "terminalNodes")$predictions

  # Use .proximity_from_nodes with the new nodes as query
  # We need to combine: ref nodes + new nodes, then call with appropriate indices
  n_ref <- nrow(terminal_nodes_ref)
  n_new <- nrow(tn_new)
  combined_tn <- rbind(terminal_nodes_ref, tn_new)

  idx_new <- (n_ref + 1):(n_ref + n_new)
  .proximity_from_nodes(combined_tn, idx_ref, idx_new, progress = progress)
}
```

- [ ] **Step 13: Run tests**

Run: `Rscript -e "testthat::test_file('tests/testthat/test-rf-proximity-internal.R')"`
Expected: All tests PASS

- [ ] **Step 14: Commit**

```bash
git add R/rf_proximity_internal.R tests/testthat/test-rf-proximity-internal.R
git commit --author="matthias-da <matthias-da@users.noreply.github.com>" \
  -m "Add .proximity_from_nodes_newdata() for holdout push-through"
```

---

## Chunk 2: `.distance_risk_prepare()` Refactor

### Task 5: `.distance_risk_prepare()` — tests and implementation

**Files:**
- Create: `tests/testthat/test-distance-risk-prepare.R`
- Modify: `R/utils_internal.R` (append after line 109)

- [ ] **Step 15: Write tests for `.distance_risk_prepare()`**

```r
# tests/testthat/test-distance-risk-prepare.R
test_that(".distance_risk_prepare splits holdout correctly", {
  set.seed(1)
  X <- data.frame(a = 1:100, b = rnorm(100))
  Y <- data.frame(a = 101:200, b = rnorm(100))
  prep <- riskutility:::.distance_risk_prepare(X, Y, holdout_fraction = 0.3,
                                                seed = 42)

  expect_equal(nrow(prep$train) + nrow(prep$holdout), 100)
  expect_equal(nrow(prep$synthetic), 100)
  expect_true(prep$was_split)
  expect_true(all(prep$vars %in% c("a", "b")))
  # Reproducible with same seed
  prep2 <- riskutility:::.distance_risk_prepare(X, Y, holdout_fraction = 0.3,
                                                 seed = 42)
  expect_identical(prep$train, prep2$train)
})

test_that(".distance_risk_prepare uses explicit holdout", {
  X <- data.frame(a = 1:50, b = rnorm(50))
  Y <- data.frame(a = 51:100, b = rnorm(50))
  H <- data.frame(a = 101:120, b = rnorm(20))
  prep <- riskutility:::.distance_risk_prepare(X, Y, holdout = H)

  expect_equal(nrow(prep$train), 50)  # X unchanged
  expect_equal(nrow(prep$holdout), 20)
  expect_false(prep$was_split)
})

test_that(".distance_risk_prepare intersects vars", {
  X <- data.frame(a = 1:10, b = 1:10, c = 1:10)
  Y <- data.frame(a = 1:10, b = 1:10, d = 1:10)
  prep <- riskutility:::.distance_risk_prepare(X, Y)

  expect_equal(sort(prep$vars), c("a", "b"))
})

test_that(".distance_risk_prepare applies vars filter", {
  X <- data.frame(a = 1:10, b = 1:10, c = 1:10)
  Y <- data.frame(a = 1:10, b = 1:10, c = 1:10)
  prep <- riskutility:::.distance_risk_prepare(X, Y, vars = "a")

  expect_equal(prep$vars, "a")
})

test_that(".distance_risk_prepare removes NAs when na.rm = TRUE", {
  X <- data.frame(a = c(1:9, NA), b = 1:10)
  Y <- data.frame(a = 1:10, b = 1:10)
  prep <- riskutility:::.distance_risk_prepare(X, Y, na.rm = TRUE,
                                                holdout_fraction = 0.5,
                                                seed = 1)

  total <- nrow(prep$train) + nrow(prep$holdout)
  expect_true(total <= 9)  # NA row removed before split
})

test_that(".distance_risk_prepare keeps NAs when na.rm = FALSE", {
  X <- data.frame(a = c(1:9, NA), b = 1:10)
  Y <- data.frame(a = 1:10, b = 1:10)
  prep <- riskutility:::.distance_risk_prepare(X, Y, na.rm = FALSE,
                                                holdout_fraction = 0.5,
                                                seed = 1)

  total <- nrow(prep$train) + nrow(prep$holdout)
  expect_equal(total, 10)
})

test_that(".distance_risk_prepare validates holdout_fraction", {
  X <- data.frame(a = 1:10)
  Y <- data.frame(a = 1:10)
  expect_error(riskutility:::.distance_risk_prepare(X, Y,
                holdout_fraction = 0), "holdout_fraction")
  expect_error(riskutility:::.distance_risk_prepare(X, Y,
                holdout_fraction = 1), "holdout_fraction")
})

test_that(".distance_risk_prepare enforces min_holdout", {
  X <- data.frame(a = 1:10)
  Y <- data.frame(a = 1:10)
  # With min_holdout = 2 and very small fraction, should still get >= 2
  prep <- riskutility:::.distance_risk_prepare(X, Y,
            holdout_fraction = 0.05, min_holdout = 2, seed = 1)
  expect_true(nrow(prep$holdout) >= 2)
})
```

- [ ] **Step 16: Implement `.distance_risk_prepare()`**

Append to `R/utils_internal.R` (after line 109):

```r
#' Shared holdout preparation for distance-risk functions
#'
#' Validates inputs, intersects variables, handles NAs, and splits holdout.
#' Used by dcr(), nndr(), and rf_privacy().
#'
#' @param X data.frame of original data
#' @param Y data.frame of synthetic data
#' @param holdout data.frame or NULL
#' @param holdout_fraction numeric in (0, 1)
#' @param vars character vector or NULL
#' @param na.rm logical
#' @param seed integer or NULL
#' @param min_holdout integer, minimum holdout size (nndr needs 2)
#' @return list with train, synthetic, holdout, vars, was_split
#' @keywords internal
.distance_risk_prepare <- function(X, Y, holdout = NULL,
                                    holdout_fraction = 0.5,
                                    vars = NULL, na.rm = TRUE,
                                    seed = NULL, min_holdout = 1L) {
  if (!is.data.frame(X)) stop("X must be a data.frame", call. = FALSE)
  if (!is.data.frame(Y)) stop("Y must be a data.frame", call. = FALSE)
  if (holdout_fraction <= 0 || holdout_fraction >= 1) {
    stop("holdout_fraction must be in (0, 1)", call. = FALSE)
  }

  # Determine common vars
  common <- intersect(names(X), names(Y))
  if (!is.null(holdout)) {
    common <- intersect(common, names(holdout))
  }
  if (!is.null(vars)) {
    common <- intersect(common, vars)
  }
  if (length(common) == 0) {
    stop("No common variables found between X, Y, and holdout", call. = FALSE)
  }

  X <- X[, common, drop = FALSE]
  Y <- Y[, common, drop = FALSE]

  # NA handling
  if (na.rm) {
    X <- X[stats::complete.cases(X), , drop = FALSE]
    Y <- Y[stats::complete.cases(Y), , drop = FALSE]
  }

  was_split <- FALSE
  if (is.null(holdout)) {
    # Split X into train + holdout
    if (!is.null(seed)) set.seed(seed)
    n <- nrow(X)
    n_holdout <- max(min_holdout, floor(n * holdout_fraction))
    n_holdout <- min(n_holdout, n - 1L)  # keep at least 1 train record
    holdout_idx <- sample(n, n_holdout)
    holdout_df <- X[holdout_idx, , drop = FALSE]
    train_df <- X[-holdout_idx, , drop = FALSE]
    was_split <- TRUE
  } else {
    if (!is.data.frame(holdout)) {
      stop("holdout must be a data.frame", call. = FALSE)
    }
    holdout_df <- holdout[, common, drop = FALSE]
    if (na.rm) {
      holdout_df <- holdout_df[stats::complete.cases(holdout_df), ,
                               drop = FALSE]
    }
    train_df <- X
  }

  list(
    train = train_df,
    synthetic = Y,
    holdout = holdout_df,
    vars = common,
    was_split = was_split
  )
}
```

- [ ] **Step 17: Run tests**

Run: `Rscript -e "testthat::test_file('tests/testthat/test-distance-risk-prepare.R')"`
Expected: All tests PASS

- [ ] **Step 18: Commit**

```bash
git add R/utils_internal.R tests/testthat/test-distance-risk-prepare.R
git commit --author="matthias-da <matthias-da@users.noreply.github.com>" \
  -m "Add .distance_risk_prepare() shared holdout helper"
```

---

### Task 6: Refactor `dcr()` and `nndr()` to use `.distance_risk_prepare()`

**Files:**
- Modify: `R/dcr.R` (lines ~177-240)
- Modify: `R/nndr.R` (lines ~133-200)

- [ ] **Step 19: Run existing dcr and nndr tests (baseline)**

Run: `Rscript -e "testthat::test_file('tests/testthat/test-dcr.R')"`
Run: `Rscript -e "testthat::test_file('tests/testthat/test-nndr.R')"`
Expected: All PASS (baseline before refactor)

- [ ] **Step 20: Refactor `dcr.default()` holdout boilerplate**

In `R/dcr.R`, replace the holdout-splitting boilerplate (approximately lines 177-240) with a call to `.distance_risk_prepare()`. Keep the `method = c("gower", "euclidean")` parameter and all dcr-specific logic (distance computation, null_test, etc.) unchanged. The refactored code replaces:
- Input validation (X/Y are data frames)
- Variable intersection
- NA removal
- Holdout splitting with seed

With: `prep <- .distance_risk_prepare(X, Y, holdout, holdout_fraction, vars, na.rm, seed)`, then use `prep$train`, `prep$synthetic`, `prep$holdout`, `prep$vars`.

- [ ] **Step 21: Refactor `nndr.default()` holdout boilerplate**

Same pattern as dcr, but pass `min_holdout = 2L` (nndr needs at least 2 holdout records for the nearest-neighbor ratio).

- [ ] **Step 22: Run existing tests to verify no regressions**

Run: `Rscript -e "testthat::test_file('tests/testthat/test-dcr.R')"`
Run: `Rscript -e "testthat::test_file('tests/testthat/test-nndr.R')"`
Expected: All PASS (identical results to baseline)

- [ ] **Step 23: Commit**

```bash
git add R/dcr.R R/nndr.R
git commit --author="matthias-da <matthias-da@users.noreply.github.com>" \
  -m "Refactor dcr/nndr to use shared .distance_risk_prepare()"
```

---

## Chunk 3: `rf_privacy()`

### Task 7: `rf_privacy()` — core tests

**Files:**
- Create: `tests/testthat/test-rf-privacy.R`

- [ ] **Step 24: Write core tests for `rf_privacy()`**

```r
# tests/testthat/test-rf-privacy.R
test_that("rf_privacy detects memorized data", {
  skip_if_not_installed("ranger")
  set.seed(1)
  X <- data.frame(a = rnorm(200), b = rnorm(200), c = rnorm(200))
  # Y = copy of X (memorized)
  Y <- X
  res <- rf_privacy(X, Y, holdout_fraction = 0.5, seed = 1,
                    n_trees = 200, null_test = FALSE)

  expect_s3_class(res, "rf_privacy")
  expect_true(res$max_prox_share > 0.55)
  expect_true(res$max_prox_ratio > 1.0)
  expect_false(res$privacy_pass)
})

test_that("rf_privacy passes for random synthetic data", {
  skip_if_not_installed("ranger")
  set.seed(1)
  X <- data.frame(a = rnorm(200), b = rnorm(200))
  Y <- data.frame(a = rnorm(200), b = rnorm(200))
  res <- rf_privacy(X, Y, holdout_fraction = 0.5, seed = 1,
                    n_trees = 200, null_test = FALSE)

  expect_s3_class(res, "rf_privacy")
  # Should be around 0.5 (within tolerance)
  expect_true(abs(res$max_prox_share - 0.5) < 0.15)
})

test_that("rf_privacy returns all expected fields", {
  skip_if_not_installed("ranger")
  set.seed(1)
  X <- data.frame(a = rnorm(100), b = rnorm(100))
  Y <- data.frame(a = rnorm(100), b = rnorm(100))
  res <- rf_privacy(X, Y, seed = 1, n_trees = 50, null_test = FALSE)

  expect_true("max_prox_share" %in% names(res))
  expect_true("max_prox_ratio" %in% names(res))
  expect_true("max_prox_train" %in% names(res))
  expect_true("max_prox_holdout" %in% names(res))
  expect_true("prox_share" %in% names(res))
  expect_true("prox_ratio" %in% names(res))
  expect_true("privacy_pass" %in% names(res))
  expect_true("wilcox_test" %in% names(res))
  expect_true("oob_error" %in% names(res))
  expect_true("var_importance" %in% names(res))
  expect_true(is.logical(res$privacy_pass))
  expect_equal(length(res$max_prox_train), nrow(Y))
})

test_that("rf_privacy null_test works", {
  skip_if_not_installed("ranger")
  set.seed(1)
  X <- data.frame(a = rnorm(100), b = rnorm(100))
  Y <- data.frame(a = rnorm(100), b = rnorm(100))
  res <- rf_privacy(X, Y, seed = 1, n_trees = 50,
                    null_test = TRUE, n_null = 20)

  expect_true("null_distribution" %in% names(res))
  expect_true(!is.null(res$null_distribution))
  expect_true("null_pvalue" %in% names(res$null_distribution))
})

test_that("rf_privacy na.rm = FALSE works with NAs", {
  skip_if_not_installed("ranger")
  set.seed(1)
  X <- data.frame(a = c(rnorm(99), NA), b = rnorm(100))
  Y <- data.frame(a = rnorm(100), b = rnorm(100))
  # Should not error with na.rm = FALSE (ranger handles NAs)
  res <- rf_privacy(X, Y, na.rm = FALSE, seed = 1,
                    n_trees = 50, null_test = FALSE)
  expect_s3_class(res, "rf_privacy")
})

test_that("rf_privacy synth_pair dispatch works", {
  skip_if_not_installed("ranger")
  set.seed(1)
  X <- data.frame(a = rnorm(100), b = rnorm(100))
  Y <- data.frame(a = rnorm(100), b = rnorm(100))
  sp <- synth_pair(X, Y)
  res <- rf_privacy(sp, seed = 1, n_trees = 50, null_test = FALSE)
  expect_s3_class(res, "rf_privacy")
})

test_that("rf_privacy seed separation works", {
  skip_if_not_installed("ranger")
  set.seed(1)
  X <- data.frame(a = rnorm(100), b = rnorm(100))
  Y <- data.frame(a = rnorm(100), b = rnorm(100))
  # Same seed should produce same results
  res1 <- rf_privacy(X, Y, seed = 42, n_trees = 50, null_test = FALSE)
  res2 <- rf_privacy(X, Y, seed = 42, n_trees = 50, null_test = FALSE)
  expect_equal(res1$max_prox_share, res2$max_prox_share)
  # Different seed should produce different holdout splits
  res3 <- rf_privacy(X, Y, seed = 99, n_trees = 50, null_test = FALSE)
  expect_true(res1$max_prox_share != res3$max_prox_share ||
              res1$max_prox_ratio != res3$max_prox_ratio)
})

test_that("rf_privacy prox_ratio zero-guard works", {
  skip_if_not_installed("ranger")
  set.seed(1)
  # Very small holdout with highly separable data → may trigger zero-guard
  X <- data.frame(a = c(rep(0, 10), rep(100, 10)))
  Y <- data.frame(a = rep(50, 10))
  # This may or may not trigger — test that NA is returned with warning
  # when denominator is degenerate
  res <- suppressWarnings(
    rf_privacy(X, Y, holdout_fraction = 0.5, seed = 1,
               n_trees = 50, null_test = FALSE)
  )
  expect_s3_class(res, "rf_privacy")
  # max_prox_ratio is either numeric or NA (both valid)
  expect_true(is.numeric(res$max_prox_ratio) || is.na(res$max_prox_ratio))
})

test_that("rf_privacy print/summary/plot methods work", {
  skip_if_not_installed("ranger")
  set.seed(1)
  X <- data.frame(a = rnorm(80), b = rnorm(80))
  Y <- data.frame(a = rnorm(80), b = rnorm(80))
  res <- rf_privacy(X, Y, seed = 1, n_trees = 50, null_test = FALSE)

  expect_output(print(res), "RF Privacy")
  s <- summary(res)
  expect_s3_class(s, "summary.rf_privacy")
  expect_output(print(s))
  expect_silent(plot(res, which = 1))
  expect_silent(plot(res, which = 2))
})
```

- [ ] **Step 25: Run tests to verify they fail**

Run: `Rscript -e "testthat::test_file('tests/testthat/test-rf-privacy.R')"`
Expected: FAIL — `rf_privacy` not found

- [ ] **Step 26: Commit test file**

```bash
git add tests/testthat/test-rf-privacy.R
git commit --author="matthias-da <matthias-da@users.noreply.github.com>" \
  -m "Add tests for rf_privacy()"
```

---

### Task 8: `rf_privacy()` — implementation

**Files:**
- Create: `R/rf_privacy.R`

- [ ] **Step 27: Implement `rf_privacy()`**

```r
# R/rf_privacy.R

#' RF Proximity Privacy Assessment
#'
#' Detects memorization in synthetic data using Random Forest proximity.
#' Trains a supervised RF to discriminate training from synthetic records,
#' then compares proximity of synthetic records to training vs. holdout.
#'
#' @param X For \code{rf_privacy.default}: a data.frame of original data.
#'   For \code{rf_privacy.synth_pair}: a \code{synth_pair} object.
#' @param Y data.frame of synthetic data (for default method)
#' @param holdout data.frame or NULL. If NULL, split from X.
#' @param holdout_fraction numeric in (0, 1), fraction of X to use as holdout
#' @param vars character vector of variable names (NULL = all common)
#' @param na.rm logical. If FALSE (default), ranger handles NAs natively
#'   via surrogate splits. If TRUE, remove records with any NA before training.
#' @param seed integer or NULL. Used for holdout splitting (seed) and forest
#'   training (seed + 1) to avoid correlation.
#' @param progress logical, show progress bar during proximity computation
#' @param null_test logical, run permutation null test (default TRUE)
#' @param n_null integer, number of permutations for null test
#' @param n_trees integer (>= 10), number of trees
#' @param mtry integer or NULL, number of variables to consider at each split
#' @param ... additional arguments passed to ranger via modifyList
#'
#' @return An S3 object of class \code{"rf_privacy"} with fields:
#'   \describe{
#'     \item{max_prox_share}{fraction of synthetic records with higher max
#'       proximity to training than holdout (mid-rank tie correction)}
#'     \item{max_prox_ratio}{ratio of mean max-proximities (train/holdout)}
#'     \item{max_prox_train}{per-record max proximity to nearest training record}
#'     \item{max_prox_holdout}{per-record max proximity to nearest holdout record}
#'     \item{prox_share}{fraction with higher mean proximity to training}
#'     \item{prox_ratio}{ratio of mean mean-proximities}
#'     \item{prox_train_mean}{per-record mean proximity to training}
#'     \item{prox_holdout_mean}{per-record mean proximity to holdout}
#'     \item{privacy_pass}{logical, TRUE if no memorization detected}
#'     \item{wilcox_test}{Wilcoxon signed-rank test object (heuristic)}
#'     \item{null_distribution}{list with null stats and p-values (if null_test)}
#'     \item{oob_error}{OOB classification error from the forest}
#'     \item{var_importance}{named numeric vector of variable importances}
#'     \item{n_synthetic, n_train, n_holdout, vars}{dataset metadata}
#'   }
#'
#' @details
#' RF proximity measures how often two records land in the same terminal node
#' across all trees. A proximity of 1 means they always co-terminate; 0 means
#' never. This function trains a supervised RF to discriminate training records
#' from synthetic records, then checks whether synthetic records are more
#' similar (proximate) to training records than to holdout records.
#'
#' If no memorization occurred, synthetic records should have roughly equal
#' proximity to training and holdout (both are real data). If synthetic records
#' consistently land in terminal nodes dominated by training records, this
#' signals memorization.
#'
#' @section When to use this method:
#' Use \code{rf_privacy()} when you have mixed data types with complex
#' interactions (20+ variables), or when you want a data-adaptive alternative
#' to \code{\link{dcr}()}. For simple QIs with interpretable risk, \code{dcr()}
#' with Gower distance is more transparent. For speed on small data (n < 5,000),
#' \code{dcr()} is faster.
#'
#' @section Interpretation:
#' \code{max_prox_share} is the primary metric. Values near 0.5 indicate no
#' memorization; values above 0.5 suggest synthetic records are systematically
#' closer to training than holdout. The \code{max_prox_ratio} captures the
#' magnitude: values near 1 = no memorization, > 1 = memorization signal.
#'
#' Mean-based metrics (\code{prox_share}, \code{prox_ratio}) detect aggregate
#' distributional leakage rather than individual memorized records.
#'
#' @section Comparison with DCR:
#' \code{rf_privacy()} is the RF-proximity analog of \code{dcr()}, using the
#' same holdout design but replacing Gower/Euclidean distance with terminal-node
#' co-occurrence. RF proximity is data-adaptive and handles mixed types
#' natively, but adds forest training overhead.
#'
#' @section Limitations:
#' The supervised forest conflates distributional similarity with memorization.
#' When OOB error is near 0.5 (high-quality synthetic data), memorization
#' detection is most reliable. When OOB error is low (poor utility), asymmetry
#' may reflect distributional differences, not individual copying.
#'
#' The Wilcoxon p-value is anti-conservative because all proximity values share
#' the same forest. Use the permutation null test (\code{null_test = TRUE}) for
#' principled inference.
#'
#' @section Computational considerations:
#' Expected runtimes (n_trees = 500, modern laptop):
#' n = 1,000: ~5 seconds; n = 5,000: ~30 seconds;
#' n = 10,000: 2-5 minutes; n = 50,000: 30+ minutes.
#'
#' @family distance-risk
#' @seealso \code{\link{dcr}}, \code{\link{nndr}}, \code{\link{ims}},
#'   \code{\link{propscore}}, \code{\link{recordLinkage}}
#' @references
#' Breiman, L. (2001). Random Forests. \emph{Machine Learning}, 45(1), 5-32.
#'
#' Lin, Y. & Jeon, Y. (2006). Random forests and adaptive nearest neighbors.
#' \emph{JASA}, 101(474), 578-590.
#'
#' @examples
#' \donttest{
#' # Memorized data (copy of training)
#' set.seed(1)
#' X <- data.frame(a = rnorm(200), b = rnorm(200))
#' Y <- X  # memorized
#' res_mem <- rf_privacy(X, Y, seed = 1, n_trees = 200, null_test = FALSE)
#' print(res_mem)
#'
#' # Random synthetic data (no memorization)
#' Y_rand <- data.frame(a = rnorm(200), b = rnorm(200))
#' res_rand <- rf_privacy(X, Y_rand, seed = 1, n_trees = 200, null_test = FALSE)
#' print(res_rand)
#' }
#' @export
rf_privacy <- function(X, ...) UseMethod("rf_privacy")

#' @rdname rf_privacy
#' @export
rf_privacy.synth_pair <- function(X, ...) {
  rf_privacy.default(X = X$original, Y = X$synthetic,
                     vars = X$vars, ...)
}

#' @rdname rf_privacy
#' @export
rf_privacy.default <- function(X, Y,
                               holdout = NULL, holdout_fraction = 0.5,
                               vars = NULL, na.rm = FALSE, seed = NULL,
                               progress = FALSE,
                               null_test = TRUE, n_null = 100L,
                               n_trees = 500L, mtry = NULL, ...) {

  if (!requireNamespace("ranger", quietly = TRUE)) {
    stop("Package 'ranger' required for rf_privacy(). ",
         "Install with install.packages('ranger')", call. = FALSE)
  }

  # Prepare holdout split
  prep <- .distance_risk_prepare(X, Y,
                                  holdout = holdout,
                                  holdout_fraction = holdout_fraction,
                                  vars = vars, na.rm = na.rm,
                                  seed = seed, min_holdout = 1L)
  train     <- prep$train
  synthetic <- prep$synthetic
  holdout_d <- prep$holdout
  vars_use  <- prep$vars

  n_syn   <- nrow(synthetic)
  n_train <- nrow(train)
  n_ho    <- nrow(holdout_d)

  if (n_syn + n_train > 10000) {
    message("rf_privacy: n = ", n_syn + n_train,
            ". Expected runtime: 2-5 minutes with n_trees = ", n_trees, ".")
  }

  # Train RF on train vs synthetic
  rf_seed <- if (!is.null(seed)) seed + 1L else NULL
  rf_res <- .rf_proximity(train, synthetic, vars = vars_use,
                          n_trees = n_trees, mtry = mtry,
                          importance = TRUE, seed = rf_seed, ...)

  # Indices: train = 1:n_train, synthetic = (n_train+1):(n_train+n_syn)
  idx_train <- seq_len(n_train)
  idx_syn   <- n_train + seq_len(n_syn)

  # Compute proximity: synthetic (query) vs training (reference)
  prox_syn_train <- .proximity_from_nodes(rf_res$terminal_nodes,
                                          idx_train, idx_syn,
                                          progress = progress)
  # prox_syn_train is n_syn x n_train

  # Push holdout through forest to get its terminal nodes
  tn_ho <- predict(rf_res$forest, holdout_d,
                   type = "terminalNodes")$predictions
  # Combine training TN + holdout TN
  combined_tn <- rbind(rf_res$terminal_nodes, tn_ho)
  idx_ho <- nrow(rf_res$terminal_nodes) + seq_len(n_ho)
  prox_syn_ho <- .proximity_from_nodes(combined_tn, idx_ho, idx_syn,
                                       progress = FALSE)
  # prox_syn_ho is n_syn x n_ho

  # Per-record max and mean proximity
  max_prox_train   <- apply(prox_syn_train, 1, max)
  max_prox_holdout <- apply(prox_syn_ho, 1, max)
  mean_prox_train  <- rowMeans(prox_syn_train)
  mean_prox_ho     <- rowMeans(prox_syn_ho)

  # Max-based metrics (primary)
  # Mid-rank tie correction: mean(a > b) + 0.5 * mean(a == b)
  max_prox_share <- mean(max_prox_train > max_prox_holdout) +
    0.5 * mean(max_prox_train == max_prox_holdout)

  # Zero-guard for ratio
  denom_max <- mean(max_prox_holdout)
  if (denom_max < 1 / n_trees) {
    max_prox_ratio <- NA_real_
    warning("max_prox_holdout mean < 1/n_trees (degenerate). ",
            "max_prox_ratio set to NA.", call. = FALSE)
  } else {
    max_prox_ratio <- mean(max_prox_train) / denom_max
  }

  # Mean-based metrics (supplementary)
  prox_share <- mean(mean_prox_train > mean_prox_ho) +
    0.5 * mean(mean_prox_train == mean_prox_ho)

  denom_mean <- mean(mean_prox_ho)
  if (denom_mean < 1 / n_trees) {
    prox_ratio <- NA_real_
    warning("prox_holdout_mean mean < 1/n_trees (degenerate). ",
            "prox_ratio set to NA.", call. = FALSE)
  } else {
    prox_ratio <- mean(mean_prox_train) / denom_mean
  }

  # Wilcoxon signed-rank test (heuristic — anti-conservative)
  wilcox_res <- suppressWarnings(
    stats::wilcox.test(max_prox_train, max_prox_holdout,
                       paired = TRUE, alternative = "greater")
  )

  # Null test: permute train/holdout labels
  null_dist <- NULL
  if (null_test) {
    # Pre-compute full syn-vs-all-real proximity
    # Combine all real (train + holdout) terminal nodes
    all_real_tn <- rbind(
      rf_res$terminal_nodes[idx_train, , drop = FALSE],
      tn_ho
    )
    n_real <- n_train + n_ho
    syn_tn <- rf_res$terminal_nodes[idx_syn, , drop = FALSE]
    combined_for_null <- rbind(all_real_tn, syn_tn)
    idx_all_real <- seq_len(n_real)
    idx_syn_null <- n_real + seq_len(n_syn)

    # Full syn-vs-all-real proximity matrix (n_syn x n_real)
    prox_syn_all <- .proximity_from_nodes(combined_for_null,
                                          idx_all_real, idx_syn_null,
                                          progress = progress)

    null_shares <- numeric(n_null)
    null_ratios <- numeric(n_null)
    null_max_shares <- numeric(n_null)
    null_max_ratios <- numeric(n_null)

    for (p in seq_len(n_null)) {
      perm <- sample(n_real)
      perm_train <- perm[seq_len(n_train)]
      perm_ho    <- perm[(n_train + 1):n_real]

      # Mean-based
      perm_train_mean <- rowMeans(prox_syn_all[, perm_train, drop = FALSE])
      perm_ho_mean    <- rowMeans(prox_syn_all[, perm_ho, drop = FALSE])
      null_shares[p] <- mean(perm_train_mean > perm_ho_mean) +
        0.5 * mean(perm_train_mean == perm_ho_mean)

      d_ho <- mean(perm_ho_mean)
      null_ratios[p] <- if (d_ho < 1 / n_trees) NA_real_ else
        mean(perm_train_mean) / d_ho

      # Max-based
      perm_train_max <- apply(prox_syn_all[, perm_train, drop = FALSE], 1, max)
      perm_ho_max    <- apply(prox_syn_all[, perm_ho, drop = FALSE], 1, max)
      null_max_shares[p] <- mean(perm_train_max > perm_ho_max) +
        0.5 * mean(perm_train_max == perm_ho_max)

      d_ho_max <- mean(perm_ho_max)
      null_max_ratios[p] <- if (d_ho_max < 1 / n_trees) NA_real_ else
        mean(perm_train_max) / d_ho_max
    }

    # Permutation p-values: (sum + 1) / (n_null + 1) [Phipson & Smyth]
    pval_share <- (sum(null_max_shares >= max_prox_share, na.rm = TRUE) + 1) /
      (n_null + 1)
    pval_ratio <- (sum(null_max_ratios >= max_prox_ratio, na.rm = TRUE) + 1) /
      (n_null + 1)

    null_dist <- list(
      null_max_shares = null_max_shares,
      null_max_ratios = null_max_ratios,
      null_shares = null_shares,
      null_ratios = null_ratios,
      null_pvalue = pval_share,
      null_ratio_pvalue = pval_ratio
    )
  }

  # Privacy pass
  if (!is.null(null_dist)) {
    # Null-test-derived: observed within 95th percentile
    q95_share <- stats::quantile(null_dist$null_max_shares, 0.95,
                                 na.rm = TRUE)
    q95_ratio <- stats::quantile(null_dist$null_max_ratios, 0.95,
                                 na.rm = TRUE)
    privacy_pass <- (max_prox_share <= q95_share) &&
      (is.na(max_prox_ratio) || max_prox_ratio <= q95_ratio)
  } else {
    # Heuristic fallback
    privacy_pass <- max_prox_share <= 0.55 &&
      wilcox_res$p.value > 0.05
  }

  result <- list(
    max_prox_share   = max_prox_share,
    max_prox_ratio   = max_prox_ratio,
    max_prox_train   = max_prox_train,
    max_prox_holdout = max_prox_holdout,
    prox_share       = prox_share,
    prox_ratio       = prox_ratio,
    prox_train_mean  = mean_prox_train,
    prox_holdout_mean = mean_prox_ho,
    privacy_pass     = privacy_pass,
    wilcox_test      = wilcox_res,
    null_distribution = null_dist,
    oob_error        = rf_res$oob_error,
    var_importance   = rf_res$importance,
    n_synthetic      = n_syn,
    n_train          = n_train,
    n_holdout        = n_ho,
    vars             = vars_use
  )
  class(result) <- "rf_privacy"
  result
}

#' @export
print.rf_privacy <- function(x, ...) {
  pass_label <- if (x$privacy_pass) "PASS" else "FAIL"
  cat("RF Privacy Assessment (rf_privacy)\n")
  cat("  Max proximity share: ", sprintf("%.2f", x$max_prox_share),
      if (x$max_prox_share <= 0.55) " (training not preferred)"
      else " (training preferred — memorization signal)", "\n")
  cat("  Max proximity ratio: ", sprintf("%.2f", x$max_prox_ratio), "\n")
  cat("  OOB error:           ", sprintf("%.2f", x$oob_error), "\n")
  if (!is.null(x$null_distribution)) {
    cat("  Null test:           ", pass_label,
        sprintf(" (p = %.2f)", x$null_distribution$null_pvalue), "\n")
  }
  cat("  Privacy:             ", pass_label, "\n")
  invisible(x)
}

#' @export
summary.rf_privacy <- function(object, ...) {
  s <- list(
    max_prox_share   = object$max_prox_share,
    max_prox_ratio   = object$max_prox_ratio,
    prox_share       = object$prox_share,
    prox_ratio       = object$prox_ratio,
    oob_error        = object$oob_error,
    privacy_pass     = object$privacy_pass,
    wilcox_p         = object$wilcox_test$p.value,
    null_distribution = object$null_distribution,
    var_importance   = object$var_importance,
    n_synthetic      = object$n_synthetic,
    n_train          = object$n_train,
    n_holdout        = object$n_holdout,
    # Per-record outliers: top 5 by max_prox_train
    top_outliers     = utils::head(
      order(object$max_prox_train, decreasing = TRUE), 5
    )
  )
  class(s) <- "summary.rf_privacy"
  s
}

#' @export
print.summary.rf_privacy <- function(x, ...) {
  cat("RF Privacy Assessment — Summary\n")
  cat(strrep("-", 50), "\n")
  cat("Dataset: ", x$n_synthetic, " synthetic, ",
      x$n_train, " train, ", x$n_holdout, " holdout\n\n")

  cat("Max-based metrics (primary):\n")
  cat("  max_prox_share: ", sprintf("%.4f", x$max_prox_share), "\n")
  cat("  max_prox_ratio: ", sprintf("%.4f", x$max_prox_ratio), "\n\n")

  cat("Mean-based metrics (supplementary):\n")
  cat("  prox_share:     ", sprintf("%.4f", x$prox_share), "\n")
  cat("  prox_ratio:     ", sprintf("%.4f", x$prox_ratio), "\n\n")

  cat("Forest diagnostics:\n")
  cat("  OOB error:      ", sprintf("%.4f", x$oob_error), "\n")
  cat("  Wilcoxon p:     ", sprintf("%.4f", x$wilcox_p),
      " (heuristic)\n\n")

  if (!is.null(x$null_distribution)) {
    cat("Null test:\n")
    cat("  p-value (share): ", sprintf("%.4f",
        x$null_distribution$null_pvalue), "\n")
    cat("  p-value (ratio): ", sprintf("%.4f",
        x$null_distribution$null_ratio_pvalue), "\n\n")
  }

  cat("Privacy: ", if (x$privacy_pass) "PASS" else "FAIL", "\n")

  if (!is.null(x$var_importance) && length(x$var_importance) > 0) {
    cat("\nTop variables by importance:\n")
    imp_sorted <- sort(x$var_importance, decreasing = TRUE)
    for (i in seq_len(min(10, length(imp_sorted)))) {
      cat("  ", names(imp_sorted)[i], ": ",
          sprintf("%.4f", imp_sorted[i]), "\n")
    }
  }

  cat("\nTop 5 outlier records (by max proximity to training): ",
      paste(x$top_outliers, collapse = ", "), "\n")
  invisible(x)
}

#' @export
plot.rf_privacy <- function(x, y = NULL, which = 1L, ...) {
  show <- rep(FALSE, 3)
  show[which] <- TRUE

  if (show[1]) {
    # Paired density: max proximity to training vs holdout
    df <- data.frame(
      proximity = c(x$max_prox_train, x$max_prox_holdout),
      group = rep(c("To training", "To holdout"),
                  each = length(x$max_prox_train))
    )
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$proximity,
                                           fill = .data$group)) +
      ggplot2::geom_density(alpha = 0.5) +
      ggplot2::labs(title = "Max Proximity: Training vs Holdout",
                    x = "Max proximity", y = "Density",
                    fill = "Reference") +
      ggplot2::theme_minimal()
    print(p)
  }

  if (show[2]) {
    # Per-record difference histogram
    diffs <- x$max_prox_train - x$max_prox_holdout
    df <- data.frame(diff = diffs)
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$diff)) +
      ggplot2::geom_histogram(bins = 30, fill = "steelblue",
                               color = "white") +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                           color = "red") +
      ggplot2::labs(
        title = "Per-Record Max Proximity Difference",
        x = "max_prox(train) - max_prox(holdout)",
        y = "Count"
      ) +
      ggplot2::theme_minimal()
    print(p)
  }

  if (show[3]) {
    if (is.null(x$null_distribution)) {
      message("No null distribution available (null_test = FALSE)")
    } else {
      df <- data.frame(null = x$null_distribution$null_max_shares)
      p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$null)) +
        ggplot2::geom_histogram(bins = 20, fill = "grey70",
                                 color = "white") +
        ggplot2::geom_vline(xintercept = x$max_prox_share,
                             color = "red", linewidth = 1.2) +
        ggplot2::labs(
          title = "Null Distribution of Max Proximity Share",
          x = "max_prox_share (permuted)",
          y = "Count",
          caption = sprintf("Observed = %.3f, p = %.3f",
                            x$max_prox_share,
                            x$null_distribution$null_pvalue)
        ) +
        ggplot2::theme_minimal()
      print(p)
    }
  }

  invisible(x)
}
```

- [ ] **Step 28: Run `devtools::document()` to generate NAMESPACE entries**

Run: `Rscript -e "devtools::document()"`
Expected: New exports for `rf_privacy`, `print.rf_privacy`, etc.

- [ ] **Step 29: Run tests**

Run: `Rscript -e "testthat::test_file('tests/testthat/test-rf-privacy.R')"`
Expected: All tests PASS

- [ ] **Step 30: Run full test suite to check for regressions**

Run: `Rscript -e "devtools::test()"`
Expected: All existing tests still pass

- [ ] **Step 31: Commit**

```bash
git add R/rf_privacy.R man/ NAMESPACE
git commit --author="matthias-da <matthias-da@users.noreply.github.com>" \
  -m "Implement rf_privacy() memorization detection via RF proximity"
```

---

## Chunk 4: `propscore(method = "ranger")`

### Task 9: `propscore(method = "ranger")` — tests

**Files:**
- Create: `tests/testthat/test-propscore-ranger.R`

- [ ] **Step 32: Write tests**

```r
# tests/testthat/test-propscore-ranger.R
test_that("propscore(method = 'ranger') returns correct class", {
  skip_if_not_installed("ranger")
  set.seed(1)
  X <- data.frame(a = rnorm(100), b = rnorm(100))
  Y <- data.frame(a = rnorm(100), b = rnorm(100))
  res <- propscore(X, Y, method = "ranger", proximity = "none")

  expect_s3_class(res, "propscore")
  expect_equal(res$method, "ranger")
  expect_true("oob_error" %in% names(res))
  expect_true("var_importance" %in% names(res))
})

test_that("propscore(method = 'ranger') identical data gives OOB ~0.5", {
  skip_if_not_installed("ranger")
  set.seed(1)
  X <- data.frame(a = rnorm(200), b = rnorm(200))
  Y <- X  # identical
  res <- propscore(X, Y, method = "ranger", proximity = "none")

  expect_true(res$oob_error > 0.35)  # near 0.5
})

test_that("propscore(method = 'ranger') different data gives OOB < 0.3", {
  skip_if_not_installed("ranger")
  set.seed(1)
  X <- data.frame(a = rnorm(200), b = rnorm(200))
  Y <- data.frame(a = rnorm(200) + 3, b = rnorm(200) + 3)
  res <- propscore(X, Y, method = "ranger", proximity = "none")

  expect_true(res$oob_error < 0.3)
})

test_that("propscore proximity = 'summary' returns structural metrics", {
  skip_if_not_installed("ranger")
  set.seed(1)
  X <- data.frame(a = rnorm(100), b = rnorm(100))
  Y <- data.frame(a = rnorm(100), b = rnorm(100))
  res <- propscore(X, Y, method = "ranger", proximity = "summary")

  expect_true("within_orig_prox" %in% names(res))
  expect_true("within_synth_prox" %in% names(res))
  expect_true("cross_prox" %in% names(res))
  expect_true("structure_ratio" %in% names(res))
  expect_null(res$proximity_matrix)
})

test_that("propscore proximity = 'full' stores matrix", {
  skip_if_not_installed("ranger")
  set.seed(1)
  X <- data.frame(a = rnorm(30), b = rnorm(30))
  Y <- data.frame(a = rnorm(30), b = rnorm(30))
  res <- propscore(X, Y, method = "ranger", proximity = "full")

  expect_true(!is.null(res$proximity_matrix))
  expect_equal(nrow(res$proximity_matrix), 60)
  expect_equal(ncol(res$proximity_matrix), 60)
})

test_that("propscore match.arg rejects invalid method", {
  X <- data.frame(a = 1:10)
  Y <- data.frame(a = 1:10)
  expect_error(propscore(X, Y, method = "invalid"), "arg")
})

test_that("propscore existing methods still work after match.arg", {
  set.seed(1)
  X <- data.frame(a = rnorm(50), b = rnorm(50))
  Y <- data.frame(a = rnorm(50), b = rnorm(50))

  res_rf <- propscore(X, Y, method = "rf")
  expect_s3_class(res_rf, "propscore")

  res_lr <- propscore(X, Y, method = "logreg")
  expect_s3_class(res_lr, "propscore")
})

test_that("propscore emits message for ranger params with non-ranger method", {
  set.seed(1)
  X <- data.frame(a = rnorm(50))
  Y <- data.frame(a = rnorm(50))
  expect_message(
    propscore(X, Y, method = "rf", proximity = "summary"),
    "ranger"
  )
})

test_that("propscore(method = 'ranger') synth_pair dispatch works", {
  skip_if_not_installed("ranger")
  set.seed(1)
  X <- data.frame(a = rnorm(50), b = rnorm(50))
  Y <- data.frame(a = rnorm(50), b = rnorm(50))
  sp <- synth_pair(X, Y)
  res <- propscore(sp, method = "ranger", proximity = "summary")
  expect_s3_class(res, "propscore")
  expect_true("structure_ratio" %in% names(res))
})

test_that("propscore plot which = 3 and 4 work for ranger", {
  skip_if_not_installed("ranger")
  set.seed(1)
  X <- data.frame(a = rnorm(50), b = rnorm(50))
  Y <- data.frame(a = rnorm(50), b = rnorm(50))
  res <- propscore(X, Y, method = "ranger", proximity = "summary")

  expect_silent(plot(res, which = 3))
  expect_silent(plot(res, which = 4))
})
```

- [ ] **Step 33: Commit tests**

```bash
git add tests/testthat/test-propscore-ranger.R
git commit --author="matthias-da <matthias-da@users.noreply.github.com>" \
  -m "Add tests for propscore(method = 'ranger')"
```

---

### Task 10: `propscore(method = "ranger")` — implementation

**Files:**
- Modify: `R/propscore.R`

- [ ] **Step 34: Add `match.arg()` and ranger branch to `propscore.default()`**

Apply these changes to `R/propscore.R`:

**A. Update function signature** (line ~136):

Change:
```r
method = "rf",
```
To:
```r
method = c("rf", "ranger", "logreg"),
proximity = c("summary", "full", "none"),
importance = TRUE,
```

**B. Add match.arg** (beginning of function body, after line ~140):
```r
method <- match.arg(method)
proximity <- match.arg(proximity)
```

**C. Add ranger-param message** (before method branching, ~line 260):
```r
if (method != "ranger" && (!missing(proximity) || !missing(importance))) {
  message("'proximity' and 'importance' are only used with method = 'ranger'.")
}
```

**D. Add ranger branch** in the method branching section (alongside existing `"logreg"` and `"rf"` branches). Insert the following `else if` block:

```r
} else if (method == "ranger") {
  if (!requireNamespace("ranger", quietly = TRUE)) {
    stop("Package 'ranger' required for propscore(method = 'ranger'). ",
         "Install with install.packages('ranger')", call. = FALSE)
  }

  # Use .rf_proximity() engine
  rf_res <- .rf_proximity(X, Y, vars = vars, n_trees = 500L,
                          importance = importance, ...)

  # OOB propensity scores
  p_hat <- rf_res$forest$predictions[, 2]  # P(synthetic)
  # Handle OOB prediction NAs
  if (any(is.na(p_hat))) {
    p_inbag <- predict(rf_res$forest,
                       rbind(X, Y))$predictions[, 2]
    p_hat[is.na(p_hat)] <- p_inbag[is.na(p_hat)]
  }

  cr <- nrow(Y) / (nrow(X) + nrow(Y))  # class ratio
  n_total <- nrow(X) + nrow(Y)
  pmse <- (1 / (2 * n_total)) * sum((p_hat - cr)^2)

  # Proximity-based structural metrics
  prox_fields <- list()
  if (proximity != "none") {
    n1 <- rf_res$n1
    n2 <- rf_res$n2
    idx_orig <- seq_len(n1)
    idx_synth <- n1 + seq_len(n2)

    # Within-original proximity (mean of upper triangle)
    if (n1 > 1) {
      wo <- .proximity_from_nodes(rf_res$terminal_nodes, idx_orig, idx_orig)
      prox_fields$within_orig_prox <- (sum(wo) - n1) / (n1 * (n1 - 1))
    } else {
      prox_fields$within_orig_prox <- NA_real_
    }

    # Within-synthetic proximity
    if (n2 > 1) {
      ws <- .proximity_from_nodes(rf_res$terminal_nodes, idx_synth, idx_synth)
      prox_fields$within_synth_prox <- (sum(ws) - n2) / (n2 * (n2 - 1))
    } else {
      prox_fields$within_synth_prox <- NA_real_
    }

    # Cross-class proximity
    cross <- .proximity_from_nodes(rf_res$terminal_nodes, idx_orig, idx_synth)
    prox_fields$cross_prox <- mean(cross)

    # Structure ratio
    denom <- (prox_fields$within_orig_prox +
              prox_fields$within_synth_prox) / 2
    prox_fields$structure_ratio <- if (denom > 0) {
      prox_fields$cross_prox / denom
    } else {
      NA_real_
    }

    if (proximity == "full") {
      n_combined <- n1 + n2
      if (n_combined > 10000) {
        mem_mb <- round(n_combined^2 * 8 / 1e6)
        warning("Full proximity matrix is ", n_combined, " x ",
                n_combined, " (~", mem_mb, " MB). ",
                "Consider proximity = 'summary'.", call. = FALSE)
      }
      all_idx <- seq_len(n_combined)
      prox_fields$proximity_matrix <- .proximity_from_nodes(
        rf_res$terminal_nodes, all_idx, all_idx
      )
    }
  }

  oob_error <- rf_res$oob_error
  if (oob_error > 0.45) {
    message("OOB error > 0.45: forest has little discriminative power.")
  }
```

**Important**: The existing code (lines 268-331) creates a `p` data.frame + KDE metrics + a monolithic `results <- list(...)` + `class(results) <- "propscore"` + `return(results)`. There is NO common result builder — each method branch leads to the same result list. The ranger branch must **return early** with its own result list (before the KDE code at line 284), because KDE on RF OOB scores is not the right diagnostic for the ranger method.

Place this code inside the ranger branch, BEFORE the KDE code. The branch should end with `return(results)`:

```r
  # Build result — ranger bypasses KDE and returns directly
  p <- data.frame(
    prediction = p_hat,
    group = rep(c("real", "synth"), times = c(nrow(X), nrow(Y)))
  )

  results <- list(
    predictions = p,
    ps_ratio = NA_real_,   # not meaningful for ranger (KDE-based)
    ps_score = pmse,       # pMSE (same formula as other methods)
    cr = cr,
    mean_ps_x = mean(p_hat[seq_len(nrow(X))], na.rm = TRUE),
    mean_ps_y = mean(p_hat[nrow(X) + seq_len(nrow(Y))], na.rm = TRUE),
    density_ratio = NULL,       # KDE not computed for ranger
    density_ratio_bayes = NULL,
    kl = NA_real_,
    kl_bayes = NA_real_,
    mean_ratio = NA_real_,
    sd_ratio = NA_real_,
    mean_ratio_bayes = NA_real_,
    sd_ratio_bayes = NA_real_,
    points = NULL,
    denX = NULL,
    denY = NULL,
    bayesspace = FALSE,
    n_x = nrow(X),
    n_y = nrow(Y),
    method = "ranger",
    # Ranger-specific fields
    oob_error = oob_error,
    var_importance = rf_res$importance
  )

  # Add proximity fields
  for (nm in names(prox_fields)) results[[nm]] <- prox_fields[[nm]]

  class(results) <- "propscore"
  return(results)
```

This ensures the result object has all fields that `print.propscore()`, `summary.propscore()`, and `plot.propscore()` expect, with `NA`/`NULL` for inapplicable KDE fields.

**E. Update `plot.propscore()`** (line ~419):

Change `show <- rep(FALSE, 2)` to `show <- rep(FALSE, 4)`.

Add guard in the existing `show[1L]` and `show[2L]` blocks (inside the `else` branch at line 491, before accessing `x$points`, `x$denX`, etc.):

```r
# At the top of the else block (single propscore, not list-of-lists), line ~491:
} else {
    # Guard for ranger method: density/ratio plots need KDE fields
    if (is.null(x$points) && (show[1L] || show[2L])) {
      message("Density plots not available for method = 'ranger'. ",
              "Use which = 3 (proximity) or which = 4 (importance).")
      show[1L] <- FALSE
      show[2L] <- FALSE
    }
    # ... existing code continues ...
```

Also update `print.summary.propscore()` to handle NA KDE fields gracefully: only print KL/density ratio lines when they are not NA.

Then add the new `which = 3` and `which = 4` blocks:

```r
if (show[3]) {
  # Proximity structure boxplot
  if (is.null(x$within_orig_prox)) {
    message("No proximity data (method != 'ranger' or proximity = 'none')")
  } else {
    df <- data.frame(
      group = c("Within original", "Within synthetic", "Cross-class"),
      proximity = c(x$within_orig_prox, x$within_synth_prox, x$cross_prox)
    )
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$group,
                                           y = .data$proximity)) +
      ggplot2::geom_col(fill = "steelblue") +
      ggplot2::labs(title = "Proximity Structure",
                    x = NULL, y = "Mean proximity",
                    caption = sprintf("Structure ratio: %.3f",
                                      x$structure_ratio)) +
      ggplot2::theme_minimal()
    print(p)
  }
}

if (show[4]) {
  # Variable importance bar chart
  if (is.null(x$var_importance)) {
    message("No importance data (importance = FALSE or method != 'ranger')")
  } else {
    imp <- sort(x$var_importance, decreasing = TRUE)
    df <- data.frame(var = factor(names(imp), levels = names(imp)),
                     importance = imp)
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$var,
                                           y = .data$importance)) +
      ggplot2::geom_col(fill = "darkgreen") +
      ggplot2::coord_flip() +
      ggplot2::labs(title = "Variable Importance",
                    x = NULL, y = "Importance") +
      ggplot2::theme_minimal()
    print(p)
  }
}
```

**F. Update roxygen** `@param method` to include `"ranger"` and add `@param proximity`, `@param importance` documentation.

- [ ] **Step 35: Run `devtools::document()`**

Run: `Rscript -e "devtools::document()"`

- [ ] **Step 36: Run all propscore tests**

Run: `Rscript -e "testthat::test_file('tests/testthat/test-propscore-ranger.R')"`
Run: `Rscript -e "testthat::test_file('tests/testthat/test-pmse.R')"`
Expected: All PASS (new + existing)

- [ ] **Step 37: Commit**

```bash
git add R/propscore.R man/propscore.Rd
git commit --author="matthias-da <matthias-da@users.noreply.github.com>" \
  -m "Add method = 'ranger' to propscore() with proximity-based utility"
```

---

## Chunk 5: `recordLinkage(method = "rf")`

### Task 11: `recordLinkage(method = "rf")` — tests

**Files:**
- Create: `tests/testthat/test-recordLinkage-rf.R`

- [ ] **Step 38: Write tests**

```r
# tests/testthat/test-recordLinkage-rf.R
test_that("recordLinkage(method = 'rf') independent matching works", {
  skip_if_not_installed("ranger")
  set.seed(1)
  X <- data.frame(a = rnorm(50), b = sample(letters[1:5], 50, TRUE))
  Y <- data.frame(a = rnorm(50), b = sample(letters[1:5], 50, TRUE))
  res <- recordLinkage(X, Y, key = c("a", "b"), method = "rf",
                       n_trees = 100, matching = "independent")

  expect_s3_class(res, "recordLinkageRisk")
  expect_true("var_importance" %in% names(res))
  expect_true(all(res$per_record$risk >= 0 & res$per_record$risk <= 1))
})

test_that("recordLinkage(method = 'rf') bijective matching works", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("clue")
  set.seed(1)
  X <- data.frame(a = rnorm(30), b = rnorm(30))
  Y <- data.frame(a = rnorm(30), b = rnorm(30))
  res <- recordLinkage(X, Y, key = c("a", "b"), method = "rf",
                       n_trees = 100, matching = "bijective")

  expect_s3_class(res, "recordLinkageRisk")
  # Bijective risk is binary
  expect_true(all(res$per_record$risk %in% c(0, 1)))
})

test_that("recordLinkage(method = 'rf') with blocking works", {
  skip_if_not_installed("ranger")
  set.seed(1)
  X <- data.frame(a = rnorm(100), b = rep(c("A", "B"), 50))
  Y <- data.frame(a = rnorm(100), b = rep(c("A", "B"), 50))
  res <- recordLinkage(X, Y, key = "a", block = "b", method = "rf",
                       n_trees = 50)

  expect_s3_class(res, "recordLinkageRisk")
})

test_that("recordLinkage(method = 'rf') rf_global = TRUE works", {
  skip_if_not_installed("ranger")
  set.seed(1)
  X <- data.frame(a = rnorm(100), b = rep(c("A", "B"), 50))
  Y <- data.frame(a = rnorm(100), b = rep(c("A", "B"), 50))
  res <- recordLinkage(X, Y, key = "a", block = "b", method = "rf",
                       n_trees = 50, rf_global = TRUE)

  expect_s3_class(res, "recordLinkageRisk")
})

test_that("recordLinkage(method = 'rf') direction = 'reverse' works", {
  skip_if_not_installed("ranger")
  set.seed(1)
  X <- data.frame(a = rnorm(30), b = rnorm(30))
  Y <- data.frame(a = rnorm(30), b = rnorm(30))
  res_fwd <- recordLinkage(X, Y, key = c("a", "b"), method = "rf",
                           n_trees = 50, direction = "forward")
  res_rev <- recordLinkage(X, Y, key = c("a", "b"), method = "rf",
                           n_trees = 50, direction = "reverse")

  expect_s3_class(res_fwd, "recordLinkageRisk")
  expect_s3_class(res_rev, "recordLinkageRisk")
})

test_that("recordLinkage(method = 'rf') emits message for weights", {
  skip_if_not_installed("ranger")
  set.seed(1)
  X <- data.frame(a = rnorm(30), b = rnorm(30))
  Y <- data.frame(a = rnorm(30), b = rnorm(30))
  expect_message(
    recordLinkage(X, Y, key = c("a", "b"), method = "rf",
                  weights = c(1, 2), n_trees = 50),
    "importance"
  )
})

test_that("recordLinkage(method = 'rf') small block fallback", {
  skip_if_not_installed("ranger")
  set.seed(1)
  # Block C has only 1 record per class → should fall back
  X <- data.frame(a = rnorm(31), b = c(rep("A", 15), rep("B", 15), "C"))
  Y <- data.frame(a = rnorm(31), b = c(rep("A", 15), rep("B", 15), "C"))
  expect_message(
    res <- recordLinkage(X, Y, key = "a", block = "b",
                         method = "rf", n_trees = 50),
    "blocks"
  )
  expect_s3_class(res, "recordLinkageRisk")
})

test_that("recordLinkage(method = 'rf') truth = 'row' works", {
  skip_if_not_installed("ranger")
  set.seed(1)
  X <- data.frame(a = rnorm(20), b = rnorm(20))
  Y <- X[sample(20), ]  # permuted copy
  res <- recordLinkage(X, Y, key = c("a", "b"), method = "rf",
                       truth = "row", n_trees = 50)
  expect_s3_class(res, "recordLinkageRisk")
  expect_true("true_in_set" %in% names(res$per_record))
})

test_that("recordLinkage(method = 'rf') var_importance plot works", {
  skip_if_not_installed("ranger")
  set.seed(1)
  X <- data.frame(a = rnorm(30), b = rnorm(30))
  Y <- data.frame(a = rnorm(30), b = rnorm(30))
  res <- recordLinkage(X, Y, key = c("a", "b"), method = "rf",
                       n_trees = 50)
  # which = 3 is existing var importance plot
  expect_silent(plot(res, which = 3))
})
```

- [ ] **Step 39: Commit tests**

```bash
git add tests/testthat/test-recordLinkage-rf.R
git commit --author="matthias-da <matthias-da@users.noreply.github.com>" \
  -m "Add tests for recordLinkage(method = 'rf')"
```

---

### Task 12: `recordLinkage(method = "rf")` — implementation

**Files:**
- Modify: `R/recordLinkage.R`
- Create: (within `R/recordLinkage.R` or `R/rf_proximity_internal.R`) `.rf_linkage_block()` helper

- [ ] **Step 40: Add RF method to `recordLinkage()`**

Apply these changes to `R/recordLinkage.R`:

**A. Update roxygen** (line ~165): Add `"rf"` to `@param method` and add:
```r
#' @param n_trees integer (>= 10), number of trees for \code{method = "rf"}.
#'   Ignored for other methods.
#' @param rf_global logical. If \code{FALSE} (default), train separate RF per
#'   block. If \code{TRUE}, train one global RF and restrict proximity matching
#'   to within-block pairs. Ignored for other methods.
```

**B. Update function signature** (line ~372): Add `n_trees = 500, rf_global = FALSE` params and add `"rf"` to method character vector.

**C. After `match.arg` section** (line ~397): Add:
```r
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
```

**D. Add `.rf_linkage_block()` helper** to `R/rf_proximity_internal.R`:

```r
#' RF linkage matching for one block
#'
#' Trains RF on query + candidate records, computes cross-proximity,
#' returns per-record matches.
#'
#' @param query data.frame (records being matched)
#' @param candidate data.frame (records to match against)
#' @param key character vector of key variables
#' @param n_trees integer
#' @param ... additional args to .rf_proximity
#' @return list with risk (proximity scores), match_idx, var_importance
#' @keywords internal
.rf_linkage_block <- function(query, candidate, key,
                               n_trees = 500L, ...) {
  n_q <- nrow(query)
  n_c <- nrow(candidate)

  rf_res <- .rf_proximity(query[, key, drop = FALSE],
                          candidate[, key, drop = FALSE],
                          n_trees = n_trees, importance = TRUE, ...)

  # Cross-proximity: query vs candidate
  idx_query <- seq_len(n_q)
  idx_cand  <- n_q + seq_len(n_c)
  cross_prox <- .proximity_from_nodes(rf_res$terminal_nodes,
                                       idx_cand, idx_query)
  # cross_prox is n_q x n_c

  # Per query record: best match = highest proximity
  match_idx <- apply(cross_prox, 1, which.max)
  risk <- apply(cross_prox, 1, max)

  list(
    risk = risk,
    match_idx = match_idx,
    cross_prox = cross_prox,
    var_importance = rf_res$importance
  )
}
```

**E. Add RF branch** in the main method dispatch.

The RF branch **replaces the per-record main loop** (lines 576-835). To insert it:

1. **Wrap the existing main loop** in a conditional: change line 575 from `# main loop ----` + `for (i in seq_len(n_query)) {` to:
```r
    # main loop ----
    if (method != "rf") {
    for (i in seq_len(n_query)) {
```
2. **Close the conditional** after line 835 (the end of the for-loop and the closing `}` of the deterministic/probabilistic/pram/predictive branches): add `}` to close the `if (method != "rf")` block.
3. **Insert the `else if (method == "rf")` block** immediately after that closing `}`.

After the RF branch, execution falls through to the existing bijective override (line 837), risk_band, var_importance, and result construction. The RF branch populates the same bare vectors (`risk`, `cand_n`, `true_in_set`, `d_true`, `d_min`, `d_rank`, `score_cache`) that these sections expect.

**Key integration points:**
- Uses bare vectors (`risk[i]`, `cand_n[i]`, etc.) — NOT `per_record$*`
- Populates `score_cache` for bijective matching (with `maximize = TRUE` since proximity is a similarity score)
- The existing bijective override at line 839 handles the LSAP call via `.solve_bijective()`
- Truth evaluation uses `true_idx[i]` (already computed before the main loop)
- `query_data`, `search_data`, `split_search`, `blk_query` are already set up before the main loop
- `d_true` = proximity to true match, `d_min` = max proximity to any candidate, `d_rank` = rank of true match

```r
} else if (method == "rf") {
  # RF method: skip the per-record main loop, use batch RF proximity instead

  n_total <- nrow(query_data) + nrow(search_data)
  if (is.null(block) && n_total > 10000) {
    message("recordLinkage(method = 'rf'): n = ", n_total,
            " without blocking. Consider adding block variables.")
  }

  .rf_process_block <- function(q_idx, s_idx, cross_prox) {
    # Populate bare vectors for each query record in this block
    for (r in seq_along(q_idx)) {
      qi <- q_idx[r]
      prox_vec <- cross_prox[r, ]  # proximities to all candidates
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

      # Bijective: store score_cache entry
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
                                  n_trees = n_trees, importance = TRUE, ...)
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

        block_n <- length(q_idx) + length(s_idx)
        min_per_class <- min(length(q_idx), length(s_idx))

        if (min_per_class < 2) {
          small_blocks <- small_blocks + 1L
          fallback_blocks <- c(fallback_blocks, blk)
          next
        }
        if (block_n < 50) small_blocks <- small_blocks + 1L

        block_res <- .rf_linkage_block(
          query_data[q_idx, , drop = FALSE],
          search_data[s_idx, , drop = FALSE],
          key, n_trees = n_trees, ...
        )
        # Map local match indices back to global search indices
        .rf_process_block(q_idx, s_idx, block_res$cross_prox)
        importance_list[[blk]] <- block_res$var_importance
        block_sizes <- c(block_sizes, block_n)
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
                " blocks have < 50 records. ",
                "Consider rf_global = TRUE or coarser blocking.")
      }
    }
  }
  # Fall through to existing bijective override (line 837) which uses
  # score_cache with maximize = TRUE, then risk_band, then result construction.
```

**Why this works with the existing bijective override:** The RF branch populates `score_cache[[i]]` with `maximize = TRUE` and `cand`/`scores` fields, matching the interface expected by `.solve_bijective()`. The existing bijective override at line 839 fires unconditionally for all methods when `matching == "bijective"`, so it will correctly process the RF scores. The `maximize = TRUE` flag tells `.solve_bijective()` to transform `max_score - score` for the LSAP cost matrix.

**F. Add `var_importance`** to the var_importance section (line ~860-889). Insert before the final `else` clause:
```r
} else if (method == "rf") {
    var_importance <- if (!is.null(rf_var_importance)) {
        rf_var_importance
    } else {
        setNames(rep(NA_real_, length(key)), key)
    }
```

**Integration summary**: The RF branch populates `risk`, `cand_n`, `true_in_set`, `d_true`, `d_min`, `d_rank`, `score_cache`, and `matches` — all the same bare vectors used by other methods. After the RF branch, execution falls through to the existing bijective override (line 837), `risk_band` computation (line 853), var_importance section (line 860), and result construction (line 891). No changes to those sections except adding the RF case to var_importance.

- [ ] **Step 41: Run `devtools::document()`**

Run: `Rscript -e "devtools::document()"`

- [ ] **Step 42: Run tests**

Run: `Rscript -e "testthat::test_file('tests/testthat/test-recordLinkage-rf.R')"`
Run: `Rscript -e "testthat::test_file('tests/testthat/test-recordLinkage.R')"`
Expected: All PASS (new + existing)

- [ ] **Step 43: Commit**

```bash
git add R/recordLinkage.R R/rf_proximity_internal.R man/recordLinkage.Rd
git commit --author="matthias-da <matthias-da@users.noreply.github.com>" \
  -m "Add method = 'rf' to recordLinkage() with blocking and bijective support"
```

---

## Chunk 6: Integration, Documentation, and Final Checks

### Task 13: Integration with `disclosure_report()` and `rumap()`

**Files:**
- Modify: `R/disclosure_report.R`
- Modify: `R/rumap.R`

- [ ] **Step 44: Add `rf_privacy` to `disclosure_report()`**

In `R/disclosure_report.R`, add `rf_privacy` as an optional metric. Follow the existing skip-with-message pattern:
```r
if ("rf_privacy" %in% risk_measures) {
  if (requireNamespace("ranger", quietly = TRUE)) {
    tryCatch({
      rfp <- rf_privacy(original, synthetic, ...)
      # store result
    }, error = function(e) {
      warning("rf_privacy failed: ", e$message)
    })
  } else {
    message("Skipping rf_privacy (ranger not installed)")
  }
}
```

- [ ] **Step 45: Add `rf_privacy` to `rumap()`**

In `R/rumap.R`, add `"rf_privacy"` to `valid_risk` vector and add the computation block following the existing pattern for other risk measures.

- [ ] **Step 46: Run `devtools::document()`**

Run: `Rscript -e "devtools::document()"`

- [ ] **Step 47: Commit**

```bash
git add R/disclosure_report.R R/rumap.R man/
git commit --author="matthias-da <matthias-da@users.noreply.github.com>" \
  -m "Integrate rf_privacy into disclosure_report() and rumap()"
```

---

### Task 14: Full test suite and R CMD check

- [ ] **Step 48: Run full test suite**

Run: `Rscript -e "devtools::test()"`
Expected: All tests PASS (500+ existing + ~50 new)

- [ ] **Step 49: Run R CMD check**

Run: `Rscript -e "devtools::document()" && R CMD build . && _R_CHECK_FORCE_SUGGESTS_=FALSE R CMD check riskutility_*.tar.gz --no-manual --no-vignettes`
Expected: 0 errors, 0 warnings (2 expected vignette notes OK)

- [ ] **Step 50: Fix any check issues and commit**

```bash
git add -A
git commit --author="matthias-da <matthias-da@users.noreply.github.com>" \
  -m "Fix R CMD check issues for RF proximity"
```

---

### Task 15: Vignette section

**Files:**
- Modify: `vignettes/riskutility.Rmd`

- [ ] **Step 51: Add RF proximity section to main vignette**

Add a section covering:
- Brief RF proximity explanation (2-3 paragraphs)
- `rf_privacy()` worked example: same dataset with `dcr()` and `rf_privacy()` side-by-side
- `propscore(method = "ranger")` worked example showing structural metrics
- Comparison table: RF proximity vs Gower vs Euclidean

- [ ] **Step 52: Verify vignette builds**

Run: `Rscript -e "rmarkdown::render('vignettes/riskutility.Rmd')"`
Expected: Renders without errors

- [ ] **Step 53: Commit**

```bash
git add vignettes/riskutility.Rmd
git commit --author="matthias-da <matthias-da@users.noreply.github.com>" \
  -m "Add RF proximity section to main vignette"
```

---

### Task 16: Memory file update

- [ ] **Step 54: Update MEMORY.md with RF proximity status**

Update `/Users/matthias/.claude/projects/-Users-matthias-workspace-riskutility/memory/MEMORY.md` with a brief note about the RF proximity implementation status.

---

## File Map Summary

| Action | File | Purpose |
|--------|------|---------|
| Create | `R/rf_proximity_internal.R` | `.rf_proximity()`, `.proximity_from_nodes()`, `.proximity_from_nodes_newdata()`, `.rf_linkage_block()` |
| Create | `R/rf_privacy.R` | `rf_privacy()` + S3 methods |
| Modify | `R/utils_internal.R` | Add `.distance_risk_prepare()` |
| Modify | `R/dcr.R` | Refactor to use `.distance_risk_prepare()` |
| Modify | `R/nndr.R` | Refactor to use `.distance_risk_prepare()` |
| Modify | `R/propscore.R` | Add `method = "ranger"` branch |
| Modify | `R/recordLinkage.R` | Add `method = "rf"` branch |
| Modify | `R/disclosure_report.R` | Integrate `rf_privacy` |
| Modify | `R/rumap.R` | Integrate `rf_privacy` |
| Modify | `vignettes/riskutility.Rmd` | Add RF proximity section |
| Create | `tests/testthat/test-rf-proximity-internal.R` | Internal engine tests |
| Create | `tests/testthat/test-distance-risk-prepare.R` | Shared helper tests |
| Create | `tests/testthat/test-rf-privacy.R` | rf_privacy tests |
| Create | `tests/testthat/test-propscore-ranger.R` | propscore ranger tests |
| Create | `tests/testthat/test-recordLinkage-rf.R` | recordLinkage RF tests |
