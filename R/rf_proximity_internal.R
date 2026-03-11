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

  if (length(vars) == 0L) {
    stop("No common variables found between data1 and data2.", call. = FALSE)
  }

  if (nrow(data1) < 1L || nrow(data2) < 1L) {
    stop("data1 and data2 must have at least 1 row each.", call. = FALSE)
  }

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
    default_args$respect.unordered.factors <- "order"
    message("High-cardinality unordered factor detected (> 53 levels). ",
            "Using respect.unordered.factors = 'order'.")
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

  # prox[i, j] = fraction of trees where idx2[i] and idx1[j] share a node
  prox <- matrix(0L, nrow = n2, ncol = n1)

  if (progress && n_trees > 10) {
    pb <- utils::txtProgressBar(min = 0, max = n_trees, style = 3)
  } else {
    pb <- NULL
  }

  # Tree-by-tree accumulator
  for (t in seq_len(n_trees)) {
    tn1 <- nodes1[, t]  # terminal node IDs for idx1 in tree t
    tn2 <- nodes2[, t]  # terminal node IDs for idx2 in tree t

    # Group idx1 positions by their terminal node in this tree
    node_ids <- unique(tn1)
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

  # Combine: ref nodes + new nodes, then call with appropriate indices
  n_ref <- nrow(terminal_nodes_ref)
  n_new <- nrow(tn_new)
  combined_tn <- rbind(terminal_nodes_ref, tn_new)

  idx_new <- (n_ref + 1):(n_ref + n_new)
  .proximity_from_nodes(combined_tn, idx_ref, idx_new, progress = progress)
}
