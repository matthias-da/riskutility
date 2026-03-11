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
