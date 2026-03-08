#' Compare Distributions using the Wasserstein Distance
#'
#' This function computes the one-dimensional Wasserstein distance between a numeric variable in two datasets,
#' X (original) and Y (anonymized or synthetic). For continuous variables, the Wasserstein distance is approximated
#' by integrating the absolute differences between the quantile functions computed on a grid of probabilities.
#' If sampling weights are provided, weighted quantiles are computed using \code{Hmisc::wtd.quantile}. For nominal
#' (categorical) variables, the function computes the total variation distance (half the L1 distance between the probability
#' distributions) assuming a unit cost for mismatches.
#'
#' When grouping variables are provided (via \code{cat_vars}), the Wasserstein distance is computed within each group.
#'
#' @param X A data.frame or data.table containing the original dataset.
#' @param Y A data.frame or data.table containing the anonymized/synthetic dataset.
#' @param num_var A character string specifying the numeric (or nominal) variable to compare.
#' @param cat_vars Optional. A character vector of categorical variables for grouping. Default is NULL.
#' @param weight_X Optional. A character string specifying the sampling weight variable in X.
#' @param weight_Y Optional. A character string specifying the sampling weight variable in Y.
#' @param var_type A character string specifying the type of variable: "continuous" or "nominal".
#'        Default is "auto", which infers "continuous" if the variable is numeric, else "nominal".
#' @param n_grid Number of grid points for quantile approximation (for continuous variables). Default is 1000.
#'
#' @param ... additional arguments passed to methods
#'
#' @return A data.table with the computed Wasserstein distance for each group (or overall if no grouping is provided).
#'         For continuous variables, the Wasserstein distance is approximated by the mean absolute difference between the
#'         quantile functions. For nominal variables, the result is the total variation distance.
#'
#' @importFrom data.table as.data.table rbindlist
#' @family comparison
#' @author Matthias Templ
#' @export
#'
#' @examples
#' set.seed(123)
#' X <- data.frame(
#'   income = rnorm(500, mean = 50000, sd = 10000),
#'   gender = sample(c("Male", "Female"), 500, replace = TRUE),
#'   region = sample(c("North", "South", "East", "West"), 500, replace = TRUE),
#'   weight = runif(500, 0.5, 1.5)
#' )
#'
#' Y <- data.frame(
#'   income = rnorm(1000, mean = 48000, sd = 12000),
#'   gender = sample(c("Male", "Female"), 1000, replace = TRUE),
#'   region = sample(c("North", "South", "East", "West"), 1000, replace = TRUE),
#'   weight = runif(1000, 0.5, 1.5)
#' )
#'
#' # Compare overall continuous distributions (weighted)
#' compare_wasserstein(X, Y, num_var = "income", var_type = "continuous",
#'                     weight_X = "weight", weight_Y = "weight")
#'
#' # Compare distributions within groups defined by gender
#' compare_wasserstein(X, Y, num_var = "income", cat_vars = c("gender"),
#'                     var_type = "continuous", weight_X = "weight", weight_Y = "weight")
#' compare_wasserstein(X, Y, num_var = "income", cat_vars = c("gender","region"),
#'                     var_type = "continuous", weight_X = "weight", weight_Y = "weight")
#'
#' # Example for nominal variables
#'
#' # Create a synthetic dataset X with a nominal variable "gender"
#' set.seed(123)
#' X_nom <- data.frame(
#'   gender = sample(c("Male", "Female"), 500, replace = TRUE),
#'   weight = runif(500, 0.5, 1.5)
#' )
#'
#' # Create a synthetic dataset Y with a nominal variable "gender"
#' Y_nom <- data.frame(
#'   gender = sample(c("Male", "Female"), 1000, replace = TRUE),
#'   weight = runif(1000, 0.5, 1.5)
#' )
#'
#' # Compare the distributions of the nominal variable using the Wasserstein distance.
#' # Here, var_type is explicitly set to "nominal" so that the function uses frequency tables.
#' result_nominal <- compare_wasserstein(X_nom, Y_nom, num_var = "gender",
#'                                       var_type = "nominal",
#'                                       weight_X = "weight", weight_Y = "weight")
#'
#' print(result_nominal)
#'
compare_wasserstein <- function(X, ...) {
  UseMethod("compare_wasserstein")
}

#' @rdname compare_wasserstein
#' @export
compare_wasserstein.synth_pair <- function(X, ...) {
  compare_wasserstein.default(X = X$original, Y = X$synthetic, ...)
}

#' @rdname compare_wasserstein
#' @export
compare_wasserstein.default <- function(X, Y, num_var, cat_vars = NULL,
                                weight_X = NULL, weight_Y = NULL,
                                var_type = "auto", n_grid = 1000, ...) {
  # Convert datasets to data.table (use copy to avoid side effects)
  X <- copy(as.data.table(X))
  Y <- copy(as.data.table(Y))

  # Check if num_var exists in both datasets
  if (!(num_var %in% names(X) && num_var %in% names(Y))) {
    stop("The specified variable must be present in both X and Y.")
  }

  # If var_type is "auto", determine type based on X[[num_var]]
  if (var_type == "auto") {
    if (is.numeric(X[[num_var]])) {
      var_type <- "continuous"
    } else {
      var_type <- "nominal"
    }
  }

  # If grouping variables are provided, check they exist
  if (!is.null(cat_vars)) {
    if (!all(cat_vars %in% names(X)) || !all(cat_vars %in% names(Y))) {
      stop("All grouping variables must be present in both X and Y.")
    }
    # Merge grouping variables into one key for each dataset
    X[, group_key := do.call(paste, c(.SD, sep = "_")), .SDcols = cat_vars]
    Y[, group_key := do.call(paste, c(.SD, sep = "_")), .SDcols = cat_vars]
    groups <- unique(rbind(X[, .(group_key)], Y[, .(group_key)]))$group_key
  } else {
    groups <- NA  # single overall group
  }

  # Function to compute Wasserstein distance for continuous variables
  compute_continuous_wasserstein <- function(x, y, wx = NULL, wy = NULL) {
    p_grid <- seq(0, 1, length.out = n_grid)
    # If weights provided, use Hmisc::wtd.quantile
    if (!is.null(wx)) {
      if (!requireNamespace("Hmisc", quietly = TRUE)) {
        stop("Package 'Hmisc' is required for weighted quantiles. Please install it.")
      }
      qx <- as.numeric(Hmisc::wtd.quantile(x, weights = wx, probs = p_grid, na.rm = TRUE))
    } else {
      qx <- as.numeric(quantile(x, probs = p_grid, na.rm = TRUE))
    }
    if (!is.null(wy)) {
      if (!requireNamespace("Hmisc", quietly = TRUE)) {
        stop("Package 'Hmisc' is required for weighted quantiles. Please install it.")
      }
      qy <- as.numeric(Hmisc::wtd.quantile(y, weights = wy, probs = p_grid, na.rm = TRUE))
    } else {
      qy <- as.numeric(quantile(y, probs = p_grid, na.rm = TRUE))
    }
    # Approximate the 1-Wasserstein distance (L1 norm between quantile functions)
    W <- mean(abs(qx - qy))
    return(W)
  }

  # Function to compute Wasserstein (total variation) for nominal variables
  compute_nominal_wasserstein <- function(x, y, wx = NULL, wy = NULL) {
    # Create frequency tables (weighted if provided)
    if (!is.null(wx)) {
      tab_x <- xtabs(wx ~ x)
    } else {
      tab_x <- table(x)
    }
    if (!is.null(wy)) {
      tab_y <- xtabs(wy ~ y)
    } else {
      tab_y <- table(y)
    }
    # Convert counts to proportions
    p_x <- tab_x / sum(tab_x)
    p_y <- tab_y / sum(tab_y)
    # Align categories
    all_levels <- union(names(p_x), names(p_y))
    p_x <- p_x[all_levels]; p_x[is.na(p_x)] <- 0
    p_y <- p_y[all_levels]; p_y[is.na(p_y)] <- 0
    # For nominal variables, assume cost=1 for mismatches, so Wasserstein distance equals half L1 distance
    W <- sum(abs(p_x - p_y)) / 2
    return(W)
  }

  # Initialize results list
  results_list <- list()

  # Loop over groups (if groups is NA, then do overall)
  if (all(is.na(groups))) {
    groups <- NA
  }

  for (g in groups) {
    if (is.na(g)) {
      X_sub <- X
      Y_sub <- Y
      group_label <- "Overall"
    } else {
      X_sub <- X[group_key == g]
      Y_sub <- Y[group_key == g]
      group_label <- g
    }

    if (nrow(X_sub) == 0 || nrow(Y_sub) == 0) next

    if (var_type == "continuous") {
      W <- compute_continuous_wasserstein(X_sub[[num_var]], Y_sub[[num_var]],
                                          if (!is.null(weight_X)) X_sub[[weight_X]] else NULL,
                                          if (!is.null(weight_Y)) Y_sub[[weight_Y]] else NULL)
    } else if (var_type == "nominal") {
      W <- compute_nominal_wasserstein(X_sub[[num_var]], Y_sub[[num_var]],
                                       if (!is.null(weight_X)) X_sub[[weight_X]] else NULL,
                                       if (!is.null(weight_Y)) Y_sub[[weight_Y]] else NULL)
    } else {
      stop("var_type must be either 'continuous' or 'nominal'.")
    }

    results_list[[group_label]] <- data.table(group = group_label, wasserstein = W, var_type = var_type)
  }

  result_dt <- rbindlist(results_list, fill = TRUE)
  return(result_dt)
}
