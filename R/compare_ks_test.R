#' Compare Distributions using the Kolmogorov-Smirnov Test
#'
#' This function performs a non-parametric Kolmogorov-Smirnov (KS) test to compare the cumulative distribution
#' functions of a numeric variable between an original dataset (X) and an anonymized/synthetic dataset (Y). The test
#' measures the maximum deviation between the two empirical cumulative distribution functions. The function can also
#' perform the test separately for groups defined by one or more categorical variables.
#'
#' @param X A data.frame or data.table containing the original dataset.
#' @param Y A data.frame or data.table containing the anonymized/synthetic dataset.
#' @param num_var A character string specifying the numeric variable to compare.
#' @param cat_vars Optional. A character vector of categorical variables used for grouping. Default is NULL, in which
#'        case the KS test is performed on the entire datasets.
#'
#' @param ... additional arguments passed to methods
#'
#' @return A data.table with columns for the grouping variables (if provided), the KS test statistic, and the corresponding
#'         p-value. If no grouping is provided, a data.table with one row is returned.
#' @family comparison
#' @export
#'
#' @examples
#' set.seed(123)
#' X <- data.frame(
#'   age = sample(20:80, 500, replace = TRUE),
#'   gender = sample(c("Male", "Female"), 500, replace = TRUE),
#'   region = sample(c("North", "South", "East", "West"), 500, replace = TRUE),
#'   income = rnorm(500, mean = 50000, sd = 10000)
#' )
#'
#' Y <- data.frame(
#'   age = sample(20:80, 1000, replace = TRUE),
#'   gender = sample(c("Male", "Female"), 1000, replace = TRUE),
#'   region = sample(c("North", "South", "East", "West"), 500, replace = TRUE),
#'   income = rnorm(1000, mean = 48000, sd = 12000)
#' )
#'
#' # Perform KS test on the entire datasets
#' compare_ks_test(X, Y, num_var = "income")
#'
#' # Perform KS test within groups defined by gender
#' compare_ks_test(X, Y, num_var = "income", cat_vars = c("gender"))
#' compare_ks_test(X, Y, num_var = "income", cat_vars = c("gender","region"))
compare_ks_test <- function(X, ...) {
  UseMethod("compare_ks_test")
}

#' @rdname compare_ks_test
#' @export
compare_ks_test.synth_pair <- function(X, ...) {
  compare_ks_test.default(X = X$original, Y = X$synthetic, ...)
}

#' @rdname compare_ks_test
#' @export
compare_ks_test.default <- function(X, Y, num_var, cat_vars = NULL, ...) {
  # Convert X and Y to data.table if they are not already
  X <- as.data.table(X)
  Y <- as.data.table(Y)

  # Check that the numeric variable exists in both datasets
  if (!(num_var %in% names(X) && num_var %in% names(Y))) {
    stop("The specified numeric variable must be present in both X and Y.")
  }

  # If grouping variables are provided, check they exist in both datasets
  if (!is.null(cat_vars)) {
    if (!all(cat_vars %in% names(X)) || !all(cat_vars %in% names(Y))) {
      stop("All grouping categorical variables must be present in both X and Y.")
    }
  }

  # If no grouping is specified, perform a KS test on the full datasets
  if (is.null(cat_vars)) {
    ks_result <- ks.test(X[[num_var]], Y[[num_var]])
    result <- data.table(ks_statistic = ks_result$statistic, p_value = ks_result$p.value)
    return(result)
  }

  # If grouping is specified, perform KS test for each group
  # Get the unique combinations of grouping variables from the union of both datasets
  groups_X <- unique(X[, ..cat_vars])
  groups_Y <- unique(Y[, ..cat_vars])
  all_groups <- unique(rbind(groups_X, groups_Y))

  result_list <- list()

  # Iterate over each group combination
  for (i in seq_len(nrow(all_groups))) {
    group_vals <- all_groups[i, ]
    # Subset each dataset for the current group using a join
    X_grp <- X[group_vals, on = cat_vars]
    Y_grp <- Y[group_vals, on = cat_vars]

    # Only perform the test if both subsets have data
    if (nrow(X_grp) > 0 && nrow(Y_grp) > 0) {
      ks_result <- ks.test(X_grp[[num_var]], Y_grp[[num_var]])
      # Create a data.table row with the group values, KS statistic and p-value
      result_list[[i]] <- data.table(all_groups[i, ],
                                     ks_statistic = ks_result$statistic,
                                     p_value = ks_result$p.value)
    }
  }

  result <- rbindlist(result_list, fill = TRUE)
  return(result)
}
