#' Compare Frequencies using Chi-square Goodness-of-Fit Test with Optional Grouping, Weights, and Simulated p-values
#'
#' This function performs a Chi-square goodness-of-fit test to compare the observed frequencies in an
#' anonymized/synthetic dataset (Y) with the expected frequencies from an original dataset (X), based on one or
#' more categorical variables (cat_vars). If sampling weights are provided via weight_X and weight_Y, weighted
#' frequency tables are computed using xtabs with a weight formula; otherwise, unweighted counts are computed.
#' Additionally, if group_vars is provided, the Chi-square test is computed separately for each unique group defined by the
#' grouping variables. For each group (or overall if no grouping is provided), the expected frequencies from X are scaled
#' to match the total count in Y before performing the test.
#'
#' @param X A data.frame or data.table containing the original dataset.
#' @param Y A data.frame or data.table containing the anonymized/synthetic dataset.
#' @param cat_vars A character vector specifying the categorical variables to form the contingency table.
#' @param group_vars Optional. A character vector specifying the grouping variables for which the test should be performed separately.
#'        Default is NULL, in which case the test is applied to the overall datasets.
#' @param weight_X Optional. A character string specifying the sampling weight variable in X.
#' @param weight_Y Optional. A character string specifying the sampling weight variable in Y.
#' @param simulate_p Logical. Whether to simulate p-values in the Chi-square test. Default is TRUE.
#' @param B Numeric. The number of simulations to perform if simulate_p is TRUE. Default is 2000.
#'
#' @return A data.table containing the grouping variables (if provided), the Chi-square statistic, degrees of freedom, and p-value
#'         for each group (or one overall row if group_vars is NULL).
#'
#' @details The function computes frequency tables for the variables in cat_vars for both datasets X and Y.
#' If sampling weights are provided, xtabs is used with a formula of the form "weight ~ var1 + var2 + ...";
#' otherwise, unweighted counts are computed using xtabs. A complete grid of all possible combinations of factor levels is
#' created to ensure that missing cells are filled with zeros. The expected frequencies from X are scaled to match the total count in Y
#' before performing the Chi-square goodness-of-fit test. When group_vars is provided, the function iterates over each unique combination
#' of grouping variables and performs the test on each subset. The simulation of p-values (when simulate_p is TRUE) helps to handle cases
#' with small expected frequencies.
#'
#' @importFrom data.table as.data.table rbindlist
#' @importFrom stats chisq.test xtabs
#' @export
#'
#' @examples
#' set.seed(123)
#' X <- data.frame(
#'   gender = factor(sample(c("Male", "Female"), 500, replace = TRUE)),
#'   region = factor(sample(c("North", "South", "East", "West"), 500, replace = TRUE)),
#'   occupation = factor(sample(c("Engineer", "Doctor", "Artist", "Teacher"), 500, replace = TRUE)),
#'   weight = runif(500, 0.5, 1.5)
#' )
#'
#' Y <- data.frame(
#'   gender = factor(sample(c("Male", "Female"), 1000, replace = TRUE)),
#'   region = factor(sample(c("North", "South", "East", "West"), 1000, replace = TRUE)),
#'   occupation = factor(sample(c("Engineer", "Doctor", "Artist", "Teacher"), 500, replace = TRUE)),
#'   weight = runif(1000, 0.5, 1.5)
#' )
#'
#' # Overall test (weighted) with simulated p-values
#' result_overall <- compare_chisq_gof(X, Y, cat_vars = c("gender", "region"),
#'                                     weight_X = "weight", weight_Y = "weight",
#'                                     simulate_p = FALSE, B = 2000)
#' print(result_overall)
#'
#' # Test conditionally by grouping variables (e.g., "region" and "gender")
#' result_by_group <- compare_chisq_gof(X, Y, cat_vars = c("gender"), group_vars = c("region", "occupation"),
#'                                      weight_X = "weight", weight_Y = "weight",
#'                                      simulate_p = TRUE, B = 2000)
#' print(result_by_group)
compare_chisq_gof <- function(X, Y, cat_vars, group_vars = NULL, weight_X = NULL, weight_Y = NULL,
                              simulate_p = TRUE, B = 2000) {
  library(data.table)

  # Convert X and Y to data.table if not already
  X <- as.data.table(X)
  Y <- as.data.table(Y)

  # Check that all specified categorical variables exist in both datasets
  if (!all(cat_vars %in% names(X)) || !all(cat_vars %in% names(Y))) {
    stop("All categorical variables specified in cat_vars must be present in both X and Y.")
  }

  # If grouping variables are provided, check they exist in both datasets
  if (!is.null(group_vars)) {
    if (!all(group_vars %in% names(X)) || !all(group_vars %in% names(Y))) {
      stop("All grouping variables specified in group_vars must be present in both X and Y.")
    }
    # Create a group key in both datasets based on group_vars
    X[, group_key := do.call(paste, c(.SD, sep = "_")), .SDcols = group_vars]
    Y[, group_key := do.call(paste, c(.SD, sep = "_")), .SDcols = group_vars]
    groups <- unique(rbind(X[, .(group_key)], Y[, .(group_key)]))$group_key
  } else {
    # No grouping: assign a default group key
    X[, group_key := "Overall"]
    Y[, group_key := "Overall"]
    groups <- "Overall"
  }

  results_list <- list()

  # Loop over each group and perform Chi-square GOF test
  for (g in groups) {
    if (g == "Overall") {
      X_sub <- X
      Y_sub <- Y
      group_label <- "Overall"
    } else {
      X_sub <- X[group_key == g]
      Y_sub <- Y[group_key == g]
      group_label <- g
    }

    # Only proceed if both subsets have data
    if (nrow(X_sub) == 0 || nrow(Y_sub) == 0) next

    # Compute weighted frequency table for X using xtabs with weights
    if (!is.null(weight_X) && weight_X %in% names(X_sub)) {
      expected_table <- xtabs(as.formula(paste(weight_X, "~", paste(cat_vars, collapse = "+"))),
                              data = X_sub)
    } else {
      expected_table <- xtabs(~ ., data = X_sub[, cat_vars, with = FALSE])
    }

    # Compute weighted frequency table for Y
    if (!is.null(weight_Y) && weight_Y %in% names(Y_sub)) {
      observed_table <- xtabs(as.formula(paste(weight_Y, "~", paste(cat_vars, collapse = "+"))),
                              data = Y_sub)
    } else {
      observed_table <- xtabs(~ ., data = Y_sub[, cat_vars, with = FALSE])
    }

    # Create a complete grid of all possible factor combinations for cat_vars
    complete_levels <- lapply(X_sub[, cat_vars, with = FALSE], function(x) levels(as.factor(x)))
    complete_grid <- do.call(expand.grid, complete_levels)

    # Merge complete grid with frequency tables to fill missing combinations with zeros
    expected_df <- as.data.frame(expected_table)
    observed_df <- as.data.frame(observed_table)

    complete_expected <- merge(complete_grid, expected_df, by = cat_vars, all.x = TRUE)
    complete_observed <- merge(complete_grid, observed_df, by = cat_vars, all.x = TRUE)

    complete_expected[is.na(complete_expected$Freq), "Freq"] <- 0
    complete_observed[is.na(complete_observed$Freq), "Freq"] <- 0

    # Recreate full contingency tables from the complete data frames
    expected_full <- xtabs(Freq ~ ., data = complete_expected)
    observed_full <- xtabs(Freq ~ ., data = complete_observed)

    # Scale expected frequencies to match the total count in Y_sub
    expected_scaled <- expected_full * (sum(observed_full) / sum(expected_full))

    # Perform Chi-square goodness-of-fit test
    chi_res <- chisq.test(x = observed_full,
                          p = as.vector(expected_scaled) / sum(expected_scaled),
                          simulate.p.value = simulate_p, B = B)

    results_list[[group_label]] <- data.table(group = group_label,
                                              chi_sq_statistic = chi_res$statistic,
                                              df = chi_res$parameter,
                                              p_value = chi_res$p.value)
  }

  result_dt <- rbindlist(results_list, fill = TRUE)
  return(result_dt)
}

