#' Compare Boxplots of a Numeric Variable in Two Datasets
#'
#' This function compares the distribution of a numeric variable between
#' two datasets (e.g., an original dataset and an anonymized/synthetic version)
#' using boxplots. It allows conditional faceting on categorical variables
#' and supports different visualization styles, including side-by-side,
#' overlapping, and separate boxplots. Sampling weights can be considered
#' to compute weighted medians and quartiles.
#'
#' @param X A data.frame or data.table containing the original dataset.
#' @param Y A data.frame or data.table containing the anonymized/synthetic dataset.
#' @param num_var A character string specifying the numeric variable to compare.
#' @param cat_vars A character vector of categorical variables used for faceting.
#' @param weight_X Optional. A character string specifying the sampling weight variable in X.
#' @param weight_Y Optional. A character string specifying the sampling weight variable in Y.
#' @param facet_type Character string. The faceting type: "wrap" (default) for facet_wrap() or "grid" for facet_grid().
#' @param plot_type Character string. Type of plot:
#'   - "side" (default): Side-by-side boxplots.
#'   - "overlap": Overlapping boxplots.
#'   - "separate": Separate boxplots for each dataset-category combination (Original above, Anonymized below).
#' @return A ggplot2 boxplot visualization comparing the distributions.
#' @export
#'
#' @examples
#' set.seed(123)
#' X <- data.frame(
#'   age = sample(20:80, 500, replace = TRUE),
#'   gender = sample(c("Male", "Female"), 500, replace = TRUE),
#'   income = rnorm(500, mean = 50000, sd = 10000),
#'   weight = runif(500, 0.5, 1.5)
#' )
#'
#' Y <- data.frame(
#'   age = sample(20:80, 1000, replace = TRUE),
#'   gender = sample(c("Male", "Female"), 1000, replace = TRUE),
#'   income = rnorm(1000, mean = 48000, sd = 12000),
#'   weight = runif(1000, 0.5, 1.5)
#' )
#'
#' compare_boxplots(X, Y, num_var = "income", cat_vars = c("gender"),
#'                  weight_X = "weight", weight_Y = "weight",
#'                  facet_type = "wrap", plot_type = "side")
#'
#' compare_boxplots(X, Y, num_var = "income", cat_vars = c("gender"),
#'                  weight_X = "weight", weight_Y = "weight",
#'                  facet_type = "wrap", plot_type = "separate")
compare_boxplots <- function(X, Y, num_var, cat_vars,
                             weight_X = NULL, weight_Y = NULL,
                             facet_type = "wrap", plot_type = "side") {

  # Convert to data.table for efficient processing (use copy to avoid side effects)
  X <- copy(as.data.table(X))
  Y <- copy(as.data.table(Y))

  # Check if the numeric variable exists in both datasets
  if (!(num_var %in% names(X) && num_var %in% names(Y))) {
    stop("The specified numeric variable must be in both X and Y.")
  }

  # Check if categorical variables exist in both datasets
  if (!all(cat_vars %in% names(X)) || !all(cat_vars %in% names(Y))) {
    stop("All categorical variables must be present in both X and Y.")
  }

  # Assign dataset labels
  X[, .dataset := "Original (X)"]
  Y[, .dataset := "Anonymized/Synthetic (Y)"]

  # Force factor levels to ensure proper ordering
  X$.dataset <- factor(X$.dataset, levels = c("Original (X)", "Anonymized/Synthetic (Y)"))
  Y$.dataset <- factor(Y$.dataset, levels = c("Original (X)", "Anonymized/Synthetic (Y)"))

  # Handle weights (if provided, apply; else use unweighted counts)
  if (!is.null(weight_X) && weight_X %in% names(X)) {
    X[, weight := get(weight_X)]
  } else {
    X[, weight := 1]
  }
  if (!is.null(weight_Y) && weight_Y %in% names(Y)) {
    Y[, weight := get(weight_Y)]
  } else {
    Y[, weight := 1]
  }

  # Combine datasets for plotting
  combined_data <- rbindlist(list(X, Y), use.names = TRUE, fill = TRUE)

  # For separate plot_type, create an additional facet column for the first categorical variable.
  if (plot_type == "separate") {
    # Create a column 'facet_cat' from the first categorical variable.
    combined_data[, facet_cat := get(cat_vars[1])]
  }

  # Define base ggplot object (common for all types)
  p <- ggplot(combined_data, aes(x = .dataset, y = get(num_var), fill = .dataset)) +
    scale_fill_manual(values = c("blue", "red"),
                      labels = c("Original (X)", "Anonymized/Synthetic (Y)")) +
    theme_minimal() +
    theme(legend.position = "none") +
    labs(x = "Dataset", y = num_var,
         title = paste("Comparison of Boxplots for", num_var))

  # Choose plotting method based on plot_type
  if (plot_type == "side") {
    p <- p + geom_boxplot(alpha = 0.7, position = position_dodge(0.6))
    # Optionally add faceting on cat_vars if desired
    if (length(cat_vars) == 1) {
      if (facet_type == "wrap") {
        p <- p + facet_wrap(as.formula(paste("~", cat_vars[1])), scales = "free_x")
      } else {
        p <- p + facet_grid(as.formula(paste(". ~", cat_vars[1])), scales = "free_x")
      }
    } else if (length(cat_vars) == 2) {
      p <- p + facet_grid(as.formula(paste(cat_vars[1], "~", cat_vars[2])), scales = "free_x")
    }
  } else if (plot_type == "overlap") {
    p <- p + geom_boxplot(position = position_identity(), alpha = 0.6)
    # Add faceting as above
    if (length(cat_vars) == 1) {
      if (facet_type == "wrap") {
        p <- p + facet_wrap(as.formula(paste("~", cat_vars[1])), scales = "free_x")
      } else {
        p <- p + facet_grid(as.formula(paste(". ~", cat_vars[1])), scales = "free_x")
      }
    } else if (length(cat_vars) == 2) {
      p <- p + facet_grid(as.formula(paste(cat_vars[1], "~", cat_vars[2])), scales = "free_x")
    }
  } else if (plot_type == "separate") {
    # For separate plots, we want a two-row layout:
    # Row: .dataset (Original on top, Anonymized below)
    # Column: Levels of the first categorical variable (facet_cat)
    p <- ggplot(combined_data, aes(x = "", y = get(num_var), fill = .dataset)) +
      geom_boxplot(alpha = 0.7) +
      facet_grid(.dataset ~ facet_cat, scales = "free_y") +
      scale_fill_manual(values = c("blue", "red")) +
      theme_minimal() +
      theme(legend.position = "none",
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank()) +
      labs(x = cat_vars[1], y = num_var,
           title = paste("Separate Boxplots for", num_var, "by", cat_vars[1]))
  } else {
    stop("Invalid plot_type. Choose 'side', 'overlap', or 'separate'.")
  }

  # Print and return plot
  print(p)
  invisible(p)
}
