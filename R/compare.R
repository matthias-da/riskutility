#' Compare Two Data Sets
#'
#' This function calculates and compares the empirical cumulative distribution
#' functions (ECDF) of a sample (\code{X}) and a population (\code{Y}) for
#' specified variables. It can handle weighted samples and provides options for
#' conditional CDFs, and approximation.
#'
#' @param X A data frame or data table.
#' @param Y A data frame or data table to compare with \code{X}.
#' @param variables A character vector specifying the variables for which the
#' comparisons are made. Must be present in both \code{X} and \code{Y}. The
#' class of variables determines the default kind of calculation. See details.
#' @param kind To be set when changing the default kind of calculation. See
#' details.
#' @param conditional A character vector specifying conditioning variables, or
#' \code{NULL} if no conditioning is desired. Must be of length 1 if specified.
#' @param weights A character vector of length 1 specifying the weights variable
#' in \code{X} and/or \code{Y}, or \code{NULL} if no weights are used.
#' @param approx A logical vector of length 2 indicating whether to approximate
#' the CDFs in \code{X} and \code{Y}, respectively. \code{TRUE} is recommended
#' for large data.
#' @param n_approx A numeric value specifying the number of points to use for
#' the approximation.
#' @param bounds A logical value indicating whether to use bounds in the CDF
#' calculation.
#' @export
#' @return A list with class \code{"compare"} containing the following
#' components:
#' \describe{
#'   \item{\code{formula}}{A formula representing the CDF comparison.}
#'   \item{\code{ecdf}}{A data frame with the ECDF values.}
#'   \item{\code{kind}}{A character string indicating the type of comparison
#'   ("numeric").}
#' }
#'
#' @details
#' The following default calculations are taken:
#'
#' If \code{variables} links to one or two numeric variables: the (weighted)
#' ECDF is calculated for both data sets. With \code{kind} one can change to
#' \code{"density"}, \code{"density_bayes"}, \code{"density_ratio"} and
#' \code{"density_ratio_bayes"}.
#'
#' If \code{variables} links to one, two or three categorical variables:
#' barcharts and mosaic plots are made.
#'
#' If \code{variables} links to one numeric and one categorical variable:
#' boxplots statistics are calculated.
#'
#' If weights are provided, the weighted versions are
#' computed. The function supports conditional estimates, which are computed by
#' conditioning on the specified variable.
#'
#' If the \code{approx} argument is set to \code{TRUE}, the ECDFs are
#' approximated using the specified number of points (\code{n_approx}).
#'
#' @importFrom data.table data.table as.data.table
#'
#' @examples
#' # Example data
#' S <- data.frame(
#'   age = sample(20:80, 100, replace = TRUE),
#'   income = rnorm(100, mean = 50000, sd = 10000),
#'   weights = runif(100, 0.5, 2)
#' )
#' U <- data.frame(
#'   age = sample(20:80, 10000, replace = TRUE),
#'   income = rnorm(10000, mean = 50000, sd = 10000)
#' )
#'
#' # Compare ECDFs for age and income
#' result <- compare(S, U, variables = c("age", "income"), weights = "weights")
#'
#' # View the result
#' print(result$formula)
#' head(result$ecdf)
#' plot(result)
#'
compare <- function(X, Y, variables = NULL, kind = NULL, conditional = NULL,
                          weights = NULL, approx = c(FALSE, TRUE),
                          n_approx = 10000, bounds = TRUE)
{
  # Function to check and convert to data.table if necessary
  convert_to_data_table <- function(X) {
    if (!inherits(X, "data.table")) {
      X <- as.data.table(X)
    }
    return(X)
  }

  # Convert to data.table if necessary
  X <- convert_to_data_table(X)
  Y <- convert_to_data_table(Y)

  n <- nrow(X)
  N <- nrow(Y)

  validate_weights <- function(weights) {
    # Check if weights is NULL
    if (is.null(weights)) {
      return(TRUE)
    }

    # Check if weights is a character vector of length 1
    if (is.character(weights) && length(weights) == 1) {
      return(TRUE)
    }

    # If neither condition is met, throw an error
    stop("weights must be either NULL or a character vector of length 1.")
  }

  validate_weights(weights)

  if(is.null(weights)){
    X$weights <- rep(1, n)
    weights.samp <- weights <- "weights"
  } else{
    weights.samp <- weights
  }

  if (!(weights %in% colnames(X))) {
    stop("The variables names specified in argument 'weight' must be available
         in the sample.\n")
  }
  if (weights %in% colnames(Y)) {
    weights.pop <- "weights"
  } else{
    weights.pop <- NULL
  }

  if (!is.character(variables) || length(variables) == 0) {
    stop("'variables' must be a character vector of positive length!\n")
  }
  if (!(all(variables %in% colnames(Y)) & (all(variables %in% colnames(X))))) {
    stop("The variables names specified in argument 'variables' must be available
         both in the population and the sample.\n")
  }
  if (!is.null(conditional) && !is.character(conditional)) {
    stop("'conditional' must be a character vector or NULL!\n")
    if (length(conditional) != 1) {
      stop("argument 'conditional' must have length 1!\n")
    }
  }
  if (!(all(conditional %in% colnames(Y)) & (all(conditional %in% colnames(X))))) {
    stop("The variables names specified in argument 'conditional' must be
         available both in the population and the sample.")
  }
  if (!is.logical(approx) || length(approx) == 0)
    approx <- formals()$approx
  else approx <- sapply(rep(approx, length.out = 2), isTRUE)
  if (any(approx) && (!is.numeric(n_approx) || length(n_approx) == 0)) {
    stop("'n_approx' is not numeric or does not have positive length.")
  }
  else n_approx <- ifelse(approx, n_approx[1], NA)
  bounds <- isTRUE(bounds)
  lab <- c("Sample", "Population")
  tmp <- simPop::getCdf(variables, weights.samp, conditional, X, approx = approx[1],
                n = n_approx[1], name = lab[1])
  values <- tmp$values
  values$.x <- as.numeric(as.character(values$.x))
  app <- t(tmp$approx)
  tmp <- simPop::getCdf(variables, weights.pop, conditional, Y, approx = approx[2],
                n = n_approx[2], name = lab[2])
  values$.x <- as.numeric(as.character(values$.x))
  values <- rbind(values, tmp$values)

  app <- rbind(app, tmp$approx)
  form <- ".y~.x"
  if (length(variables) > 1)
    conditional <- c(".var", conditional)
  if (!is.null(conditional)) {
    conditional <- paste(conditional, collapse = " + ")
    form <- paste(form, conditional, sep = " | ")
  }

  results <- list("formula" = form,
               "ecdf" = values,
               "kind" = "numeric")
  class(results) <- "compare"
  return(results)
}
NULL




#' Plot Method for Objects of Class "compare"
#'
#' This function provides visualizations for objects of class \code{"compare"}
#' produced by the \code{\link{compare}} function. It generates plots of the
#' empirical cumulative distribution functions (ECDF) and can optionally
#' produce interactive plots.
#'
#' @param x An object of class \code{"compare"}.
#' @param y Not used.
#' @param ... Additional arguments passed to the plotting function.
#' @param which A character string specifying the type of plot to generate.
#' Options include \code{"ecdf"}, \code{"density"}, \code{"density_bayes"},
#' \code{"density_ratio"}, and \code{"density_ratio_bayes"}.
#' Default is \code{"ecdf"}.
#' @param interactive A logical value indicating whether to produce an
#' interactive plot using \code{plotly}. Default is \code{FALSE}.
#'
#' @return Invisibly returns the ggplot2 object or the plotly object if
#' \code{interactive = TRUE}.
#'
#' @details
#' The \code{plot.compare} function generates various types of plots for
#' objects of class \code{"compare"}. The primary plot type is the ECDF plot,
#' which can be generated by setting \code{which} to \code{"ecdf"}.
#' Additional plot types (\code{"density"}, \code{"density_bayes"},
#' \code{"density_ratio"}, and \code{"density_ratio_bayes"}) are currently not
#' implemented and will raise an error if requested.
#'
#' The function supports both static and interactive plots.
#' If \code{interactive = TRUE}, the function will use \code{plotly} to
#' create an interactive version of the plot.
#'
#' @importFrom plotly ggplotly
#' @importFrom ggh4x facet_grid2
#'
#' @examples
#' S <- data.frame(
#'   age = sample(20:80, 100, replace = TRUE),
#'   income = rnorm(100, mean = 50000, sd = 10000),
#'   weights = runif(100, 0.5, 2)
#' )
#' U <- data.frame(
#'   age = sample(20:80, 10000, replace = TRUE),
#'   income = rnorm(10000, mean = 50000, sd = 10000)
#' )
#'
#' # Compare ECDFs for age and income
#' result <- compare(S, U, variables = c("age", "income"), weights = "weights")
#'
#' # Plot the ECDF
#' plot(result)
#'
#' \dontrun{
#' # Interactive plot
#' plot(result, interactive = TRUE)
#'
#' ## approx. 20 seconds computation time
#' data("eusilc13puf", package = "simPop")
#' # Function to replace NAs in factor columns with a new level
#' replace_na_in_factor <- function(factor_col, new_level = "not unique") {
#'   # Convert the factor to character
#'   char_col <- as.character(factor_col)
#'   # Replace NA with the new level
#'   char_col[is.na(char_col)] <- new_level
#'   # Convert back to factor and include the new level
#'   factor_col <- factor(char_col, levels = unique(c(levels(factor_col), new_level)))
#'   return(factor_col)
#' }
#'
#' # Apply the function to the relevant columns
#' eusilc13puf$pb220a <- replace_na_in_factor(eusilc13puf$pb220a)
#' eusilc13puf$pl031 <- replace_na_in_factor(eusilc13puf$pl031, new_level = "child")
#' eusilc13puf[is.na(eusilc13puf)] <- 0
#' eusilc13puf$age <- as.numeric(as.character(eusilc13puf$age))
#'
#' inp <- simPop::specifyInput(data=eusilc13puf, hhid="db030",
#'                             hhsize="hsize", strata="db040", weight="rb050")
#' simPop <- simPop::simStructure(data = inp, method = "direct",
#'                               basicHHvars=c("age", "rb090", "hsize", "db040", "pb220a"))
#' simPop <- simPop::simCategorical(simPop, additional=c("pl031"),
#'                                  method = "multinom", nr_cpus = 1)
#' # multinomial model with random draws
#' simPop <- simPop::simContinuous(simPop, additional="pgrossIncome",
#'                                 regModel = ~rb090+hsize+pl031+pb220a+age,
#'                                 upper=200000, equidist=FALSE, nr_cpus=1)
#'
#' eusilc13puf_synth <- data.frame(simPop::pop(simPop))
#' eusilc13puf_synth$age <- as.numeric(as.character(eusilc13puf_synth$age))
#'
#'
#' c1 <- compare(eusilc13puf, eusilc13puf_synth,
#'               variables = c("age", "pgrossIncome"),
#'               weights = "rb050",
#'               n_approx = 10000)
#' p1 <- plot(c1)
#' p1
#' p1$ggplot_object + theme_dark()
#' c2 <- compare(eusilc13puf, eusilc13puf_synth,
#'               variables = c("age", "pgrossIncome"),
#'               weights = "rb050",
#'               conditional = "pb220a",
#'               n_approx = 10000)
#' plot(c2)
#' c3 <- compare(eusilc13puf, eusilc13puf_synth,
#'               variables = c("pgrossIncome"),
#'               weights = "rb050",
#'               conditional = "rb090",
#'               n_approx = 10000)
#' plot(c3)
#' c4 <- compareSU_cdf(eusilc13puf, eusilc13puf_synth,
#'                     variables = c("age"),
#'                     weights = "rb050",
#'                     conditional = "rb090",
#'                    n_approx = 10000)
#' plot(c4)
#' c5 <- compareSU_cdf(eusilc13puf, eusilc13puf_synth,
#'                     variables = c("age"),
#'                     weights = "rb050",
#'                     conditional = "db040",
#'                     n_approx = 10000)
#' plot(c5)
#' }
#'
#' @export
plot.compare <- function(x, y, ..., which = "ecdf", interactive = FALSE){

  if(x$kind == "numeric" && (which == 1 || which == "ecdf")){
    # Function to convert formula to ggplot2 syntax with color and appropriate
    # faceting
    convert_formula_to_ggplot <- function(formula, data) {
      # Parse the formula
      formula_parts <- strsplit(as.character(formula), " \\| ")[[1]]
      yx_part <- formula_parts[1]

      # Extract x and y
      yx_split <- strsplit(yx_part, "~")[[1]]
      y_var <- trimws(yx_split[1])
      x_var <- trimws(yx_split[2])

      # Check if there are facets
      if (length(formula_parts) > 1) {
        facet_part <- formula_parts[2]
        facets <- strsplit(facet_part, "\\+")[[1]]
        facets <- trimws(facets)

        # Determine the appropriate faceting
        if (length(facets) == 1) {
          # Use facet_wrap with conditional scales = 'free_x'
          if (".var" %in% names(data) && length(unique(data$.var)) > 1) {
            facet_cmd <- paste0("facet_wrap(~ ", facets[1], ", scales = 'free_x')")
          } else {
            facet_cmd <- paste0("facet_wrap(~ ", facets[1], ")")
          }
        } else {
          # Use facet_grid with conditional scales = 'free_x'
          if (".var" %in% names(data) && length(unique(data$.var)) > 1) {
            facet_cmd <- paste0("ggh4x::facet_grid2(", facets[1], " ~ ", facets[2], ", scales = 'free_x', independent = 'x')")
          } else {
            facet_cmd <- paste0("ggh4x::facet_grid2(", facets[1], " ~ ", facets[2], ")")
          }
        }
      } else {
        # No faceting
        facet_cmd <- ""
      }

      # Create ggplot2 command with color
      ggplot_cmd <- paste0(
        "ggplot(data, aes(x = ", x_var, ", y = ", y_var, ",
        colour = .name)) + ",
        "geom_line() + theme_bw() + theme(legend.title=element_blank()) +
        xlab('x') + ylab('F(x)')",
        if (facet_cmd != "") paste0(" + ", facet_cmd) else ""
      )

      # Evaluate and return the ggplot2 command
      plot <- eval(parse(text = ggplot_cmd))
      plot <- plot
      return(list("ggplot_object" = plot, "cmd" = ggplot_cmd))
    }

    # Convert and plot using ggplot2
    ecdfplot <- convert_formula_to_ggplot(x$formula, x$ecdf)
    if(interactive){
      p <- plotly::ggplotly(ecdfplot$ggplot_object)
      print(p)
    } else{
      print(ecdfplot$ggplot_object)
    }
    cmd <- ecdfplot$cmd
  }
  ###
  if(x$kind == "numeric" && (which == 1 || which == "density")){
    stop("not implemented yet")
  }
  if(x$kind == "numeric" && (which == 1 || which == "density_bayes")){
    stop("not implemented yet")
  }
  if(x$kind == "numeric" && (which == 1 || which == "density_ratio")){
    stop("not implemented yet")
  }
  if(x$kind == "numeric" && (which == 1 || which == "density_ratio_bayes")){
    stop("not implemented yet")
  }
  ###

  invisible(ecdfplot)
}


