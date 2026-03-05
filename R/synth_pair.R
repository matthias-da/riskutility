#' Create a Synthetic Data Comparison Pair
#'
#' Creates a container object that holds both the original and synthetic datasets
#' along with metadata for consistent analysis across multiple risk and utility measures.
#'
#' @param original Original/training data (data.frame or data.table)
#' @param synthetic Synthetic/anonymized data (data.frame or data.table)
#' @param key_vars Character vector of quasi-identifier variables (for attribution-based risk).
#'   Required for \code{\link{dcap}}, \code{\link{tcap}}, \code{\link{disco}}, \code{\link{weap}}.
#' @param target_var Character string naming the sensitive target variable (for CAP metrics).
#'   Required for \code{\link{dcap}}, \code{\link{tcap}}, \code{\link{disco}}.
#' @param weight_original Character string naming the weight variable in original data, or NULL.
#' @param weight_synthetic Character string naming the weight variable in synthetic data, or NULL.
#' @param vars Character vector of variables to use for utility comparisons.
#'   If NULL (default), all common variables between original and synthetic are used.
#' @param holdout Optional data.frame for distance-based risk measures (\code{\link{dcr}}, \code{\link{nndr}}).
#'   If NULL, these functions will use internal splitting.
#' @param source Character string describing data source (e.g., "synthpop", "simPop", "custom").
#' @param metadata Named list of additional metadata to store with the pair.
#'
#' @return An object of class "synth_pair" containing:
#' \itemize{
#'   \item original: the original dataset
#'   \item synthetic: the synthetic dataset
#'   \item key_vars, target_var: for attribution-based risk measures
#'   \item weight_original, weight_synthetic: weight variable names
#'   \item vars: variables for utility comparison
#'   \item cat_vars, num_vars: auto-detected categorical and numeric variables
#'   \item holdout: optional holdout data for distance-based measures
#'   \item source: data source identifier
#'   \item metadata: additional metadata
#'   \item n_original, n_synthetic: dataset sizes
#' }
#'
#' @details
#' The \code{synth_pair} class provides a convenient container for synthetic data
#' evaluation workflows. Instead of repeatedly passing the same parameters to
#' multiple functions, you can create a single object that carries all the
#' context needed for analysis.
#'
#' Many functions in this package have methods for \code{synth_pair} objects,
#' allowing simplified calls like \code{dcap(pair)} instead of
#' \code{dcap(X, Y, key_vars = ..., target_var = ...)}.
#'
#' For users of synthpop, simPop, or sdcMicro packages, convenience constructors
#' are available:
#' \itemize{
#'   \item \code{\link{from_synthpop}}: Create from synthpop's synds objects
#'   \item \code{\link{from_simPop}}: Create from simPop's simPopObj objects
#'   \item \code{\link{from_sdcMicro}}: Create from sdcMicro's sdcMicroObj objects
#' }
#'
#' @seealso
#' \code{\link{from_synthpop}}, \code{\link{from_simPop}}, \code{\link{from_sdcMicro}}
#' for package-specific constructors,
#' \code{\link{dcap}}, \code{\link{tcap}}, \code{\link{ims}}, \code{\link{hellinger}}
#' for functions with synth_pair methods
#'
#' @importFrom utils head
#' @family containers
#' @export
#' @examples
#' # Create example data
#' set.seed(123)
#' original <- data.frame(
#'   age = sample(20:70, 200, replace = TRUE),
#'   gender = sample(c("M", "F"), 200, replace = TRUE),
#'   region = sample(c("N", "S", "E", "W"), 200, replace = TRUE),
#'   income = sample(c("low", "medium", "high"), 200, replace = TRUE)
#' )
#'
#' synthetic <- data.frame(
#'   age = sample(20:70, 200, replace = TRUE),
#'   gender = sample(c("M", "F"), 200, replace = TRUE),
#'   region = sample(c("N", "S", "E", "W"), 200, replace = TRUE),
#'   income = sample(c("low", "medium", "high"), 200, replace = TRUE)
#' )
#'
#' # Create a synth_pair object
#' pair <- synth_pair(
#'   original = original,
#'   synthetic = synthetic,
#'   key_vars = c("age", "gender", "region"),
#'   target_var = "income"
#' )
#'
#' print(pair)
#'
#' # Now functions can be called directly on the pair
#' # dcap(pair)
#' # ims(pair)
#' # hellinger(pair)
synth_pair <- function(original,
                       synthetic,
                       key_vars = NULL,
                       target_var = NULL,
                       weight_original = NULL,
                       weight_synthetic = NULL,
                       vars = NULL,
                       holdout = NULL,
                       source = "custom",
                       metadata = list()) {

 # Input validation
  if (!is.data.frame(original)) {
    stop("'original' must be a data.frame or data.table")
  }
  if (!is.data.frame(synthetic)) {
    stop("'synthetic' must be a data.frame or data.table")
 }

  # Convert to data.frame for consistency
  original <- as.data.frame(original)
  synthetic <- as.data.frame(synthetic)

  # Find common variables
  common_vars <- intersect(names(original), names(synthetic))
  if (length(common_vars) == 0) {
    stop("No common variables found between original and synthetic datasets")
  }

  # Set vars to common variables if not specified
  if (is.null(vars)) {
    vars <- common_vars
    # Exclude weight variables from default vars
    if (!is.null(weight_original)) vars <- setdiff(vars, weight_original)
    if (!is.null(weight_synthetic)) vars <- setdiff(vars, weight_synthetic)
  }

  # Validate key_vars and target_var if provided
 if (!is.null(key_vars)) {
    missing_key <- setdiff(key_vars, common_vars)
    if (length(missing_key) > 0) {
      stop(paste("key_vars not found in both datasets:", paste(missing_key, collapse = ", ")))
    }
  }

  if (!is.null(target_var)) {
    if (!target_var %in% common_vars) {
      stop(paste("target_var", target_var, "not found in both datasets"))
    }
  }

  # Validate weight variables
  if (!is.null(weight_original) && !weight_original %in% names(original)) {
    stop(paste("weight_original", weight_original, "not found in original data"))
  }
  if (!is.null(weight_synthetic) && !weight_synthetic %in% names(synthetic)) {
    stop(paste("weight_synthetic", weight_synthetic, "not found in synthetic data"))
  }

  # Auto-detect variable types from vars
  cat_vars <- names(original)[sapply(original[, vars, drop = FALSE],
                                      function(x) is.factor(x) || is.character(x))]
  num_vars <- names(original)[sapply(original[, vars, drop = FALSE], is.numeric)]

  # Exclude weight variables from cat/num vars
  if (!is.null(weight_original)) {
    cat_vars <- setdiff(cat_vars, weight_original)
    num_vars <- setdiff(num_vars, weight_original)
  }

  structure(
    list(
      original = original,
      synthetic = synthetic,
      key_vars = key_vars,
      target_var = target_var,
      weight_original = weight_original,
      weight_synthetic = weight_synthetic,
      vars = vars,
      cat_vars = cat_vars,
      num_vars = num_vars,
      holdout = holdout,
      source = source,
      metadata = metadata,
      n_original = nrow(original),
      n_synthetic = nrow(synthetic)
    ),
    class = "synth_pair"
  )
}


#' Print method for synth_pair objects
#'
#' @param x An object of class "synth_pair"
#' @param ... Additional arguments (ignored)
#' @export
print.synth_pair <- function(x, ...) {
  cat("Synthetic Data Comparison Pair\n")
  cat("==============================\n\n")

  cat("Source:", x$source, "\n")
  cat("Original records:", x$n_original, "\n")
  cat("Synthetic records:", x$n_synthetic, "\n\n")

  cat("Variables:", length(x$vars), "total\n")
  cat("  Categorical:", length(x$cat_vars),
      if (length(x$cat_vars) > 0) paste0("(", paste(head(x$cat_vars, 3), collapse = ", "),
                                          if (length(x$cat_vars) > 3) ", ..." else "", ")") else "",
      "\n")
  cat("  Numeric:", length(x$num_vars),
      if (length(x$num_vars) > 0) paste0("(", paste(head(x$num_vars, 3), collapse = ", "),
                                          if (length(x$num_vars) > 3) ", ..." else "", ")") else "",
      "\n\n")

  if (!is.null(x$key_vars)) {
    cat("Key variables:", paste(x$key_vars, collapse = ", "), "\n")
  } else {
    cat("Key variables: (not set - required for dcap, tcap, disco)\n")
  }

  if (!is.null(x$target_var)) {
    cat("Target variable:", x$target_var, "\n")
  } else {
    cat("Target variable: (not set - required for dcap, tcap, disco)\n")
  }

  if (!is.null(x$weight_original) || !is.null(x$weight_synthetic)) {
    cat("Weights:",
        if (!is.null(x$weight_original)) paste("original =", x$weight_original) else "original = (none)",
        "|",
        if (!is.null(x$weight_synthetic)) paste("synthetic =", x$weight_synthetic) else "synthetic = (none)",
        "\n")
  }

  if (!is.null(x$holdout)) {
    cat("Holdout data:", nrow(x$holdout), "records\n")
  }

  invisible(x)
}


#' Summary method for synth_pair objects
#'
#' @param object An object of class "synth_pair"
#' @param ... Additional arguments (ignored)
#' @export
summary.synth_pair <- function(object, ...) {
  # Compute quick data quality checks
  na_original <- sum(sapply(object$original[, object$vars, drop = FALSE], function(x) sum(is.na(x))))
  na_synthetic <- sum(sapply(object$synthetic[, object$vars, drop = FALSE], function(x) sum(is.na(x))))

  summ <- list(
    source = object$source,
    n_original = object$n_original,
    n_synthetic = object$n_synthetic,
    n_vars = length(object$vars),
    n_cat_vars = length(object$cat_vars),
    n_num_vars = length(object$num_vars),
    vars = object$vars,
    cat_vars = object$cat_vars,
    num_vars = object$num_vars,
    key_vars = object$key_vars,
    target_var = object$target_var,
    has_weights = c(original = !is.null(object$weight_original),
                    synthetic = !is.null(object$weight_synthetic)),
    has_holdout = !is.null(object$holdout),
    na_count = c(original = na_original, synthetic = na_synthetic),
    metadata = object$metadata
  )

  class(summ) <- "summary.synth_pair"
  return(summ)
}


#' Print method for summary.synth_pair objects
#'
#' @param x An object of class "summary.synth_pair"
#' @param ... Additional arguments (ignored)
#' @export
print.summary.synth_pair <- function(x, ...) {
  cat("Summary: Synthetic Data Comparison Pair\n")
  cat("=======================================\n\n")

  cat("Data Dimensions:\n")
  cat("  Original:", x$n_original, "records\n")
  cat("  Synthetic:", x$n_synthetic, "records\n")
  cat("  Size ratio:", sprintf("%.2f", x$n_synthetic / x$n_original), "\n\n")

  cat("Variables (", x$n_vars, " total):\n", sep = "")
  cat("  Categorical:", x$n_cat_vars, "\n")
  if (x$n_cat_vars > 0) cat("   ", paste(x$cat_vars, collapse = ", "), "\n")
  cat("  Numeric:", x$n_num_vars, "\n")
  if (x$n_num_vars > 0) cat("   ", paste(x$num_vars, collapse = ", "), "\n")
  cat("\n")

  cat("Risk Measure Configuration:\n")
  cat("  Key variables:", if (!is.null(x$key_vars)) paste(x$key_vars, collapse = ", ") else "(not set)", "\n")
  cat("  Target variable:", if (!is.null(x$target_var)) x$target_var else "(not set)", "\n")
  cat("  Holdout data:", if (x$has_holdout) "available" else "not provided", "\n\n")

  cat("Data Quality:\n")
  cat("  Missing values (original):", x$na_count["original"], "\n")
  cat("  Missing values (synthetic):", x$na_count["synthetic"], "\n")
  cat("  Weights:", if (x$has_weights["original"]) "original" else "",
      if (x$has_weights["original"] && x$has_weights["synthetic"]) "+" else "",
      if (x$has_weights["synthetic"]) "synthetic" else "",
      if (!any(x$has_weights)) "(none)" else "", "\n")

  invisible(x)
}


#' Extract comparison data from synthpop object
#'
#' Creates a \code{\link{synth_pair}} object from a synthpop \code{synds} or
#' \code{synds.list} object.
#'
#' @param x A synthpop object of class "synds" or "synds.list"
#' @param original The original data used to create the synthetic data.
#'   This is required as synthpop does not store the original data.
#' @param key_vars Character vector of quasi-identifier variables (optional)
#' @param target_var Character string naming the sensitive target variable (optional)
#' @param m For synds.list objects with multiple syntheses, which one to use (default: 1)
#' @param ... Additional arguments passed to \code{\link{synth_pair}}
#'
#' @return An object of class "synth_pair"
#'
#' @details
#' The synthpop package creates synthetic data using sequential modeling.
#' It returns objects of class \code{synds} (single synthesis) or \code{synds.list}
#' (multiple syntheses). This function extracts the synthetic data and combines
#' it with the original data into a \code{synth_pair} for analysis.
#'
#' Note: synthpop does not store the original data in the output object, so you
#' must provide it explicitly.
#'
#' @seealso \code{\link{synth_pair}}, \code{\link{from_simPop}}
#'
#' @export
#' @examples
#' \donttest{
#' if (requireNamespace("synthpop", quietly = TRUE)) {
#'   library(synthpop)
#'
#'   # Create some example data
#'   original <- data.frame(
#'     age = sample(20:60, 100, replace = TRUE),
#'     sex = factor(sample(c("M", "F"), 100, replace = TRUE)),
#'     income = rnorm(100, 50000, 10000)
#'   )
#'
#'   # Generate synthetic data
#'   syn_obj <- syn(original, seed = 123)
#'
#'   # Create synth_pair
#'   pair <- from_synthpop(syn_obj, original,
#'                         key_vars = c("age", "sex"),
#'                         target_var = "income")
#'   print(pair)
#' }
#' }
from_synthpop <- function(x, original, key_vars = NULL, target_var = NULL, m = 1, ...) {
  # Check for synthpop classes
 if (!inherits(x, c("synds", "synds.list"))) {
    stop("x must be a 'synds' or 'synds.list' object from the synthpop package")
  }

  if (missing(original) || is.null(original)) {
    stop("'original' data must be provided (synthpop does not store it)")
  }

  # Extract synthetic data
  if (inherits(x, "synds.list")) {
    if (m > length(x$syn)) {
      stop(paste("m =", m, "but only", length(x$syn), "synthetic datasets available"))
    }
    synthetic <- x$syn[[m]]
    n_syntheses <- length(x$syn)
  } else {
    synthetic <- x$syn
    n_syntheses <- 1
  }

  # Build metadata
  meta <- list(
    n_syntheses = n_syntheses,
    selected_m = m,
    method = if (!is.null(x$method)) x$method else NA,
    seed = if (!is.null(x$seed)) x$seed else NA
  )

  synth_pair(
    original = original,
    synthetic = synthetic,
    key_vars = key_vars,
    target_var = target_var,
    source = "synthpop",
    metadata = meta,
    ...
  )
}


#' Extract comparison data from simPop object
#'
#' Creates a \code{\link{synth_pair}} object from a simPop \code{simPopObj}.
#'
#' @param x A simPop object of class "simPopObj"
#' @param key_vars Character vector of quasi-identifier variables (optional)
#' @param target_var Character string naming the sensitive target variable (optional)
#' @param use_sample_weights Logical, whether to extract and use sample weights
#'   from the simPop object (default: TRUE)
#' @param ... Additional arguments passed to \code{\link{synth_pair}}
#'
#' @return An object of class "synth_pair"
#'
#' @details
#' The simPop package generates synthetic populations from survey samples.
#' Unlike synthpop which creates synthetic samples, simPop creates a full
#' synthetic population. This has implications for interpretation:
#'
#' \itemize{
#'   \item \strong{Original}: The survey sample (typically has sampling weights)
#'   \item \strong{Synthetic}: The generated population (typically no weights, or equal weights)
#' }
#'
#' The simPop object stores both the original sample and the synthetic population,
#' so unlike \code{\link{from_synthpop}}, you don't need to provide the original data.
#'
#' @seealso \code{\link{synth_pair}}, \code{\link{from_synthpop}}
#'
#' @export
#' @examples
#' \donttest{
#' if (requireNamespace("simPop", quietly = TRUE)) {
#'   # simPop example would go here
#'   # Note: simPop has complex dependencies, example kept minimal
#' }
#' }
from_simPop <- function(x, key_vars = NULL, target_var = NULL,
                        use_sample_weights = TRUE, ...) {

  if (!inherits(x, "simPopObj")) {
    stop("x must be a 'simPopObj' from the simPop package")
  }

  # Check if simPop is available
  if (!requireNamespace("simPop", quietly = TRUE)) {
    stop("Package 'simPop' is required for this function. Please install it.",
         call. = FALSE)
  }

  # Extract sample and population using accessors (not @ slot access)
  original <- as.data.frame(simPop::samp(x))
  synthetic <- as.data.frame(simPop::pop(x))

  # Get weight variable if requested
  weight_original <- NULL
  if (use_sample_weights) {
    sample_obj <- simPop::sampleObj(x)
    weight_var <- sample_obj@weight
    if (length(weight_var) > 0 && nchar(weight_var) > 0 && weight_var %in% names(original)) {
      weight_original <- weight_var
    }
  }

  # Build metadata from the simPopObj directly
  sample_obj <- simPop::sampleObj(x)
  meta <- list(
    strata = if (length(sample_obj@strata) > 0) sample_obj@strata else NA,
    hhid = if (length(sample_obj@hhid) > 0) sample_obj@hhid else NA,
    hhsize = if (length(sample_obj@hhsize) > 0) sample_obj@hhsize else NA
  )

  synth_pair(
    original = original,
    synthetic = synthetic,
    key_vars = key_vars,
    target_var = target_var,
    weight_original = weight_original,
    weight_synthetic = NULL,  # Population typically has no weights
    source = "simPop",
    metadata = meta,
    ...
  )
}


#' Extract comparison data from sdcMicro object
#'
#' Creates a \code{\link{synth_pair}} object from an sdcMicro \code{sdcMicroObj}.
#'
#' @param x An sdcMicro object of class "sdcMicroObj"
#' @param target_var Character string naming the sensitive target variable (optional).
#'   If NULL and sensible variables are defined in the sdcMicro object, the first
#'   sensible variable is used.
#' @param use_weights Logical, whether to extract and use sampling weights
#'   from the sdcMicro object (default: TRUE)
#' @param ... Additional arguments passed to \code{\link{synth_pair}}
#'
#' @return An object of class "synth_pair"
#'
#' @details
#' The sdcMicro package provides statistical disclosure control methods for
#' microdata (local suppression, recoding, microaggregation, PRAM, noise addition,
#' etc.). The \code{sdcMicroObj} S4 class stores both the original and anonymized
#' data. This function extracts both into a \code{synth_pair} for evaluation
#' with riskutility measures.
#'
#' Note that sdcMicro produces \emph{anonymized} data (perturbation of real records),
#' not truly synthetic data. Nonetheless, many risk and utility measures apply
#' equally. In particular, attribution-based measures (CAP family), distance-based
#' measures, and all utility measures can be used to compare original vs.
#' anonymized microdata.
#'
#' Variable roles are automatically extracted from the sdcMicro object:
#' \itemize{
#'   \item \strong{key_vars}: Categorical quasi-identifiers (\code{keyVars} slot)
#'   \item \strong{target_var}: Sensitive variable (\code{sensibleVar} slot, first entry)
#'   \item \strong{weight_original}: Sampling weight (\code{weightVar} slot)
#'   \item \strong{num_vars}: Numeric key variables (\code{numVars} slot)
#' }
#'
#' The anonymized data is reconstructed via \code{sdcMicro::extractManipData()},
#' which combines manipulated key, numeric, and PRAM variables with unmodified
#' columns from the original data.
#'
#' @seealso \code{\link{synth_pair}}, \code{\link{from_synthpop}}, \code{\link{from_simPop}}
#'
#' @export
#' @examples
#' \donttest{
#' if (requireNamespace("sdcMicro", quietly = TRUE)) {
#'   # Create a simple sdcMicro object
#'   data("testdata2", package = "sdcMicro")
#'   sdc <- sdcMicro::createSdcObj(
#'     dat = testdata2,
#'     keyVars = c("urbrur", "roof", "walls", "water", "sex"),
#'     numVars = c("expend", "income", "savings"),
#'     weightVar = "sampling_weight"
#'   )
#'   # Apply some anonymization
#'   sdc <- sdcMicro::localSuppression(sdc)
#'   # Create synth_pair
#'   pair <- from_sdcMicro(sdc)
#'   print(pair)
#' }
#' }
from_sdcMicro <- function(x, target_var = NULL, use_weights = TRUE, ...) {

  if (!inherits(x, "sdcMicroObj")) {
    stop("x must be an 'sdcMicroObj' from the sdcMicro package")
  }

  if (!requireNamespace("sdcMicro", quietly = TRUE)) {
    stop("Package 'sdcMicro' is required for this function. Please install it.",
         call. = FALSE)
  }

  # Extract original data
  original <- as.data.frame(x@origData)

  # Extract anonymized data (combines manipulated + unmodified columns)
  synthetic <- as.data.frame(sdcMicro::extractManipData(x))

  # Convert key variable indices to names
  key_vars <- NULL
  if (!is.null(x@keyVars) && length(x@keyVars) > 0) {
    key_vars <- names(original)[x@keyVars]
  }

  # Convert numeric variable indices to names (stored as metadata)
  num_var_names <- NULL
  if (!is.null(x@numVars) && length(x@numVars) > 0) {
    num_var_names <- names(original)[x@numVars]
  }

  # Convert sensitive variable indices to names
  if (is.null(target_var) && !is.null(x@sensibleVar) && length(x@sensibleVar) > 0) {
    target_var <- names(original)[x@sensibleVar[1]]
  }

  # Extract weight variable name
  weight_original <- NULL
  if (use_weights && !is.null(x@weightVar) && length(x@weightVar) > 0) {
    wname <- names(original)[x@weightVar]
    if (length(wname) > 0 && nchar(wname) > 0 && wname %in% names(original)) {
      weight_original <- wname
    }
  }

  # Build metadata
  meta <- list(
    num_key_vars = num_var_names,
    hhId = if (!is.null(x@hhId) && length(x@hhId) > 0) names(original)[x@hhId] else NA,
    strataVar = if (!is.null(x@strataVar) && length(x@strataVar) > 0) names(original)[x@strataVar] else NA,
    has_local_suppression = !is.null(x@localSuppression),
    has_pram = !is.null(x@pram),
    deleted_vars = if (!is.null(x@deletedVars)) x@deletedVars else character(0)
  )

  # Combine key_vars (categorical) and num_var_names for the vars parameter
  all_key_vars <- unique(c(key_vars, num_var_names))

  synth_pair(
    original = original,
    synthetic = synthetic,
    key_vars = if (length(all_key_vars) > 0) all_key_vars else NULL,
    target_var = target_var,
    weight_original = weight_original,
    weight_synthetic = weight_original,  # sdcMicro preserves records, same weights
    source = "sdcMicro",
    metadata = meta,
    ...
  )
}
