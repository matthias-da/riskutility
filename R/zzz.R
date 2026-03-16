.onLoad <- function(libname, pkgname) {
  # Suppress notes about global variables used in data.table/ggplot NSE
}

# Global variables used in non-standard evaluation contexts
utils::globalVariables(c(
  # sdcMicro functions
  "createSdcObj", "pram", "microaggregation", "extractManipData", "pop",
  # General
  "value", "what", ".",
  # data.table and ggplot2 NSE
  ".dataset", "weight", "facet_cat", "facet_group", "group_key", "..vars",
  "..cat_vars", "..cont_vars",
  "V1", "V2", "dataset", "Dim1", "Dim2",
  ".var", ".name", ".label", "x_val", "y_val",
  # Additional variables used in package
  "variable", "freq", "pct", "cum_freq", "relative_freq", "rel_freq",
  "stat", "X_value", "Y_value",
  "PC1", "PC2", "xend", "yend", "label",
  "train_percentile", "holdout_percentile",
  "i", "j", "n", "count", "Category",
  # density ratio variables
  "density_ratios", "density_ratios_bayes",
  # dcap/tcap data.table variables
  "key_sig", "target", "n_total", "orig_idx", "n_correct", "unique_syn", "unique_orig",
  # rapid plot variables
  "at_risk",
  # torch nn_module self reference
  "self"
))

#' @importFrom grDevices rgb
#' @importFrom stats IQR aov cor cov cov2cor ecdf ks.test mad mahalanobis
#' @importFrom stats qchisq setNames var weighted.mean
NULL
