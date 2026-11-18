# RAPID: conditions, argument handling and data preparation

# ---- Conditions --------------------------------------------------------------
# Every error and warning raised by RAPID code goes through these two functions,
# so that callers can catch them by class.

#' Signal a classed RAPID error
#'
#' @param message Character; pasted together without separator.
#' @param class Suffix of the condition class (`"input"`, `"model"`,
#'   `"levels"`, `"tolerance"`, `"dependency"`, `"object"`).
#' @noRd
.rapid_abort <- function(message, class = NULL) {
  stop(structure(
    class = c(if (!is.null(class)) paste0("rapid_error_", class),
              "rapid_error", "error", "condition"),
    list(message = paste(message, collapse = ""), call = NULL)
  ))
}

#' Signal a classed RAPID warning
#'
#' @param message Character; pasted together without separator.
#' @param class Suffix of the condition class (`"deprecated"`,
#'   `"ignored_argument"`, `"legacy_tolerance"`, `"skipped_model"`).
#' @noRd
.rapid_warn <- function(message, class = NULL) {
  warning(structure(
    class = c(if (!is.null(class)) paste0("rapid_warning_", class),
              "rapid_warning", "warning", "condition"),
    list(message = paste(message, collapse = ""), call = NULL)
  ))
}
