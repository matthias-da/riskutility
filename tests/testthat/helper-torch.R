# Skip a test unless the torch C++ backend (libtorch / Lantern) is available.
#
# `skip_if_not_installed("torch")` only checks that the *R package* is
# installed. On CRAN / win-builder the torch package is present but its
# backend is not (it is downloaded separately via torch::install_torch()),
# so any tensor operation fails with "Lantern is not loaded". Backend-dependent
# tests must therefore gate on the backend itself, not on the package.
skip_if_no_torch <- function() {
  testthat::skip_if_not_installed("torch")
  backend_ok <- tryCatch(isTRUE(torch::torch_is_installed()),
                         error = function(e) FALSE)
  if (!backend_ok) {
    testthat::skip("torch backend (libtorch / Lantern) not installed")
  }
}
