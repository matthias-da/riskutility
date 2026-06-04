## Submission

This is the initial CRAN submission of riskutility (version 0.1.0).

The package provides disclosure-risk and data-utility metrics for anonymized
and synthetic data (attribution-based CAP/TCAP/WEAP/DiSCO and RAPID risk
measures, distance-based DCR/NNDR/IMS, classical SDC privacy models, and a
range of utility measures), together with a multivariate Risk-Utility map.

## R CMD check results

0 errors | 0 warnings | 0 notes

## Test environments

* local macOS, R 4.5.x: `R CMD check --as-cran` (with manual)
* win-builder (devel and release)  <!-- run before submitting -->
* R-hub: Windows, Linux, macOS      <!-- run before submitting -->

## Notes for the CRAN team

* All examples that rely on optional packages or are slow are wrapped in
  `\donttest{}` and seeded with `set.seed()`.
* The package suggests several optional backends (e.g. 'torch', 'xgboost',
  'sdcMicro', 'synthpop'); their use is guarded with `requireNamespace()` and
  example/vignette chunks are conditionally evaluated.

## Downstream dependencies

There are currently no reverse dependencies (new package).
