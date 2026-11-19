## Submission

This is a feature release of riskutility (version 0.3.0; version 0.2.0 was
published on 2026-09-02).

The RAPID functions (`rapid()` and its companion functions) are rebuilt on the
measure of the paper that introduces RAPID (Templ, Thees and Müller, 2026,
under review). Several defaults change reported values. NEWS.md lists every
such change together with the call that reproduces the 0.2.0 result, and the
help page `?rapid` has a section "Changes in 0.3.0". A second vignette,
`rapid`, documents the measure and the result object.

## R CMD check results

0 errors | 0 warnings | 0 notes

## Test environments

* local macOS, R 4.5.x: `R CMD check --as-cran` (with manual)
* win-builder (devel and release)  <!-- run before submitting -->
* R-hub: Windows, Linux, macOS      <!-- run before submitting -->

## Notes for the CRAN team

* Examples that rely on optional packages or are slow are wrapped in
  `\donttest{}`.
* The package suggests several optional backends (e.g. 'torch', 'xgboost',
  'rpart', 'sdcMicro', 'synthpop'); their use is guarded with
  `requireNamespace()` and example/vignette chunks are conditionally evaluated.

## Downstream dependencies

There are currently no reverse dependencies.
