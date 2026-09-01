## Submission

This is a maintenance and feature release of riskutility (version 0.2.0;
first release 0.1.0 published 2026-06-22).

It contains bug fixes surfaced while preparing a software paper (see NEWS.md;
three fixes change reported values), one API addition to `tcap()`, and one
change made in coordination with the synthpop maintainers and the CRAN team:
the forthcoming synthpop 1.9-3 redefines its TCAP disclosure measure, which a
riskutility test compared against. `tcap()` now implements both definitions
explicitly (`tcap_matched` for synthpop <= 1.9-2, `tcap_conditional` for the
Little et al. 2025 definition used from synthpop 1.9-3), and all tests that
compare against synthpop output are version-aware and skipped on CRAN, so a
future synthpop release cannot be blocked by riskutility's checks again.

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

There are currently no reverse dependencies.
