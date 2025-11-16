## powerbrmsINLA 1.1.0

### Summary of changes

* Added `brms_inla_power_parallel()` to provide a parallel wrapper around existing simulation engines (fixed and sequential designs), preserving current behaviour while reducing wall-clock time.
* Extended plotting functionality:
  * `plot_bf_assurance_curve_smooth()` for Bayes factor assurance with Wilson intervals.
  * `plot_bf_expected_evidence()` and `plot_bf_heatmap()` for expected log10 BF visualisation across effect grids and sample sizes.
  * `plot_precision_fan_chart()` as a convenience wrapper for robustness/precision plots.
  * `add_decision_overlay()` to overlay sample-size decisions on assurance / power plots.
* Refined `decide_sample_size()` to allow flexible Bayes factor cutoffs and clearer rationales.
* Minor internal clean-up of plotting helpers and decision helpers; no breaking changes to existing exported interfaces.

### R CMD check results

* Local (macOS Sequoia 15.4.1, R 4.5.0):

  * `R CMD check --as-cran` via `devtools::check()`:  
    `0 errors | 0 warnings | 0 notes`

* Windows (R-devel) via win-builder:

  * R Under development (unstable) (2025-11-14 r89021 ucrt):  
    `0 errors | 0 warnings | 0 notes`

### R-hub

* R-hub v2 GitHub Actions workflow was run on:
  * linux (R-devel), macos-arm64 (R-devel), windows (R-devel).

* These jobs failed during dependency resolution because `INLA (>= 22.05.07)` is not on CRAN and the default pak configuration could not see the additional INLA repository:

  * “Could not solve package dependencies: Can't install dependency INLA (>= 22.05.07); INLA: Can't find package called INLA.”

* On win-builder and local checks, `INLA` is correctly found via the `Additional_repositories` field in `DESCRIPTION`, and R CMD check completes without issues.
