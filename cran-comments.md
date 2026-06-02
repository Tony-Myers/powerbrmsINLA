## powerbrmsINLA 1.2.0

### Summary of changes

* New `compute_assurance()` function for unconditional Bayesian assurance
  (O'Hagan & Stevens, 2001) computed as a weighted average of conditional
  power over a design prior on the effect size.
* New `assurance_prior_weights()` convenience wrapper for constructing
  normalised design-prior weights.
* New `decide_sample_size()` with both assurance-mode (design prior) and
  conditional mode for recommending sample sizes from simulation output.
* New `validate_inla_vs_brms()` for spot-checking INLA posterior estimates
  against brms/Stan.
* `brms_inla_power()` now supports multi-effect grids, brms-to-INLA prior
  translation with full audit trail, marginal-likelihood Bayes factors, and
  automatic INLA thread detection.
* `error_sd` and `group_sd` accept distributional specifications
  (`halfnormal`, `lognormal`, `uniform`) for variance-uncertainty integration.
* 15 new plotting functions for assurance, Bayes factor, decision-rule,
  precision, and multi-effect visualisation.
* Print methods for `brms_inla_power`, `powerbrmsINLA_assurance`, and
  `powerbrmsINLA_sample_size` objects.

### Downstream dependencies

None (checked via `revdepcheck`).

### R CMD check results

* Local (macOS, R 4.5.x):

  * `R CMD check --as-cran` via `devtools::check()`:
    `0 errors | 0 warnings | 0 notes`

* Windows (R-devel) via win-builder:

  * R Under development (unstable) (2026-06-01 r90092 ucrt):
    `0 errors | 0 warnings | 1 note`
  * The single NOTE is the expected "Suggests or Enhances not in mainstream
    repositories: INLA, bayesassurance" (see below).

* INLA is listed in `Suggests` (not `Imports`) and is available via the
  `Additional_repositories` field in DESCRIPTION.  All INLA-dependent code
  paths are guarded by `requireNamespace(“INLA”, quietly = TRUE)`.

### Note on bayesassurance

`bayesassurance` is listed in `Suggests` for an optional validation vignette
only.  It is not required by any exported function.  The package was not
available on the check platform but this does not affect results.
