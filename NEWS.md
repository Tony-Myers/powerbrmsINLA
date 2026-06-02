# powerbrmsINLA 1.2.0

## Major additions

* New `compute_assurance()` function for unconditional Bayesian assurance
  (O'Hagan & Stevens, 2001) computed as a weighted average of conditional
  power over a design prior on the effect size.
* New `assurance_prior_weights()` convenience wrapper for constructing
  normalised design-prior weights (normal, uniform, beta) over an effect grid.
* New `decide_sample_size()` function with both assurance mode (design prior)
  and conditional mode for recommending sample sizes from simulation output.
* New `validate_inla_vs_brms()` function for spot-checking INLA posterior
  estimates against brms/Stan.
* Print methods for `brms_inla_power`, `powerbrmsINLA_assurance`, and
  `powerbrmsINLA_sample_size` objects.

## Plotting

* Added `plot_assurance_curve()` and `plot_assurance_with_robustness()` for
  unconditional assurance visualisation.
* Added `plot_bf_assurance_curve_smooth()`, `plot_bf_assurance_curve()`,
  `plot_bf_expected_evidence()`, and `plot_bf_heatmap()` for Bayes factor
  visualisation.
* Added `plot_decision_assurance_curve()`, `plot_decision_threshold_contour()`,
  and `add_decision_overlay()` for decision-rule visualisation.
* Added `plot_design_prior()` for visualising design priors.
* Added `plot_interaction_surface()` for multi-effect grid visualisation.
* Added `plot_power_contour()`, `plot_power_heatmap()`, and
  `plot_power_assurance_overlay()` for conditional power visualisation.
* Added `plot_precision_assurance_curve()` and `plot_precision_fan_chart()`.

## Engine improvements

* `brms_inla_power()` now supports multi-effect grids (data.frame
  `effect_grid`), brms-to-INLA prior translation with full audit trail,
  marginal-likelihood Bayes factors (`bf_method = "marglik"`), and automatic
  INLA thread detection.
* `brms_inla_power_sequential()` rewritten with multi-effect support and
  prior translation.
* `brms_inla_power_two_stage()` now uses the modernised engine internally.

## Bug fixes and housekeeping

* Removed duplicate internal function definitions (`.to_inla_family()`,
  `.scale_fill_viridis_discrete()`).
* Added `requireNamespace("MASS")` guard for negative binomial data
  generation.
* Expanded `.Rbuildignore` to exclude `.claude/`, `.DS_Store`, `.Rcheck/`,
  and `.tar.gz` artefacts.
* New test suites for assurance computation, decision helpers, prior bridge,
  and plotting functions.

# powerbrmsINLA 1.1.0

* Added `brms_inla_power_parallel()` for parallel simulations.
* Added `decide_sample_size()` and `add_decision_overlay()` helpers.
* Added new Bayes factor and precision assurance plotting functions.

# powerbrmsINLA 1.1.1

* `error_sd` and `group_sd` now accept distributional specifications (`halfnormal`, `lognormal`, `uniform`) for variance-uncertainty integration; new `validate_sd_spec()` helper exported.
* Added validation test suite (`test-validation-classical.R`, `test-validation-bayesassurance.R`) and accompanying vignette benchmarking against `power.t.test()` and `bayesassurance::assurance_nd_na()`.
* CRAN housekeeping: exclude `.github`, `LICENSE.md`, and `cran-comments.md` from the source tarball via `.Rbuildignore`.
