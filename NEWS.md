# powerbrmsINLA 1.3.0

## Bug fixes and validation

* `brms_inla_power()` now raises a hard error (rather than a warning) when
  `effect_name` does not match a formula-level fixed-effect term and the
  automatic data generator is in use. Previously such a name was silently
  ignored when building the linear predictor, so the requested effect was never
  applied to the simulated data even though its grid value was recorded. The
  error message clarifies that the formula term should be used (e.g.
  `effect_name = "Condition"`), not a fitted coefficient level (e.g.
  `"Condition1"`). When a custom `data_generator` is supplied the check is
  downgraded to a warning, since the naming convention is then under user
  control. This validation now runs before the INLA dependency check.
* `decide_sample_size()` in conditional mode now requires at least one decision
  target (`direction`, `threshold`, `rope_in`, or `bf10`); supplying none
  previously produced an arbitrary smallest-n recommendation. Assurance mode
  already enforced this.
* `decide_sample_size()` no longer mistakes the per-cell SD-moment summary
  columns (`mean_sampled_error_sd`, `sd_sampled_error_sd`,
  `mean_sampled_group_sd`, `sd_sampled_group_sd`) for effect-grid/design
  columns. Only `sampled_*` was excluded previously, which left the
  `mean_sampled_*` / `sd_sampled_*` columns to fragment the effect grouping.
* Added INLA-independent regression tests covering all three fixes.

## Sequential Bayesian analysis (new)

* New `sequential_design()` for prespecifying a sequential Bayesian analysis:
  model, analysis priors, monitored effect, decision metric (direction /
  threshold / ROPE), per-look success and futility thresholds, and the
  planned look schedule. The object carries an MD5 fingerprint of all
  decision-relevant fields for quotation in preregistrations; any change to
  the rules changes the fingerprint.
* New `sequential_analysis()` for monitoring a real study as data accumulate:
  fits the prespecified model with INLA at each interim look, applies the
  prespecified stopping rule, and records an auditable trajectory (estimates,
  credible intervals, monitored probabilities, decisions, and any protocol
  deviations). Refuses to continue after a recorded stop unless explicitly
  overridden, in which case the continuation is logged as a deviation.
* New `plot_sequential_monitor()` for trajectory plots of the monitored
  posterior probability against the prespecified stopping boundaries, or of
  the effect estimate with its credible interval across looks.
* New `brms_inla_sequential_trial()` simulates the operating characteristics
  of a sequential design before data collection: probabilities of stopping
  for success / futility / equivalence, expected sample size, the per-look
  stopping distribution, and the exaggeration of effect estimates at early
  stops. True effects may be fixed scenarios or drawn from a design prior
  (sequential assurance). Accepts `sequential_design()` objects so the
  simulated design is exactly the design that is later monitored.

## Breaking changes

* In `brms_inla_power_sequential()` summaries, the column previously named
  `assurance` is now `conditional_power`. The old name was misleading: the
  quantity is conditional power at a fixed effect-grid value, not
  unconditional assurance (which `compute_assurance()` provides).

## Bug fixes

* `brms_inla_power_two_stage()` no longer errors when called with default
  arguments (`error_sd` and `obs_per_group` previously defaulted to `NULL`,
  which the input validator rejects); defaults are now 1 and 10, matching
  `brms_inla_power()`.
* `brms_inla_power_sequential()` now returns an object of class
  `"brms_inla_power"`, so the print method applies; it also validates
  `effect_name` and `target` on entry.
* Two test files (`sequencial-test.R`, `set-up.R`) did not match testthat's
  file-name conventions and were silently never executed; they have been
  renamed (`test-sequential-stopping.R`, `setup.R`), updated, and now run.

## Other improvements

* All simulation engines now fail early with an informative message when
  INLA is not installed, instead of failing inside the simulation loop.
* `.plot_decision_assurance_curve_from_summary()` is no longer exported (it
  is an internal helper).
* The sequential engine's default seed is now 123, matching the other
  engines.
* Fixed the internal `.geom_point_lw()` helper, which passed `linewidth`
  to `ggplot2::geom_point()` on ggplot2 >= 3.4; the parameter was silently
  ignored ("Ignoring unknown parameters" warning) and points rendered at
  their default size. Point sizing now uses `size`, as geom_point expects.
* Removed duplicated internal plotting helpers (single copies retained in
  `utils-helpers.R`); pruned stray `globalVariables` entries; large PNG/PDF
  artefacts at the package root are excluded from builds via `.Rbuildignore`.

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
