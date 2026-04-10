# powerbrmsINLA 1.1.0

* Added `brms_inla_power_parallel()` for parallel simulations.
* Added `decide_sample_size()` and `add_decision_overlay()` helpers.
* Added new Bayes factor and precision assurance plotting functions.

# powerbrmsINLA 1.1.1

* `error_sd` and `group_sd` now accept distributional specifications (`halfnormal`, `lognormal`, `uniform`) for variance-uncertainty integration; new `validate_sd_spec()` helper exported.
* Added validation test suite (`test-validation-classical.R`, `test-validation-bayesassurance.R`) and accompanying vignette benchmarking against `power.t.test()` and `bayesassurance::assurance_nd_na()`.
* CRAN housekeeping: exclude `.github`, `LICENSE.md`, and `cran-comments.md` from the source tarball via `.Rbuildignore`.
