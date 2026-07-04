# powerbrmsINLA

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/powerbrmsINLA)](https://cran.r-project.org/package=powerbrmsINLA)
<!-- badges: end -->

## Overview

**powerbrmsINLA** provides tools for **Bayesian power analysis** and
**assurance calculations** using the statistical frameworks of
[`brms`](https://cran.r-project.org/package=brms) and
[`INLA`](https://www.r-inla.org/).

It includes simulation-based approaches, support for multiple decision
rules (direction, threshold, ROPE, Bayes factors, precision), sequential
and two-stage adaptive designs, and a comprehensive suite of
visualisation functions.

## What's New in 1.3.0

- **Sequential Bayesian analysis module**: `sequential_design()` for
  prespecifying a sequential analysis (with an MD5 fingerprint of all
  decision-relevant fields for preregistration), `sequential_analysis()`
  for interim monitoring of real accumulating data with an auditable
  decision trail, `plot_sequential_monitor()` for trajectory plots, and
  `brms_inla_sequential_trial()` for simulating a sequential design's
  operating characteristics (stopping probabilities, expected sample size,
  early-stop exaggeration).
- **Validation hardening**: `brms_inla_power()` now raises an error when
  `effect_name` does not match a formula-level fixed-effect term and the
  built-in data generator is in use (previously such a name was silently
  ignored).
- `decide_sample_size()` in conditional mode now requires at least one
  decision target, and no longer mistakes the per-cell SD-moment summary
  columns for effect-grid columns.
- **Breaking change**: in `brms_inla_power_sequential()` summaries, the
  column previously named `assurance` is now `conditional_power` (the old
  name was statistically misleading).
- Bug fix: `brms_inla_power_two_stage()` no longer errors when called with
  default `error_sd` / `obs_per_group`.
- All engines now fail early with an informative message when INLA is
  absent.

## What's New in 1.2.0

- **Unconditional Bayesian assurance** via `compute_assurance()` — averages
  conditional power over a design prior on the effect size (O'Hagan &
  Stevens, 2001).
- **`assurance_prior_weights()`** for constructing normalised design-prior
  weights (normal, uniform, beta) over an effect grid.
- **`decide_sample_size()`** with both assurance mode (design prior) and
  conditional mode for recommending sample sizes from simulation output.
- **`validate_inla_vs_brms()`** for spot-checking INLA posterior estimates
  against brms/Stan.
- **brms-to-INLA prior translation** with full audit trail — specify
  analysis priors using `brms::prior()` syntax.
- **Variance uncertainty integration** — `error_sd` and `group_sd` now
  accept distributional specifications (`halfnormal`, `lognormal`,
  `uniform`) so that power is integrated over variance uncertainty.
- **Marginal-likelihood Bayes factors** (`bf_method = "marglik"`) alongside
  the existing Savage-Dickey method.
- **Automatic INLA thread detection** when `inla_num_threads = NULL`.
- **15 new plotting functions** for assurance, Bayes factor, decision-rule,
  precision, and multi-effect visualisation.
- Print methods for `brms_inla_power`, `powerbrmsINLA_assurance`, and
  `powerbrmsINLA_sample_size` objects.

See `NEWS.md` for the full changelog.

## Installation

Install from CRAN:

``` r
install.packages("powerbrmsINLA")
```

INLA is listed under `Suggests` and must be installed separately:

``` r
if (!requireNamespace("INLA", quietly = TRUE)) {
  install.packages(
    "INLA",
    repos = c(getOption("repos"),
              INLA = "https://inla.r-inla-download.org/R/stable"),
    dep = TRUE
  )
}
```

To install the development version from GitHub:

``` r
# install.packages("remotes")
remotes::install_github("Tony-Myers/powerbrmsINLA")
```

## Quick Example

``` r
library(powerbrmsINLA)

# Step 1: Conditional power simulation
results <- brms_inla_power(
  formula      = y ~ treatment,
  effect_name  = "treatment",
  effect_grid  = c(0.2, 0.5, 0.8),
  sample_sizes = c(50, 100),
  nsims        = 50,
  seed         = 123
)
results$summary

# Step 2: Unconditional assurance (new in 1.2.0)
assurance <- compute_assurance(
  results,
  prior_weights = list(dist = "normal", mean = 0.5, sd = 0.2),
  metric = "direction"
)
print(assurance)

# Step 3: Sample size recommendation
decide_sample_size(
  results,
  direction = 0.80,
  prior_weights = list(dist = "normal", mean = 0.5, sd = 0.2)
)
```

## Model Complexity Considerations

For optimal performance:

- **Simple to moderate models**: All sample sizes supported.
- **Complex random effects** (e.g., `(1 + time | subject)`): Recommend
  n >= 50 subjects.
- **Large effect grids**: Consider starting with fewer simulations
  (nsims = 50-100) for initial exploration, or use the sequential/two-stage
  engines.

## Citation

If you use powerbrmsINLA in published work, please cite:

> Myers, T. (2026). powerbrmsINLA: Bayesian Power Analysis Using 'brms'
> and 'INLA'. R package version 1.2.0.
> https://cran.r-project.org/package=powerbrmsINLA

## License

This package is released under the MIT License.
See the [LICENSE](LICENSE) file for details.
