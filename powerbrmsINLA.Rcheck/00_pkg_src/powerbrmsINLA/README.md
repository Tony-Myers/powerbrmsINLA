
<!-- README.md is generated from README.Rmd. Please edit README.Rmd -->

# powerbrmsINLA

## Overview

**powerbrmsINLA** provides tools for **Bayesian power analysis** and
**assurance calculations** using the statistical frameworks of
[`brms`](https://cran.r-project.org/package=brms) and
[`INLA`](https://www.r-inla.org/).

It includes simulation-based and analytical approaches, support for
multiple decision rules (`direction`, `threshold`, `rope`), sequential
and two-stage designs, and visualisation helpers for power curves,
precision, Bayes factors, and robustness.

<!-- badges: start -->
  [![R-CMD-check](https://github.com/Tony-Myers/powerbrmsINLA/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Tony-Myers/powerbrmsINLA/actions/workflows/R-CMD-check.yaml)
  <!-- badges: end -->

## Installation

You can install the development version from GitHub:

``` r
# install.packages("remotes")
remotes::install_github("https://github.com/Tony-Myers/powerbrmsINLA")
```

## Example

Here is a minimal example to get started. For speed in a README, the
code is not evaluated on knit.

``` r
library(brmsINLApower)

set.seed(1)

# Simple dataset generator using .auto_data_generator()
gen_fun <- brmsINLApower:::.auto_data_generator(
  formula     = y ~ x,
  effect_name = "x",
  family      = gaussian()
)

# Run Bayesian power analysis with reduced settings for speed
res <- brms_inla_power(
  formula        = y ~ x,
  effect_name    = "x",
  effect_grid    = 0.5,
  sample_sizes   = c(20, 40),
  nsims          = 5,
  data_generator = gen_fun
)

# Inspect summary results
res$summary

# Plot a power heatmap
plot_power_heatmap(res)
```

## Model Complexity Considerations

For optimal performance:

- **Simple to moderate models**: All sample sizes supported
- **Complex random effects** (e.g., `(1 + time | subject)`): Recommend n ≥ 50 subjects
- **Large effect grids**: Consider starting with fewer simulations (nsims = 50-100) for initial exploration

The package handles the vast majority of Bayesian power analysis scenarios. For computationally demanding models, standard Bayesian modeling best practices apply (adequate sample sizes, model complexity appropriate to data).

## Package documentation

If you use [`pkgdown`](https://pkgdown.r-lib.org/) you can build a
website:

``` r
usethis::use_pkgdown()           # once, to set up pkgdown
pkgdown::build_site()            # build the site locally
# usethis::use_pkgdown_github_pages()  # set up GitHub Pages
```

## License

This package is released under the MIT License.  
See the [LICENSE](LICENSE) file for details.
