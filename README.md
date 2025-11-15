
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
library(powerbrmsINLA)

# Run Bayesian power analysis
results <- brms_inla_power(
  formula = outcome ~ treatment,
  effect_name = "treatment", 
  effect_grid = c(0.2, 0.5, 0.8),
  sample_sizes = c(50, 100),
  nsims = 5  # Reduced for speed
)

# Inspect summary results
results$summary

# Plot power heatmap  
plot_power_heatmap(results)
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
