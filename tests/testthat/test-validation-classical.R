# tests/testthat/test-validation-classical.R
#
# Validates that brms_inla_power() with a direction-based Bayesian decision
# rule (prob_threshold = 0.975) and vague priors reproduces classical
# one-sided power (alpha = 0.025) from power.t.test() to within ±0.05.
#
# nsims is kept at 500 (rather than the 2000 recommended for publication-
# quality validation) so that each INLA fit on this machine (~0.25-1 s/fit)
# keeps the test suite in a reasonable wall-clock window.  The Monte Carlo
# SE at nsims = 500 is ≈ 0.018 for p ≈ 0.8, giving a ±3σ window of ±0.054,
# which is tight enough for the ±0.05 tolerance.  Use the validation
# vignette with nsims = 2000 for publication-quality comparisons.
#
# Interpretation of sample sizes
# --------------------------------
# The auto data generator creates a binary 'treatment' predictor sampled
# uniformly from {0, 1}, so n observations give approximately n/2 per group.
# To match power.t.test(n = n_per_group), brms_inla_power() must receive
# sample_sizes = 2 * n_per_group (total N).  This is explicit via the
# `n_per_group` and `n_total` variables.
#
# All tests are skipped on CRAN and when INLA is not installed.

library(testthat)


# ============================================================
# Helper: custom data generator for two-sample Gaussian design
# ============================================================

.make_two_sample_gen <- function(error_sd) {
  function(n, effect) {
    delta <- as.numeric(effect[[1L]])
    treat <- rbinom(n, 1L, 0.5)
    y     <- delta * treat + rnorm(n, 0, error_sd)
    data.frame(treatment = treat, y = y, stringsAsFactors = FALSE)
  }
}


# ============================================================
# Test case 1: Gaussian two-sample, delta = 0.5, SD = 1
# n = 64 per group (n_total = 128)
# Classical benchmark: power.t.test(n=64, ...) = 0.8015
# ============================================================

test_that("TC1: direction power matches classical power.t.test within +/-0.05 (n=64/grp, delta=0.5)", {
  skip_on_cran()
  skip_if_not_installed("INLA")

  n_per_group <- 64L
  n_total     <- 2L * n_per_group  # 128 total observations
  delta       <- 0.5
  sigma       <- 1.0

  # Classical one-sided power at alpha = 0.025
  classical <- power.t.test(
    n           = n_per_group,
    delta       = delta,
    sd          = sigma,
    sig.level   = 0.025,
    type        = "two.sample",
    alternative = "one.sided"
  )$power

  suppressMessages(
    res <- brms_inla_power(
      formula        = y ~ treatment,
      effect_name    = "treatment",
      effect_grid    = delta,
      sample_sizes   = n_total,
      nsims          = 500L,
      prob_threshold = 0.975,
      error_sd       = sigma,
      data_generator = .make_two_sample_gen(sigma),
      seed           = 101L,
      progress       = "none"
    )
  )

  bayesian <- res$summary$power_direction

  expect_length(bayesian, 1L)
  expect_true(
    abs(bayesian - classical) <= 0.05,
    label = sprintf(
      "Bayesian power (%.4f) vs classical (%.4f): diff = %.4f, exceeds 0.05",
      bayesian, classical, abs(bayesian - classical)
    )
  )
})


# ============================================================
# Test case 2: Chen et al. (2018) motivational interviewing
# delta = 2.26, pooled SD = 6.536
# n per group = c(60, 100, 133, 200) → n_total = c(120, 200, 266, 400)
# ============================================================

test_that("TC2: direction power matches Chen et al. classical benchmark within +/-0.05", {
  skip_on_cran()
  skip_if_not_installed("INLA")

  n_per_group <- c(60L, 100L, 133L, 200L)
  n_total     <- 2L * n_per_group
  delta       <- 2.26
  sigma       <- 6.536

  # Classical one-sided powers at alpha = 0.025 for each per-group n
  classical_powers <- vapply(n_per_group, function(ng) {
    power.t.test(
      n           = ng,
      delta       = delta,
      sd          = sigma,
      sig.level   = 0.025,
      type        = "two.sample",
      alternative = "one.sided"
    )$power
  }, numeric(1L))

  suppressMessages(
    res <- brms_inla_power(
      formula        = y ~ treatment,
      effect_name    = "treatment",
      effect_grid    = delta,
      sample_sizes   = n_total,
      nsims          = 500L,
      prob_threshold = 0.975,
      error_sd       = sigma,
      data_generator = .make_two_sample_gen(sigma),
      seed           = 102L,
      progress       = "none"
    )
  )

  summ <- res$summary[order(res$summary$n), ]

  expect_equal(nrow(summ), length(n_per_group))

  diffs <- abs(summ$power_direction - classical_powers)

  for (i in seq_along(n_per_group)) {
    expect_true(
      diffs[i] <= 0.05,
      label = sprintf(
        "n_per_group=%d: Bayesian=%.4f, classical=%.4f, diff=%.4f (must be <=0.05)",
        n_per_group[i], summ$power_direction[i], classical_powers[i], diffs[i]
      )
    )
  }
})


# ============================================================
# Structural / sanity checks (no INLA required)
# ============================================================

test_that("Classical power for TC1 matches known value", {
  p <- power.t.test(
    n = 64, delta = 0.5, sd = 1, sig.level = 0.025,
    type = "two.sample", alternative = "one.sided"
  )$power
  expect_true(p > 0 & p < 1)
  expect_equal(round(p, 2L), 0.80)
})

test_that("Classical powers for TC2 match Chen et al. (2018) Table 1", {
  published <- c(n20  = 0.186,
                 n100 = 0.682,
                 n133 = 0.800,
                 n200 = 0.932)
  for (nm in names(published)) {
    ng <- as.integer(sub("n", "", nm))
    cp <- power.t.test(
      n = ng, delta = 2.26, sd = 6.536, sig.level = 0.025,
      type = "two.sample", alternative = "one.sided"
    )$power
    expect_equal(round(cp, 3L), unname(published[nm]),
                 tolerance = 0.003,
                 label = paste("Published power for n =", ng))
  }
})
