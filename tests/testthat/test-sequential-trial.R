# tests/testthat/test-sequential-trial.R
# Tests for brms_inla_sequential_trial() and .sample_design_prior()

library(testthat)

# ==========================================================================
# (a) Design-prior sampler (no INLA required)
# ==========================================================================

test_that(".sample_design_prior draws from the requested distribution", {
  set.seed(1)
  x <- powerbrmsINLA:::.sample_design_prior(
    list(dist = "normal", mean = 0.5, sd = 0.1), n = 5000L
  )
  expect_length(x, 5000L)
  expect_equal(mean(x), 0.5, tolerance = 0.01)
  expect_equal(sd(x), 0.1, tolerance = 0.01)

  u <- powerbrmsINLA:::.sample_design_prior(
    list(dist = "uniform", min = 0.2, max = 0.8), n = 1000L
  )
  expect_true(all(u >= 0.2 & u <= 0.8))

  b <- powerbrmsINLA:::.sample_design_prior(
    list(dist = "beta", shape1 = 2, shape2 = 5), n = 1000L
  )
  expect_true(all(b > 0 & b < 1))
})

test_that(".sample_design_prior rejects bad specifications", {
  expect_error(powerbrmsINLA:::.sample_design_prior(list(mean = 0.5)),
               "dist")
  expect_error(powerbrmsINLA:::.sample_design_prior(
    list(dist = "normal", mean = 0, sd = -1)), "sd > 0")
  expect_error(powerbrmsINLA:::.sample_design_prior(
    list(dist = "uniform", min = 1, max = 0)), "min < max")
  expect_error(powerbrmsINLA:::.sample_design_prior(
    list(dist = "cauchy", location = 0)), "Unsupported")
})

# ==========================================================================
# (b) Argument validation (runs before the INLA dependency check)
# ==========================================================================

test_that("invalid looks are rejected", {
  expect_error(
    brms_inla_sequential_trial(
      formula = y ~ treatment, effect_name = "treatment",
      looks = c(100, 50), nsims = 2
    ),
    "strictly increasing"
  )
  expect_error(
    brms_inla_sequential_trial(
      formula = y ~ treatment, effect_name = "treatment",
      looks = c(-10, 50), nsims = 2
    ),
    "strictly increasing|positive"
  )
})

test_that("multiple effect names are rejected", {
  expect_error(
    brms_inla_sequential_trial(
      formula = y ~ a + b, effect_name = c("a", "b"),
      looks = c(50, 100), nsims = 2
    ),
    "length 1"
  )
})

test_that("threshold vectors are validated", {
  expect_error(
    brms_inla_sequential_trial(
      formula = y ~ treatment, effect_name = "treatment",
      looks = c(50, 100), prob_success = 0.4, nsims = 2
    ),
    "prob_success"
  )
  expect_error(
    brms_inla_sequential_trial(
      formula = y ~ treatment, effect_name = "treatment",
      looks = c(50, 100), prob_futility = 0.6, nsims = 2
    ),
    "prob_futility"
  )
  expect_error(
    brms_inla_sequential_trial(
      formula = y ~ treatment, effect_name = "treatment",
      looks = c(50, 100), prob_success = c(0.99, 0.97, 0.95), nsims = 2
    ),
    "prob_success"
  )
})

test_that("rope metric requires rope_bounds", {
  expect_error(
    brms_inla_sequential_trial(
      formula = y ~ treatment, effect_name = "treatment",
      looks = c(50, 100), metric = "rope", nsims = 2
    ),
    "rope_bounds"
  )
})

test_that("bad true_effect specifications are rejected", {
  expect_error(
    brms_inla_sequential_trial(
      formula = y ~ treatment, effect_name = "treatment",
      looks = c(50, 100), true_effect = "big", nsims = 2
    ),
    "true_effect"
  )
  expect_error(
    brms_inla_sequential_trial(
      formula = y ~ treatment, effect_name = "treatment",
      looks = c(50, 100),
      true_effect = list(dist = "cauchy", location = 0), nsims = 2
    ),
    "Unsupported"
  )
})

# ==========================================================================
# (c) Integration tests (INLA required)
# ==========================================================================

test_that("sequential trial returns a well-formed object (fixed effects)", {
  skip_if_not_installed("INLA")
  skip_on_cran()

  res <- brms_inla_sequential_trial(
    formula          = y ~ treatment,
    effect_name      = "treatment",
    true_effect      = c(0, 0.8),
    looks            = c(30, 60),
    nsims            = 20,
    metric           = "direction",
    prob_success     = 0.95,
    prob_futility    = 0.05,
    seed             = 42,
    progress         = "none"
  )

  expect_s3_class(res, "powerbrmsINLA_seq_trial")
  expect_named(res, c("results", "summary", "look_summary",
                      "diagnostics", "settings"))
  expect_equal(nrow(res$results), 40L)
  expect_setequal(unique(res$results$scenario), c("0", "0.8"))
  expect_true(all(res$results$decision %in%
                    c("success", "futility", "inconclusive", "failed")))

  s <- res$summary
  expect_true(all(c("p_success", "p_futility", "p_inconclusive",
                    "expected_n", "p_stop_early") %in% names(s)))
  # Probabilities sum to one within each scenario
  expect_equal(s$p_success + s$p_futility + s$p_equivalence + s$p_inconclusive,
               rep(1, nrow(s)), tolerance = 1e-12)
  # Expected n bounded by the look schedule
  expect_true(all(s$expected_n >= 30 & s$expected_n <= 60))
  # A large true effect should succeed more often than a null effect
  expect_gt(s$p_success[s$scenario == "0.8"],
            s$p_success[s$scenario == "0"])
  # Regression guard: mean_true_at_success must reflect the data column,
  # not a masked summary scalar (was NA for all scenarios before the fix)
  s8 <- s[s$scenario == "0.8", ]
  if (s8$p_success > 0) {
    expect_equal(s8$mean_true_at_success, 0.8, tolerance = 1e-12)
  }
  expect_equal(s$true_effect[s$scenario == "0.8"], 0.8, tolerance = 1e-12)

  expect_output(print(res), "Sequential Bayesian trial")
})

test_that("design-prior mode reports sequential assurance", {
  skip_if_not_installed("INLA")
  skip_on_cran()

  res <- brms_inla_sequential_trial(
    formula      = y ~ treatment,
    effect_name  = "treatment",
    true_effect  = list(dist = "normal", mean = 0.5, sd = 0.15),
    looks        = c(30, 60),
    nsims        = 10,
    metric       = "direction",
    seed         = 42,
    progress     = "none"
  )

  expect_s3_class(res, "powerbrmsINLA_seq_trial")
  expect_equal(unique(res$results$scenario), "design_prior")
  # True effects vary across trials in design-prior mode
  expect_gt(stats::sd(res$results$true_effect), 0)
  expect_output(print(res), "ASSURANCE")
})

test_that("rope metric can stop for equivalence", {
  skip_if_not_installed("INLA")
  skip_on_cran()

  res <- brms_inla_sequential_trial(
    formula      = y ~ treatment,
    effect_name  = "treatment",
    true_effect  = 0,
    looks        = c(60, 120),
    nsims        = 15,
    metric       = "rope",
    rope_bounds  = c(-0.3, 0.3),
    prob_success = 0.90,
    seed         = 42,
    progress     = "none"
  )

  expect_true(all(res$results$decision %in%
                    c("success", "equivalence", "inconclusive", "failed")))
})
