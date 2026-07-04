# tests/testthat/test-sequential-analysis.R
# Tests for sequential_design(), sequential_analysis(), and
# plot_sequential_monitor() (new in 1.3.0).

library(testthat)

make_design <- function(...) {
  sequential_design(
    formula      = y ~ treatment,
    effect_name  = "treatment",
    metric       = "direction",
    looks        = c(20, 40, 60),
    prob_success = 0.95,
    prob_futility = 0.05,
    ...
  )
}

# ==========================================================================
# (a) Design constructor and fingerprint (no INLA required)
# ==========================================================================

test_that("sequential_design validates its inputs", {
  expect_error(
    sequential_design(formula = y ~ a + b, effect_name = c("a", "b"),
                      looks = c(20, 40)),
    "length 1"
  )
  expect_error(
    sequential_design(formula = y ~ x, effect_name = "x",
                      looks = c(40, 20)),
    "strictly increasing"
  )
  expect_error(
    sequential_design(formula = y ~ x, effect_name = "x",
                      looks = c(20, 40), prob_success = 0.4),
    "prob_success"
  )
  expect_error(
    sequential_design(formula = y ~ x, effect_name = "x",
                      looks = c(20, 40), prob_futility = 0.7),
    "prob_futility"
  )
  expect_error(
    sequential_design(formula = y ~ x, effect_name = "x",
                      looks = c(20, 40), metric = "rope"),
    "rope_bounds"
  )
})

test_that("identical designs share a fingerprint; changed rules do not", {
  d1 <- make_design()
  d2 <- make_design()
  expect_identical(d1$fingerprint, d2$fingerprint)

  d3 <- sequential_design(
    formula = y ~ treatment, effect_name = "treatment",
    metric = "direction", looks = c(20, 40, 60),
    prob_success = 0.975,        # changed rule
    prob_futility = 0.05
  )
  expect_false(identical(d1$fingerprint, d3$fingerprint))

  # Per-look thresholds are part of the fingerprint
  d4 <- sequential_design(
    formula = y ~ treatment, effect_name = "treatment",
    metric = "direction", looks = c(20, 40, 60),
    prob_success = c(0.99, 0.97, 0.95), prob_futility = 0.05
  )
  expect_false(identical(d1$fingerprint, d4$fingerprint))
})

test_that("design print method reports the fingerprint", {
  d <- make_design(label = "Test study")
  expect_s3_class(d, "powerbrmsINLA_seq_design")
  expect_output(print(d), "Fingerprint")
  expect_output(print(d), "Test study")
})

# ==========================================================================
# (b) sequential_analysis() guards (no INLA required)
# ==========================================================================

test_that("sequential_analysis rejects invalid inputs", {
  d <- make_design()

  expect_error(sequential_analysis(list(), data.frame(y = 1)),
               "sequential_design")

  # Missing model variables are caught before any fitting
  expect_error(
    sequential_analysis(d, data.frame(z = rnorm(20))),
    "missing variables"
  )
})

test_that("a tampered design is refused", {
  d <- make_design()
  d$prob_success <- rep(0.6, 3)   # tampering after creation
  expect_error(
    sequential_analysis(d, data.frame(y = rnorm(20),
                                      treatment = rbinom(20, 1, 0.5))),
    "fingerprint mismatch"
  )
})

test_that("analysis after a recorded stop is refused without override", {
  d <- make_design()
  fake_monitor <- structure(
    list(
      design      = d,
      fingerprint = d$fingerprint,
      looks       = tibble::tibble(
        look = 1L, time = "t", n = 20L, planned_n = 20L,
        est_mean = 0.6, est_sd = 0.2, ci_lower = 0.2, ci_upper = 1.0,
        pr = 0.99, pr_in_rope = NA_real_,
        threshold_success = 0.95, threshold_futility = 0.05,
        decision = "stop_success", deviations = "", note = ""
      ),
      status = "stop_success"
    ),
    class = c("powerbrmsINLA_seq_monitor", "list")
  )
  dat <- data.frame(y = rnorm(40), treatment = rbinom(40, 1, 0.5))
  expect_error(
    sequential_analysis(fake_monitor, dat),
    "recorded the decision"
  )
})

# ==========================================================================
# (c) Integration: monitoring simulated 'real' data (INLA required)
# ==========================================================================

test_that("a full monitoring cycle works on accumulating data", {
  skip_if_not_installed("INLA")
  skip_on_cran()

  set.seed(7)
  d <- make_design(label = "Integration test")

  # Simulate a 'real' study with a large true effect
  n_max <- 60L
  treatment <- rbinom(n_max, 1, 0.5)
  y <- 1.2 * treatment + rnorm(n_max)
  full <- data.frame(y = y, treatment = treatment)

  mon1 <- sequential_analysis(d, full[1:20, ])
  expect_s3_class(mon1, "powerbrmsINLA_seq_monitor")
  expect_equal(nrow(mon1$looks), 1L)
  expect_true(mon1$looks$pr[1] >= 0 && mon1$looks$pr[1] <= 1)
  expect_true(mon1$status %in%
                c("ongoing", "stop_success", "stop_futility"))

  if (identical(mon1$status, "ongoing")) {
    mon2 <- sequential_analysis(mon1, full[1:40, ])
    expect_equal(nrow(mon2$looks), 2L)
    # Data must accumulate: re-analysing fewer rows is an error
    expect_error(sequential_analysis(mon2, full[1:30, ]), "ALL accumulated")
    mon_final <- mon2
  } else {
    mon_final <- mon1
  }

  expect_output(print(mon_final), "monitor")
  p1 <- plot_sequential_monitor(mon_final, type = "probability")
  p2 <- plot_sequential_monitor(mon_final, type = "estimate")
  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
})

test_that("the trial simulator accepts a design object", {
  skip_if_not_installed("INLA")
  skip_on_cran()

  d <- make_design()
  res <- brms_inla_sequential_trial(
    design      = d,
    true_effect = 0.8,
    nsims       = 10,
    seed        = 11,
    progress    = "none"
  )
  expect_s3_class(res, "powerbrmsINLA_seq_trial")
  expect_identical(res$settings$looks, d$looks)
  expect_identical(res$settings$metric, d$metric)
})

test_that("the trial simulator requires either a design or core arguments", {
  expect_error(brms_inla_sequential_trial(nsims = 5), "Supply either")
})
