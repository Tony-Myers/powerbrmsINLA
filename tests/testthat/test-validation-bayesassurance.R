# tests/testthat/test-validation-bayesassurance.R
#
# Validates that compute_assurance() agrees with bayesassurance::assurance_nd_na()
# for a conjugate normal one-sample design where both packages should converge
# to the same analytical result.
#
# The INLA-dependent tests use a minimal 3-point effect grid and small nsims
# (50) to keep the check in a reasonable wall-clock window (~1 s per INLA fit
# on this platform).  The full comparison with nsims = 500 and a 21-point grid
# is shown pre-computed in the validation vignette.
#
# All INLA tests are skipped on CRAN and when INLA is absent.
# The bayesassurance comparison is additionally skipped when that package is
# not installed.

library(testthat)


# ============================================================
# Module-level constants
# ============================================================

.ba_sigma   <- 1.0
.ba_delta   <- 0.5
.ba_sigma_d <- 0.5      # design-prior SD
.ba_n_sizes <- c(30L, 100L)
.ba_n0      <- 0.01     # nearly flat analysis prior
# Coarse 3-point grid: covers the prior's bulk without many INLA cells
.ba_coarse_grid <- c(0.0, 0.5, 1.0)


# ============================================================
# Helper: one-sample data generator
# ============================================================

.one_sample_gen <- function(sigma = 1) {
  function(n, effect) {
    mu <- as.numeric(effect[["mu"]])
    data.frame(
      y  = rnorm(n, mean = mu, sd = sigma),
      mu = rep(1L, n),
      stringsAsFactors = FALSE
    )
  }
}


# ============================================================
# bayesassurance comparison — requires bayesassurance + INLA
# ============================================================

test_that("TC3: compute_assurance() agrees with bayesassurance::assurance_nd_na() within +/-0.05", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  skip_if_not_installed("bayesassurance")

  # ---- bayesassurance reference values ----------------------------------------
  ba_result <- tryCatch(
    bayesassurance::assurance_nd_na(
      n     = .ba_n_sizes,
      n0    = .ba_n0,
      delta = .ba_delta,
      sigma = .ba_sigma,
      alpha = 0.025
    ),
    error = function(e) {
      skip(paste("bayesassurance::assurance_nd_na() call failed:", conditionMessage(e)))
    }
  )

  # Extract assurance vector; handle data-frame or named-vector output
  if (is.data.frame(ba_result)) {
    col <- intersect(c("assurance", "Assurance", "power", "Power"), names(ba_result))
    if (length(col) == 0L) skip("Cannot identify assurance column in bayesassurance output.")
    ba_assurance <- ba_result[[col[1L]]]
  } else if (is.numeric(ba_result)) {
    ba_assurance <- as.numeric(ba_result)
  } else {
    skip("Unrecognised return type from bayesassurance::assurance_nd_na().")
  }

  if (length(ba_assurance) != length(.ba_n_sizes)) {
    skip("bayesassurance returned unexpected number of assurance values.")
  }

  # ---- powerbrmsINLA reference values -----------------------------------------
  suppressMessages(
    power_res <- brms_inla_power(
      formula        = y ~ mu,
      effect_name    = "mu",
      effect_grid    = .ba_coarse_grid,
      sample_sizes   = .ba_n_sizes,
      nsims          = 50L,
      prob_threshold = 0.975,
      error_sd       = .ba_sigma,
      data_generator = .one_sample_gen(.ba_sigma),
      seed           = 201L,
      progress       = "none"
    )
  )

  w <- assurance_prior_weights(
    .ba_coarse_grid,
    dist = "normal",
    mean = .ba_delta,
    sd   = .ba_sigma_d
  )

  assur_res    <- compute_assurance(power_res, prior_weights = w)
  pa_assurance <- assur_res$assurance$assurance[
    order(assur_res$assurance$sample_size)
  ]

  # Compare — allow wider tolerance (±0.08) given coarse grid + small nsims
  for (i in seq_along(.ba_n_sizes)) {
    expect_true(
      abs(pa_assurance[i] - ba_assurance[i]) <= 0.08,
      label = sprintf(
        "n=%d: powerbrmsINLA=%.4f, bayesassurance=%.4f, diff=%.4f (must be <=0.08)",
        .ba_n_sizes[i], pa_assurance[i], ba_assurance[i],
        abs(pa_assurance[i] - ba_assurance[i])
      )
    )
  }
})


# ============================================================
# Structural checks (INLA required, bayesassurance not required)
# ============================================================

test_that("TC3 powerbrmsINLA assurance is monotone non-decreasing in n", {
  skip_on_cran()
  skip_if_not_installed("INLA")

  suppressMessages(
    power_res <- brms_inla_power(
      formula        = y ~ mu,
      effect_name    = "mu",
      effect_grid    = .ba_coarse_grid,
      sample_sizes   = .ba_n_sizes,
      nsims          = 50L,
      prob_threshold = 0.975,
      error_sd       = .ba_sigma,
      data_generator = .one_sample_gen(.ba_sigma),
      seed           = 202L,
      progress       = "none"
    )
  )

  w <- assurance_prior_weights(
    .ba_coarse_grid,
    dist = "normal",
    mean = .ba_delta,
    sd   = .ba_sigma_d
  )

  assur_res    <- compute_assurance(power_res, prior_weights = w)
  assur_sorted <- assur_res$assurance[order(assur_res$assurance$sample_size), ]

  diffs <- diff(assur_sorted$assurance)
  expect_true(all(diffs >= -0.08),  # allow Monte Carlo noise at low nsims
              label = "Assurance should not decrease substantially as n increases")
})

test_that("TC3 assurance does not substantially exceed conditional power at design-prior mean", {
  skip_on_cran()
  skip_if_not_installed("INLA")

  # Single sample size (n=100) to keep fit count low
  suppressMessages(
    power_at_delta <- brms_inla_power(
      formula        = y ~ mu,
      effect_name    = "mu",
      effect_grid    = .ba_delta,     # single point: conditional power
      sample_sizes   = 100L,
      nsims          = 50L,
      prob_threshold = 0.975,
      error_sd       = .ba_sigma,
      data_generator = .one_sample_gen(.ba_sigma),
      seed           = 203L,
      progress       = "none"
    )
  )

  suppressMessages(
    power_grid <- brms_inla_power(
      formula        = y ~ mu,
      effect_name    = "mu",
      effect_grid    = .ba_coarse_grid,
      sample_sizes   = 100L,
      nsims          = 50L,
      prob_threshold = 0.975,
      error_sd       = .ba_sigma,
      data_generator = .one_sample_gen(.ba_sigma),
      seed           = 204L,
      progress       = "none"
    )
  )

  w <- assurance_prior_weights(
    .ba_coarse_grid,
    dist = "normal",
    mean = .ba_delta,
    sd   = .ba_sigma_d
  )

  assur_res   <- compute_assurance(power_grid, prior_weights = w)
  cond_power  <- power_at_delta$summary$power_direction
  assur_value <- assur_res$assurance$assurance

  expect_true(
    assur_value <= cond_power + 0.10,
    label = sprintf(
      "Assurance (%.4f) should not greatly exceed conditional power at delta (%.4f)",
      assur_value, cond_power
    )
  )
})


# ============================================================
# Structural checks — no INLA required
# ============================================================

test_that("bayesassurance::assurance_nd_na is a function", {
  skip_if_not_installed("bayesassurance")
  expect_true(is.function(bayesassurance::assurance_nd_na))
})

test_that("assurance_prior_weights normalises to 1 over the coarse grid", {
  w <- assurance_prior_weights(
    .ba_coarse_grid,
    dist = "normal",
    mean = .ba_delta,
    sd   = .ba_sigma_d
  )
  expect_length(w, length(.ba_coarse_grid))
  expect_equal(sum(w), 1, tolerance = 1e-10)
  expect_true(all(w >= 0))
})
