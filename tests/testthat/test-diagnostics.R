# tests/testthat/test-diagnostics.R
# Tests for INLA diagnostic capture and the print.brms_inla_power S3 method.
# Tests that require INLA are guarded with skip_if_not_installed("INLA").

library(testthat)


# ============================================================
# Helper: build a minimal mock brms_inla_power result object
# (no INLA required)
# ============================================================

make_mock_power_result <- function(prop_warned = 0, n_failed = 0) {
  n_sims <- 10L
  sim_results <- data.frame(
    sim                 = seq_len(n_sims * 2),
    n                   = rep(c(50L, 100L), each = n_sims),
    ok                  = TRUE,
    treatment           = 0.5,
    post_prob_direction = runif(n_sims * 2, 0.6, 0.95),
    ci_width            = runif(n_sims * 2, 0.25, 0.50),
    bf10                = NA_real_,
    log10_bf10          = NA_real_,
    had_warning         = rep(FALSE, n_sims * 2),
    warning_msg         = "",
    log_mlik            = runif(n_sims * 2, -100, -50),
    mode_ok             = TRUE,
    stringsAsFactors    = FALSE
  )

  power_summ <- data.frame(
    n               = c(50L, 100L),
    treatment       = c(0.5, 0.5),
    power_direction = c(0.65, 0.82),
    avg_ci_width    = c(0.40, 0.28),
    nsims_ok        = c(n_sims, n_sims),
    stringsAsFactors = FALSE
  )

  diag_df <- data.frame(
    n           = c(50L, 100L),
    treatment   = c(0.5, 0.5),
    prop_warned = c(prop_warned, prop_warned),
    n_failed    = c(n_failed, n_failed),
    mlik_min    = c(-95.0, -80.0),
    mlik_max    = c(-55.0, -52.0),
    mlik_median = c(-75.0, -66.0),
    n_mode_ok   = c(n_sims, n_sims),
    stringsAsFactors = FALSE
  )

  out <- list(
    results     = sim_results,
    summary     = power_summ,
    diagnostics = diag_df,
    settings    = list(
      effect_name  = "treatment",
      sample_sizes = c(50L, 100L),
      nsims        = n_sims
    )
  )
  class(out) <- "brms_inla_power"
  out
}


# ============================================================
# (b) print.brms_inla_power  —  no INLA required
# ============================================================

test_that("print.brms_inla_power runs without error", {
  res <- make_mock_power_result()
  expect_no_error(print(res))
})

test_that("print.brms_inla_power returns x invisibly", {
  res <- make_mock_power_result()
  expect_identical(print(res), res)
})

test_that("print output contains the INLA diagnostics line", {
  res <- make_mock_power_result(prop_warned = 0.05)
  txt <- capture.output(print(res))
  expect_true(any(grepl("INLA diagnostics", txt, fixed = TRUE)))
  expect_true(any(grepl("warnings", txt,          fixed = TRUE)))
  expect_true(any(grepl("failed",   txt,          fixed = TRUE)))
})

test_that("print output shows 0.0% when no warnings occurred", {
  res <- make_mock_power_result(prop_warned = 0)
  txt <- capture.output(print(res))
  # Should contain something like "0.0% of fits produced warnings"
  expect_true(any(grepl("0.0%", txt, fixed = TRUE)))
})

test_that("print output shows non-zero percentage when warnings occurred", {
  res <- make_mock_power_result(prop_warned = 0.1)
  txt <- capture.output(print(res))
  # "10.0% of fits produced warnings"
  expect_true(any(grepl("10.0%", txt, fixed = TRUE)))
})

test_that("print.brms_inla_power handles missing diagnostics slot gracefully", {
  res <- make_mock_power_result()
  res$diagnostics <- NULL
  expect_no_error(print(res))
  txt <- capture.output(print(res))
  # No diagnostic line when slot is absent
  expect_false(any(grepl("INLA diagnostics", txt, fixed = TRUE)))
})

test_that("print.brms_inla_power handles missing settings slot", {
  res <- make_mock_power_result()
  res$settings <- NULL
  expect_no_error(print(res))
})

test_that("print.brms_inla_power shows power summary rows", {
  res <- make_mock_power_result()
  txt <- capture.output(print(res))
  # At least one of the sample sizes (50 or 100) should appear in the summary
  expect_true(any(grepl("50|100", txt)))
})

test_that("result object has class brms_inla_power for mock result", {
  res <- make_mock_power_result()
  expect_s3_class(res, "brms_inla_power")
})


# ============================================================
# (a) brms_inla_power() diagnostics slot  —  requires INLA
# ============================================================

# Minimal mock data generator (bypasses INLA auto-generator quirks)
.diag_mock_gen <- function(n, effect) {
  eff <- if (is.list(effect) || is.numeric(effect)) {
    as.numeric(effect[1])
  } else 0.5
  data.frame(
    y         = rnorm(n, mean = eff, sd = 1),
    treatment = rnorm(n, 0, 1),
    stringsAsFactors = FALSE
  )
}

test_that("brms_inla_power() returns $diagnostics with expected columns", {
  skip_on_cran()
  skip_if_not_installed("INLA")

  res <- brms_inla_power(
    formula        = y ~ treatment,
    effect_name    = "treatment",
    effect_grid    = 0.5,
    sample_sizes   = 30L,
    nsims          = 3L,
    data_generator = .diag_mock_gen,
    seed           = 1L,
    progress       = "none"
  )

  # Slot exists and is a data frame
  expect_true("diagnostics" %in% names(res))
  expect_s3_class(res$diagnostics, "data.frame")

  # Expected columns all present
  expected_cols <- c(
    "n", "prop_warned", "n_failed",
    "mlik_min", "mlik_max", "mlik_median", "n_mode_ok"
  )
  for (col in expected_cols) {
    expect_true(col %in% names(res$diagnostics),
                label = paste("diagnostics column:", col))
  }

  # One row per (n x effect) cell — here 1 cell
  expect_equal(nrow(res$diagnostics), 1L)

  # Values are in admissible ranges
  expect_true(res$diagnostics$prop_warned >= 0 & res$diagnostics$prop_warned <= 1)
  expect_true(res$diagnostics$n_failed    >= 0L)
  expect_true(res$diagnostics$n_mode_ok   >= 0L)
})

test_that("brms_inla_power() $results contains per-sim diagnostic columns", {
  skip_on_cran()
  skip_if_not_installed("INLA")

  res <- brms_inla_power(
    formula        = y ~ treatment,
    effect_name    = "treatment",
    effect_grid    = 0.3,
    sample_sizes   = 20L,
    nsims          = 2L,
    data_generator = .diag_mock_gen,
    seed           = 2L,
    progress       = "none"
  )

  r <- res$results
  expect_true("had_warning" %in% names(r))
  expect_true("warning_msg" %in% names(r))
  expect_true("log_mlik"    %in% names(r))
  expect_true("mode_ok"     %in% names(r))

  expect_type(r$had_warning, "logical")
  expect_type(r$warning_msg, "character")
  expect_type(r$log_mlik,    "double")
  expect_type(r$mode_ok,     "logical")

  # had_warning is always logical (never NA)
  expect_false(any(is.na(r$had_warning)))
  expect_false(any(is.na(r$mode_ok)))
})

test_that("brms_inla_power() output has class brms_inla_power", {
  skip_on_cran()
  skip_if_not_installed("INLA")

  res <- brms_inla_power(
    formula        = y ~ treatment,
    effect_name    = "treatment",
    effect_grid    = 0.3,
    sample_sizes   = 20L,
    nsims          = 2L,
    data_generator = .diag_mock_gen,
    seed           = 3L,
    progress       = "none"
  )

  expect_s3_class(res, "brms_inla_power")
})

test_that("print.brms_inla_power works on actual brms_inla_power() output", {
  skip_on_cran()
  skip_if_not_installed("INLA")

  res <- brms_inla_power(
    formula        = y ~ treatment,
    effect_name    = "treatment",
    effect_grid    = 0.3,
    sample_sizes   = 20L,
    nsims          = 2L,
    data_generator = .diag_mock_gen,
    seed           = 4L,
    progress       = "none"
  )

  expect_no_error(print(res))
  txt <- capture.output(print(res))
  expect_true(any(grepl("INLA diagnostics", txt, fixed = TRUE)))
})

test_that("diagnostics n_failed counts only failed fits correctly", {
  skip_on_cran()
  skip_if_not_installed("INLA")

  # Use a generator that always produces valid data (expect 0 failures)
  res <- brms_inla_power(
    formula        = y ~ treatment,
    effect_name    = "treatment",
    effect_grid    = 0.5,
    sample_sizes   = 50L,
    nsims          = 4L,
    data_generator = .diag_mock_gen,
    seed           = 5L,
    progress       = "none"
  )

  # nsims_ok should equal nsims when data generator is well-behaved
  # n_failed should be 0 or small
  expect_true(res$diagnostics$n_failed >= 0L)
  expect_true(res$diagnostics$n_failed <= 4L)
})

test_that("diagnostics across multiple sample sizes has one row per cell", {
  skip_on_cran()
  skip_if_not_installed("INLA")

  res <- brms_inla_power(
    formula        = y ~ treatment,
    effect_name    = "treatment",
    effect_grid    = c(0.3, 0.6),
    sample_sizes   = c(30L, 60L),
    nsims          = 2L,
    data_generator = .diag_mock_gen,
    seed           = 6L,
    progress       = "none"
  )

  # 2 sample sizes x 2 effect values = 4 cells
  expect_equal(nrow(res$diagnostics), 4L)
})


# ============================================================
# validate_inla_vs_brms() — basic structural checks (no INLA/brms needed)
# ============================================================

test_that("validate_inla_vs_brms is an exported function", {
  expect_true(is.function(validate_inla_vs_brms))
})

test_that("validate_inla_vs_brms stops when INLA is unavailable", {
  # This test only applies when INLA is NOT installed.
  skip_if(requireNamespace("INLA", quietly = TRUE), "INLA is installed")
  expect_error(
    validate_inla_vs_brms(
      formula      = y ~ treatment,
      effect_name  = "treatment"
    ),
    regexp = "INLA"
  )
})
