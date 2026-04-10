# tests/testthat/test-variance-uncertainty.R
# Tests for distributional error_sd / group_sd specifications and the
# validate_sd_spec() helper.  Tests that require INLA are guarded with
# skip_on_cran() + skip_if_not_installed("INLA").

library(testthat)


# ============================================================
# validate_sd_spec() — no INLA required
# ============================================================

test_that("validate_sd_spec accepts a positive numeric scalar", {
  expect_no_error(validate_sd_spec(1.0))
  expect_no_error(validate_sd_spec(0.01))
  expect_no_error(validate_sd_spec(100))
})

test_that("validate_sd_spec rejects non-positive or non-finite scalars", {
  expect_error(validate_sd_spec(0),     regexp = "positive")
  expect_error(validate_sd_spec(-1),    regexp = "positive")
  expect_error(validate_sd_spec(Inf),   regexp = "positive|finite")
  expect_error(validate_sd_spec(NA_real_), regexp = "positive|finite")
})

test_that("validate_sd_spec rejects a character string", {
  expect_error(validate_sd_spec("1.0"), regexp = "scalar or a list")
})

test_that("validate_sd_spec accepts a valid halfnormal spec", {
  spec <- list(dist = "halfnormal", sd = 0.3, location = 1.0)
  expect_no_error(validate_sd_spec(spec))
})

test_that("validate_sd_spec accepts halfnormal without location (uses default 0)", {
  expect_no_error(validate_sd_spec(list(dist = "halfnormal", sd = 0.5)))
})

test_that("validate_sd_spec rejects halfnormal with missing or non-positive sd", {
  expect_error(validate_sd_spec(list(dist = "halfnormal")),       regexp = "sd > 0")
  expect_error(validate_sd_spec(list(dist = "halfnormal", sd = 0)),  regexp = "sd > 0")
  expect_error(validate_sd_spec(list(dist = "halfnormal", sd = -1)), regexp = "sd > 0")
})

test_that("validate_sd_spec accepts a valid lognormal spec", {
  expect_no_error(validate_sd_spec(list(dist = "lognormal", meanlog = 0, sdlog = 0.3)))
})

test_that("validate_sd_spec rejects lognormal with missing or bad sdlog", {
  expect_error(validate_sd_spec(list(dist = "lognormal", meanlog = 0)),        regexp = "sdlog")
  expect_error(validate_sd_spec(list(dist = "lognormal", meanlog = 0, sdlog = 0)),  regexp = "sdlog")
  expect_error(validate_sd_spec(list(dist = "lognormal", meanlog = 0, sdlog = -1)), regexp = "sdlog")
})

test_that("validate_sd_spec accepts a valid uniform spec", {
  expect_no_error(validate_sd_spec(list(dist = "uniform", min = 0.5, max = 2.0)))
})

test_that("validate_sd_spec rejects uniform with min < 0 or max <= min", {
  expect_error(validate_sd_spec(list(dist = "uniform", min = -1, max = 1)),  regexp = "min >= 0")
  expect_error(validate_sd_spec(list(dist = "uniform", min = 1,  max = 0.5)), regexp = "max > min")
  expect_error(validate_sd_spec(list(dist = "uniform", min = 1,  max = 1)),   regexp = "max > min")
})

test_that("validate_sd_spec rejects an unsupported dist name", {
  expect_error(
    validate_sd_spec(list(dist = "gamma", shape = 2, rate = 1)),
    regexp = "unsupported"
  )
})

test_that("validate_sd_spec rejects a list without a dist element", {
  expect_error(validate_sd_spec(list(sd = 1)), regexp = "dist")
})

test_that("validate_sd_spec returns x invisibly on success", {
  spec <- list(dist = "halfnormal", sd = 0.5)
  out  <- validate_sd_spec(spec)
  expect_identical(out, spec)
})

test_that("validate_sd_spec arg_name appears in error messages", {
  expect_error(validate_sd_spec(-1, arg_name = "error_sd"), regexp = "error_sd")
  expect_error(validate_sd_spec(list(dist = "gamma"), arg_name = "group_sd"), regexp = "group_sd")
})


# ============================================================
# .sample_sd_spec() draws — no INLA required
# ============================================================

test_that("halfnormal draws are always positive", {
  set.seed(1L)
  spec <- list(dist = "halfnormal", sd = 1.0, location = 0)
  draws <- replicate(200L, powerbrmsINLA:::.sample_sd_spec(spec))
  expect_true(all(draws > 0))
})

test_that("halfnormal with non-zero location still returns positive values", {
  set.seed(2L)
  spec <- list(dist = "halfnormal", sd = 0.3, location = 1.0)
  draws <- replicate(200L, powerbrmsINLA:::.sample_sd_spec(spec))
  expect_true(all(draws >= 0))
  # With location 1 and sd 0.3, virtually all values should be > 0.5
  expect_true(mean(draws > 0.5) > 0.9)
})

test_that("lognormal draws are always positive", {
  set.seed(3L)
  spec <- list(dist = "lognormal", meanlog = 0, sdlog = 0.3)
  draws <- replicate(200L, powerbrmsINLA:::.sample_sd_spec(spec))
  expect_true(all(draws > 0))
})

test_that("uniform draws stay within [min, max]", {
  set.seed(4L)
  spec <- list(dist = "uniform", min = 0.5, max = 2.0)
  draws <- replicate(200L, powerbrmsINLA:::.sample_sd_spec(spec))
  expect_true(all(draws >= 0.5 & draws <= 2.0))
})


# ============================================================
# brms_inla_power() integration — requires INLA
# ============================================================

# Minimal mock data generator (bypasses auto-generator quirks)
.var_mock_gen <- function(n, effect) {
  eff <- as.numeric(effect[[1L]])
  data.frame(
    y         = rnorm(n, mean = eff, sd = 1),
    treatment = rnorm(n, 0, 1),
    stringsAsFactors = FALSE
  )
}

test_that("(a) scalar error_sd still works identically", {
  skip_on_cran()
  skip_if_not_installed("INLA")

  res <- brms_inla_power(
    formula        = y ~ treatment,
    effect_name    = "treatment",
    effect_grid    = 0.5,
    sample_sizes   = 30L,
    nsims          = 3L,
    data_generator = .var_mock_gen,
    error_sd       = 1.0,
    seed           = 10L,
    progress       = "none"
  )

  expect_s3_class(res, "brms_inla_power")
  expect_true(nrow(res$results) > 0L)
  # sampled_error_sd should be all NA for scalar input
  expect_true(all(is.na(res$results$sampled_error_sd)))
})

test_that("(b) distributional error_sd produces sampled_error_sd column", {
  skip_on_cran()
  skip_if_not_installed("INLA")

  res <- brms_inla_power(
    formula        = y ~ treatment,
    effect_name    = "treatment",
    effect_grid    = 0.5,
    sample_sizes   = 30L,
    nsims          = 4L,
    data_generator = .var_mock_gen,
    error_sd       = list(dist = "halfnormal", sd = 0.3, location = 1.0),
    seed           = 20L,
    progress       = "none"
  )

  expect_true("sampled_error_sd" %in% names(res$results))
  ok_rows <- res$results[res$results$ok, ]
  expect_false(any(is.na(ok_rows$sampled_error_sd)))
  expect_true(all(ok_rows$sampled_error_sd > 0))
})

test_that("(c) halfnormal sampled_error_sd values are always positive", {
  skip_on_cran()
  skip_if_not_installed("INLA")

  res <- brms_inla_power(
    formula        = y ~ treatment,
    effect_name    = "treatment",
    effect_grid    = 0.5,
    sample_sizes   = 30L,
    nsims          = 6L,
    data_generator = .var_mock_gen,
    error_sd       = list(dist = "halfnormal", sd = 0.5),
    seed           = 30L,
    progress       = "none"
  )

  drawn <- res$results$sampled_error_sd
  drawn <- drawn[!is.na(drawn)]
  expect_true(length(drawn) > 0L)
  expect_true(all(drawn > 0))
})

test_that("(d) invalid error_sd spec produces an informative error before any simulation", {
  # No INLA needed — validate_sd_spec() runs before any INLA calls
  expect_error(
    brms_inla_power(
      formula      = y ~ treatment,
      effect_name  = "treatment",
      effect_grid  = 0.5,
      sample_sizes = 30L,
      nsims        = 2L,
      error_sd     = list(dist = "gamma", shape = 2),
      seed         = 40L,
      progress     = "none"
    ),
    regexp = "unsupported"
  )
})

test_that("(d) missing sd in halfnormal spec errors informatively", {
  expect_error(
    brms_inla_power(
      formula      = y ~ treatment,
      effect_name  = "treatment",
      effect_grid  = 0.5,
      sample_sizes = 30L,
      nsims        = 2L,
      error_sd     = list(dist = "halfnormal"),
      seed         = 41L,
      progress     = "none"
    ),
    regexp = "sd > 0"
  )
})

test_that("(d) negative scalar error_sd errors", {
  expect_error(
    brms_inla_power(
      formula      = y ~ treatment,
      effect_name  = "treatment",
      effect_grid  = 0.5,
      sample_sizes = 30L,
      nsims        = 2L,
      error_sd     = -1,
      seed         = 42L,
      progress     = "none"
    ),
    regexp = "positive"
  )
})

test_that("(e) lognormal error_sd spec works without error", {
  skip_on_cran()
  skip_if_not_installed("INLA")

  expect_no_error(
    brms_inla_power(
      formula        = y ~ treatment,
      effect_name    = "treatment",
      effect_grid    = 0.5,
      sample_sizes   = 30L,
      nsims          = 3L,
      data_generator = .var_mock_gen,
      error_sd       = list(dist = "lognormal", meanlog = 0, sdlog = 0.2),
      seed           = 50L,
      progress       = "none"
    )
  )
})

test_that("(e) uniform error_sd spec works without error", {
  skip_on_cran()
  skip_if_not_installed("INLA")

  expect_no_error(
    brms_inla_power(
      formula        = y ~ treatment,
      effect_name    = "treatment",
      effect_grid    = 0.5,
      sample_sizes   = 30L,
      nsims          = 3L,
      data_generator = .var_mock_gen,
      error_sd       = list(dist = "uniform", min = 0.5, max = 1.5),
      seed           = 60L,
      progress       = "none"
    )
  )
})

test_that("summary reports mean_sampled_error_sd when distributional", {
  skip_on_cran()
  skip_if_not_installed("INLA")

  res <- brms_inla_power(
    formula        = y ~ treatment,
    effect_name    = "treatment",
    effect_grid    = 0.5,
    sample_sizes   = 30L,
    nsims          = 4L,
    data_generator = .var_mock_gen,
    error_sd       = list(dist = "halfnormal", sd = 0.3, location = 1.0),
    seed           = 70L,
    progress       = "none"
  )

  expect_true("mean_sampled_error_sd" %in% names(res$summary))
  expect_true("sd_sampled_error_sd"   %in% names(res$summary))
  expect_true(is.finite(res$summary$mean_sampled_error_sd))
  expect_true(res$summary$mean_sampled_error_sd > 0)
})

test_that("summary mean_sampled_error_sd is NA when scalar", {
  skip_on_cran()
  skip_if_not_installed("INLA")

  res <- brms_inla_power(
    formula        = y ~ treatment,
    effect_name    = "treatment",
    effect_grid    = 0.5,
    sample_sizes   = 30L,
    nsims          = 3L,
    data_generator = .var_mock_gen,
    error_sd       = 1.0,
    seed           = 80L,
    progress       = "none"
  )

  expect_true(is.na(res$summary$mean_sampled_error_sd))
  expect_true(is.na(res$summary$sd_sampled_error_sd))
})

test_that("settings stores original error_sd and group_sd specs", {
  skip_on_cran()
  skip_if_not_installed("INLA")

  spec <- list(dist = "halfnormal", sd = 0.3, location = 1.0)
  res <- brms_inla_power(
    formula        = y ~ treatment,
    effect_name    = "treatment",
    effect_grid    = 0.5,
    sample_sizes   = 30L,
    nsims          = 2L,
    data_generator = .var_mock_gen,
    error_sd       = spec,
    seed           = 90L,
    progress       = "none"
  )

  expect_identical(res$settings$error_sd, spec)
})
