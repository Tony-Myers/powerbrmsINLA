test_that("global normal fixed-effect prior maps to INLA precision", {
  pri <- brms::prior(normal(0, 2), class = "b")
  mapped <- .map_brms_priors_to_inla(pri)

  expect_equal(mapped$control_fixed$mean, 0)
  expect_equal(mapped$control_fixed$prec, 1 / 2^2)
  expect_true(any(mapped$prior_audit$applied))
})

test_that("coefficient-specific normal fixed-effect prior behaviour is preserved", {
  pri <- brms::prior(normal(0, 0.5), class = "b", coef = "treatment")
  mapped <- .map_brms_priors_to_inla(pri)

  expect_equal(mapped$control_fixed$mean$treatment, 0)
  expect_equal(mapped$control_fixed$prec$treatment, 4)
})

test_that("exponential residual SD prior maps to INLA pc.prec", {
  pri <- brms::prior(exponential(1), class = "sigma")
  mapped <- .map_brms_priors_to_inla(pri, inla_family = "gaussian")

  expect_equal(mapped$control_family$hyper$prec$prior, "pc.prec")
  expect_equal(mapped$control_family$hyper$prec$param, c(-log(0.05), 0.05))
  expect_true(mapped$prior_audit$applied[1])
})

test_that("family_control precedence prevents translated sigma prior overwrite", {
  pri <- brms::prior(exponential(1), class = "sigma")
  mapped <- .map_brms_priors_to_inla(
    pri,
    inla_family = "gaussian",
    family_control_supplied = TRUE
  )

  expect_equal(mapped$control_family, list())
  expect_false(mapped$prior_audit$applied[1])
  expect_match(mapped$prior_audit$note[1], "family_control")
})

test_that("exponential group-level SD prior maps to iid random-effect pc.prec", {
  pri <- brms::prior(exponential(2), class = "sd", group = "id")
  mapped <- .map_brms_priors_to_inla(pri)

  expect_equal(mapped$hyper_by_re$id$prec$prior, "pc.prec")
  expect_equal(mapped$hyper_by_re$id$prec$param, c(-log(0.05) / 2, 0.05))

  tf <- .brms_to_inla_formula2(y ~ x + (1 | id), hyper_by_re = mapped$hyper_by_re)
  f_txt <- paste(deparse(tf$inla_formula), collapse = "")
  expect_true(grepl("pc.prec", f_txt, fixed = TRUE))
  expect_true("id" %in% tf$re_hyper_groups)
})

test_that("unsupported priors are audited as not applied", {
  pri <- brms::prior(cauchy(0, 1), class = "b", coef = "x")
  mapped <- .map_brms_priors_to_inla(pri)

  expect_true(nrow(mapped$prior_audit) >= 1)
  expect_false(mapped$prior_audit$applied[1])
  expect_equal(mapped$prior_audit$method[1], "unsupported")
})

test_that("coefficient-specific group-level SD priors are not over-applied", {
  pri <- brms::prior(exponential(2), class = "sd", group = "id", coef = "x")
  mapped <- .map_brms_priors_to_inla(pri)

  expect_equal(mapped$hyper_by_re, list())
  expect_false(mapped$prior_audit$applied[1])
  expect_equal(mapped$prior_audit$method[1], "unsupported")
})

test_that("correlated random-effect terms are audited as unsupported correlation priors", {
  mapped <- .map_brms_priors_to_inla(NULL)
  tf <- .brms_to_inla_formula2(y ~ x + (1 + x | id))
  mapped <- .audit_re_correlation_terms(mapped, tf$re_specs)

  expect_true(any(mapped$prior_audit$class == "cor"))
  expect_false(mapped$prior_audit$applied[mapped$prior_audit$class == "cor"][1])
  expect_match(mapped$prior_audit$note[mapped$prior_audit$class == "cor"][1], "correlated")
})

test_that("assurance and plotting helpers tolerate prior audit metadata", {
  syn_summary <- data.frame(
    n = rep(c(50, 100), each = 3),
    treatment = rep(c(0.2, 0.5, 0.8), 2),
    power_direction = c(0.40, 0.65, 0.85, 0.60, 0.82, 0.95),
    stringsAsFactors = FALSE
  )
  syn_result <- list(
    summary = syn_summary,
    settings = list(
      effect_name = "treatment",
      prior_translation = .empty_prior_audit()
    )
  )

  assur <- compute_assurance(
    syn_result,
    prior_weights = c("0.2" = 1 / 3, "0.5" = 1 / 3, "0.8" = 1 / 3)
  )

  expect_s3_class(assur, "powerbrmsINLA_assurance")
  expect_s3_class(plot_assurance_curve(assur), "ggplot")
  expect_s3_class(plot_power_assurance_overlay(syn_result, assur), "ggplot")
})
