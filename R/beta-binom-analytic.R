#' Highest Density Interval from an Inverse CDF
#'
#' Computes an HDI of given mass from any distribution for which you have a
#' quantile function (inverse CDF).
#'
#' @param qfun Quantile function, e.g., \code{qbeta}, \code{qnorm}, ...
#' @param width Desired HDI mass (e.g., 0.95).
#' @param tol Optimizer tolerance.
#' @param ... Additional arguments passed to \code{qfun}.
#'
#' @return Named numeric vector with elements \code{ll} and \code{ul}.
#' @export
hdi_of_icdf <- function(qfun, width = 0.95, tol = 1e-8, ...) {
  stopifnot(is.function(qfun), width > 0, width < 1)
  incredible <- 1 - width
  wfun <- function(low, qfun, width, ...) qfun(low + width, ...) - qfun(low, ...)
  opt <- optimize(wfun, c(0, incredible), qfun = qfun, width = width, tol = tol, ...)
  low <- opt$minimum
  c(ll = qfun(low, ...), ul = qfun(low + width, ...))
}

#' Analytic Assurance for Beta-Binomial Designs
#'
#' Computes assurance (power) using generating and audience Beta priors for
#' a binomial count via a Beta-Binomial predictive distribution.
#'
#' @param n Sample size (number of trials).
#' @param gen_prior_a,gen_prior_b Generating Beta prior parameters.
#' @param aud_prior_a,aud_prior_b Audience Beta prior parameters.
#' @param hdi_mass HDI mass (e.g., 0.95).
#' @param rope Length-2 numeric vector for ROPE bounds, or NULL for max-width rule.
#' @param hdi_max_width Positive width threshold for the HDI (used if \code{rope=NULL}).
#'
#' @return Assurance value between 0 and 1.
#' @export
beta_binom_power <- function(n,
                             gen_prior_a, gen_prior_b,
                             aud_prior_a, aud_prior_b,
                             hdi_mass = 0.95,
                             rope = NULL,
                             hdi_max_width = NULL) {
  stopifnot(n >= 0, xor(is.null(rope), is.null(hdi_max_width)))
  stopifnot(all(c(gen_prior_a, gen_prior_b, aud_prior_a, aud_prior_b) > 0))

  z <- 0:n
  # Predictive mass for z under generating prior (Beta-Binomial)
  log_pz <- lchoose(n, z) +
    lbeta(z + gen_prior_a, n - z + gen_prior_b) -
    lbeta(gen_prior_a, gen_prior_b)
  log_pz <- log_pz - max(log_pz)
  pz <- exp(log_pz)
  pz <- pz / sum(pz)

  # Posterior HDI for each z under audience prior
  hdi <- t(vapply(
    z,
    function(zz)
      hdi_of_icdf(qbeta, width = hdi_mass,
                  shape1 = zz + aud_prior_a,
                  shape2 = n - zz + aud_prior_b),
    numeric(2L)
  ))

  pass <- if (!is.null(rope)) {
    (hdi[, 1] > rope[2]) | (hdi[, 2] < rope[1])
  } else {
    (hdi[, 2] - hdi[, 1]) < hdi_max_width
  }

  sum(pz[pass])
}

#' Minimum n for Target Assurance (Beta-Binomial)
#'
#' @inheritParams beta_binom_power
#' @param gen_prior_mode Generating prior mode in (0,1).
#' @param gen_prior_n Generating prior concentration (>= 2).
#' @param desired_power Target assurance value in (0,1).
#' @param aud_prior_mode Audience prior mode in (0,1).
#' @param aud_prior_n Audience prior concentration (>= 2).
#' @param n_start Starting sample size for search.
#' @param n_max Maximum sample size to try.
#' @param verbose If TRUE, prints progress.
#'
#' @return Smallest n meeting the target assurance.
#' @export
min_n_beta_binom <- function(gen_prior_mode, gen_prior_n,
                             desired_power,
                             aud_prior_mode = 0.5, aud_prior_n = 2,
                             hdi_mass = 0.95,
                             rope = NULL, hdi_max_width = NULL,
                             n_start = 20, n_max = 1e5,
                             verbose = TRUE) {
  stopifnot(xor(is.null(rope), is.null(hdi_max_width)))
  stopifnot(gen_prior_mode > 0, gen_prior_mode < 1, gen_prior_n >= 2)
  stopifnot(aud_prior_mode > 0, aud_prior_mode < 1, aud_prior_n >= 2)
  stopifnot(desired_power > 0, desired_power < 1)

  ag <- gen_prior_mode * (gen_prior_n - 2) + 1
  bg <- (1 - gen_prior_mode) * (gen_prior_n - 2) + 1
  aa <- aud_prior_mode * (aud_prior_n - 2) + 1
  ba <- (1 - aud_prior_mode) * (aud_prior_n - 2) + 1

  n <- as.integer(n_start)
  repeat {
    pwr <- beta_binom_power(n, ag, bg, aa, ba,
                            hdi_mass = hdi_mass,
                            rope = rope,
                            hdi_max_width = hdi_max_width)
    if (verbose) cat("n =", n, ", assurance =", round(pwr, 6), "\n")
    if (pwr >= desired_power || n >= n_max) return(n)
    n <- n + 1L
  }
}

#' Beta-Prior Weights Over an Effect Grid
#'
#' Computes prior weights over a grid of true effect values by evaluating
#' a Beta(mode, n) prior. If the grid is not in (0,1), it is rescaled linearly.
#'
#' @param effects Numeric vector of effect values (grid).
#' @param mode Prior mode in (0,1).
#' @param n Prior concentration (> 2).
#'
#' @return Normalised numeric weights over the grid (sum to 1).
#' @export
beta_weights_on_grid <- function(effects, mode, n) {
  stopifnot(is.numeric(effects), length(effects) >= 2, mode > 0, mode < 1, n > 2)
  a <- mode * (n - 2) + 1
  b <- (1 - mode) * (n - 2) + 1
  rng <- range(effects)
  x <- effects
  if (rng[1] < 0 || rng[2] > 1) {
    if (diff(rng) == 0) stop("All 'effects' are identical; cannot rescale.")
    message("beta_weights_on_grid(): effects not in [0,1]; rescaling linearly.")
    x <- (effects - rng[1]) / (rng[2] - rng[1])
  }
  w <- dbeta(x, a, b)
  w <- w / sum(w)
  names(w) <- as.character(effects)
  w
}
