#' Print method for brms_inla_power result objects
#'
#' Displays a concise summary of a simulation result, including key settings,
#' the first rows of the power summary table, and — when present — a one-line
#' INLA diagnostic summary stating the proportion of fits that produced
#' warnings and the number of failed fits.
#'
#' @param x A list of class `"brms_inla_power"` as returned by
#'   [brms_inla_power()], [brms_inla_power_parallel()], or
#'   [brms_inla_power_two_stage()].
#' @param n_rows Maximum number of summary rows to print (default 10).
#' @param ... Unused; present for S3 compatibility.
#'
#' @return `x`, invisibly.
#' @export
print.brms_inla_power <- function(x, n_rows = 10L, ...) {
  # Settings: present in main/parallel engine as $settings, in two-stage as $params
  cfg <- x$settings %||% x$params %||% list()

  cat("Bayesian power / assurance simulation (powerbrmsINLA)\n")
  cat("======================================================\n")

  if (!is.null(cfg$effect_name)) {
    cat("Effect(s)   :", paste(cfg$effect_name, collapse = ", "), "\n")
  }
  if (!is.null(x$summary) && "n" %in% names(x$summary)) {
    cat("Sample sizes:", paste(sort(unique(x$summary$n)), collapse = ", "), "\n")
  }
  if (!is.null(cfg$nsims)) {
    cat("Simulations :", cfg$nsims, "per cell\n")
  }

  # ---- Diagnostic one-liner -----------------------------------------------
  if (!is.null(x$diagnostics) && nrow(x$diagnostics) > 0L) {
    pct_warned <- mean(x$diagnostics$prop_warned, na.rm = TRUE) * 100
    n_fail     <- sum(x$diagnostics$n_failed,    na.rm = TRUE)
    cat(sprintf(
      "INLA diagnostics: %.1f%% of fits produced warnings; %d fit(s) failed.\n",
      pct_warned, as.integer(n_fail)
    ))
  }

  # ---- Power summary -------------------------------------------------------
  if (!is.null(x$summary) && nrow(x$summary) > 0L) {
    cat("\nPower summary")
    if (nrow(x$summary) > n_rows) {
      cat(" (first", n_rows, "of", nrow(x$summary), "rows)")
    }
    cat(":\n")
    print(head(as.data.frame(x$summary), n_rows), row.names = FALSE)
  }

  invisible(x)
}


#' Spot-check INLA posterior estimates against brms/Stan
#'
#' @description
#' Generates a small number of synthetic datasets using the same internal
#' data-generator as [brms_inla_power()], fits each dataset with both INLA
#' (as the package does) and `brms::brm()` (via Stan), and compares the
#' posterior mean, posterior SD, and 95 % credible-interval bounds for the
#' primary effect parameter.
#'
#' This is a **qualitative spot-check**, not a formal equivalence test.
#' Differences are expected due to the Laplace approximation used by INLA
#' versus full MCMC in Stan, and will vary across models and priors.
#'
#' @section Performance warning:
#' This function runs `n_check` full brms/Stan fits.  Even with
#' `brms_chains = 2` and `brms_iter = 1000`, this can take several minutes
#' for non-trivial models.  Use `n_check = 3` and `brms_chains = 2` for
#' quick sanity checks.
#'
#' @param formula Model formula (same as [brms_inla_power()]).
#' @param family GLM family object (default `gaussian()`).
#' @param priors A `brms::prior()` specification or `NULL` for defaults.
#'   Supported priors are translated to INLA controls where possible and audited
#'   in the returned settings.
#' @param data_generator Optional function `f(n, effect)` returning a
#'   data frame.  If `NULL`, the automatic generator is used.
#' @param effect_name Character vector of fixed-effect names (same as
#'   [brms_inla_power()]).
#' @param effect_value Numeric value (or named vector) of the true effect
#'   size used when generating data for the comparison.  Default `0.5`.
#' @param sample_size Integer sample size for each comparison dataset.
#'   Default `100`.
#' @param n_check Number of independent datasets to compare.  Default `5`.
#' @param brms_iter Total MCMC iterations per chain passed to
#'   `brms::brm()`.  Default `2000`.
#' @param brms_chains Number of MCMC chains.  Default `4`.
#' @param tolerance Numeric threshold: the flag is set to `TRUE` when the
#'   maximum absolute difference in posterior means exceeds `tolerance`.
#'   If `NULL` (default), `0.1 * mean(INLA posterior SD)` is used.
#' @param seed Integer random seed.  Default `42`.
#' @param error_sd,group_sd,obs_per_group,predictor_means,predictor_sds,family_args
#'   Passed to the automatic data generator when `data_generator = NULL`
#'   (see [brms_inla_power()] for details).
#' @param inla_num_threads INLA threading string (e.g. `"4:1"`).  Auto-
#'   detected when `NULL`.
#'
#' @return A list with components:
#'   \describe{
#'     \item{`comparisons`}{Data frame with one row per dataset and columns
#'       `sim`, `inla_ok`, `brms_ok`, `inla_mean`, `brms_mean`, `diff_mean`,
#'       `inla_sd`, `brms_sd`, `inla_ci_lower`, `inla_ci_upper`,
#'       `brms_ci_lower`, `brms_ci_upper`.}
#'     \item{`flag`}{`TRUE` if any `|diff_mean|` exceeds `tolerance`.}
#'     \item{`max_abs_diff`}{Maximum observed `|diff_mean|`.}
#'     \item{`tolerance`}{The tolerance value actually used.}
#'     \item{`settings`}{List of key settings for reproducibility.}
#'   }
#'
#' @export
validate_inla_vs_brms <- function(
    formula,
    family           = gaussian(),
    priors           = NULL,
    data_generator   = NULL,
    effect_name,
    effect_value     = 0.5,
    sample_size      = 100L,
    n_check          = 5L,
    brms_iter        = 2000L,
    brms_chains      = 4L,
    tolerance        = NULL,
    seed             = 42L,
    error_sd         = 1,
    group_sd         = 0.5,
    obs_per_group    = 10,
    predictor_means  = NULL,
    predictor_sds    = NULL,
    family_args      = list(),
    inla_num_threads = NULL
) {
  # ---- Dependency checks ---------------------------------------------------
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop(
      "Package 'brms' is required for validate_inla_vs_brms(). ",
      "Install it with: install.packages('brms')",
      call. = FALSE
    )
  }
  if (!requireNamespace("INLA", quietly = TRUE)) {
    stop(
      "Package 'INLA' is required for validate_inla_vs_brms(). ",
      "See https://www.r-inla.org for installation instructions.",
      call. = FALSE
    )
  }

  # ---- Input validation ----------------------------------------------------
  stopifnot(
    is.character(effect_name), length(effect_name) >= 1L,
    is.numeric(effect_value),  length(effect_value) >= 1L,
    is.numeric(sample_size),   sample_size >= 1L,
    is.numeric(n_check),       n_check >= 1L
  )

  set.seed(seed)

  # ---- Setup ---------------------------------------------------------------
  if (is.null(inla_num_threads)) {
    n_cores <- parallel::detectCores()
    if (!is.numeric(n_cores) || length(n_cores) != 1L || !is.finite(n_cores)) n_cores <- 1L
    inla_num_threads <- if (n_cores >= 4L) "4:1" else if (n_cores >= 2L) "2:1" else "1:1"
  }

  if (is.null(data_generator)) {
    data_generator <- .auto_data_generator(
      formula         = formula,
      effect_name     = effect_name,
      family          = family,
      family_args     = family_args,
      error_sd        = error_sd,
      group_sd        = group_sd,
      obs_per_group   = obs_per_group,
      predictor_means = predictor_means,
      predictor_sds   = predictor_sds
    )
  } else {
    stopifnot(is.function(data_generator))
  }

  fam_inla <- .to_inla_family(family)$inla
  needs_N  <- fam_inla %in% c("binomial", "betabinomial")
  needs_E  <- fam_inla %in% c("poisson")
  pmap     <- .map_brms_priors_to_inla(priors, inla_family = fam_inla)
  tf       <- .brms_to_inla_formula2(formula, hyper_by_re = pmap$hyper_by_re)
  inla_f   <- tf$inla_formula
  re_specs <- tf$re_specs
  pmap     <- .mark_unmatched_re_priors(pmap, tf$re_hyper_groups)
  pmap     <- .audit_re_correlation_terms(pmap, re_specs)

  # Build the named effect vector for the data generator
  eff_val <- if (length(effect_value) >= length(effect_name)) {
    setNames(as.numeric(effect_value[seq_along(effect_name)]), effect_name)
  } else {
    setNames(rep(as.numeric(effect_value[1L]), length(effect_name)), effect_name)
  }
  primary_eff <- effect_name[1L]
  brms_prior  <- if (!is.null(priors)) priors else NULL

  # ---- Run n_check comparisons ---------------------------------------------
  rows <- vector("list", n_check)

  for (s in seq_len(n_check)) {
    message(sprintf("validate_inla_vs_brms(): comparison %d / %d", s, n_check))

    dat <- data_generator(as.integer(sample_size), eff_val)

    # Add INLA-specific RE index columns (ignored by brms)
    if (length(re_specs) > 0L) {
      for (re in re_specs) {
        gid <- as.integer(as.factor(dat[[re$group]]))
        if (isTRUE(re$has_intercept) && is.null(dat[[re$id_intercept]]))
          dat[[re$id_intercept]] <- gid
        if (!is.null(re$slope) && is.null(dat[[re$id_slope]]))
          dat[[re$id_slope]] <- gid
      }
    }

    # ---- INLA fit ----------------------------------------------------------
    inla_args <- list(
      formula           = inla_f,
      data              = dat,
      family            = fam_inla,
      control.fixed     = pmap$control_fixed %||% list(),
      control.family    = pmap$control_family %||% list(),
      control.predictor = list(link = 1),
      num.threads       = inla_num_threads,
      verbose           = FALSE
    )
    if (needs_N && !is.null(dat$.Ntrials)) inla_args$Ntrials <- dat$.Ntrials
    if (needs_E && !is.null(dat$.E))       inla_args$E       <- dat$.E
    if (!is.null(dat$.scale))              inla_args$scale   <- dat$.scale

    fit_inla <- tryCatch(
      suppressWarnings(suppressMessages(do.call(INLA::inla, inla_args))),
      error = function(e) e
    )

    inla_ok <- !inherits(fit_inla, "error") && !is.null(fit_inla$summary.fixed)
    if (inla_ok) {
      fnms <- rownames(fit_inla$summary.fixed)
      nm   <- if (primary_eff %in% fnms) primary_eff else {
        cands <- grep(paste0("^", primary_eff), fnms, value = TRUE)
        if (length(cands)) cands[1L] else NA_character_
      }
      if (!is.na(nm)) {
        inla_mean  <- as.numeric(fit_inla$summary.fixed[nm, "mean"])
        inla_sd    <- as.numeric(fit_inla$summary.fixed[nm, "sd"])
        inla_ci_lo <- as.numeric(fit_inla$summary.fixed[nm, "0.025quant"])
        inla_ci_hi <- as.numeric(fit_inla$summary.fixed[nm, "0.975quant"])
      } else {
        inla_ok <- FALSE
        inla_mean <- inla_sd <- inla_ci_lo <- inla_ci_hi <- NA_real_
      }
    } else {
      inla_mean <- inla_sd <- inla_ci_lo <- inla_ci_hi <- NA_real_
    }

    # ---- brms / Stan fit ---------------------------------------------------
    fit_brms <- tryCatch(
      suppressWarnings(suppressMessages(
        brms::brm(
          formula = formula,
          data    = dat,
          family  = family,
          prior   = brms_prior,
          iter    = as.integer(brms_iter),
          chains  = as.integer(brms_chains),
          refresh = 0L,
          seed    = seed + s
        )
      )),
      error = function(e) e
    )

    brms_ok <- !inherits(fit_brms, "error")
    if (brms_ok) {
      fe <- tryCatch(brms::fixef(fit_brms), error = function(e) NULL)
      if (!is.null(fe) && is.matrix(fe)) {
        bnms <- rownames(fe)
        bnm  <- if (primary_eff %in% bnms) primary_eff else {
          cands <- grep(paste0("^", primary_eff), bnms, value = TRUE)
          if (length(cands)) cands[1L] else NA_character_
        }
        if (!is.na(bnm)) {
          brms_mean  <- as.numeric(fe[bnm, "Estimate"])
          brms_sd    <- as.numeric(fe[bnm, "Est.Error"])
          brms_ci_lo <- as.numeric(fe[bnm, "Q2.5"])
          brms_ci_hi <- as.numeric(fe[bnm, "Q97.5"])
        } else {
          brms_ok <- FALSE
          brms_mean <- brms_sd <- brms_ci_lo <- brms_ci_hi <- NA_real_
        }
      } else {
        brms_ok <- FALSE
        brms_mean <- brms_sd <- brms_ci_lo <- brms_ci_hi <- NA_real_
      }
    } else {
      brms_mean <- brms_sd <- brms_ci_lo <- brms_ci_hi <- NA_real_
    }

    rows[[s]] <- data.frame(
      sim           = s,
      inla_ok       = inla_ok,
      brms_ok       = brms_ok,
      inla_mean     = inla_mean,
      brms_mean     = brms_mean,
      diff_mean     = inla_mean - brms_mean,
      inla_sd       = inla_sd,
      brms_sd       = brms_sd,
      inla_ci_lower = inla_ci_lo,
      inla_ci_upper = inla_ci_hi,
      brms_ci_lower = brms_ci_lo,
      brms_ci_upper = brms_ci_hi,
      stringsAsFactors = FALSE
    )
  }

  comparisons <- do.call(rbind, rows)
  rownames(comparisons) <- NULL

  # ---- Compute flag --------------------------------------------------------
  if (is.null(tolerance)) {
    avg_sd <- mean(comparisons$inla_sd, na.rm = TRUE)
    tol    <- if (is.finite(avg_sd) && avg_sd > 0) 0.1 * avg_sd else 0.1
  } else {
    stopifnot(is.numeric(tolerance), length(tolerance) == 1L, tolerance > 0)
    tol <- as.numeric(tolerance)
  }

  max_diff <- if (any(is.finite(comparisons$diff_mean)))
    max(abs(comparisons$diff_mean), na.rm = TRUE) else NA_real_
  flag <- is.finite(max_diff) && max_diff > tol

  list(
    comparisons  = comparisons,
    flag         = flag,
    max_abs_diff = max_diff,
    tolerance    = tol,
    settings     = list(
      formula      = formula,
      effect_name  = effect_name,
      effect_value = effect_value,
      sample_size  = sample_size,
      n_check      = n_check,
      brms_iter    = brms_iter,
      brms_chains  = brms_chains,
      seed         = seed,
      prior_translation = pmap$prior_audit
    )
  )
}
