#' Compute weights from a parametric prior distribution over an effect grid
#'
#' Evaluates a named distribution at each value in `effect_values`, normalises
#' the result to sum to 1, and returns a named numeric vector compatible with
#' [compute_assurance()] and [beta_weights_on_grid()].
#'
#' @param dist_spec A list with at minimum a `$dist` element naming the
#'   distribution.  Remaining elements are distribution parameters:
#'   * `"normal"` — `mean`, `sd`
#'   * `"uniform"` — `min`, `max` (defaults: range of `effect_values`)
#'   * `"beta"` — either `shape1`, `shape2` *or* `mode`, `n` (concentration;
#'     must be > 2).  Effects are rescaled to \eqn{[0, 1]} automatically when
#'     outside that interval.
#' @param effect_values Numeric vector of effect-grid values at which to
#'   evaluate the density.
#' @return Named numeric vector of normalised weights (sums to 1).
#' @keywords internal
.compute_weights_from_dist <- function(dist_spec, effect_values) {
  dist <- dist_spec$dist
  if (!is.character(dist) || !nzchar(dist)) {
    stop("dist_spec$dist must be a non-empty character string.", call. = FALSE)
  }

  if (identical(dist, "normal")) {
    mu    <- dist_spec$mean %||% 0
    sigma <- dist_spec$sd   %||% 1
    if (!is.numeric(sigma) || length(sigma) != 1L || sigma <= 0) {
      stop("Normal prior requires sd > 0.", call. = FALSE)
    }
    w <- stats::dnorm(effect_values, mean = mu, sd = sigma)

  } else if (identical(dist, "uniform")) {
    rng   <- range(effect_values)
    lo    <- dist_spec$min %||% rng[1L]
    hi    <- dist_spec$max %||% rng[2L]
    if (lo >= hi) stop("Uniform prior requires min < max.", call. = FALSE)
    w <- stats::dunif(effect_values, min = lo, max = hi)

  } else if (identical(dist, "beta")) {
    # Support two parameterisations: shape1/shape2 or mode/n (concentration)
    if (!is.null(dist_spec$shape1) && !is.null(dist_spec$shape2)) {
      a <- dist_spec$shape1
      b <- dist_spec$shape2
      if (a <= 0 || b <= 0) stop("Beta prior requires shape1 > 0 and shape2 > 0.", call. = FALSE)
      x <- effect_values
    } else if (!is.null(dist_spec$mode) && !is.null(dist_spec$n)) {
      mode_val <- dist_spec$mode
      n_val    <- dist_spec$n
      if (mode_val <= 0 || mode_val >= 1) {
        stop("Beta prior (mode/n parameterisation) requires 0 < mode < 1.", call. = FALSE)
      }
      if (n_val <= 2) {
        stop("Beta prior (mode/n parameterisation) requires n > 2.", call. = FALSE)
      }
      a <- mode_val * (n_val - 2) + 1
      b <- (1 - mode_val) * (n_val - 2) + 1
      # Rescale effects to [0,1] when outside that range
      rng <- range(effect_values)
      if (rng[1L] < 0 || rng[2L] > 1) {
        if (diff(rng) < .Machine$double.eps * 100) {
          stop("All effect values are identical; cannot rescale for beta prior.", call. = FALSE)
        }
        message("compute_assurance / .compute_weights_from_dist: ",
                "effects not in [0,1]; rescaling linearly.")
        x <- (effect_values - rng[1L]) / (rng[2L] - rng[1L])
      } else {
        x <- effect_values
      }
    } else {
      stop("Beta prior requires either {shape1, shape2} or {mode, n} parameters.",
           call. = FALSE)
    }
    w <- stats::dbeta(x, shape1 = a, shape2 = b)

  } else {
    stop("Unsupported distribution '", dist,
         "'. Supported: 'normal', 'uniform', 'beta'.",
         call. = FALSE)
  }

  if (!any(is.finite(w) & w > 0)) {
    stop("All computed weights are zero or non-finite. ",
         "Check your prior specification and effect grid.", call. = FALSE)
  }
  w / sum(w, na.rm = TRUE)
}


#' Create prior weights over an effect grid for use with compute_assurance()
#'
#' A convenience wrapper that evaluates a parametric prior distribution at
#' each value in `effects`, normalises the densities to sum to 1, and returns
#' a named numeric vector that can be passed directly to [compute_assurance()]
#' as `prior_weights`.
#'
#' The output format is identical to that of [beta_weights_on_grid()], which
#' can also be used directly when a beta prior is appropriate.
#'
#' @param effects Numeric vector of effect-grid values (the same vector
#'   passed to [brms_inla_power()] as `effect_grid`).
#' @param dist Character string naming the distribution.  One of
#'   `"normal"`, `"uniform"`, or `"beta"`.
#' @param ... Named distribution parameters forwarded to the internal
#'   density calculator.  See [compute_assurance()] `@details` for the
#'   full list per distribution.
#'
#' @return Named numeric vector of normalised weights (sums to 1), with
#'   names equal to `as.character(effects)`.
#' @export
#'
#' @examples
#' effects <- c(0.1, 0.3, 0.5, 0.7, 0.9)
#' # Normal prior centred on 0.5
#' assurance_prior_weights(effects, dist = "normal", mean = 0.5, sd = 0.2)
#' # Uniform prior (equal weight)
#' assurance_prior_weights(effects, dist = "uniform")
assurance_prior_weights <- function(effects, dist = c("normal", "uniform", "beta"), ...) {
  stopifnot(is.numeric(effects), length(effects) >= 1L)
  dist <- match.arg(dist)
  spec <- c(list(dist = dist), list(...))
  w    <- .compute_weights_from_dist(spec, effects)
  names(w) <- as.character(effects)
  w
}


#' Compute unconditional Bayesian assurance from simulation results
#'
#' @description
#' Computes **unconditional Bayesian assurance** — the probability of a
#' successful trial outcome averaged over prior uncertainty about the true
#' effect size — from the output of [brms_inla_power()] or related engines.
#'
#' @details
#' ## Assurance vs. conditional Bayesian power
#'
#' The simulations run by [brms_inla_power()] are *conditional*: for each
#' point on the effect grid the engine estimates the probability that the
#' chosen decision rule is satisfied.  This is Bayesian design power
#' (a function of the unknown true effect).
#'
#' **Assurance**, in the sense of O'Hagan & Stevens (2001) and
#' O'Hagan, Stevens & Campbell (2005), is the *unconditional* version:
#'
#' \deqn{
#'   A(n) = \int \mathrm{Power}(n, \delta)\, \pi(\delta)\, d\delta
#'         \approx \sum_j w_j \cdot \mathrm{Power}(n, \delta_j)
#' }
#'
#' where \eqn{\pi(\delta)} is a **design prior** on the effect size and
#' \eqn{w_j} are the normalised prior weights over the discrete effect grid
#' \eqn{\{\delta_j\}}.  Assurance therefore accounts for the investigator's
#' genuine uncertainty about the effect, not just a single "assumed" value
#' (Ristl et al., 2019; Kunzmann et al., 2021).
#'
#' ## Variance uncertainty
#'
#' If the simulation was run with multiple sampled variance parameters (stored
#' in columns such as `sampled_error_sd` or `sampled_group_sd` in the results),
#' the averaging over those values is already implicit in each per-cell power
#' estimate, so no additional action is required here.
#'
#' ## Supplying prior weights
#'
#' Two forms are accepted for `prior_weights`:
#'
#' * **Named numeric vector** — names must match the effect-grid values used
#'   in `power_result` (as produced by `as.character()`, which is the format
#'   used by [beta_weights_on_grid()] and [assurance_prior_weights()]).
#'   Unnamed vectors are accepted only when their length equals the number of
#'   unique effect values, in which case they are applied in ascending order.
#' * **Distribution list** — a list with at minimum `$dist` naming one of
#'   `"normal"` (`mean`, `sd`), `"uniform"` (`min`, `max`), or `"beta"`
#'   (`shape1`/`shape2` **or** `mode`/`n`).  Weights are computed by
#'   evaluating the density at each grid point and normalising.
#'   *Supported for single-effect results only.*
#'
#' For **multi-effect** grids the `prior_weights` argument must be a numeric
#' vector of length equal to the number of unique effect combinations in the
#' summary (sorted lexicographically by effect columns).  Use
#' [assurance_prior_weights()] or [beta_weights_on_grid()] to construct
#' compatible weights for single-effect cases.
#'
#' @param power_result A list returned by [brms_inla_power()],
#'   [brms_inla_power_parallel()], [brms_inla_power_two_stage()], or
#'   [brms_inla_power_sequential()]; or a plain data frame with at least
#'   columns `n` and the relevant `power_*` metric column.
#' @param prior_weights Either (a) a named numeric vector of weights over
#'   effect-size values (must sum to 1 within tolerance `0.01`), or (b) a
#'   list specifying a distribution — see **Details**.  The output of
#'   [beta_weights_on_grid()] and [assurance_prior_weights()] is directly
#'   compatible.
#' @param metric Character string selecting the decision metric.  Must match
#'   a column present in the summary.  One of:
#'   * `"direction"` → `power_direction`
#'   * `"threshold"` → `power_threshold`
#'   * `"rope"`      → `power_rope`
#'   * `"bf"`        → `bf_hit_10`
#' @param weight_tol Numeric tolerance for the weights-sum-to-1 check
#'   (default `0.01`).
#'
#' @return A list of class `"powerbrmsINLA_assurance"` containing:
#'   \describe{
#'     \item{`assurance`}{Data frame with columns `sample_size` and
#'       `assurance`.}
#'     \item{`metric`}{The decision metric used.}
#'     \item{`power_col`}{Name of the summary column used for power.}
#'     \item{`prior_spec`}{The `prior_weights` argument as supplied (useful
#'       for reproducibility).}
#'     \item{`weights`}{Named numeric vector of the normalised weights
#'       actually applied.}
#'     \item{`eff_cols`}{Character vector of effect-grid column names
#'       identified in the summary.}
#'   }
#'
#' @references
#' O'Hagan, A., & Stevens, J. W. (2001). Bayesian assessment of sample size
#' for clinical trials of cost-effectiveness. *Medical Decision Making*,
#' 21(3), 219–230. \doi{10.1177/0272989X0102100307}
#'
#' O'Hagan, A., Stevens, J. W., & Campbell, M. J. (2005). Assurance in
#' clinical trial design. *Pharmaceutical Statistics*, 4(3), 187–201.
#' \doi{10.1002/pst.175}
#'
#' Ristl, R., Glimm, E., Stallard, N., & Posch, M. (2019). Optimal design
#' and analysis of two-stage adaptive enrichment trials. *Biometrical Journal*,
#' 61(6), 1461–1481.
#'
#' Kunzmann, K., Grayling, M. J., Lee, K. M., Robertson, D. S., Rufibach,
#' K., & Wason, J. M. S. (2021). A review of Bayesian perspectives on sample
#' size derivation for confirmatory trials. *The American Statistician*, 75(4),
#' 424–432. \doi{10.1080/00031305.2021.1901782}
#'
#' @export
#'
#' @examples
#' # Build a small synthetic power_result without running INLA
#' syn_summary <- data.frame(
#'   n               = rep(c(50, 100, 200), each = 3),
#'   treatment       = rep(c(0.2, 0.5, 0.8), 3),
#'   power_direction = c(0.40, 0.65, 0.85,
#'                       0.60, 0.82, 0.95,
#'                       0.72, 0.90, 0.98),
#'   stringsAsFactors = FALSE
#' )
#' syn_result <- list(
#'   summary  = syn_summary,
#'   settings = list(effect_name = "treatment")
#' )
#'
#' # (a) Uniform weights — assurance is the simple mean of per-cell powers
#' w_uniform <- c("0.2" = 1/3, "0.5" = 1/3, "0.8" = 1/3)
#' out <- compute_assurance(syn_result, prior_weights = w_uniform)
#' print(out)
#'
#' # (b) Normal design prior centred on a medium-sized effect
#' out2 <- compute_assurance(
#'   syn_result,
#'   prior_weights = list(dist = "normal", mean = 0.5, sd = 0.2)
#' )
#' print(out2)
#'
#' # (c) Using assurance_prior_weights() to build the weight vector explicitly
#' w_norm <- assurance_prior_weights(c(0.2, 0.5, 0.8), dist = "normal",
#'                                   mean = 0.5, sd = 0.2)
#' out3 <- compute_assurance(syn_result, prior_weights = w_norm)
compute_assurance <- function(
    power_result,
    prior_weights,
    metric     = c("direction", "threshold", "rope", "bf"),
    weight_tol = 0.01
) {
  metric <- match.arg(metric)

  # ---- 1. Extract summary --------------------------------------------------
  s <- if (is.list(power_result) && !is.null(power_result$summary)) {
    as.data.frame(power_result$summary)
  } else if (is.data.frame(power_result)) {
    power_result
  } else {
    stop("'power_result' must be a list with a $summary element, or a data.frame.",
         call. = FALSE)
  }

  if (!"n" %in% names(s)) {
    stop("The summary must contain a column 'n'.", call. = FALSE)
  }
  if (nrow(s) == 0L) {
    stop("The summary has 0 rows; nothing to compute.", call. = FALSE)
  }

  # ---- 2. Identify the power column ----------------------------------------
  power_col <- switch(
    metric,
    direction = "power_direction",
    threshold = "power_threshold",
    rope      = "power_rope",
    bf        = "bf_hit_10"
  )

  if (!power_col %in% names(s)) {
    stop("Column '", power_col, "' not found in the summary. ",
         "Ensure the corresponding metric was computed in brms_inla_power().",
         call. = FALSE)
  }

  # ---- 3. Identify effect-grid columns -------------------------------------
  non_eff_cols <- c(
    "n",
    grep("^power_",     names(s), value = TRUE),
    grep("^bf_hit_",    names(s), value = TRUE),
    grep("^mean_log10", names(s), value = TRUE),
    "avg_ci_width", "ci_coverage", "ci_width",
    "ciw_q05", "ciw_q25", "ciw_q50", "ciw_q75", "ciw_q95",
    "bf_hit_3", "bf_hit_10", "mean_log10_bf", "nsims_ok",
    "sims_used", "assurance", "conditional_power", "effect_val"
  )
  eff_cols <- setdiff(names(s), non_eff_cols)

  # Prefer the authoritative effect_name from settings when available
  settings   <- if (is.list(power_result) && !is.null(power_result$settings))
    power_result$settings else NULL
  if (!is.null(settings$effect_name)) {
    known_effs <- intersect(settings$effect_name, names(s))
    if (length(known_effs) > 0L) eff_cols <- known_effs
  }

  if (length(eff_cols) == 0L) {
    stop("No effect-grid columns found in the summary. ",
         "Check that brms_inla_power() returned the expected structure.",
         call. = FALSE)
  }

  is_multi <- length(eff_cols) > 1L

  # ---- 4. Process prior_weights --------------------------------------------
  if (is.list(prior_weights) && !is.null(prior_weights$dist)) {
    # Distribution specification
    if (is_multi) {
      stop("Distribution-based prior_weights are only supported for single-effect results. ",
           "For multi-effect grids, supply a named numeric vector of weights.",
           call. = FALSE)
    }
    eff_col        <- eff_cols[1L]
    unique_effects <- sort(unique(s[[eff_col]]))
    w              <- .compute_weights_from_dist(prior_weights, unique_effects)
    names(w)       <- as.character(unique_effects)
    prior_spec     <- prior_weights

  } else if (is.numeric(prior_weights)) {
    w          <- prior_weights
    prior_spec <- list(type = "user-supplied weights", weights = w)

  } else {
    stop("'prior_weights' must be a named numeric vector or a list with a $dist element.",
         call. = FALSE)
  }

  # ---- 5. Validate weights sum to 1 ----------------------------------------
  w_sum <- sum(w, na.rm = TRUE)
  if (!is.finite(w_sum) || abs(w_sum - 1) > weight_tol) {
    stop("'prior_weights' must sum to 1 (tolerance ", weight_tol,
         "). Current sum: ", round(w_sum, 6L), call. = FALSE)
  }
  w <- w / sum(w)  # normalise for floating-point safety

  # ---- 6. Map weights to effect combinations in the summary ----------------
  sample_sizes <- sort(unique(s$n))

  if (is_multi) {
    # Multi-effect: weights are positional over unique combinations sorted by
    # effect columns (ascending)
    unique_combos <- unique(s[eff_cols])
    # Sort lexicographically
    ord <- do.call(order, lapply(eff_cols, function(ec) unique_combos[[ec]]))
    unique_combos <- unique_combos[ord, , drop = FALSE]
    n_combos <- nrow(unique_combos)

    if (length(w) != n_combos) {
      stop("For multi-effect results, 'prior_weights' must have length equal to the ",
           "number of unique effect combinations (", n_combos, "). ",
           "Got ", length(w), ".", call. = FALSE)
    }
    names(w) <- seq_len(n_combos)

    # Build a key for each row in summary and in unique_combos
    row_key   <- do.call(paste, c(as.list(s[eff_cols]),          list(sep = "\x1f")))
    combo_key <- do.call(paste, c(as.list(unique_combos[eff_cols]), list(sep = "\x1f")))
    combo_idx <- match(row_key, combo_key)
    if (any(is.na(combo_idx))) {
      stop("Some rows in the summary could not be matched to a unique effect combination. ",
           "This is an internal error; please report.", call. = FALSE)
    }
    s$.w <- w[combo_idx]

  } else {
    # Single effect
    eff_col        <- eff_cols[1L]
    unique_effects <- sort(unique(s[[eff_col]]))

    if (is.null(names(w))) {
      # Unnamed: apply positionally over sorted unique effects
      if (length(w) != length(unique_effects)) {
        stop("Unnamed 'prior_weights' must have length equal to the number of unique ",
             "effect values (", length(unique_effects), "). Got ", length(w), ".",
             call. = FALSE)
      }
      names(w) <- as.character(unique_effects)
    }

    # Match by character representation (same convention as beta_weights_on_grid)
    eff_char <- as.character(s[[eff_col]])
    w_matched <- w[eff_char]

    # Fallback: numeric tolerance matching for floating-point mismatches
    if (any(is.na(w_matched))) {
      for (i in which(is.na(w_matched))) {
        ev    <- s[[eff_col]][i]
        diffs <- abs(as.numeric(names(w)) - ev)
        best  <- which.min(diffs)
        if (is.finite(diffs[best]) && diffs[best] < 1e-9) {
          w_matched[i] <- w[best]
        }
      }
    }

    still_na <- which(is.na(w_matched))
    if (length(still_na) > 0L) {
      missing_vals <- unique(s[[eff_col]][still_na])
      stop("'prior_weights' names do not cover effect values: ",
           paste(missing_vals, collapse = ", "),
           ". Check that weight names match the effect grid.",
           call. = FALSE)
    }
    s$.w <- w_matched
  }

  # ---- 7. Compute assurance per sample size --------------------------------
  assurance_vals <- vapply(sample_sizes, function(ni) {
    s_n <- s[s$n == ni, , drop = FALSE]
    pwr <- s_n[[power_col]]
    wts <- s_n$.w
    # Only include finite power/weight pairs
    ok  <- is.finite(pwr) & is.finite(wts)
    if (!any(ok)) return(NA_real_)
    wts_ok  <- wts[ok]
    wts_ok  <- wts_ok / sum(wts_ok)  # renormalise within n if any rows dropped
    sum(pwr[ok] * wts_ok)
  }, numeric(1L))

  # ---- 8. Build and return output ------------------------------------------
  result <- list(
    assurance  = data.frame(sample_size = sample_sizes, assurance = assurance_vals,
                             stringsAsFactors = FALSE),
    metric     = metric,
    power_col  = power_col,
    prior_spec = prior_spec,
    weights    = w,
    eff_cols   = eff_cols
  )
  class(result) <- "powerbrmsINLA_assurance"
  result
}


#' Print method for powerbrmsINLA_assurance objects
#'
#' Displays the unconditional Bayesian assurance by sample size together with
#' a plain-language description of the prior used.
#'
#' @param x An object of class `"powerbrmsINLA_assurance"` returned by
#'   [compute_assurance()].
#' @param digits Number of significant digits for assurance values (default 4).
#' @param ... Unused; present for S3 compatibility.
#'
#' @return `x`, invisibly.
#' @export
print.powerbrmsINLA_assurance <- function(x, digits = 4L, ...) {
  cat("Bayesian assurance (unconditional probability)\n")
  cat("==============================================\n")
  cat("Decision metric :", x$metric, "(column:", x$power_col, ")\n")
  cat("Effect column(s):", paste(x$eff_cols, collapse = ", "), "\n")

  # Describe the prior
  ps <- x$prior_spec
  if (!is.null(ps$dist)) {
    param_str <- paste(
      vapply(setdiff(names(ps), "dist"), function(nm)
        paste0(nm, " = ", ps[[nm]]), character(1L)),
      collapse = ", "
    )
    cat("Design prior    :", ps$dist, "(", param_str, ")\n")
  } else {
    cat("Design prior    : user-supplied weights over",
        length(x$weights), "effect value(s)\n")
    if (!is.null(names(x$weights)) && length(x$weights) <= 10L) {
      w_str <- paste(
        paste0(names(x$weights), " = ", round(x$weights, 4L)),
        collapse = ", "
      )
      cat("  Weights       :", w_str, "\n")
    }
  }

  cat("\nAssurance by sample size:\n")
  df          <- x$assurance
  df$assurance <- round(df$assurance, digits)
  print(df, row.names = FALSE)
  invisible(x)
}
