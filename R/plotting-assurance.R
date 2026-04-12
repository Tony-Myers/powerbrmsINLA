#' Plot Assurance Curve(s) vs Sample Size
#'
#' Plots unconditional Bayesian assurance against sample size for one or more
#' design priors, with an optional horizontal reference line at a target
#' assurance level.  Each prior produces one line, making it easy to compare
#' how the choice of design prior changes the required sample size.
#'
#' @param assurance_results A single `powerbrmsINLA_assurance` object (from
#'   [compute_assurance()]) **or** a list of such objects.
#' @param labels Character vector naming each design prior, e.g.
#'   `c("Sceptical", "Primary", "Optimistic")`.  Defaults to
#'   `"Design prior"` for a single object or `"Prior 1"`, `"Prior 2"`, … for
#'   multiple objects.
#' @param target Numeric in \eqn{(0, 1)}: horizontal reference line added at
#'   this assurance level (default 0.80).  Set to `NULL` to suppress the line.
#' @param title,subtitle Optional plot title and subtitle strings.
#' @return A ggplot object.
#'
#' @examples
#' syn_summary <- data.frame(
#'   n               = rep(c(50, 100, 200), each = 3),
#'   treatment       = rep(c(0.2, 0.5, 0.8), 3),
#'   power_direction = c(0.40, 0.65, 0.85, 0.60, 0.82, 0.95, 0.72, 0.90, 0.98)
#' )
#' pr <- list(summary = syn_summary, settings = list(effect_name = "treatment"))
#' a1 <- compute_assurance(pr, list(dist = "normal", mean = 0.5, sd = 0.15))
#' a2 <- compute_assurance(pr, list(dist = "normal", mean = 0.5, sd = 0.30))
#' plot_assurance_curve(list(a1, a2), labels = c("Informative", "Diffuse"))
#'
#' @export
plot_assurance_curve <- function(assurance_results,
                                  labels   = NULL,
                                  target   = 0.80,
                                  title    = NULL,
                                  subtitle = NULL) {
  # Accept a single object or a list
  if (inherits(assurance_results, "powerbrmsINLA_assurance")) {
    assurance_results <- list(assurance_results)
  }

  n_priors <- length(assurance_results)
  if (n_priors == 0L) stop("'assurance_results' must contain at least one object.", call. = FALSE)

  for (i in seq_len(n_priors)) {
    if (!inherits(assurance_results[[i]], "powerbrmsINLA_assurance")) {
      stop("Each element of 'assurance_results' must be a 'powerbrmsINLA_assurance' object.",
           call. = FALSE)
    }
  }

  if (is.null(labels)) {
    labels <- if (n_priors == 1L) "Design prior" else paste0("Prior ", seq_len(n_priors))
  }
  if (length(labels) != n_priors) {
    stop("'labels' must have the same length as 'assurance_results' (", n_priors, ").",
         call. = FALSE)
  }

  # Build combined data frame
  plot_dat <- do.call(rbind, lapply(seq_len(n_priors), function(i) {
    df        <- assurance_results[[i]]$assurance
    df$prior  <- labels[i]
    df
  }))
  plot_dat$prior <- factor(plot_dat$prior, levels = labels)

  p <- ggplot2::ggplot(
    plot_dat,
    ggplot2::aes(x = sample_size, y = assurance, colour = prior, group = prior)
  ) +
    .geom_line_lw(width = 1) +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_colour_viridis_d(name = "Design prior") +
    ggplot2::scale_y_continuous(
      limits = c(0, 1),
      labels = scales::percent_format(accuracy = 1)
    ) +
    ggplot2::labs(
      x        = "Sample size",
      y        = "Assurance",
      title    = title    %||% "Assurance by sample size",
      subtitle = subtitle
    ) +
    ggplot2::theme_minimal()

  if (!is.null(target)) {
    p <- p + ggplot2::geom_hline(
      yintercept = target,
      linetype   = "dashed",
      colour     = "grey40",
      alpha      = 0.8
    )
  }

  p
}


#' Plot Conditional Power Curves with Assurance Overlay
#'
#' Shows the conditional power curve for each effect-size value in the design
#' grid (thin lines), overlaid with the unconditional Bayesian assurance curve
#' (thick line).  This communicates that assurance is the prior-weighted
#' average of the conditional power curves.
#'
#' @param power_result Output from [brms_inla_power()] (or compatible wrappers)
#'   containing a `$summary` element.
#' @param assurance_result A `powerbrmsINLA_assurance` object from
#'   [compute_assurance()].
#' @param metric Decision metric to plot: `"direction"` (default),
#'   `"threshold"`, `"rope"`, or `"bf"`.
#' @param title,subtitle Optional plot title and subtitle strings.
#' @return A ggplot object.
#'
#' @examples
#' syn_summary <- data.frame(
#'   n               = rep(c(50, 100, 200), each = 3),
#'   treatment       = rep(c(0.2, 0.5, 0.8), 3),
#'   power_direction = c(0.40, 0.65, 0.85, 0.60, 0.82, 0.95, 0.72, 0.90, 0.98)
#' )
#' pr <- list(summary = syn_summary, settings = list(effect_name = "treatment"))
#' a  <- compute_assurance(pr, list(dist = "normal", mean = 0.5, sd = 0.15))
#' plot_power_assurance_overlay(pr, a)
#'
#' @export
plot_power_assurance_overlay <- function(power_result,
                                          assurance_result,
                                          metric   = c("direction", "threshold", "rope", "bf"),
                                          title    = NULL,
                                          subtitle = NULL) {
  stopifnot("summary" %in% names(power_result))
  if (!inherits(assurance_result, "powerbrmsINLA_assurance")) {
    stop("'assurance_result' must be a 'powerbrmsINLA_assurance' object.", call. = FALSE)
  }

  metric    <- match.arg(metric)
  power_col <- switch(
    metric,
    direction = "power_direction",
    threshold = "power_threshold",
    rope      = "power_rope",
    bf        = "bf_hit_10"
  )

  dat      <- as.data.frame(power_result$summary)
  eff_cols <- assurance_result$eff_cols

  if (!power_col %in% names(dat)) {
    stop("Column '", power_col, "' not found in the power result summary. ",
         "Ensure this metric was computed in brms_inla_power().", call. = FALSE)
  }

  # Build a human-readable label per effect combination
  if (length(eff_cols) == 1L) {
    dat$eff_label <- as.character(dat[[eff_cols[1L]]])
    legend_name   <- eff_cols[1L]
  } else {
    dat$eff_label <- do.call(
      paste,
      c(lapply(eff_cols, function(ec) paste0(ec, "=", dat[[ec]])), list(sep = ", "))
    )
    legend_name <- "Effect"
  }

  assurance_dat <- assurance_result$assurance

  p <- ggplot2::ggplot() +
    # Thin lines: conditional power per effect value
    .geom_line_lw(
      data    = dat,
      mapping = ggplot2::aes(
        x      = n,
        y      = !!rlang::sym(power_col),
        colour = eff_label,
        group  = eff_label
      ),
      width = 0.7,
      alpha = 0.7
    ) +
    ggplot2::geom_point(
      data    = dat,
      mapping = ggplot2::aes(
        x      = n,
        y      = !!rlang::sym(power_col),
        colour = eff_label
      ),
      size  = 1.5,
      alpha = 0.7
    ) +
    ggplot2::scale_colour_viridis_d(name = legend_name) +
    # Thick line: assurance (weighted average of the conditional curves)
    .geom_line_lw(
      data         = assurance_dat,
      mapping      = ggplot2::aes(x = sample_size, y = assurance),
      colour       = "#C0392B",
      width        = 1.8,
      inherit.aes  = FALSE
    ) +
    ggplot2::geom_point(
      data        = assurance_dat,
      mapping     = ggplot2::aes(x = sample_size, y = assurance),
      colour      = "#C0392B",
      size        = 3,
      inherit.aes = FALSE
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1),
      labels = scales::percent_format(accuracy = 1)
    ) +
    ggplot2::labs(
      x        = "Sample size",
      y        = "Power / Assurance",
      title    = title    %||% paste0("Conditional power & assurance (", metric, ")"),
      subtitle = subtitle %||%
        "Thin lines: conditional power per effect size; thick red line: assurance (weighted average)"
    ) +
    ggplot2::theme_minimal()

  p
}


#' Plot Design Prior Density over the Effect Grid
#'
#' Visualises one or more design prior distributions over the effect-size range
#' used in a power analysis.  Vertical dashed lines mark the grid points at
#' which power was simulated, helping assess prior-grid alignment.
#'
#' @param prior_spec A single prior specification list (e.g.
#'   `list(dist = "normal", mean = 0.5, sd = 0.15)`) **or** a list of such
#'   lists.  Supported distributions: `"normal"` (`mean`, `sd`),
#'   `"uniform"` (`min`, `max`), `"beta"` (`shape1`, `shape2` or `mode`, `n`).
#' @param effect_grid Numeric vector of effect-grid values passed to
#'   [brms_inla_power()].  Vertical dashed lines are drawn at each value.
#' @param labels Character vector naming each prior.  Defaults to
#'   `"Design prior"` for a single spec or `"Prior 1"`, `"Prior 2"`, … for
#'   multiple.
#' @param n_points Number of x points at which to evaluate the density curve
#'   (default 300).
#' @param title,subtitle Optional plot title and subtitle strings.
#' @return A ggplot object.
#'
#' @examples
#' plot_design_prior(
#'   prior_spec  = list(dist = "normal", mean = 0.5, sd = 0.15),
#'   effect_grid = c(0.2, 0.5, 0.8)
#' )
#'
#' # Multiple priors
#' plot_design_prior(
#'   prior_spec  = list(
#'     list(dist = "normal", mean = 0.5, sd = 0.10),
#'     list(dist = "normal", mean = 0.5, sd = 0.30)
#'   ),
#'   effect_grid = c(0.2, 0.5, 0.8),
#'   labels      = c("Informative", "Diffuse")
#' )
#'
#' @export
plot_design_prior <- function(prior_spec,
                               effect_grid,
                               labels   = NULL,
                               n_points = 300L,
                               title    = NULL,
                               subtitle = NULL) {
  if (!is.numeric(effect_grid) || length(effect_grid) == 0L) {
    stop("'effect_grid' must be a non-empty numeric vector.", call. = FALSE)
  }

  # Accept a single spec or a list of specs
  if (is.list(prior_spec) && !is.null(prior_spec$dist)) {
    prior_spec <- list(prior_spec)
  }
  n_priors <- length(prior_spec)
  if (n_priors == 0L) stop("'prior_spec' must contain at least one specification.", call. = FALSE)

  if (is.null(labels)) {
    labels <- if (n_priors == 1L) "Design prior" else paste0("Prior ", seq_len(n_priors))
  }
  if (length(labels) != n_priors) {
    stop("'labels' must have the same length as 'prior_spec' (", n_priors, ").", call. = FALSE)
  }

  # Build x range: effect_grid range padded by 20 % (min 0.2 each side)
  rng <- range(effect_grid)
  pad <- max(0.2 * diff(rng), 0.2)
  x_seq <- seq(rng[1L] - pad, rng[2L] + pad, length.out = n_points)

  # Compute density for each prior specification
  dens_dat <- do.call(rbind, lapply(seq_len(n_priors), function(i) {
    y         <- .eval_prior_density(prior_spec[[i]], x_seq, effect_grid)
    data.frame(x = x_seq, density = y, prior = labels[i],
               stringsAsFactors = FALSE)
  }))
  dens_dat$prior <- factor(dens_dat$prior, levels = labels)

  p <- ggplot2::ggplot(dens_dat, ggplot2::aes(x = x, y = density, colour = prior)) +
    .geom_line_lw(width = 1) +
    ggplot2::geom_vline(
      xintercept = effect_grid,
      linetype   = "dashed",
      colour     = "grey55",
      alpha      = 0.75
    ) +
    ggplot2::scale_colour_viridis_d(name = "Design prior") +
    ggplot2::labs(
      x        = "Effect size",
      y        = "Density",
      title    = title    %||% "Design prior density",
      subtitle = subtitle %||% "Dashed lines: effect-grid values used in simulation"
    ) +
    ggplot2::theme_minimal()

  p
}


# ---- Internal helpers --------------------------------------------------------

#' Evaluate a prior density over a continuous x grid
#'
#' @param spec  List with `$dist` and distribution parameters.
#' @param x     Numeric vector of evaluation points.
#' @param effect_grid Numeric effect-grid values (used for beta rescaling range).
#' @return Numeric vector of density values (same length as `x`).
#' @keywords internal
.eval_prior_density <- function(spec, x, effect_grid) {
  dist <- spec$dist
  if (!is.character(dist) || !nzchar(dist)) {
    stop("Prior spec must contain a non-empty 'dist' element.", call. = FALSE)
  }

  if (identical(dist, "normal")) {
    mu    <- spec$mean %||% 0
    sigma <- spec$sd   %||% 1
    if (!is.numeric(sigma) || sigma <= 0) {
      stop("Normal prior requires sd > 0.", call. = FALSE)
    }
    return(stats::dnorm(x, mean = mu, sd = sigma))
  }

  if (identical(dist, "uniform")) {
    rng <- range(effect_grid)
    lo  <- spec$min %||% rng[1L]
    hi  <- spec$max %||% rng[2L]
    if (lo >= hi) stop("Uniform prior requires min < max.", call. = FALSE)
    return(stats::dunif(x, min = lo, max = hi))
  }

  if (identical(dist, "beta")) {
    if (!is.null(spec$shape1) && !is.null(spec$shape2)) {
      a <- spec$shape1
      b <- spec$shape2
      if (a <= 0 || b <= 0) stop("Beta prior requires shape1 > 0 and shape2 > 0.", call. = FALSE)
    } else if (!is.null(spec$mode) && !is.null(spec$n)) {
      mode_val <- spec$mode
      n_val    <- spec$n
      if (mode_val <= 0 || mode_val >= 1) {
        stop("Beta prior (mode/n parameterisation) requires 0 < mode < 1.", call. = FALSE)
      }
      if (n_val <= 2) {
        stop("Beta prior (mode/n parameterisation) requires n > 2.", call. = FALSE)
      }
      a <- mode_val * (n_val - 2) + 1
      b <- (1 - mode_val) * (n_val - 2) + 1
    } else {
      stop("Beta prior requires either {shape1, shape2} or {mode, n} parameters.",
           call. = FALSE)
    }
    # Rescale x to [0, 1] using the effect_grid range
    rng <- range(effect_grid)
    if (rng[1L] < 0 || rng[2L] > 1) {
      if (diff(rng) < .Machine$double.eps * 100) {
        stop("All effect values are identical; cannot rescale for beta prior.", call. = FALSE)
      }
      x_scaled <- (x - rng[1L]) / (rng[2L] - rng[1L])
    } else {
      x_scaled <- x
    }
    # Density is zero outside [0, 1] on the rescaled scale; dbeta handles this
    return(stats::dbeta(x_scaled, shape1 = a, shape2 = b))
  }

  stop("Unsupported distribution '", dist,
       "'. Use 'normal', 'uniform', or 'beta'.", call. = FALSE)
}
