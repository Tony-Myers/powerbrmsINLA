# File: R/utils-helpers.R
# Combined internal helpers for brmsINLApower package.

#' Internal Coalesce Operator
#' Returns the left-hand side if it is not NULL, otherwise the right-hand side.
#' @param x Left-hand value.
#' @param y Right-hand value.
#' @return If x is not NULL, returns x; else y.
#' @keywords internal
#' @name or_or
#' @export
`%||%` <- rlang::`%||%`


#' Compute Mean Assurance for a Given Metric (Multi-Effect Compatible)
#' Summarises simulation results and computes proportion passing for decision rule metric.
#' @param df Data frame containing simulation results with columns n, effect cols, ok, and metric columns.
#' @param metric One of "direction", "threshold", "rope".
#' @param prob_threshold Numeric, threshold for "threshold" and "rope" metrics.
#' @param rope_rule Reserved for future use.
#' @param direction_p Numeric cutoff for "direction" metric (default 0.5).
#' @param effect_cols Optional character vector of effect columns to group by. If NULL, auto-detects "true_effect" or others.
#' @return Tibble grouped by n and effects with assurance column.
#' @keywords internal
.compute_assurance <- function(
    df,
    metric,
    prob_threshold,
    rope_rule = c(">="),
    direction_p = 0.5,
    effect_cols = NULL
) {
  metric <- match.arg(metric, c("direction", "threshold", "rope"))
  colname <- switch(metric,
                    direction = "post_prob_direction",
                    threshold = "post_prob_threshold",
                    rope = "post_prob_rope")

  if (is.null(effect_cols)) {
    nonmeta_cols <- setdiff(names(df), c(
      "n", "ok", "sim",
      "post_prob_direction", "post_prob_threshold", "post_prob_rope",
      "ci_width", "ci_lower", "ci_upper",
      "bf10", "log10_bf10", "mean_log10_bf",
      "power_direction", "power_threshold", "power_rope",
      "avg_post_prob_direction", "avg_post_prob_threshold", "avg_post_prob_rope",
      "bf_hit_3", "bf_hit_10", "true_effect"
    ))
    effect_cols <- if ("true_effect" %in% names(df)) "true_effect"
    else if (length(nonmeta_cols) > 0) nonmeta_cols
    else "true_effect"
  }

  if (!all(c("n", colname, "ok") %in% names(df))) {
    stop("Expected columns 'n', effect columns, 'ok', and ", colname, " in results.")
  }

  df <- df[df$ok, , drop = FALSE]

  if (metric == "direction") {
    df$pass <- df[[colname]] >= direction_p
  } else {
    df$pass <- df[[colname]] >= prob_threshold
  }

  dplyr::summarise(
    dplyr::group_by(df, n, !!!rlang::syms(effect_cols)),
    assurance = mean(pass, na.rm = TRUE),
    .groups = "drop"
  )
}


#' Wilson Confidence Interval Early Stopping Rule
#' Determines whether to stop early based on Wilson binomial confidence interval.
#' @param hits Number of successes observed
#' @param trials Total number of trials
#' @param target Target proportion
#' @param margin Margin around target for stopping
#' @param conf Confidence level for Wilson CI
#' @return List with stop (logical) and ci (numeric vector)
#' @keywords internal
.should_stop_binom <- function(hits, trials, target, margin = 0.02, conf = 0.95) {
  if (trials == 0) return(list(stop = FALSE, ci = c(0, 1)))

  p_hat <- hits / trials
  z <- qnorm(1 - (1 - conf) / 2)

  # Wilson confidence interval
  denom <- 1 + z^2 / trials
  center <- (p_hat + z^2 / (2 * trials)) / denom
  half_width <- z * sqrt((p_hat * (1 - p_hat) + z^2 / (4 * trials)) / trials) / denom

  ci_lower <- center - half_width
  ci_upper <- center + half_width

  # Stop if CI is entirely above target + margin or below target - margin
  stop_high <- ci_lower > target + margin
  stop_low <- ci_upper < target - margin

  list(stop = stop_high || stop_low, ci = c(ci_lower, ci_upper))
}


#' Determine ggplot2 Line Width Argument Name by Version
#' Returns the correct argument name for line width in ggplot2,
#' depending on package version ("linewidth" for >= 3.4.0, else "size").
#' @return Character string of argument name.
#' @keywords internal
.gg_line_arg <- function() {
  if (utils::packageVersion("ggplot2") >= "3.4.0") "linewidth" else "size"
}


#' Add Contour Lines to a ggplot2 Plot
#' Wrapper around `geom_contour` with preset defaults for colour, alpha, width.
#' Uses the correct linewidth/size argument depending on ggplot2 version.
#' @param mapping Mapping aesthetic.
#' @param data Data frame.
#' @param breaks Break points for contours.
#' @param colour Colour of contour lines.
#' @param alpha Transparency level.
#' @param width Line width.
#' @param bins Number of bins for contour fill.
#' @return A ggplot2 layer adding contour lines.
#' @keywords internal
.add_contour_lines <- function(mapping = NULL, data = NULL,
                               breaks = NULL, colour = "white",
                               alpha = 0.3, width = 0.2, bins = NULL) {
  arg <- .gg_line_arg()
  args <- list(mapping = mapping,
               data = data,
               breaks = breaks,
               colour = colour,
               alpha = alpha,
               bins = bins)
  args[[arg]] <- width
  do.call(ggplot2::geom_contour, args)
}


#' Create a ggplot2 Point Layer with Version-Compatible Width
#' Creates a `geom_point` with a width argument adapted to ggplot2 version.
#' @param mapping Mapping aesthetic.
#' @param data Data frame.
#' @param ... Additional parameters passed to `geom_point`.
#' @param width Numeric line width for points, default 1.5.
#' @return ggplot2 layer for points.
#' @keywords internal
.geom_point_lw <- function(mapping = NULL, data = NULL, ..., width = 1.5) {
  arg <- .gg_line_arg()
  args <- c(list(mapping = mapping, data = data, ...), setNames(list(width), arg))
  do.call(ggplot2::geom_point, args)
}


#' Create a ggplot2 Line Layer with Version-Compatible Width
#' Creates a `geom_line` with a width argument adapted to ggplot2 version.
#' @param mapping Mapping aesthetic.
#' @param data Data frame.
#' @param ... Additional parameters passed to `geom_line`.
#' @param width Numeric line width for lines, default 1.
#' @return ggplot2 layer for lines.
#' @keywords internal
.geom_line_lw <- function(mapping = NULL, data = NULL, ..., width = 1) {
  arg <- .gg_line_arg()
  args <- c(list(mapping = mapping, data = data, ...), setNames(list(width), arg))
  do.call(ggplot2::geom_line, args)
}


#' Scale Fill for Viridis Discrete Data
#' @param name Character legend title (default "Assurance")
#' @return ggplot2 fill scale object
#' @keywords internal
.scale_fill_viridis_discrete <- function(name = "Assurance") {
  if ("scale_fill_viridis_d" %in% getNamespaceExports("ggplot2")) {
    ggplot2::scale_fill_viridis_d(name = name)
  } else {
    ggplot2::scale_fill_stepsn(colours = viridisLite::viridis(12), name = name)
  }
}


#' Scale Fill for Viridis Continuous Data
#' @param name Legend title
#' @param limits Numeric vector length 2 for limits
#' @param breaks Numeric vector for breaks
#' @param labels Function or vector for labels
#' @return ggplot2 fill scale object
#' @keywords internal
.scale_fill_viridis_continuous <- function(
    name = "Assurance",
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    labels = scales::percent_format(accuracy = 1)
) {
  if ("scale_fill_viridis_c" %in% getNamespaceExports("ggplot2")) {
    ggplot2::scale_fill_viridis_c(name = name, limits = limits,
                                  breaks = breaks, labels = labels)
  } else {
    ggplot2::scale_fill_gradientn(
      colours = viridisLite::viridis(256), name = name,
      limits = limits, breaks = breaks, labels = labels
    )
  }
}

#' Validate an SD Specification for error_sd or group_sd
#'
#' Checks whether the input is a valid positive numeric scalar or one of the
#' supported distributional list specifications.  Called automatically by
#' [brms_inla_power()] before the simulation loop; can also be called
#' directly for interactive validation.
#'
#' Supported distributional formats:
#' * `list(dist = "halfnormal", sd = X, location = Y)` — draws
#'   `|Normal(location, sd)|`; `location` defaults to 0.
#' * `list(dist = "lognormal", meanlog = X, sdlog = Y)` — draws from a
#'   log-normal distribution.
#' * `list(dist = "uniform", min = X, max = Y)` — draws from Uniform(min, max);
#'   requires `min >= 0`.
#'
#' @param x A positive numeric scalar **or** a named list with element `dist`.
#' @param arg_name Character string used in error messages (default `"x"`).
#'
#' @return `x`, invisibly.  Called for its side effects (stopping on invalid
#'   input).
#' @export
validate_sd_spec <- function(x, arg_name = "x") {
  # ---- Scalar path ----------------------------------------------------------
  if (is.numeric(x)) {
    if (length(x) == 1L && is.finite(x) && x > 0) return(invisible(x))
    stop(
      sprintf(
        "`%s` must be a single positive finite number or a distributional list. Got: numeric of length %d.",
        arg_name, length(x)
      ),
      call. = FALSE
    )
  }

  # ---- List path ------------------------------------------------------------
  if (!is.list(x)) {
    stop(
      sprintf(
        "`%s` must be a positive numeric scalar or a list with a `dist` element. Got class: %s.",
        arg_name, paste(class(x), collapse = "/")
      ),
      call. = FALSE
    )
  }

  dist <- x[["dist"]]
  if (is.null(dist) || !is.character(dist) || length(dist) != 1L) {
    stop(
      sprintf(
        '`%s` list must have a character `dist` element. Supported: "halfnormal", "lognormal", "uniform".',
        arg_name
      ),
      call. = FALSE
    )
  }

  if (dist == "halfnormal") {
    if (is.null(x[["sd"]]) || !is.numeric(x[["sd"]]) || length(x[["sd"]]) != 1L || x[["sd"]] <= 0) {
      stop(
        sprintf('`%s = list(dist = "halfnormal", ...)` requires a single numeric `sd > 0`.', arg_name),
        call. = FALSE
      )
    }
    loc <- x[["location"]] %||% 0
    if (!is.numeric(loc) || length(loc) != 1L || !is.finite(loc)) {
      stop(
        sprintf('`%s = list(dist = "halfnormal", ...)` requires `location` to be a finite numeric scalar.', arg_name),
        call. = FALSE
      )
    }
  } else if (dist == "lognormal") {
    if (is.null(x[["meanlog"]]) || !is.numeric(x[["meanlog"]]) || length(x[["meanlog"]]) != 1L ||
        !is.finite(x[["meanlog"]])) {
      stop(
        sprintf('`%s = list(dist = "lognormal", ...)` requires a finite numeric `meanlog`.', arg_name),
        call. = FALSE
      )
    }
    if (is.null(x[["sdlog"]]) || !is.numeric(x[["sdlog"]]) || length(x[["sdlog"]]) != 1L ||
        x[["sdlog"]] <= 0) {
      stop(
        sprintf('`%s = list(dist = "lognormal", ...)` requires a single numeric `sdlog > 0`.', arg_name),
        call. = FALSE
      )
    }
  } else if (dist == "uniform") {
    if (is.null(x[["min"]]) || !is.numeric(x[["min"]]) || length(x[["min"]]) != 1L) {
      stop(
        sprintf('`%s = list(dist = "uniform", ...)` requires a numeric `min`.', arg_name),
        call. = FALSE
      )
    }
    if (is.null(x[["max"]]) || !is.numeric(x[["max"]]) || length(x[["max"]]) != 1L) {
      stop(
        sprintf('`%s = list(dist = "uniform", ...)` requires a numeric `max`.', arg_name),
        call. = FALSE
      )
    }
    if (x[["min"]] < 0) {
      stop(
        sprintf('`%s = list(dist = "uniform", ...)` requires `min >= 0`.', arg_name),
        call. = FALSE
      )
    }
    if (x[["max"]] <= x[["min"]]) {
      stop(
        sprintf('`%s = list(dist = "uniform", ...)` requires `max > min`.', arg_name),
        call. = FALSE
      )
    }
  } else {
    stop(
      sprintf(
        '`%s` has unsupported `dist = "%s"`. Supported distributions: "halfnormal", "lognormal", "uniform".',
        arg_name, dist
      ),
      call. = FALSE
    )
  }

  invisible(x)
}


#' Draw One Sample from an SD Specification
#'
#' Internal helper used by [brms_inla_power()] to draw a per-iteration value
#' from a distributional SD specification produced by [validate_sd_spec()].
#' Only called when `is.list(spec)`.
#'
#' @param spec A validated distributional list (see [validate_sd_spec()]).
#' @return A single positive numeric draw.
#' @keywords internal
.sample_sd_spec <- function(spec) {
  dist <- spec[["dist"]]
  if (dist == "halfnormal") {
    loc <- spec[["location"]] %||% 0
    abs(stats::rnorm(1L, mean = loc, sd = spec[["sd"]]))
  } else if (dist == "lognormal") {
    stats::rlnorm(1L, meanlog = spec[["meanlog"]], sdlog = spec[["sdlog"]])
  } else {  # uniform
    stats::runif(1L, min = spec[["min"]], max = spec[["max"]])
  }
}
