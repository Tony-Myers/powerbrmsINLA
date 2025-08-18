
#' Plot Bayes Factor Assurance Curve (Multi-Effect Grid Friendly)
#'
#' Plots the proportion of simulations in which BF10 meets or exceeds a threshold,
#' grouped by any effect grid variable(s) and sample size.
#'
#' @param power_results List returned by `brms_inla_power*` or two-stage variant.
#' @param bf_threshold Numeric; BF10 threshold to count as a "success" (default: 3).
#' @param x_effect Name of effect grid column for x-axis (default: first detected grid column).
#' @param facet_by Optional grid column(s) for faceting.
#' @param effect_filters Optional named list to restrict/show only selected grid rows, e.g. list(treatment=0).
#' @param effect_weights Optional named numeric vector of weights for selected x_effect values.
#' @param title,subtitle Optional plot labels.
#' @return ggplot object.
#' @export
plot_bf_assurance_curve <- function(
    power_results,
    bf_threshold = 3,
    x_effect = NULL,
    facet_by = NULL,
    effect_filters = NULL,
    effect_weights = NULL,
    title = NULL,
    subtitle = NULL
) {
  df <- dplyr::filter(power_results$results, ok)
  # Detect effect grid columns
  effect_names <- setdiff(
    colnames(df),
    c("n", "ok", "bf10", "post_prob_direction", "post_prob_threshold",
      "post_prob_rope", "ci_width", "ci_lower", "ci_upper", "log10_bf10")
  )
  # (Handle user filter, e.g. restrict to treatment==0 only)
  if (!is.null(effect_filters) && length(effect_filters)) {
    for (ef in names(effect_filters)) {
      df <- df[df[[ef]] == effect_filters[[ef]], , drop = FALSE]
    }
  }
  # Main x-axis effect (default: first grid column)
  x_effect <- x_effect %||% effect_names[1]
  stopifnot(!is.null(x_effect), x_effect %in% names(df))
  facet_names <- facet_by %||% character(0)
  facet_names <- intersect(facet_names, colnames(df))
  # Group by: sample size + x_effect + facet(s)
  tmp <- df %>%
    dplyr::group_by(n, !!rlang::sym(x_effect), !!!rlang::syms(facet_names)) %>%
    dplyr::summarise(
      assurance = mean(bf10 >= bf_threshold, na.rm = TRUE),
      .groups = "drop"
    )
  # If requested: weighted means over x_effect
  if (!is.null(effect_weights) && is.numeric(effect_weights)) {
    w <- tibble::tibble(val = as.numeric(names(effect_weights)), w = as.numeric(effect_weights))
    names(w)[1] <- x_effect
    tmp <- dplyr::left_join(tmp, w, by = x_effect) %>%
      dplyr::group_by(n, !!!rlang::syms(facet_names)) %>%
      dplyr::summarise(
        assurance = stats::weighted.mean(assurance, w, na.rm = TRUE),
        .groups = "drop"
      )
  }
  # Build plot: x_effect vs assurance
  p <- ggplot2::ggplot(
    tmp,
    ggplot2::aes(x = !!rlang::sym(x_effect), y = assurance, group = n, color = factor(n))
  ) +
    .geom_line_lw(width = 1) + ggplot2::geom_point(size = 2) +
    ggplot2::scale_color_viridis_d(name = "Sample size") +
    ggplot2::scale_y_continuous(limits = c(0,1), labels = scales::percent_format(accuracy = 1)) +
    ggplot2::labs(
      x = x_effect,
      y = paste0("Assurance: P{ BF10 >= ", bf_threshold, " }"),
      title = title %||% "Bayes-factor Assurance Curve",
      subtitle = subtitle
    ) +
    ggplot2::theme_minimal()
  # Faceting support
  if (length(facet_names) == 1) {
    p <- p + ggplot2::facet_wrap(as.formula(paste("~", facet_names)))
  } else if (length(facet_names) > 1) {
    p <- p + ggplot2::facet_grid(as.formula(paste(paste(facet_names, collapse = " + "), "~ .")))
  }
  return(p)
}

#' Plot Expected Evidence (mean log10 BF10, Multi-Effect Grid Friendly)
#'
#' Plots the average log10 BF10 against any effect grid variable, grouped/faceted.
#'
#' @param power_results Simulation results from a `brms_inla_power*` function with `compute_bayes_factor = TRUE`.
#' @param x_effect Name of effect grid column for x-axis (default: first grid column).
#' @param facet_by Optional grid column(s) to facet by (default: NULL).
#' @param n Optional sample size to filter to (NULL means plot all; else one curve per grid/facet).
#' @param agg_fun Aggregation function if >1 entries per cell (default: mean).
#' @param title,subtitle Optional plot labels.
#' @return A ggplot object.
#' @export
plot_bf_expected_evidence <- function(
    power_results,
    x_effect = NULL,
    facet_by = NULL,
    n = NULL,
    agg_fun = mean,
    title = NULL,
    subtitle = NULL
) {
  s <- power_results$summary
  # Detect effect grid columns
  effect_names <- setdiff(
    colnames(s),
    c("n", "power_direction", "power_threshold", "power_rope",
      "avg_ci_width", "ci_coverage", "ciw_q05", "ciw_q25", "ciw_q50", "ciw_q75", "ciw_q95",
      "avg_post_prob_direction", "avg_post_prob_threshold", "avg_post_prob_rope",
      "bf_hit_3", "bf_hit_10", "mean_log10_bf", "nsims_ok")
  )
  x_effect <- x_effect %||% effect_names[1]
  stopifnot(!is.null(x_effect), x_effect %in% colnames(s))
  facet_names <- facet_by %||% character(0)
  facet_names <- intersect(facet_names, colnames(s))
  group_syms <- c(facet_names)
  # Optionally filter for n
  if (!is.null(n) && "n" %in% names(s)) {
    s <- s[s$n == n, , drop = FALSE]
    group_syms <- c(group_syms, "n")
  }
  # Group/aggregate by grid column(s)
  tmp <- s %>%
    dplyr::group_by(!!rlang::sym(x_effect), !!!rlang::syms(group_syms)) %>%
    dplyr::summarise(mean_log10_bf = agg_fun(mean_log10_bf, na.rm = TRUE), .groups = "drop")
  # Plot: lines for group(s)/facet(s)
  p <- ggplot2::ggplot(
    tmp,
    ggplot2::aes(x = !!rlang::sym(x_effect), y = mean_log10_bf, color = interaction(!!!rlang::syms(group_syms)))
  ) +
    .geom_line_lw(width = 1) +
    ggplot2::labs(
      x = x_effect,
      y = "Expected log10 BF10",
      color = paste(group_syms, collapse = "+"),
      title = title %||% "Expected Evidence Accumulation",
      subtitle = subtitle
    ) +
    ggplot2::theme_minimal()
  # Facet support
  if (length(facet_names) == 1) {
    p <- p + ggplot2::facet_wrap(as.formula(paste("~", facet_names)))
  } else if (length(facet_names) > 1) {
    p <- p + ggplot2::facet_grid(as.formula(paste(paste(facet_names, collapse = " + "), "~ .")))
  }
  return(p)
}

#' Plot Bayes Factor Heatmap (mean log10 BF10, Multi-Effect Grid Friendly)
#'
#' Heatmap of mean log10 BF10 as a function of two effect grid columns (x/y), with optional faceting.
#'
#' @param power_results Simulation results from a `brms_inla_power*` function with `compute_bayes_factor = TRUE`.
#' @param x_effect Name of effect grid column for x-axis (default: first grid column).
#' @param y_effect Name of effect grid column for y-axis (default: "n").
#' @param facet_by Optional column(s) to facet by.
#' @param n Optional sample size to filter to (NULL means plot all; else show only that n).
#' @param agg_fun Aggregation function (default: mean).
#' @param title,subtitle Optional plot labels.
#' @return ggplot object.
#' @export
plot_bf_heatmap <- function(
    power_results,
    x_effect = NULL,
    y_effect = "n",
    facet_by = NULL,
    n = NULL,
    agg_fun = mean,
    title = NULL,
    subtitle = NULL
) {
  s <- power_results$summary
  # Detect effect grid columns
  effect_names <- setdiff(
    colnames(s),
    c("n", "power_direction", "power_threshold", "power_rope",
      "avg_ci_width", "ci_coverage", "ciw_q05", "ciw_q25", "ciw_q50", "ciw_q75", "ciw_q95",
      "avg_post_prob_direction", "avg_post_prob_threshold", "avg_post_prob_rope",
      "bf_hit_3", "bf_hit_10", "mean_log10_bf", "nsims_ok")
  )
  x_effect <- x_effect %||% effect_names[1]
  y_effect <- y_effect %||% ("n" %in% names(s) ? "n" : effect_names[2] %||% effect_names[1])
  stopifnot(x_effect %in% names(s), y_effect %in% names(s))
  facet_names <- facet_by %||% character(0)
  facet_names <- intersect(facet_names, colnames(s))
  plot_dat <- s
  # If requested, filter for a fixed n
  if (!is.null(n) && "n" %in% names(plot_dat)) plot_dat <- plot_dat[plot_dat$n == n, , drop = FALSE]
  # Aggregate duplicates (if needed)
  group_syms <- c(x_effect, y_effect, facet_names)
  plot_dat <- plot_dat %>%
    dplyr::group_by(!!!rlang::syms(group_syms)) %>%
    dplyr::summarise(val = agg_fun(mean_log10_bf, na.rm = TRUE), .groups = "drop")
  # Build heatmap (tile)
  p <- ggplot2::ggplot(
    plot_dat,
    ggplot2::aes(x = !!rlang::sym(x_effect), y = !!rlang::sym(y_effect), fill = val)
  ) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_viridis_c(option = "C", name = "E[log10 BF10]") +
    ggplot2::labs(
      x = x_effect,
      y = y_effect,
      title = title %||% "Expected Evidence Heatmap",
      subtitle = subtitle
    ) +
    ggplot2::theme_minimal()
  # Faceting support
  if (length(facet_names) == 1) {
    p <- p + ggplot2::facet_wrap(as.formula(paste("~", facet_names)))
  } else if (length(facet_names) > 1) {
    p <- p + ggplot2::facet_grid(as.formula(paste(paste(facet_names, collapse = " + "), "~ .")))
  }
  return(p)
}

