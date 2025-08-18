#' Plot Bayesian Power / Assurance Contour (Multi-Effect Grid Friendly)
#'
#' Draw a filled contour plot of assurance for a chosen metric,
#' as a function of two effect grid columns and sample size.
#'
#' @param power_results Output from a `brms_inla_power` function.
#' @param power_metric Which metric to plot: `"direction"`, `"threshold"`, or `"rope"`.
#' @param x_effect Name of effect grid column for x-axis (default = first effect).
#' @param y_effect Name of effect grid column for y-axis (default = "n").
#' @param facet_by Optional effect grid column(s) to facet by.
#' @param power_threshold Optional contour line for assurance (default 0.8).
#' @param show_threshold_line Logical; add a red contour at \code{power_threshold}.
#' @param title,subtitle Optional plot labels.
#' @return A ggplot object.
#' @export
plot_power_contour <- function(power_results,
                               power_metric = c("direction", "threshold", "rope"),
                               x_effect = NULL,
                               y_effect = "n",
                               facet_by = NULL,
                               power_threshold = 0.8,
                               show_threshold_line = TRUE,
                               title = NULL,
                               subtitle = NULL) {
  stopifnot("summary" %in% names(power_results))
  power_metric <- match.arg(power_metric)
  power_col <- switch(power_metric,
                      direction = "power_direction",
                      threshold = "power_threshold",
                      rope = "power_rope")
  dat <- power_results$summary
  effect_names <- setdiff(colnames(dat), c("n", names(dat)[startsWith(names(dat), "power_")],
                                           "avg_ci_width", "ci_coverage", "ciw_q05", "ciw_q25",
                                           "ciw_q50", "ciw_q75", "ciw_q95",
                                           "avg_post_prob_direction", "avg_post_prob_threshold", "avg_post_prob_rope",
                                           "bf_hit_3", "bf_hit_10", "mean_log10_bf", "nsims_ok"))
  # Select default x (first effect), y ("n" by default)
  x_effect <- x_effect %||% (effect_names[1] %||% "n")
  y_effect <- y_effect %||% "n"
  stopifnot(x_effect %in% colnames(dat), y_effect %in% colnames(dat))
  plot_dat <- dat
  facet_list <- facet_by %||% character(0)
  facet_names <- intersect(facet_list, colnames(dat))
  aes_args <- list(x = rlang::sym(x_effect),
                   y = rlang::sym(y_effect),
                   z = rlang::sym(power_col))
  p <- ggplot2::ggplot(plot_dat, do.call(ggplot2::aes, aes_args)) +
    ggplot2::geom_contour_filled(bins = 12, alpha = 0.9) +
    .add_contour_lines(colour = "white", alpha = 0.3, width = 0.2) +
    .scale_fill_viridis_discrete(name = "Assurance") +  # Changed from continuous to discrete
    ggplot2::labs(
      x = x_effect,
      y = y_effect,
      title = title %||% paste("Assurance contour (", power_metric, ")", sep = ""),
      subtitle = subtitle
    ) +
    ggplot2::theme_minimal()
  if (show_threshold_line && !is.null(power_threshold)) {
    p <- p + .add_contour_lines(mapping = ggplot2::aes(z = !!rlang::sym(power_col)),
                                breaks = power_threshold,
                                colour = "red", width = 0.8)
  }
  if (length(facet_names) == 1) {
    p <- p + ggplot2::facet_wrap(as.formula(paste("~", facet_names)))
  } else if (length(facet_names) > 1) {
    p <- p + ggplot2::facet_grid(as.formula(paste(paste(facet_names, collapse = " + "), "~ .")))
  }
  p
}
#' Plot Bayesian Power / Assurance Heatmap (Multi-Effect Grid Friendly)
#'
#' Heatmap of assurance for a chosen metric across two selected effect grid variables and sample sizes.
#'
#' @inheritParams plot_power_contour
#' @return A ggplot object.
#' @export
plot_power_heatmap <- function(power_results,
                               power_metric = c("direction", "threshold", "rope"),
                               x_effect = NULL,
                               y_effect = "n",
                               facet_by = NULL,
                               title = NULL,
                               subtitle = NULL) {
  stopifnot("summary" %in% names(power_results))
  power_metric <- match.arg(power_metric)
  power_col <- switch(power_metric,
                      direction = "power_direction",
                      threshold = "power_threshold",
                      rope = "power_rope")
  dat <- power_results$summary
  effect_names <- setdiff(colnames(dat), c("n", names(dat)[startsWith(names(dat), "power_")],
                                           "avg_ci_width", "ci_coverage", "ciw_q05", "ciw_q25",
                                           "ciw_q50", "ciw_q75", "ciw_q95",
                                           "avg_post_prob_direction", "avg_post_prob_threshold", "avg_post_prob_rope",
                                           "bf_hit_3", "bf_hit_10", "mean_log10_bf", "nsims_ok"))
  x_effect <- x_effect %||% (effect_names[1] %||% "n")
  y_effect <- y_effect %||% "n"
  stopifnot(x_effect %in% colnames(dat), y_effect %in% colnames(dat))
  plot_dat <- dat
  facet_list <- facet_by %||% character(0)
  facet_names <- intersect(facet_list, colnames(dat))
  aes_args <- list(x = rlang::sym(x_effect),
                   y = rlang::sym(y_effect),
                   fill = rlang::sym(power_col))
  p <- ggplot2::ggplot(plot_dat, do.call(ggplot2::aes, aes_args)) +
    ggplot2::geom_tile() +
    .scale_fill_viridis_continuous(name = "Assurance",
                                   limits = c(0, 1),
                                   breaks = seq(0, 1, 0.2),
                                   labels = scales::percent_format(accuracy = 1)) +
    ggplot2::labs(
      x = x_effect,
      y = y_effect,
      title = title %||% paste("Assurance heatmap for", power_metric),
      subtitle = subtitle
    ) +
    ggplot2::theme_minimal()
  if (length(facet_names) == 1) {
    p <- p + ggplot2::facet_wrap(as.formula(paste("~", facet_names)))
  } else if (length(facet_names) > 1) {
    p <- p + ggplot2::facet_grid(as.formula(paste(paste(facet_names, collapse = " + "), "~ .")))
  }
  p
}
