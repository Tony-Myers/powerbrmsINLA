#' Plot Interaction Assurance Surface/Heatmap/Lines (Multi-Effect Grid Friendly)
#'
#' Visualizes a metric (e.g., assurance) as a function of two effect grid variables
#' for a fixed sample size or averaged over n. Allows line, heatmap, or contour modes.
#'
#' @param data Data frame (typically power_results$summary).
#' @param metric Name of the summary column to plot, e.g. "power_direction", "power_threshold".
#' @param effect1 Name of effect grid column for x-axis.
#' @param effect2 Name of effect grid column for y-axis or color/facets.
#' @param n Optional sample size to filter to (else averages/plots all n's).
#' @param line Logical; if TRUE, make a lineplot (effect1 on x, one line for each effect2).
#' If FALSE, make a heatmap or contour.
#' @param facet_by Optional grid column(s) to facet by.
#' @param agg_fun Aggregation function if multiple entries per cell (default = mean).
#' @param title,subtitle Optional plot labels.
#' @return A ggplot object.
#' @export

plot_interaction_surface <- function(
    data, metric,
    effect1, effect2,
    n = NULL, line = FALSE,
    facet_by = NULL,
    agg_fun = mean,
    title = NULL,
    subtitle = NULL
) {
  stopifnot(all(c(effect1, effect2, metric) %in% names(data)))
  df <- data
  # Filter by n if provided
  if (!is.null(n) && "n" %in% names(df)) df <- df[df$n == n, , drop = FALSE]
  # Aggregate if needed
  group_syms <- c(effect1, effect2, facet_by)
  if (!is.null(n) && "n" %in% names(df)) group_syms <- c(group_syms, "n")
  # Handle duplicate grid points
  df <- df %>%
    dplyr::group_by(!!!rlang::syms(group_syms)) %>%
    dplyr::summarise(val = agg_fun(.data[[metric]], na.rm = TRUE), .groups = "drop")
  p <- if (line) {
    # Lines for each effect2, colored by effect2/facets
    base <- ggplot2::ggplot(
      df, ggplot2::aes_string(x = effect1, y = "val", group = effect2, color = effect2)
    ) +
      ggplot2::geom_line() +
      ggplot2::geom_point() +
      ggplot2::labs(
        x = effect1, y = metric,
        color = effect2,
        title = title, subtitle = subtitle
      )
    if (!is.null(facet_by)) {
      base <- base + ggplot2::facet_wrap(facet_by)
    }
    base + ggplot2::theme_minimal()
  } else {
    # Heatmap (tile) + optional contour for smooth surfaces
    base <- ggplot2::ggplot(
      df, ggplot2::aes_string(x = effect1, y = effect2, fill = "val")
    ) +
      ggplot2::geom_tile() +
      .scale_fill_viridis_continuous(name = metric, limits = c(0, 1),
                                     breaks = seq(0, 1, 0.2), labels = scales::percent_format(accuracy = 1)) +
      ggplot2::labs(
        x = effect1, y = effect2,
        title = title, subtitle = subtitle
      )
    if (!is.null(facet_by)) {
      base <- base + ggplot2::facet_wrap(facet_by)
    }
    base + ggplot2::theme_minimal()
  }
  return(p)
}
