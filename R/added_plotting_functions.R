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
    data,
    metric,
    effect1,
    effect2,
    n        = NULL,
    line     = FALSE,
    facet_by = NULL,
    agg_fun  = mean,
    title    = NULL,
    subtitle = NULL
) {
  stopifnot(is.data.frame(data))
  if (!metric %in% names(data))  stop("Metric '", metric, "' not found.",  call. = FALSE)
  if (!effect1 %in% names(data)) stop("effect1 '", effect1, "' not found.", call. = FALSE)
  if (!effect2 %in% names(data)) stop("effect2 '", effect2, "' not found.", call. = FALSE)
  
  # optional filter on n
  if (!is.null(n) && "n" %in% names(data)) {
    data <- data[data$n %in% n, , drop = FALSE]
  }
  
  facet_names <- facet_by %||% character(0)
  facet_names <- intersect(facet_names, names(data))
  
  group_cols <- c(effect1, effect2, facet_names)
  
  data_sum <- data |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
    dplyr::summarise(
      value = agg_fun(.data[[metric]], na.rm = TRUE),
      .groups = "drop"
    )
  
  eff1_sym <- rlang::sym(effect1)
  eff2_sym <- rlang::sym(effect2)
  
  if (line) {
    # effect1 on x, value on y, colour by effect2
    p <- ggplot2::ggplot(
      data_sum,
      ggplot2::aes(
        x      = !!eff1_sym,
        y      = value,
        colour = !!eff2_sym,
        group  = !!eff2_sym
      )
    ) +
      ggplot2::geom_line(linewidth = 0.9) +
      ggplot2::geom_point(size = 2) +
      ggplot2::scale_colour_viridis_d(name = effect2)
  } else {
    # heatmap: effect1 x effect2, fill = value
    p <- ggplot2::ggplot(
      data_sum,
      ggplot2::aes(
        x    = !!eff1_sym,
        y    = !!eff2_sym,
        fill = value
      )
    ) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_viridis_c(name = metric)
  }
  
  p <- p +
    ggplot2::labs(
      x        = effect1,
      y        = if (line) metric else effect2,
      title    = title,
      subtitle = subtitle
    ) +
    ggplot2::theme_minimal()
  
  if (length(facet_names) == 1L) {
    p <- p + ggplot2::facet_wrap(as.formula(paste("~", facet_names[1L])))
  } else if (length(facet_names) > 1L) {
    p <- p + ggplot2::facet_grid(
      as.formula(paste(paste(facet_names, collapse = " + "), "~ ."))
    )
  }
  
  p
}