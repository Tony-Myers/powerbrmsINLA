#' Plot Assurance with Robustness Ribbon (Multi-Effect Grid Friendly)
#'
#' Compares assurance results from multiple scenarios by showing the range
#' ("ribbon") of values across scenarios for each sample size and effect grid variable.
#'
#' @param power_results_list Named list of results objects from `brms_inla_power` or sequential/two-stage variants.
#' @param metric Which assurance metric to compute: "precision", "direction", "threshold", or "bf".
#' @param x_effect Name of effect grid column for x-axis (default: first detected grid column).
#' @param facet_by Optional effect grid column(s) to facet by.
#' @param precision_target CI width target if metric="precision".
#' @param p_star Posterior probability threshold for "direction"/"threshold".
#' @param bf_threshold BF10 threshold for "bf".
#' @param effect_filters Optional named list for filtering rows (e.g. list(treatment=0)).
#' @param effect_weights Optional named numeric vector for averaging over grid values.
#' @param show_individual_scenarios Logical; if TRUE, overlay each scenario's curve.
#' @param title,subtitle Optional plot labels.
#' @return A ggplot object.
#' @export
plot_assurance_with_robustness <- function(
    power_results_list,
    metric = c("precision","direction","threshold","bf"),
    x_effect = NULL,
    facet_by = NULL,
    precision_target = NULL,
    p_star = 0.95,
    bf_threshold = 10,
    effect_filters = NULL,
    effect_weights = NULL,
    show_individual_scenarios = FALSE,
    title = NULL,
    subtitle = NULL
) {
  metric <- match.arg(metric)
  scenario_names <- names(power_results_list)
  if (is.null(scenario_names) || any(scenario_names == "")) {
    scenario_names <- paste0("Scenario ", seq_along(power_results_list))
  }
  # Assemble scenario-wise assurance curves
  per_scen <- lapply(seq_along(power_results_list), function(i) {
    pr <- power_results_list[[i]]
    df <- pr$results %>% dplyr::filter(ok)
    # Detect effect grid columns
    effect_names <- setdiff(
      colnames(df),
      c("n", "ok", "post_prob_direction", "post_prob_threshold", "post_prob_rope", "bf10", "ci_width",
        "ci_lower", "ci_upper", "log10_bf10")
    )
    # Filter to user-requested effect values
    if (!is.null(effect_filters) && length(effect_filters)) {
      for (ef in names(effect_filters)) {
        df <- df[df[[ef]] == effect_filters[[ef]], , drop = FALSE]
      }
    }
    # Grouping
    x_effect_cur <- x_effect %||% effect_names[1]
    facet_names <- facet_by %||% character(0)
    facet_names <- intersect(facet_names, colnames(df))
    metric_col <- switch(metric,
                         precision = "ci_width",
                         direction = "post_prob_direction",
                         threshold = "post_prob_threshold",
                         bf = "bf10"
    )
    # Compute groupwise assurance for this scenario
    tmp <- dplyr::group_by(df, n, !!rlang::sym(x_effect_cur), !!!rlang::syms(facet_names))
    if (metric == "precision") {
      pt <- precision_target %||% pr$settings$precision_target
      if (is.null(pt)) stop("precision_target must be supplied for precision.")
      tmp <- tmp %>% dplyr::summarise(assurance = mean(.data[[metric_col]] <= pt, na.rm = TRUE), .groups = "drop")
    } else if (metric == "bf") {
      tmp <- tmp %>% dplyr::summarise(assurance = mean(.data[[metric_col]] >= bf_threshold, na.rm = TRUE), .groups = "drop")
    } else {
      tmp <- tmp %>% dplyr::summarise(assurance = mean(.data[[metric_col]] >= p_star, na.rm = TRUE), .groups = "drop")
    }
    # Optional weights over x_effect
    if (!is.null(effect_weights)) {
      w <- tibble::tibble(val = as.numeric(names(effect_weights)), w = as.numeric(effect_weights))
      names(w)[1] <- x_effect_cur
      tmp <- dplyr::left_join(tmp, w, by = x_effect_cur) %>%
        dplyr::group_by(n, !!!rlang::syms(facet_names)) %>%
        dplyr::summarise(
          assurance = stats::weighted.mean(assurance, w, na.rm = TRUE),
          .groups = "drop"
        )
    }
    tmp$scenario <- scenario_names[i]
    tmp
  }) %>% dplyr::bind_rows()
  # Ribbon: min-max assurance across scenarios per sample size/ grouping
  group_syms <- unique(
    c(
      "n",
      x_effect %||% names(per_scen)[1],
      intersect(facet_by, colnames(per_scen))
    )
  )
  ribbon <- per_scen %>%
    dplyr::group_by(!!!rlang::syms(group_syms)) %>%
    dplyr::summarise(
      ymin = min(assurance, na.rm = TRUE),
      ymax = max(assurance, na.rm = TRUE),
      ymid = mean(assurance, na.rm = TRUE),
      .groups = "drop"
    )
  # Plot
  p <- ggplot2::ggplot(ribbon, ggplot2::aes(
    x = !!rlang::sym(x_effect %||% names(per_scen)[1]),
    y = ymid, group = n, color = factor(n))
  ) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = ymin, ymax = ymax), alpha = 0.15) +
    .geom_line_lw(width = 1)
  # Optionally add scenario curves
  if (isTRUE(show_individual_scenarios)) {
    p <- p +
      .geom_line_lw(data = per_scen,
                    mapping = ggplot2::aes(y = assurance, color = scenario),
                    width = 0.8)
    if ("scale_colour_viridis_d" %in% getNamespaceExports("ggplot2")) {
      p <- p + ggplot2::scale_colour_viridis_d(name = "Scenario", option = "plasma")
    } else {
      p <- p + ggplot2::scale_colour_manual(values = viridisLite::viridis(length(unique(per_scen$scenario))))
    }
  } else {
    p <- p + ggplot2::scale_color_viridis_d(name = "Sample size")
  }
  # Y label according to metric
  y_lab <- switch(metric,
                  precision = "Assurance P(width <= target)",
                  direction = paste0("Assurance P{ post-prob >= ", p_star, " }"),
                  threshold = paste0("Assurance P{ post-prob >=", p_star, " }"),
                  bf = paste0("Assurance P{ BF10 >= ", bf_threshold, " }")
  )
  p <- p +
    ggplot2::scale_y_continuous(limits = c(0,1), labels = scales::percent_format(accuracy = 1)) +
    ggplot2::labs(
      x = x_effect %||% names(per_scen)[1],
      y = y_lab,
      title = title %||% "Assurance with robustness ribbon across scenarios",
      subtitle = subtitle
    ) +
    ggplot2::theme_minimal()
  # Facet support
  if (!is.null(facet_by) && length(facet_by) == 1) {
    p <- p + ggplot2::facet_wrap(as.formula(paste("~", facet_by)))
  } else if (!is.null(facet_by) && length(facet_by) > 1) {
    p <- p + ggplot2::facet_grid(as.formula(paste(paste(facet_by, collapse = " + "), "~ .")))
  }
  return(p)
}
