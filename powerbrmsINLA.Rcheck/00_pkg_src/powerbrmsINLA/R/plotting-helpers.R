#' Compute Mean Assurance for a Given Metric (Modern, Multi-Effect Compatible)
#'
#' Computes the mean assurance (proportion passing) for a given decision metric across grouped cells.
#'
#' Determine ggplot2 Line Width Argument Name by Version
#'
#' Returns the correct argument name for line width in ggplot2,
#' depending on package version ("linewidth" for >= 3.4.0, else "size").
#'
#' @return Character string of argument name.
#' @keywords internal
.gg_line_arg <- function() {
  if (utils::packageVersion("ggplot2") >= "3.4.0") "linewidth" else "size"
}


#' Add Contour Lines to a ggplot2 Plot
#'
#' Wrapper around `geom_contour` with preset defaults for color, alpha, width.
#' Uses the correct linewidth/size argument depending on ggplot2 version.
#'
#' @param mapping Mapping aesthetic.
#' @param data Data frame.
#' @param breaks Break points for contours.
#' @param colour Color of contour lines.
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
#'
#' Creates a `geom_point` with a width argument adapted to ggplot2 version.
#'
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
#'
#' Creates a `geom_line` with a width argument adapted to ggplot2 version.
#'
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


