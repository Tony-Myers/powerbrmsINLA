#' @importFrom dplyr filter group_by summarise across all_of bind_rows slice_min pull
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_gradient labs theme_minimal
#' @importFrom ggplot2 theme element_text geom_contour geom_point geom_line
#' @importFrom ggplot2 scale_fill_viridis_c scale_fill_viridis_d scale_fill_gradientn
#' @importFrom ggplot2 geom_ribbon geom_smooth coord_cartesian guides guide_legend
#' @importFrom rlang .data !!! syms
#' @importFrom tibble tibble
#' @importFrom scales percent_format
#' @importFrom viridisLite viridis
#' @importFrom brms prior
#' @importFrom stats gaussian binomial poisson rnorm rbinom rpois
#' @importFrom stats dnorm pnorm qnorm sd var as.formula terms
#' @importFrom stats dbeta qbeta optimize setNames
#' @importFrom utils packageVersion flush.console
#' @importFrom magrittr %>%
NULL

# Suppress R CMD check notes about global variables
if(getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "val", "ribbon", "n", "treatment", "effect", "power_direction",
    "power_threshold", "avg_ci_width", "assurance", "ci_width",
    "nsims_ok", "bf10", "log10_bf10", "sim", "ok", ".", "?"
  ))
}
