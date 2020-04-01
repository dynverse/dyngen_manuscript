#' A common theme for results in this manuscript
#'
#' @inheritParams ggplot2::theme
#' @importFrom ggplot2 theme
#'
#' @export
theme_common <- function(
  legend.position = "bottom",
  plot.title = element_text(hjust = 0.5),
  plot.subtitle = element_text(hjust= 0.5),
  ...
) {
  ggplot2::theme(
    legend.position = legend.position,
    plot.title = plot.title,
    plot.subtitle = plot.subtitle,
    ...
  )
}
