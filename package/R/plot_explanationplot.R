#' @importFrom reshape2 melt
#' @export
plot_density <- function(a1, title = "Distance matrix +\n warping path", subtitle = NULL){
  linedf <- tibble(
    x = a1$index1,
    y = a1$index2
  )
  rasterdf <- reshape2::melt(a1$costMatrix, varnames = c("x", "y"), value.name = "distance")
  ggplot() +
    geom_raster(aes(x, y, fill = distance), rasterdf) +
    geom_line(aes(x, y), linedf) +
    scale_fill_distiller(palette = "RdYlGn", breaks = range(rasterdf$distance), labels = c("min", "max")) +
    scale_x_continuous(expand = c(0, 0), breaks = range(linedf$x), labels = c("start", "end")) +
    scale_y_continuous(expand = c(0, 0), breaks = range(linedf$y), labels = c("start", "end")) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    ) +
    labs(title = title, x = "Sample 2", y = "Sample 1", fill = "Accumulated\ndistance")
}

#' @export
plot_pseudotime <- function(ds, palette = "Blues", ...){
  plot_dimred(ds, color_cells = "pseudotime", size_cells = 3, ...) +
    scale_color_distiller(
      palette = palette,
      breaks = function(x) range(x),
      labels = c("begin", "end"),
      # guide = guide_colorbar(title.vjust = 0.75, title.hjust = 0.75, barwidth = 5, draw.ulim = TRUE)
    ) +
    labs(colour = "Pseudotime")
}

#' @export
plot_trajectory_in_color <- function(ds, milestone_colors, labels, colors, plot_sd, ...){
  milestone_id <- ds$milestone_ids
  cols <- as_tibble(milestone_id)
  cols$color <- milestone_colors
  cols <- rename(cols, milestone_id=value)

  p_traj1 <- plot_dimred(ds, color_milestones = "given", milestones = cols, alpha_cells = 0.8, size_cells = 3, trajectory_projection_sd = plot_sd, ...)
  p_traj1 <- p_traj1 + labs(color = "") +
    scale_color_manual(labels = labels, values = colors) +
    theme(legend.position="bottom")
  p_traj1
}
