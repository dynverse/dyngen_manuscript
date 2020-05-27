#' @export
plot_density <- function(a1, title = "Distance matrix +\n warping path", subtitle = ""){
  x <- as_tibble(a1$index1)
  x$Y <- a1$index2

  gga_1 <- melt(a1$costMatrix)
  p_heat1 <- ggplot(gga_1, aes(Var1, Var2, fill= value)) +
    geom_raster(show.legend = FALSE) +
    scale_fill_distiller(palette = "RdYlGn") +
    scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
    geom_line(data=x, aes(x=value, y=Y)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    labs(title = title, subtitle = "Accumulated\ndistance", x = "Sample 2", y = "Sample 1")

  p_heat1
}

#' @export
plot_pseudotime <- function(ds, low = "#1b2944", high = "#3abbba", ...){
  p_traj1 <- plot_dimred(ds, color_cells = "pseudotime", size_cells = 3, ...)
  p_traj1 <- p_traj1 +
    scale_color_gradient(low = low, high = high, limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1"), guide = guide_colorbar(title = "Pseudotime", title.vjust = 0.75, title.hjust = 0.75, barwidth = 5, draw.ulim = T))
  #+
    # theme(legend.position="bottom") +
    # theme(legend.key.width = unit(2, "cm"))
  p_traj1
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
