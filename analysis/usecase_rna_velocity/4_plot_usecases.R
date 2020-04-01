library(tidyverse)
library(dyngen.manuscript)
library(dynplot2)

exp <- start_analysis("usecase_rna_velocity")

# reticulate::use_python("/usr/bin/python3", required = TRUE)
design_velocity <- read_rds(exp$result("design_velocity.rds"))

# PART A: Average scores --------------------------------------------------
mean_scores <- read_rds(exp$result("scores_aggregated.rds"))

plot_part_A <-
  ggplot(mean_scores, aes(paste0(method_id, "\n", params_id), score)) +
  ggbeeswarm::geom_quasirandom(aes(color = backbone))  +
  theme_bw() +
  scale_x_discrete("") +
  scale_y_continuous("Correlation (higher is better)") +
  theme_common()


# PART B: Illustration of velocity ----------------------------------------
dataset_id <- "bifurcating_2"

dataset <- read_rds(exp$dataset_file(dataset_id))
model <- read_rds(exp$model_file(dataset_id))

# compute dimred if dimred is missing
if (is.null(dataset$dimred)) {
  set.seed(1)
  dataset <- dataset %>% dynwrap::add_dimred(dyndimred::dimred_landmark_mds)
}

# Plot 1, trajectory
labels <- setNames(LETTERS[seq_along(dataset$milestone_ids)], dataset$milestone_ids)
plot_trajectory <- dynplot_dimred(dataset) +
  geom_cell_point(aes(color = milestone_percentages)) +
  scale_milestones_colour() +
  geom_trajectory_segments(size = 2, color = "#333333") +
  geom_milestone_label(aes(fill = milestone_id, label = labels[label]), color = "black") +
  theme_common(legend.position = "none") +
  ggtitle("Trajectory")

# Plot 2, expression of a gene that goes up and down
feature_oi <-
  model$feature_network %>%
  group_by(to) %>%
  summarise(both = any(effect == -1) && any(effect == 1)) %>%
  filter(both) %>%
  pull(to) %>%
  first()

expression_plot <- dynplot_dimred(dataset) +
  geom_cell_point(aes(color = select_feature_expression(feature_oi, .data))) +
  geom_trajectory_segments() +
  scale_expression_color(breaks = c(0, 1), labels = c("min", "max")) +
  theme_common() +
  ggtitle("Expression of a gene that goes up and down")

# Plot 3, ground truth velocity
transform_groundtruth_velocity <- function(x) {
  log2(scales::squish(x, c(0.1, 10), only.finite = FALSE))
}

gs_plot <- dynplot_dimred(dataset) +
  geom_cell_point(aes(color = transform_groundtruth_velocity(dataset$propensity_ratios[,feature_oi]))) +
  dynplot2:::scale_velocity_color() +
  ggtitle("Ground truth velocity") +
  theme_common()

# Combine plots
plot_part_B <- patchwork::wrap_plots(plot_trajectory, expression_plot, gs_plot, nrow = 1)




# PART C: RNA velocity estimates of different methods ---------------------
design_velocity_oi <- design_velocity %>% filter(dataset_id == !!dataset_id)

plot_part_C <- pmap(design_velocity_oi, function(dataset_id, method_id, params_id, ...) {
  velocity <- read_rds(exp$velocity_file(dataset_id, method_id, params_id))
  dataset2 <- dataset %>% scvelo::add_velocity(velocity = velocity)
  dynplot_dimred(dataset2) +
    geom_cell_point(aes(color = select_feature_velocity(feature_oi, .data))) +
    dynplot2:::scale_velocity_color(name = "", guide = "none") +
    ggtitle(method_id, subtitle = params_id) +
    theme_common()
}) %>% patchwork::wrap_plots(nrow = 1)


#' @examples
#' design_velocity_oi %>% dynutils::extract_row_to_list(1) %>% list2env(.GlobalEnv)

# PART D: Embedded RNA velocity estimates of different methods ------------
velocity_plots <- pmap(design_velocity_oi, function(dataset_id, method_id, params_id, ...) {
  velocity_file <- exp$velocity_file(dataset_id, method_id, params_id)
  velocity <- read_rds(velocity_file)
  if (method_id == "scvelo") {
    velocity$scvelo <- reticulate::py_load_object(paste0(dirname(velocity_file), "/scvelo.pkl"))
  }

  dataset2 <- dataset %>%
    scvelo::add_velocity(velocity = velocity)
  dataset2 <- dataset2 %>%
    scvelo::add_dimred_future()

  dynplot_dimred(dataset2) +
    geom_cell_point(aes(color = milestone_percentages)) +
    scale_milestones_colour() +
    geom_velocity_arrow(
      size = 1.2,
      color = "#333333",
      stat = stat_velocity_grid(grid_bandwidth = 1),
      arrow = arrow(length = unit(0.2, "cm"))
    ) +
    # geom_velocity_arrow(size = 1, color = "white") +
    ggtitle(method_id, subtitle = params_id) +
    theme_common() +
    theme(legend.position = "none")
})

plot_part_D <- patchwork::wrap_plots(velocity_plots, nrow = 1)



# COMBINE ALL PARTS -------------------------------------------------------
tag_first <- function(x, tag) {
  y <- x[[1]]

  if ("patchwork" %in% class(y)) {
    x[[1]] <- tag_first(x[[1]], tag = tag)
  } else {
    x[[1]] <- x[[1]] + labs(tag = tag)
  }
  x
}

plot_part_A <- plot_part_A + labs(tag = "A")
plot_part_B[[1]] <- plot_part_B[[1]] + labs(tag = "B")
plot_part_C[[1]] <- plot_part_C[[1]] + labs(tag = "C")
plot_part_D[[1]] <- plot_part_D[[1]] + labs(tag = "D")

pdf("results/velocity/usecase.pdf", height = 14, width = 10, useDingbats = FALSE)
patchwork::wrap_plots(
  plot_part_A,
  plot_part_B,
  plot_part_C,
  (plot_part_D & theme(plot.title = element_blank(), plot.subtitle = element_blank())),
  ncol = 1,
  heights = rep(1, 4)
)
graphics.off()

magick::image_read_pdf("results/velocity/usecase.pdf") %>%
  magick::image_convert("png") %>%
  magick::image_write("results/velocity/usecase.png") %>%
  invisible()

