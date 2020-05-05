library(tidyverse)
library(dyngen.manuscript)
library(dynplot2)
library(ggbeeswarm)

exp <- start_analysis("usecase_rna_velocity")

design_velocity <- read_rds(exp$result("design_velocity.rds"))

# PART A: Average scores --------------------------------------------------
scores <- read_rds(exp$result("scores.rds"))

ggplot(scores$summ %>% mutate(group = paste0(method_id, " ", params_id)), aes(cor, mean_cosine)) + scale_colour_brewer(palette = "Dark2") +
  geom_point(aes(colour = group))

ggplot(scores$summ %>% mutate(group = paste0(method_id, " ", params_id)), aes(cor, mean_cosine)) + scale_colour_brewer(palette = "Dark2") +
  geom_point(colour = "lightgray", data = function(x) x %>% select(-group)) +
  geom_point(aes(colour = group)) +
  facet_wrap(~group) +
  theme_bw()

# plot_part_A <-
#   ggplot(scores$summ, aes(paste0(method_id, "\n", params_id), cor)) +
#   ggbeeswarm::geom_quasirandom(aes(color = dataset_trajectory_type, shape = difficulty), size = 3)  +
#   theme_bw() +
#   theme_common() +
#   labs(x = "", y = "Correlation (higher is better)", colour = "Trajectory type", shape = "Difficulty") +
#   scale_colour_manual(values = dynwrap::trajectory_types %>% select(id, colour) %>% deframe) +
#   scale_y_continuous(limits = c(0, 1)) +
#   guides(colour = guide_legend(nrow = 3), shape = guide_legend(nrow = 3))

plot_part_A <-
  ggplot(scores$summ %>% gather(metric, score, cor, mean_cosine), aes(paste0(method_id, "\n", params_id), score)) +
  ggbeeswarm::geom_quasirandom(aes(color = dataset_trajectory_type, shape = difficulty), size = 3)  +
  theme_bw() +
  theme_common() +
  labs(x = "", y = "Correlation (higher is better)", colour = "Trajectory type", shape = "Difficulty") +
  scale_colour_manual(values = dynwrap::trajectory_types %>% select(id, colour) %>% deframe) +
  scale_y_continuous(limits = c(0, 1)) +
  guides(colour = guide_legend(nrow = 3), shape = guide_legend(nrow = 3)) +
  facet_wrap(~metric, ncol = 2)


plot_part_A

# PART B: Illustration of velocity ----------------------------------------
# dataset_id <- "bifurcating_hard_seed1_2500"
# dataset_id <- "trifurcating_hard_seed1_2500"
# dataset_id <- "cycle_hard_seed1_2500"
# dataset_id <- "bifurcating_cycle_medium_seed1_2500"
# dataset_id <- "bifurcating_loop_hard_seed1_2500"
dataset_id <- "bifurcating_converging_hard_seed1_2500"

dataset <- read_rds(exp$dataset_file(dataset_id))
model <- read_rds(exp$model_file(dataset_id))

# compute dimred if dimred is missing
if (is.null(dataset$dimred)) {
  set.seed(1)
  dataset <- dataset %>% dynwrap::add_dimred(dyndimred::dimred_mds)
}

# Plot 1, trajectory
labels <- setNames(LETTERS[seq_along(dataset$milestone_ids)], dataset$milestone_ids)
plot_trajectory <- dynplot_dimred(dataset) +
  geom_cell_point(aes(color = milestone_percentages), size = 1) +
  scale_milestones_colour() +
  geom_trajectory_segments(size = 1, color = "#333333") +
  geom_milestone_label(aes(label = labels[label]), color = "black", fill = "#EEEEEE") +
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
  geom_cell_point(aes(color = select_feature_expression(feature_oi, .data)), size = 1) +
  geom_trajectory_segments(size = 1, color = "#333333") +
  scale_expression_color(breaks = c(0, 1), labels = c("min", "max")) +
  theme_common() +
  ggtitle("Expression of a gene that goes up and down")

# Plot 3, ground truth velocity
transform_groundtruth_velocity <- function(x) {
  scales::squish(x, c(-1, 1), only.finite = FALSE)
}

gs_plot <- dynplot_dimred(dataset) +
  geom_cell_point(aes(color = transform_groundtruth_velocity(dataset$rna_velocity[,feature_oi])), size = 1) +
  dynplot2:::scale_velocity_color() +
  ggtitle("Ground truth velocity") +
  theme_common()


# Combine plots
plot_part_B <- patchwork::wrap_plots(plot_trajectory, expression_plot, gs_plot, nrow = 1)




# PART C: RNA velocity estimates of different methods ---------------------
design_velocity_oi <- design_velocity %>% filter(dataset_id == !!dataset_id)

plot_part_C <- pmap(design_velocity_oi, function(dataset_id, method_id, params_id, ...) {
  velocity <- read_rds(exp$velocity_file(dataset_id, method_id, params_id))
  dynplot_dimred(dataset) +
    geom_cell_point(aes(color = transform_groundtruth_velocity(velocity$velocity_vector[,feature_oi])), size = 1) +
    dynplot2:::scale_velocity_color(name = "", guide = "none") +
    ggtitle(method_id, subtitle = params_id) +
    theme_common()
}) %>% patchwork::wrap_plots(nrow = 1)


# PART D: Embedded RNA velocity estimates of different methods ------------
plot_part_D <- pmap(design_velocity_oi, function(dataset_id, method_id, params_id, ...) {
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
    geom_cell_point(aes(color = milestone_percentages), size = 1) +
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
}) %>%
  patchwork::wrap_plots(nrow = 1)



# PART E: RNA velocity estimates of different methods ---------------------
#' design_velocity_oi %>% dynutils::extract_row_to_list(1) %>% list2env(.GlobalEnv)
#' i <- 1
plot_part_E <- pmap(design_velocity_oi %>% mutate(i = row_number()), function(i, dataset_id, method_id, params_id, ...) {
  velocity <- read_rds(exp$velocity_file(dataset_id, method_id, params_id))
  wps <- scores$per_waypoint %>% filter(dataset_id == !!dataset_id, method_id == !!method_id, params_id == !!params_id)
  pl <- dynplot_dimred(dataset) +
    geom_cell_point(size = 1) +
    geom_segment(aes(x = from_comp_1, xend = to_comp_1, y = from_comp_2, yend = to_comp_2, colour = simil), wps, size = 2) +
    ggtitle(method_id, subtitle = params_id) +
    theme_common() +
    scale_color_distiller(palette = "RdBu", limits = c(-1, 1), name = "Cosine")

  if (i != 1) {
    pl <- pl + theme(legend.position = "none")
  }

  pl
}) %>% patchwork::wrap_plots(nrow = 1)


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
plot_part_E[[1]] <- plot_part_E[[1]] + labs(tag = "E")

g <- patchwork::wrap_plots(
  plot_part_A,
  plot_part_B,
  plot_part_C,
  (plot_part_D & theme(plot.title = element_blank(), plot.subtitle = element_blank())),
  (plot_part_E & theme(plot.title = element_blank(), plot.subtitle = element_blank())),
  ncol = 1,
  heights = c(.8, 1, .6, .6, .6)
)
ggsave(exp$result("usecase.pdf"), g, height = 14, width = 12, useDingbats = FALSE)
ggsave(exp$result("usecase.png"), g, height = 14, width = 10)
