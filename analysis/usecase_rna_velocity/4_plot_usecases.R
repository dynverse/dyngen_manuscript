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

metric_labels <- c(cor = "Velocity correlation", mean_cosine = "Vector cosine similarity")
plota_data <-
  scores$summ %>%
  mutate(
    contains_cycle = ifelse(dataset_trajectory_type %in% c("graph", "disconnected_graph", "cycle"), "Yes", "No")
  ) %>%
  gather(metric, score, cor, mean_cosine) %>%
  mutate(
    metric_label = factor(metric_labels[metric], levels = metric_labels)
  )

plot_part_A <-
  ggplot(plota_data, aes(paste0(method_id, "\n", params_id), score)) +
  # ggbeeswarm::geom_quasirandom(aes(color = dataset_trajectory_type, shape = difficulty), size = 3)  +
  ggbeeswarm::geom_quasirandom(aes(color = difficulty), size = 2)  +
  theme_bw() +
  theme_common() +
  labs(x = "", y = "Correlation (higher is better)", colour = "Difficulty") +
  scale_y_continuous(limits = c(0, 1)) +
  guides(colour = guide_legend(nrow = 3), shape = guide_legend(nrow = 3)) +
  facet_wrap(~metric_label, ncol = 2)

plot_part_A

# PART B: Illustration of velocity ----------------------------------------
# dataset_id <- "bifurcating_hard_seed1"
# dataset_id <- "consecutive_bifurcating_hard_seed1"
# dataset_id <- "disconnected_hard_seed1"
# dataset_id <- "trifurcating_hard_seed1"
# dataset_id <- "cycle_hard_seed1"
# dataset_id <- "bifurcating_cycle_medium_seed1"
# dataset_id <- "linear_simple_easy_seed1"
# dataset_id <- "bifurcating_cycle_easy_seed1"
# dataset_id <- "bifurcating_cycle_medium_seed1"
dataset_id <- "bifurcating_cycle_hard_seed3"
# dataset_id <- "bifurcating_loop_hard_seed1"
# dataset_id <- "bifurcating_converging_hard_seed1"

dataset <- read_rds(exp$dataset_file(dataset_id))
model <- read_rds(exp$model_file(dataset_id))

# mil_cols <- dynplot2::define_milestone_colors(milestone_colors = NULL, milestone_ids = dataset$milestone_ids)
# mil_pct_mat <- dataset$milestone_percentages %>% reshape2::acast(cell_id~milestone_id, value.var = "percentage", fill = 0) %>% .[rownames(dataset$dimred), dataset$milestone_ids]
# mil_rgb <- mil_pct_mat %*% mil_cols
# cols <- rgb(mil_rgb[,1], mil_rgb[,2], mil_rgb[,3], maxColorValue = 256)
# rgl::plot3d(dataset$dimred, col = cols, size = 10)
# userMatrix <- rgl::par3d()$userMatrix[1:3,1:3]
# userMatrix %>% as.vector %>% paste0(collapse = ", ")

if (dataset_id == "bifurcating_cycle_hard_seed1") {
  userVec <- c(0.0203961916267872, -0.740907371044159, -0.671227753162384, -0.0408603921532631, -0.671426177024841, 0.739882826805115, -0.99891185760498, 0.012335604056716, -0.0439676493406296)
} else if (dataset_id == "bifurcating_cycle_hard_seed3") {
  userVec <- c(-0.700981974601746, -0.500725984573364, -0.507759869098663, -0.630953311920166, 0.767287373542786, 0.114417120814323, 0.332326680421829, 0.400592744350433, -0.853821098804474)
}
userMatrix <- matrix(userVec, nrow = 3)

dimred <- dataset$dimred %*% t(userMatrix)
colnames(dimred) <- colnames(dataset$dimred)
dataset <- dataset %>% dynwrap::add_dimred(dimred = dimred)

# dynplot_dimred(dataset2) +
#   geom_cell_point(aes(color = milestone_percentages), size = 1) +
#   scale_milestones_colour() +
#   geom_trajectory_segments(size = 1, color = "#333333") +
#   geom_milestone_label(aes(label = labels[label]), color = "black", fill = "#EEEEEE") +
#   theme_common(legend.position = "none") +
#   ggtitle("Trajectory")

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
# feature_oi <- "C1_TF1"
feature_oi <- "D1_TF1"
  # model$feature_network %>%
  # group_by(to) %>%
  # summarise(both = any(effect == -1) && any(effect == 1)) %>%
  # filter(both) %>%
  # pull(to) %>%
  # first()

# feature_oi <-
#   model$feature_network %>%
#   group_by(to) %>%
#   summarise(both = any(effect == -1) && any(effect == 1)) %>%
#   filter(both) %>%
#   pull(to) %>%
#   first()

expression_plot <- dynplot_dimred(dataset) +
  geom_cell_point(aes(color = select_feature_expression(feature_oi, .data)), size = 1) +
  geom_trajectory_segments(size = 1, color = "#333333") +
  scale_expression_color(breaks = c(0, 1), labels = c("min", "max")) +
  theme_common() +
  ggtitle(paste0("Expression of gene ", feature_oi))

# Plot 3, ground truth velocity
transform_groundtruth_velocity <- function(x) {
  scales::squish(x, c(-1, 1), only.finite = FALSE)
}

gs_plot <- dynplot_dimred(dataset) +
  geom_cell_point(aes(color = transform_groundtruth_velocity(dataset$rna_velocity[,feature_oi])), size = 1) +
  dynplot2:::scale_velocity_color() +
  ggtitle(paste0("Ground truth velocity of ", feature_oi)) +
  theme_common()


# Combine plots
plot_part_B <- patchwork::wrap_plots(plot_trajectory, expression_plot, gs_plot, nrow = 1)




# PART C: RNA velocity estimates of different methods ---------------------
design_velocity_oi <- design_velocity %>% filter(dataset_id == !!dataset_id)

plot_part_C <- pmap(design_velocity_oi, function(dataset_id, method_id, params_id, ...) {
  velocity <- read_rds(exp$velocity_file(dataset_id, method_id, params_id))
  velo_vec <- velocity$velocity_vector[,feature_oi]
  maxv <- quantile(abs(velo_vec), .8)
  dynplot_dimred(dataset) +
    geom_cell_point(aes(color = transform_groundtruth_velocity(velo_vec / maxv)), size = 1) +
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
      stat = stat_velocity_grid(grid_bandwidth = 1, filter = rlang::quo(mass > max(mass) * 0.05)),
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
    geom_segment(aes(x = from_comp_1, xend = to_comp_1, y = from_comp_2, yend = to_comp_2, colour = scales::squish(simil, c(.9, 1))), wps, size = 2) +
    ggtitle(method_id, subtitle = params_id) +
    theme_common() +
    viridis::scale_color_viridis(name = "Cosine", limits = c(.9, 1))
    # scale_color_distiller(palette = "RdBu", limits = c(-1, 1), name = "Cosine")

  if (i != 1) {
    pl <- pl + theme(legend.position = "none")
  }

  pl
}) %>% patchwork::wrap_plots(nrow = 1)

wps <- pmap_df(design_velocity_oi %>% mutate(i = row_number()), function(i, dataset_id, method_id, params_id, ...) {
  velocity <- read_rds(exp$velocity_file(dataset_id, method_id, params_id))
  wps <- scores$per_waypoint
})

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
  plot_part_B,
  plot_part_C,
  (plot_part_D & theme(plot.title = element_blank(), plot.subtitle = element_blank())),
  plot_part_A,
  # (plot_part_E & theme(plot.title = element_blank(), plot.subtitle = element_blank())),
  ncol = 1,
  heights = c(1, .6, .6, .8)
)
ggsave(exp$result("usecase.pdf"), g, height = 12, width = 12, useDingbats = FALSE)
ggsave(exp$result("usecase.png"), g, height = 12, width = 12)
