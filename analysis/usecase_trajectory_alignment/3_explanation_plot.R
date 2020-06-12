library(tidyverse)
library(dyngen.manuscript)
library(patchwork)
library(tidyr)
library(reshape2)
library(dyngen)
library(dynplot)
library(gtools)
library(dtw)

set.seed(1)
exp <- start_analysis("usecase_trajectory_alignment")

backbone <- bblego(
  bblego_start("A", type = "simple", num_modules = 4),
  bblego_linear("A", "B", type = "simple", num_modules = 6),
  bblego_linear("B", "C", type = "simple", num_modules = 6),
  bblego_linear("C", "D", type = "simple", num_modules = 6),
  bblego_end("D", type = "simple", num_modules = 4)
)

healthy_model <- initialise_model(
  num_tfs = nrow(backbone$module_info),
  num_targets = 0,
  num_hks = 0,
  num_cells = 1000,
  backbone = backbone,
  simulation_params = simulation_default(census_interval = 10),
  verbose = TRUE,
  num_cores = 8,
  download_cache_dir = "~/.cache/dyngen/"
) %>%
  generate_tf_network() %>%
  generate_feature_network() %>%
  generate_kinetics() %>%
  generate_gold_standard()

model1_control <-
  healthy_model %>% normalise_goldstandard() %>% generate_cells()

model2_diseased <-
  healthy_model %>%
  generate_diff_abundance(
    froms = c("sC", "sD"),
    tos = c("sD", "sEndD"),
    ratios = c(0, 0)
  ) %>% generate_cells()

combined_model <-
  combine_models(model1_control, model2_diseased) %>%
  generate_experiment()

datasetD <- wrap_dataset(combined_model)
plot_dimred(datasetD)

dataset1 <- model1_control %>% generate_experiment() %>% wrap_dataset()
dataset2D <- model2_diseased %>% generate_experiment() %>% wrap_dataset()

## Save all 3 datasets
write_rds(datasetD, paste0(exp$dataset_folder("expplot"), "datasetD.rds"), compress = "gz")
write_rds(dataset1, paste0(exp$dataset_folder("expplot"), "dataset1.rds"), compress = "gz")
write_rds(dataset2D, paste0(exp$dataset_folder("expplot"), "dataset2D.rds"), compress = "gz")

milestone_colors <- c("#ff681c", "#ff681c", "#ff681c", "#ff681c", "#3abbba","#3abbba", "#3abbba", "#3abbba", "#ff681c", "#3abbba")
together <- plot_trajectory_in_color(datasetD, milestone_colors, c("sample 1", "sample 2"), c("#3abbba", "#ff681c"), plot_sd = 0.11)

t1_pt <- plot_pseudotime(dataset1, palette = "Oranges", trajectory_projection_sd = 0.11, plot_trajectory = T)
t2_pt <- plot_pseudotime(dataset2D, palette = "Blues", trajectory_projection_sd = 0.11, plot_trajectory = T)

expr1 <- get_cell_expression(dataset1, dataset1$milestone_network, "sA")$expression
expr2 <- get_cell_expression(dataset2D, dataset2D$milestone_network, "sA")$expression
a <- dtw(expr2, expr1, keep.internals = T, open.end = F)
dtwPlotAlignment(a)

plot_dens <- plot_density(a, show_legend = TRUE) + labs(title = NULL, subtitle = NULL)
plot_dens


e1 <- calculate_correct_pseudotime(dataset1, dataset1$milestone_network, "sA", normalized = F)
e2 <- calculate_correct_pseudotime(dataset2D, dataset2D$milestone_network, "sA", normalized = F)

# warped segmenten
traj1 <- data.frame(x = sort(e1))
traj2 <- data.frame(x = sort(e2))
traj1$color <- 1
traj2$color <- 2
traj1$y <- 0
traj2$y <- 1

traj1_warped <- traj1[a$index2,]
traj2_warped <- traj2[a$index1,]

combinatie <- traj1_warped
combinatie$x2 <- traj2_warped$x
combinatie$y2 <- traj2_warped$y
leftover_combinatie <- combinatie[seq(1, nrow(combinatie), 20), ]

alles <- data.frame(x1 = c(leftover_combinatie$x, leftover_combinatie$x2))
alles$y <- c(leftover_combinatie$y, leftover_combinatie$y2)
alles$color <- c(rep(1, nrow(leftover_combinatie)), rep(2, nrow(leftover_combinatie)))

cell_mappings <- ggplot() +
  geom_segment(data = leftover_combinatie, mapping = aes(x = x, y = y, xend = x2, yend = y2), alpha = 1) +
  geom_point(data = alles, mapping = aes(x = x1, y = y, colour = as.factor(color)), size = 3.5, show.legend = F) +
  scale_colour_manual(values = c("#fd8d3c", "#6baed6"), label = "") +
  scale_x_continuous(breaks = c(0, 1), labels = c("Start", "End")) +
  theme_void() +
  theme(
    axis.title = element_text(),
    axis.title.y = element_blank(),
    axis.line = element_line(),
    axis.line.y = element_blank(),
    axis.ticks = element_line(),
    axis.ticks.y = element_blank(),
    axis.text = element_text(),
    axis.text.y = element_blank()
  ) +
  labs(x = "Simulation time")

cell_mappings

write_rds(lst(leftover_combinatie, alles, combinatie), exp$result("cell_mappings_data.rds"), compress = "gz")


# For each 50th cell -> 26
# For each 30th cell -> 41
# For each 20th cell -> 63
demarcation_line <- 67

straight_part <- leftover_combinatie[1:demarcation_line,]
straight_part$x2 <- straight_part$x
partial_part <- leftover_combinatie[demarcation_line+1:nrow(leftover_combinatie),]
partial_part$x2 <- partial_part$x
partial_part$y2 <- 0.5
ground_truth <- rbind(straight_part, partial_part)

ground_truth_mappings <- ggplot() +
  geom_segment(data = ground_truth, mapping = aes(x = x, y = y, xend = x2, yend = y2), alpha = 1) +
  geom_point(data = alles, mapping = aes(x = x1, y = y, colour = as.factor(color)), size = 3.5, show.legend = F) +
  scale_colour_manual(values = c("#fd8d3c", "#6baed6"), label = "") +
  scale_x_continuous(breaks = c(0, 1), labels = c("Start", "End")) +
  theme_void() +
  theme(
    axis.title = element_text(),
    axis.title.y = element_blank(),
    axis.line = element_line(),
    axis.line.y = element_blank(),
    axis.ticks = element_line(),
    axis.ticks.y = element_blank(),
    axis.text = element_text(),
    axis.text.y = element_blank()
  ) +
  labs(x = "Simulation time")

ground_truth_mappings

write_rds(lst(ground_truth, alles, combinatie), exp$result("ground_truth_mappings_data.rds"), compress = "gz")

part1 <-
  patchwork::wrap_plots(
    t2_pt,
    t1_pt,
    plot_spacer(),
    cell_mappings,
    plot_dens,
    widths = c(1, 1, .2, 1, 1)
  )

saveRDS(part1, file = exp$result("explanation_flat.rds"))
ggsave(part1, filename = exp$result("explanation_flat.png"), bg = 'transparent', width = 20, height = 4)
