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
  verbose = TRUE
) %>%
  generate_tf_network() %>%
  generate_feature_network() %>%
  generate_kinetics() %>%
  generate_gold_standard()

model1_control <-
  healthy_model %>%
  generate_cells()

model2_diseased <-
  healthy_model %>%
  generate_diff_abundance(
    froms = c("sC", "sD"),
    tos = c("sD", "sEndD"),
    ratios = c(0, 0)
  ) %>% generate_cells()

combined_model <-
  combine_models(list("control" = model1_control, "diseased" = model2_diseased)) %>%
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

t1_pt <- plot_pseudotime(dataset1, palette = "Oranges", trajectory_projection_sd = 0.11, plot_trajectory = TRUE)
t2_pt <- plot_pseudotime(dataset2D, palette = "Blues", trajectory_projection_sd = 0.11, plot_trajectory = TRUE)

expr1 <- get_cell_expression(dataset1)$expression
expr2 <- get_cell_expression(dataset2D)$expression
a <- dtw(expr2, expr1, keep.internals = TRUE, open.end = FALSE)
dtwPlotAlignment(a)

plot_dens <- plot_density(a, title = "Alignment between\ndiseased/healthy")

e1 <- compute_pseudotime_from_root(dataset1, normalized = FALSE)
e2 <- compute_pseudotime_from_root(dataset2D, normalized = FALSE)
maxe12 <- max(c(e1, e2))
e1 <- e1 / maxe12
e2 <- e2 / maxe12

# warped segmenten
traj1 <- tibble(x = sort(e1), y = 0, traji = 1, traj = "Healthy")
traj2 <- tibble(x = sort(e2), y = 1, traji = 2, traj = "Diseased")

ix <- seq(1, length(a$index1), by = 20)
traj1_warped <- traj1[a$index2[ix],]
traj2_warped <- traj2[a$index1[ix],]

alles <- bind_rows(traj1, traj2)

leftover_combinatie <- bind_cols(
  traj1_warped %>% select(x, y),
  traj2_warped %>% select(x2 = x, y2 = y)
)

cell_mappings <-
  ggplot() +
  geom_segment(data = leftover_combinatie, mapping = aes(x = x, y = y, xend = x2, yend = y2), alpha = 1) +
  geom_point(data = alles, mapping = aes(x = x, y = y, colour = traj), size = 3.5) +
  scale_colour_manual(values = c("#fd8d3c", "#6baed6")) +
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
  labs(x = "Simulation time", colour = "Sample")

cell_mappings

# For each 50th cell -> 26
# For each 30th cell -> 41
# For each 20th cell -> 63
demarcation_line <- 67

straight_part <-
  leftover_combinatie[seq_len(demarcation_line),] %>%
  mutate(x2 = x)
partial_part <-
  leftover_combinatie[-seq_len(demarcation_line),] %>%
  mutate(x2 = x, y2 = .5)

ground_truth <- rbind(straight_part, partial_part)

ground_truth_mappings <-
  ggplot() +
  geom_segment(aes(x = x, y = y, xend = x2, yend = y2), ground_truth, alpha = 1) +
  geom_point(aes(x = x, y = y, colour = traj), alles, size = 3.5) +
  scale_colour_manual(values = c("#fd8d3c", "#6baed6")) +
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
  labs(x = "Simulation time", colour = "Sample")

ground_truth_mappings

write_rds(
  lst(t1_pt, t2_pt, ground_truth = ground_truth_mappings, prediction = cell_mappings, plot_dens),
  exp$result("usecase_separateplots.rds"),
  compress = "gz"
)

