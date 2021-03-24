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

# read dataset
design_dataset <- read_rds(exp$result("design_datasets.rds"))
file <- exp$dataset_file(design_dataset$id[[1]])
dataset <- read_rds(file) %>% dynwrap::add_dimred(dimred = dyndimred::dimred_landmark_mds)

# split up in two
dataset1 <- dataset %>% filter_cells(model == "left", keep_dimred = TRUE)
dataset2 <- dataset %>% filter_cells(model == "right", keep_dimred = TRUE)

# add healthy/diseased grouping
dataset$cell_info <- dataset$cell_info %>% mutate(
  group = c("left" = "Healthy", "right" = "Diseased")[model]
)
group_palette <- c(Diseased = "#fd8d3c", Healthy = "#6baed6")

# construct gold standard alignment
step_progressions <- crossing(
  dataset$milestone_network %>% select(from, to),
  percentage = seq(0, 1, by = .2)
) %>%
  mutate(cell_id = paste0("step", row_number()))
step_waypoints <- list(
  milestone_percentages = dynwrap::convert_progressions_to_milestone_percentages(
    cell_id = step_progressions$waypoint_id,
    milestone_ids = dataset$milestone_ids,
    milestone_network = dataset$milestone_network,
    progressions = step_progressions
  ) %>% rename(waypoint_id = cell_id)
)
step_proj <- dynwrap::project_waypoints(trajectory = dataset, space = dataset$dimred, waypoints = step_waypoints)

ds <- data.frame(step_progressions, step_proj[step_progressions$cell_id, ]) %>%
  as_tibble() %>%
  mutate(from_ = gsub("^(left|right)_", "", from), to_ = gsub("^(left|right)_", "", to))

# remove large matrices
dataset$counts <-
  dataset$counts_spliced <-
  dataset$counts_unspliced <-
  dataset$counts_protein <-
  dataset$expression <-
  NULL

# create plot
library(dynplot2)
ground_truth <- dynplot_dimred(dataset) +
  geom_path(aes(comp_1, comp_2, group = paste(from_, to_, percentage, sep = "~")), linetype = "dashed", colour = "darkgray", ds) +
  geom_cell_point(colour = "black", size = 3) +
  geom_cell_point(aes(colour = group), size = 2.5) +
  geom_trajectory_segments(size = 2, arrow_size = .2) +
  geom_milestone_point(size = 4, colour = "black") +
  labs(colour = "Sample") +
  scale_colour_manual(values = group_palette) +
  theme_common()
ground_truth

# run method
dtw_out <- ta_methods$DTW(dataset1, dataset2)
res1 <- get_cell_expression(dataset1)
res2 <- get_cell_expression(dataset2)

df <- tibble(
  i = dtw_out$index1,
  j = dtw_out$index2,
  ci = names(res1$pseudotime)[i],
  cj = names(res2$pseudotime)[j],
  pi = dtw_out$pt1_aligned,
  pj = dtw_out$pt2_aligned,
  rn = seq_along(i)
) %>%
  filter(rn %% round(n()/15) == 1)
df2 <- df %>%
  select(pos = rn, left = ci, right = cj) %>%
  gather(group, cell_id, left, right) %>%
  mutate(pred_id = paste0("pred_", cell_id))

pred_waypoints <- list(
  milestone_percentages =
    dataset$milestone_percentages %>%
    filter(cell_id %in% df2$cell_id) %>%
    mutate(cell_id = paste0("pred_", cell_id)) %>%
    rename(waypoint_id = cell_id)
)
pred_proj <- dynwrap::project_waypoints(trajectory = dataset, space = dataset$dimred, waypoints = pred_waypoints)

pred_df <- data.frame(df2, pred_proj[df2$pred_id, ]) %>%
  as_tibble()

prediction <- dynplot_dimred(dataset) +
  geom_path(aes(comp_1, comp_2, group = pos), linetype = "dashed", colour = "darkgray", pred_df) +
  geom_cell_point(colour = "black", size = 3) +
  geom_cell_point(aes(colour = group), size = 2.5) +
  geom_trajectory_segments(size = 2, arrow_size = .2) +
  geom_milestone_point(size = 4, colour = "black") +
  labs(colour = "Sample") +
  scale_colour_manual(values = group_palette) +
  theme_common()
prediction

prediction_small <- dynplot_dimred(dataset) +
  geom_path(aes(comp_1, comp_2, group = pos), linetype = "dashed", colour = "darkgray", pred_df) +
  geom_cell_point(colour = "black", size = 3) +
  geom_cell_point(aes(colour = group), size = 2.5) +
  geom_trajectory_segments(size = 1, arrow_size = .2, grid::arrow(type = "closed", length = unit(0.06, "inches"))) +
  geom_milestone_point(size = 2, colour = "black") +
  labs(colour = "Sample") +
  scale_colour_manual(values = group_palette) +
  theme_common()


# save plots
plots <- lst(
  ground_truth,
  prediction,
  prediction_small
)

write_rds(plots, exp$result("usecase_separateplots.rds"), compress = "xz")

