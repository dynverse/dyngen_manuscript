library(tidyverse)
library(dyngen.manuscript)

exp <- start_analysis("usecase_rna_velocity")

design_datasets <- read_rds(exp$result("design_datasets.rds"))
design_velocity <- read_rds(exp$result("design_velocity.rds"))

method_id <- "scvelo"
params_id <- "dynamical"
dataset_id <- "cycle_easy_seed1"
dataset_id <- "bifurcating_easy_seed1_1000"
dataset_id <- "converging_hard_seed1_1000"
dataset_id <- "linear_hard_seed1_1000"


# GROUNDTRUTH -------------------------------------------------------------
dataset <- read_rds(exp$dataset_file(dataset_id)) %>% dynwrap::add_dimred(dimred = dyndimred::dimred_landmark_mds)
traj_dimred <- dataset %>% dynwrap::project_trajectory(dataset$dimred)
dataset[names(traj_dimred)] <- traj_dimred

# model <- read_rds(exp$model_file(dataset_id))
groundtruth_velocity <- dataset$rna_velocity

# PREDICTION --------------------------------------------------------------
velocity <- read_rds(exp$velocity_file(dataset_id, method_id, params_id))
dataset <- dataset %>% scvelo::add_velocity(velocity = velocity)
predicted_velocity <- dataset$velocity_vector

# COMPARE -----------------------------------------------------------------
dimred <- dataset$dimred
dimred_future <- scvelo::embed_velocity(dataset, dimred)
dimred_diff <- dimred_future - dimred

progression_segments <- dataset$dimred_segment_progressions
dimred_segments <- dataset$dimred_segment_points

dimred_diff_segments <- bw_mean(from = dimred, to = dimred_segments, data = dimred_diff)

nexts <- progression_segments %>%
  mutate(index = row_number()) %>%
  group_by(from, to) %>%
  summarize(i = list(index[-n()]), j = list(index[-1]), diff = list(diff(percentage))) %>%
  ungroup() %>%
  unnest(c(i, j, diff))
dimred_truediff_segments <- dimred_segments[nexts$j,] - dimred_segments[nexts$i,]
dimred_velodiff_segments <- dimred_diff_segments[nexts$i, ]

simil <- dyngen.manuscript::paired_simil(dimred_truediff_segments, dimred_velodiff_segments, method = "cosine")
sum(simil * nexts$diff) / sum(nexts$diff)
mean(simil)

ggplot(mapping = aes(comp_1, comp_2)) + geom_point(data = as.data.frame(dimred), color = "lightgray") + theme_bw() +
  geom_point(aes(color = simil), data.frame(dimred_segments[nexts$i,], simil)) + scale_color_distiller(palette = "RdBu")

