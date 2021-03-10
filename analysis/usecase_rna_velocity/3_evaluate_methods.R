library(tidyverse)
library(dyngen.manuscript)

exp <- start_analysis("usecase_rna_velocity")

design_datasets <- read_rds(exp$result("design_datasets.rds"))
design_velocity <- read_rds(exp$result("design_velocity.rds"))

#' @examples
#' design_velocity %>% dynutils::extract_row_to_list(16) %>% list2env(.GlobalEnv)
# file.remove(exp$result("scores.rds"))

traj_types <- pmap_dfr(design_datasets, function(id, ...) {
  if (!file.exists(exp$dataset_file(id))) return(NULL)
  dataset <- read_rds(exp$dataset_file(id))
  tibble(
    dataset_id = id,
    dataset_trajectory_type = factor(dataset$trajectory_type, dynwrap::trajectory_types$id)
  )
})

# Calculate scores -----
scores <- exp$result("scores.rds") %cache% {
  library(furrr)
  plan(multiprocess)
  outs <- future_pmap(
    .progress = TRUE,
    design_velocity %>% mutate(i = row_number()),
    function(method_id, params_id, dataset_id, i, ...) {
      cat(i, "/", nrow(design_velocity), ": ", dataset_id, "\n", sep = "")
      if (!file.exists(exp$dataset_file(dataset_id))) return(NULL)
      dataset <- read_rds(exp$dataset_file(dataset_id))
      groundtruth_velocity <- dataset$rna_velocity
      groundtruth_velocity <- as(groundtruth_velocity + runif(length(groundtruth_velocity), -1e-8, 1e-8), "dgCMatrix")

      velocity_file <- exp$velocity_file(dataset_id, method_id, params_id)
      if (!file.exists(velocity_file)) return(NULL)
      velocity <- read_rds(velocity_file)
      dataset <- dataset %>% scvelo::add_velocity(velocity = velocity)
      predicted_velocity <- dataset$velocity_vector
      predicted_velocity <- as(predicted_velocity + runif(length(groundtruth_velocity), -1e-8, 1e-8), "dgCMatrix")


      # calculate cosine similarity
      dimred <- dataset$dimred
      dimred_future <- scvelo::embed_velocity(dataset, dimred)
      dimred_diff <- dimred_future - dimred

      waypoints <- dynwrap::select_waypoints(dataset)
      dimred_waypoints <- dynwrap::project_waypoints(dataset, waypoints = waypoints, space = dimred)#, trajectory_projection_sd = .01)
      dimred_diff_waypoints <- dynwrap::project_waypoints(dataset, waypoints = waypoints, space = dimred_diff)#, trajectory_projection_sd = .01)
      nexts <- waypoints$progressions %>%
        mutate(index = row_number()) %>%
        group_by(from, to) %>%
        summarize(i = list(index[-n()]), j = list(index[-1]), diff = list(diff(percentage))) %>%
        ungroup() %>%
        unnest(c(i, j, diff))
      dimred_truediff_segments <- dimred_waypoints[nexts$j,] - dimred_waypoints[nexts$i,]
      dimred_velodiff_segments <- dimred_diff_waypoints[nexts$i, ]

      simil <- dyngen.manuscript::paired_simil(dimred_truediff_segments, dimred_velodiff_segments, method = "cosine")

      # return output
      per_cell <- tibble(
        method_id, params_id, dataset_id,
        feature_id = dataset$feature_ids,
        cor = paired_simil(predicted_velocity, groundtruth_velocity, method = "spearman", margin = 2)[feature_id],
      ) %>%
        left_join(design_datasets, c("dataset_id" = "id")) %>%
        left_join(traj_types, by = "dataset_id")

      per_waypoint <- bind_cols(
        tibble(method_id, params_id, dataset_id, simil),
        waypoints$progressions[nexts$i, ],
        dimred_waypoints[nexts$i,] %>% {magrittr::set_colnames(., paste0("from_", colnames(.)))} %>% as.data.frame,
        dimred_waypoints[nexts$j,] %>% {magrittr::set_colnames(., paste0("to_", colnames(.)))} %>% as.data.frame,
        dimred_diff_waypoints[nexts$i, ] %>% {magrittr::set_colnames(., paste0("velo_", colnames(.)))} %>% as.data.frame
      ) %>%
        left_join(design_datasets, c("dataset_id" = "id")) %>%
        left_join(traj_types, by = "dataset_id")

      summ <- tibble(
        method_id, params_id, dataset_id,
        cor = cor(as.vector(predicted_velocity), as.vector(groundtruth_velocity)),
        mean_corr = mean(per_cell$cor),
        mean_cosine = sum(simil * nexts$diff) / sum(nexts$diff)
      ) %>%
        left_join(design_datasets, c("dataset_id" = "id")) %>%
        left_join(traj_types, by = "dataset_id")

      lst(per_cell, per_waypoint, summ)
    }
  )
  list(
    per_cell = map_dfr(outs, "per_cell"),
    per_waypoint = map_dfr(outs, "per_waypoint"),
    summ = map_dfr(outs, "summ")
  )
}
