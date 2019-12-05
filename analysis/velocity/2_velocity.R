library(tidyverse)

source("analysis/velocity/2_velocity_functions.R")

devtools::load_all('~/thesis/projects/dynverse/libraries/scvelo/')

pmap(design_velocity, run_velocity)

list2env(dynutils::extract_row_to_list(design_velocity, 1), .GlobalEnv)

dataset <- load_dataset(dataset_id)
model <- load_model(dataset_id)
velociy <- load_velocity(dataset_id, method_id, params_id)

# Plot the spliced vs unspliced changes
# feature_id <- model$feature_info$feature_id[[2]]
feature_id <- "B3_TF1"
meta <- model$simulations$meta %>%
  mutate(step_ix = row_number()) %>%
  # filter(simulation_i == 1)
  identity()
counts <- model$simulations$counts[meta$step_ix, ]
plotdata <- bind_rows(
  tibble(
    expression =  counts[,paste0("w_", feature_id)],
    step_ix = meta$step_ix,
    simulation_i = meta$simulation_i,
    sim_time = meta$sim_time,
    molecule = "unspliced"
  ),
  tibble(
    expression =  counts[,paste0("x_", feature_id)],
    step_ix = meta$step_ix,
    simulation_i = meta$simulation_i,
    sim_time = meta$sim_time,
    molecule = "spliced"
  )
)
plotdata %>%
  pivot_wider(names_from = "molecule", values_from = "expression") %>%
  ggplot(aes(spliced, unspliced)) + geom_path(aes(color = sim_time))



plotdata %>%
  filter(simulation_i == 1) %>%
  ggplot(aes(step_ix, expression)) + geom_line(aes(color = molecule))


