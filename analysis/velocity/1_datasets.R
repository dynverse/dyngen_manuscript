library(tidyverse)

source("analysis/velocity/1_datasets_functions.R")
library(dyngen)

pwalk(dataset_design, create_dataset)




# characterise datasets
# devtools::load_all('~/thesis/projects/dynverse/dynwrap/')
# devtools::load_all('~/thesis/projects/dynverse/dynplot2/')
library(dynwrap)
library(dynplot2)

dataset <- load_dataset("1")
model <- load_model("1")

dataset <- dataset %>% add_dimred(dyndimred::dimred_landmark_mds, pair_with_velocity = F)

# dimred of dataset
dynplot_dimred(dataset) +
  geom_cell_point()

# features
feature_ids <- model$feature_info %>% filter(!burn) %>% group_by(module_id) %>% slice(1) %>% pull(feature_id)

dynplot_dimred(dataset) +
  geom_cell_point(aes(color = select_feature_expression(feature_id, .data))) +
  facet_wrap_data(feature_id = feature_ids) +
  scale_expression_fillcolour()


# one feature spliced vs unspliced
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

