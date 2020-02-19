library(tidyverse)

dataset_design <- tibble(
  seed = 1:3
) %>%
  crossing(backbone = c("linear", "linear_simple", "bifurcating", "cycle", "disconnected")) %>%
  mutate(id = paste0(backbone, "_", seed))


design_velocity <- tribble(
  ~method_id, ~params, ~params_id,
  "velocyto", list(assumption = "constant_velocity"), "constant_velocity",
  "velocyto", list(assumption = "constant_unspliced"), "constant_unspliced",
  "scvelo", list(mode = "deterministic"), "deterministic",
  "scvelo", list(mode = "dynamical"), "dynamical",
  "scvelo", list(mode = "stochastic"), "stochastic",
) %>%
  crossing(
    dataset_design %>% select(dataset_id = id)
  )
