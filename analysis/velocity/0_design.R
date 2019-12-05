library(tidyverse)

dataset_design <- tibble(
  seed = 1:3
) %>%
  mutate(id = as.character(seed))


design_velocity <- tribble(
  ~method_id, ~params, ~params_id,
  "scvelo", list(mode = "deterministic"), "deterministic",
  "scvelo", list(mode = "dynamical"), "dynamical",
  "scvelo", list(mode = "stochastic"), "stochastic",
) %>%
  crossing(
    dataset_design %>% select(dataset_id = id)
  )
