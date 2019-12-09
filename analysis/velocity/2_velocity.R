library(tidyverse)

source("analysis/velocity/2_velocity_functions.R")

# load either velocyto or scvelo
source("analysis/velocity/2_velocity_velocyto.R")
devtools::load_all('~/thesis/projects/dynverse/libraries/scvelo/')

# run either velocyto or scvelo
pmap(design_velocity %>% filter(method_id == "velocyto"), run_velocity)
# pmap(design_velocity %>% filter(method_id == "scvelo"), run_velocity)

list2env(dynutils::extract_row_to_list(design_velocity %>% filter(method_id == "velocyto"), 1), .GlobalEnv)

dataset <- load_dataset(dataset_id)
model <- load_model(dataset_id)

run_velocity(dataset_id, method_id, params, params_id)
velocity <- load_velocity(dataset_id, method_id, params_id)




