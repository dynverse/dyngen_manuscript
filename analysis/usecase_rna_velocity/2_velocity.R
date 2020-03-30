library(tidyverse)
library(dyngen.manuscript)

exp <- start_analysis("usecase_rna_velocity")

design_datasets <- read_rds(exp$result("design_datasets.rds"))

design_velocity <- tribble(
  ~method_id, ~params, ~params_id,
  "velocyto", list(assumption = "constant_velocity"), "constant_velocity",
  "velocyto", list(assumption = "constant_unspliced"), "constant_unspliced",
  "scvelo", list(mode = "deterministic"), "deterministic",
  "scvelo", list(mode = "dynamical"), "dynamical",
  "scvelo", list(mode = "stochastic"), "stochastic",
) %>%
  crossing(
    design_datasets %>% select(dataset_id = id)
  )

# load either velocyto or scvelo
# source("analysis/velocity/2_velocity_velocyto.R")
# devtools::load_all('~/thesis/projects/dynverse/libraries/scvelo/')

#' @examples
#' design_velocity %>% dynutils::extract_row_to_list(1) %>% list2env(.GlobalEnv)

pmap(
  design_velocity %>% filter(method_id == "velocyto"),
  function(dataset_id, method_id, params, params_id) {
    dataset <- read_rds(exp$dataset_file(dataset_id))

    velocity_file <- exp$temporary("velocity/", dataset_id, "-", method_id, "-", params_id, "/velocity.rds")

    velocity_file %cache% function() {

    }
    reread(
      experiment_file("velocity.rds"),
      function() {
        if(method_id == "scvelo") {
          velocity <- scvelo::get_velocity(dataset$expression, dataset$expression_unspliced, mode = params$mode %||% "deterministic")
        } else if (method_id == "velocyto") {
          velocity <- get_velocity_velocyto(dataset$expression, dataset$expression_unspliced, assumption = params$assumption %||% "constant_velocity")
        }

        write_rds(velocity, experiment_file("velocity.rds"))

        if(method_id == "scvelo") {
          reticulate::py_save_object(velocity$scvelo, filename = experiment_file("scvelo.pkl"))
        }
      }
    )
  }
)





list2env(dynutils::extract_row_to_list(design_velocity %>% filter(method_id == "velocyto"), 1), .GlobalEnv)

dataset <- load_dataset(dataset_id)
model <- load_model(dataset_id)

run_velocity(dataset_id, method_id, params, params_id)
velocity <- load_velocity(dataset_id, method_id, params_id)




