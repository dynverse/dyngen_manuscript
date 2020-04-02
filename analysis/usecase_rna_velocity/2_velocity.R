library(tidyverse)
library(dyngen.manuscript)

exp <- start_analysis("usecase_rna_velocity")

# if need be: reticulate::install_miniconda()
try({
  reticulate::import("velocyto")
}, error = function(e) {
  reticulate::py_install("Cython", pip = TRUE)
  reticulate::py_install("velocyto", pip = TRUE)
})

try({
  reticulate::import("scvelo")
}, error = function(e) {
  reticulate::py_install("scvelo", pip = TRUE)
})

design_velocity <- exp$result("design_velocity.rds") %cache% {
  tribble(
    ~method_id, ~params, ~params_id,
    "velocyto", list(assumption = "constant_velocity"), "constant_velocity",
    "velocyto", list(assumption = "constant_unspliced"), "constant_unspliced",
    "scvelo", list(mode = "deterministic"), "deterministic",
    "scvelo", list(mode = "dynamical"), "dynamical",
    "scvelo", list(mode = "stochastic"), "stochastic",
  ) %>%
    crossing(
      read_rds(exp$result("design_datasets.rds")) %>% select(dataset_id = id)
    )
}

#' @examples
#' design_velocity %>% dynutils::extract_row_to_list(1) %>% list2env(.GlobalEnv)

pwalk(
  design_velocity,
  function(dataset_id, method_id, params, params_id) {
    dataset <- read_rds(exp$dataset_file(dataset_id))

    exp$velocity_file(dataset_id, method_id, params_id) %cache% {
      velocity <-
        if (method_id == "scvelo") {
          scvelo::get_velocity(
            spliced = dataset$expression,
            unspliced = dataset$expression_unspliced,
            mode = params$mode %||% "deterministic"
          )
        } else if (method_id == "velocyto") {
          rnav_run_velocyto(
            spliced = dataset$expression,
            unspliced = dataset$expression_unspliced,
            assumption = params$assumption %||% "constant_velocity"
          )
        }

      pkl <- exp$temporary("velocity/", dataset_id, "-", method_id, "-", params_id, "/scvelo.pkl")
      if (method_id == "scvelo") {
        reticulate::py_save_object(velocity$scvelo, filename = pkl)
      }

      velocity
    }
  }
)


