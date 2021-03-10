library(tidyverse)
library(dyngen.manuscript)

exp <- start_analysis("usecase_rna_velocity")

# if need be: reticulate::install_miniconda()
tryCatch({
  reticulate::import("velocyto")
}, error = function(e) {
  reticulate::py_install("Cython", pip = TRUE)
  reticulate::py_install("velocyto", pip = TRUE)
})

tryCatch({
  reticulate::import("scvelo")
}, error = function(e) {
  scvelo:::install_scvelo()
})

# file.remove(design_velocity <- exp$result("design_velocity.rds"))
design_velocity <- exp$result("design_velocity.rds") %cache% {
  tribble(
    ~method_id, ~params, ~params_id,
    "velocyto", list(assumption = "constant_velocity"), "constant_velocity",
    "velocyto", list(assumption = "constant_unspliced"), "constant_unspliced",
    "scvelo", list(mode = "deterministic"), "deterministic",
    "scvelo", list(mode = "dynamical", var_names = "all"), "dynamical",
    "scvelo", list(mode = "dynamical_residuals", var_names = "all"), "dynamical_residuals",
    "scvelo", list(mode = "stochastic"), "stochastic"
  ) %>%
    crossing(
      read_rds(exp$result("design_datasets.rds")) %>% select(dataset_id = id)
    )
}

#' @examples
#' design_velocity %>% dynutils::extract_row_to_list(40) %>% list2env(.GlobalEnv)

library(furrr)
plan(multiprocess)
future_pmap(
# pwalk(
  design_velocity %>% mutate(rn = row_number()),
  function(dataset_id, method_id, params, params_id, rn) {
    cat(rn, "/", nrow(design_velocity), ": ", method_id, " ", params_id, " on ", dataset_id, "\n", sep = "")
    if (!file.exists(exp$dataset_file(dataset_id))) return(NULL)
    dataset <- read_rds(exp$dataset_file(dataset_id))

    try({
      exp$velocity_file(dataset_id, method_id, params_id) %cache% {
        params$spliced <- dataset$counts
        params$unspliced <- dataset$counts_unspliced

        velocity_fun <-
          if (method_id == "scvelo") {
            scvelo::get_velocity
          } else if (method_id == "velocyto") {
            rnav_run_velocyto
          }

        velocity <- do.call(velocity_fun, params)

        pkl <- exp$temporary("velocity/", dataset_id, "-", method_id, "-", params_id, "/scvelo.pkl")
        if (method_id == "scvelo") {
          reticulate::py_save_object(velocity$scvelo, filename = pkl)
        }

        velocity
      }
      0
    })
  }
)


design_velocity[!pmap_lgl(
  design_velocity,
  function(dataset_id, method_id, params, params_id, rn) {
    file.exists(exp$velocity_file(dataset_id, method_id, params_id))
  }),]


