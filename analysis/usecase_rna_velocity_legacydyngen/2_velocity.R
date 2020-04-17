library(tidyverse)
library(dyngen.manuscript)

exp <- start_analysis("usecase_rna_velocity_legacydyngen")

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

design_velocity <- exp$result("design_velocity.rds") %cache% {
  tribble(
    ~method_id, ~params, ~params_id,
    "velocyto", list(assumption = "constant_velocity"), "constant_velocity",
    "velocyto", list(assumption = "constant_unspliced"), "constant_unspliced",
    "scvelo", list(mode = "deterministic"), "deterministic",
    "scvelo", list(mode = "dynamical"), "dynamical_1_biorxiv",
    "scvelo", list(mode = "dynamical", var_names = "all"), "dynamical_2_varnamesall",
    "scvelo", list(mode = "dynamical", layer = "imputed"), "dynamical_3_useimputed",
    "scvelo", list(mode = "dynamical", var_names = "all", layer = "imputed"), "dynamical_4_both",
    "scvelo", list(mode = "stochastic"), "stochastic",
  ) %>%
    crossing(
      read_rds(exp$result("design_datasets.rds")) %>% select(dataset_id = id)
    )
}

#' @examples
#' design_velocity %>% dynutils::extract_row_to_list(40) %>% list2env(.GlobalEnv)

pwalk(
  design_velocity %>% mutate(rn = row_number()),
  function(dataset_id, method_id, params, params_id, rn) {
    cat(rn, "\n", sep = "")
    dataset <- read_rds(exp$dataset_file(dataset_id))

    try({
      exp$velocity_file(dataset_id, method_id, params_id) %cache% {
        params$spliced <- dataset$expression
        params$unspliced <- dataset$expression_unspliced

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
    })
  }
)


design_velocity[!pmap_lgl(
  design_velocity,
  function(dataset_id, method_id, params, params_id, rn) {
    file.exists(exp$velocity_file(dataset_id, method_id, params_id))
  }),]


