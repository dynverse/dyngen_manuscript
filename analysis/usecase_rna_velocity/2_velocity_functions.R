run_velocity <- function(
  dataset_id = "initial",
  method_id = "scvelo",
  params = list(),
  params_id = "default"
) {
  dataset <- load_dataset(dataset_id)

  experiment_file <- velocity_file(dataset_id, method_id, params_id)

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

velocity_file <- function(dataset_id, method_id, params_id) {
  dynamic_file(derived_file("velocity", dataset_id, method_id, params_id))
}


find_updown_feature <- function(model) {
  feature_ids <- model$feature_network %>%
    group_by(to) %>%
    summarise(both = any(effect == -1) && any(effect == 1)) %>%
    filter(both) %>%
    pull(to)
  feature_ids
}

plot_velocity_updown <- function(dataset, model, feature_id = NULL) {
  if (is.null(feature_id)) {
    # select a feature that goes up and down
    feature_ids <- model$feature_network %>%
      group_by(to) %>%
      summarise(both = any(effect == -1) && any(effect == 1)) %>%
      filter(both)
    if (length(feature_ids) > 0) {
      feature_id == first(feature_id)
    } else {
      feature_id <- sample(dataset$feature_ids, 1)
    }
  }

  meta <- model$simulations$meta %>%
    mutate(step_ix = row_number()) %>%
    # filter(simulation_i == 1)
    identity()
  counts <- model$simulations$counts[meta$step_ix, ]
  plotdata <- bind_rows(
    tibble(
      expression =  counts[,paste0("w_", feature_id)],
      step_ix = meta$step_ix,
      simulaton_i = meta$simulation_i,
      sim_time = meta$sim_time,
      molecule = "unspliced"
    ),
    tibble(
      expression =  counts[,paste0("x_", feature_id)],
      step_ix = meta$step_ix,
      simulaton_i = meta$simulation_i,
      sim_time = meta$sim_time,
      molecule = "spliced"
    )
  )
  plotdata %>%
    pivot_wider(names_from = "molecule", values_from = "expression") %>%
    ggplot(aes(spliced, unspliced)) + geom_path(aes(color = sim_time))
}




load_velocity <- function(dataset_id, method_id, params_id, folder = velocity_file(dataset_id, method_id, params_id)) {
  velocity <- read_rds(folder("velocity.rds"))
  if(file.exists(folder("scvelo.pkl"))) {
    velocity$scvelo <- reticulate::py_load_object(folder("scvelo.pkl"))
  }
  velocity
}

load_dataset_velocity <- function(dataset_id, method_id, params_id) {
  dataset <- load_dataset(dataset_id)
  dataset <- dataset %>% add_velocity(velocity = load_velocity(dataset_id, method_id, params_id))
  dataset
}

