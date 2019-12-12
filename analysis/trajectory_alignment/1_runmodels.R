library(dyngen)

# Simple linear backbone, with 22 modules.
backbone <- bblego(
  bblego_start("A", type = "simple", num_modules = 4),
  bblego_linear("A", "B", type = "simple", num_modules = 6),
  bblego_linear("B", "C", type = "simple", num_modules = 6),
  bblego_end("C", type = "simple", num_modules = 6)
)


runmodel <- function(){
  model <- initialise_model(
    num_tfs = 50,
    num_targets = 200,
    num_hks = 200,
    num_cells = 1000,
    backbone = backbone,
    verbose = TRUE,
    num_cores = 8,
    distance_metric = "pearson",
    tf_network_params = tf_network_default(min_tfs_per_module = 2, sample_num_regulators = function() 1),
    simulation_params = simulation_default(census_interval = .01, experiment_params = bind_rows(simulation_type_wild_type(num_simulations = 32), simulation_type_knockdown(num_simulations = 0)))
  ) %>%
    generate_tf_network() %>%
    generate_feature_network() %>%
    generate_kinetics() %>%
    generate_gold_standard() %>%
    generate_cells() %>%
    generate_experiment()
  model
}


get_noisy_data <- function(ds, noise_perc){
  smp <- sample(length(ds$counts), size = length(ds$counts) * noise_perc, replace=FALSE)
  ds$counts[smp] <- sample(ds$counts[smp])
  ds
}

get_all_noisy <- function(ds, start, stop, step){
  lapply(seq(start, stop, step), get_noisy_data, ds=ds)
}


change_speed <- function(model, target, rate) {
  param_id <- which(target == model$feature_info$feature_id)[[1]]

  model$feature_info$wpr[param_id] <- model$feature_info$wpr[param_id] * rate
  model$feature_info$xdr[param_id] <- model$feature_info$xdr[param_id] * rate

  param_wpr <- paste0("wpr_", target)
  param_xdr <- paste0("xdr_", target)

  model$simulation_system$parameters[param_wpr] <- model$simulation_system$parameters[param_wpr] * rate
  model$simulation_system$parameters[param_xdr] <- model$simulation_system$parameters[param_xdr] * rate

  model
}


surpress_module <- function(model, module){
  fts_to_surpress <- which(module == model$feature_network$to_module)
  for(index in fts_to_surpress){
    ft <- model$feature_network$from[[index]]
    model <- change_speed(model, ft, 0)
  }
  model <- model %>% generate_gold_standard() %>% generate_cells() %>% generate_experiment()
  model
}

