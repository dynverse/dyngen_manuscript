library(dyngen)

# Simple linear backbone, with 22 modules.
backbone <- bblego(
  bblego_start("A", type = "simple", num_modules = 4),
  bblego_linear("A", "B", type = "simple", num_modules = 6),
  bblego_linear("B", "C", type = "simple", num_modules = 6),
  bblego_end("C", type = "simple", num_modules = 4)
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
  ) %>%
    generate_tf_network() %>%
    generate_feature_network() %>%
    generate_kinetics() %>%
    generate_gold_standard() %>%
    generate_cells() %>%
    generate_experiment()
}


get_all_noisy_data <- function(ds, start, stop, step){
  lapply(seq(start, stop, step), function(noise_perc, ds_=ds){
    smp <- sample(length(ds_$counts), size = length(ds_$counts) * noise_perc, replace=FALSE)
    ds_$counts[smp] <- sample(ds_$counts[smp])
    ds_
  })
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



