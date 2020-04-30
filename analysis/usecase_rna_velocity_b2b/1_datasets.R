library(tidyverse)
library(rlang)
library(dyngen)
library(dyngen.manuscript)

exp <- start_analysis("usecase_rna_velocity_b2b")

# file.remove(exp$result("design_datasets.rds"))

# setup dataset design
design_datasets <- exp$result("design_datasets.rds") %cache% {
  crossing(
    seed = 1,
    # seed = 1:3,
    backbone_name = names(list_backbones()),
    tr_rate_multiplier = c(1, 5, 25)
  ) %>%
    mutate(id = paste0(backbone_name, "_", seed, "_", tr_rate_multiplier))
}

#' @examples
#' design_datasets %>% dynutils::extract_row_to_list(1) %>% list2env(.GlobalEnv)

pwalk(design_datasets, function(id, seed, backbone_name, tr_rate_multiplier) {
  if (!file.exists(exp$dataset_file(id))) {

    cat("## Generating ", id, "\n", sep = "")
    set.seed(seed)

    kinetic_params <- kinetics_default()
    kinetic_params$sampler_tfs <- function(...) {
      x <- kinetics_default()$sampler_tfs(...)
      x$transcription_rate <- x$transcription_rate * tr_rate_multiplier
      x
    }

    backbone <- dyngen::list_backbones()[[backbone_name]]()
    model <-
      initialise_model(
        id = id,
        num_tfs = nrow(backbone$module_info),
        num_targets = 70,
        num_hks = 0,
        backbone = backbone,
        num_cells = 1000,
        kinetics_params = kinetic_params,
        simulation_params = simulation_default(
          census_interval = 10,
          compute_rna_velocity = TRUE,
          store_reaction_propensities = TRUE,
          experiment_params = simulation_type_wild_type(
            num_simulations = 100
          )
        ),
        num_cores = 7,
        download_cache_dir = "~/.cache/dyngen",
        verbose = TRUE
      )
    generate_dataset(
      model,
      output_dir = exp$dataset_folder(id),
      make_plots = TRUE
    )

    gc()
  }
})


pwalk(design_datasets, function(id, seed, backbone_name, tr_rate_multiplier) {
  cat(id, "\n", sep = "")
  model <- read_rds(exp$model_file(id))
  dataset <- read_rds(exp$dataset_file(id))

  reac_prop <- model$simulations$reaction_propensities
  transcription_prop <- reac_prop[, paste0("transcription_", dataset$feature_ids)]
  degradation_prop <- reac_prop[, paste0("premrna_degradation_", dataset$feature_ids)] + reac_prop[, paste0("mrna_degradation_", dataset$feature_ids)]
  colnames(transcription_prop) <- colnames(degradation_prop) <- dataset$feature_ids
  groundtruth_velocity <- as(log2(transcription_prop + 1) - log2(degradation_prop + 1), "dgCMatrix")

  model$simulations$velocity_vector <- groundtruth_velocity

  groundtruth_velocity_exp <- groundtruth_velocity[model$experiment$cell_info$step_ix, ]
  rownames(groundtruth_velocity_exp) <- model$experiment$cell_info$cell_id
  model$experiment$velocity_vector <- groundtruth_velocity_exp

  dataset$velocity_vector <- groundtruth_velocity_exp

  write_rds(dataset, exp$dataset_file(id), compress = "gz")
  write_rds(model, exp$model_file(id), compress = "gz")
})


