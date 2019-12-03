library(tidyverse)
library(dyngen)

derived_folder <- "derived_files/network_inference/"

set.seed(1)

backbone <- backbone_bifurcating()
wanted_genes <- 200

num_tfs <- nrow(backbone$module_info) * 3
num_targets <- round((wanted_genes - num_tfs) / 2)

model <-
  initialise_model(
    id = "small_bifurcating",
    num_cells = 500,
    num_tfs = num_tfs,
    num_targets = num_targets,
    num_hks = wanted_genes - num_tfs - num_targets,
    distance_metric = "pearson",
    backbone = backbone,
    tf_network_params = tf_network_default(min_tfs_per_module = 2, sample_num_regulators = function() 2),
    feature_network_params = feature_network_default(
      target_resampling = 5000,
      damping = 1
    ),
    kinetics_params = kinetics_default(),
    gold_standard_params = gold_standard_default(),
    simulation_params = simulation_default(
      census_interval = .1,
      burn_time = 10,
      total_time = simtime_from_backbone(backbone),
      experiment_params = bind_rows(
        simulation_type_wild_type(num_simulations = 10),
        simulation_type_knockdown(num_simulations = 10, num_genes = sample(1:10, 10, replace = TRUE))
      ),
      store_reaction_propensities = TRUE,
      store_grn = TRUE,
      perform_dimred = TRUE
    ),
    experiment_params = experiment_snapshot(),
    verbose = TRUE,
    download_cache_dir = "~/.cache/dyngen",
    num_cores = 8
  )
folder <- paste0(derived_folder, model$id, "/")
dir.create(folder, recursive = TRUE)
generate_dataset(model, output_dir = folder, make_plots = TRUE, store_grn = TRUE)
