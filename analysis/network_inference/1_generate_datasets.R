library(tidyverse)
library(dyngen)
library(rlang)

derived_folder <- "derived_files/network_inference/"

if (!dir.exists(derived_folder)) dir.create(derived_folder, recursive = TRUE)

backbones <- list_backbones()

backbones$binary_tree <- NULL
backbones$branching <- NULL

walk(names(backbones), function(backbone_name) {

  folder <- paste0(derived_folder, "datasets/", backbone_name, "/")
  if (!dir.exists(folder)) dir.create(folder, recursive = TRUE)
  dataset_file <- paste0(folder, "dataset.rds")

  if (!file.exists(dataset_file)) {
    cat("============== Running ", backbone_name, " ==============\n", sep = "")

    set.seed(1)
    backbone <- backbones[[backbone_name]]()
    wanted_genes <- 500

    num_tfs <- nrow(backbone$module_info) * 3
    num_targets <- round((wanted_genes - num_tfs) / 2)
    model <-
      initialise_model(
        num_cells = 5000,
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
            simulation_type_wild_type(num_simulations = 50),
            simulation_type_knockdown(num_simulations = 150, num_genes = sample(1:10, 150, replace = TRUE))
          ),
          store_reaction_propensities = TRUE,
          store_grn = TRUE
        ),
        experiment_params = experiment_snapshot(),
        verbose = TRUE,
        download_cache_dir = "~/.cache/dyngen",
        num_cores = 8
      )

    generate_dataset(model, output_dir = folder, make_plots = TRUE, store_grn = TRUE)

    dat <- read_rds(dataset_file)
    dat$id <- backbone_name
    write_rds(dat, dataset_file, compress = "gz")
  }
})

