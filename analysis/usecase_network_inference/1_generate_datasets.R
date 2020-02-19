library(tidyverse)
library(dyngen)
library(dyngen.manuscript)

analysis <- start_analysis("usecase_network_inference")

backbones <- list_backbones()

for (backbone_name in names(backbones)) {
  cat("=============== Generating ", backbone_name, "===============\n", sep = "")
  tryCatch({
    set.seed(1)

    folder <- analysis$temporary("datasets/", backbone_name)

    if (!file.exists(paste0(folder, "/model.rds"))) {
      backbone <- backbones[[backbone_name]]()
      wanted_genes <- 200

      num_tfs <- nrow(backbone$module_info) * 3
      num_targets <- round((wanted_genes - num_tfs) / 2)

      model <-
        initialise_model(
          id = paste0("small_", backbone_name),
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
              simulation_type_wild_type(num_simulations = 25),
              simulation_type_knockdown(num_simulations = 25, num_genes = sample(1:10, 25, replace = TRUE))
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

      generate_dataset(model, output_dir = folder, make_plots = TRUE, store_grn = TRUE)
    }
  }, error = function(e) {
    print(e)
  })
}
