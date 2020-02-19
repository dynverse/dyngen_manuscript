# returns the location of the dataset and model
create_dataset <- function(id = "test") {
  assertthat::assert_that(is.function(dataset_file))

  file_dataset <- dataset_file(id, "dataset.rds")

  reread(file_dataset, function() {
    set.seed(10)
    wsr_multiplier <- 2
    backbone <- backbone_linear()
    backbone$expression_patterns$time <- backbone$expression_patterns$time * wsr_multiplier
    model <-
      initialise_model(
        num_tfs = 12,
        num_targets = 30,
        num_hks = 15,
        backbone = backbone,
        verbose = TRUE,
        num_cells = 1000,
        download_cache_dir = "~/.cache/dyngen",
        kinetics_params = kinetics_default(sample_wsr = function(n) rnorm(n, 2, 1) %>% pmax(0.5)),
        gold_standard_params = gold_standard_default(
          tau = .001 * wsr_multiplier, census_interval = 0.01 * wsr_multiplier
        ),
        simulation_params = simulation_default(
          ssa_algorithm = GillespieSSA2::ssa_etl(tau = 0.001 * wsr_multiplier),
          census_interval = .01 * wsr_multiplier,
          burn_time = simtime_from_backbone(backbone, burn = TRUE),
          total_time = simtime_from_backbone(backbone),
          num_simulations = 1
        )
      )
    model <- generate_tf_network(model)
    model <- generate_feature_network(model)
    model <- generate_kinetics(model)
    model <- generate_gold_standard(model)
    model <- generate_cells(model)
    model <- generate_experiment(model)
    dataset <- wrap_dataset(model)
    dataset$id <- id

    write_rds(model, dataset_file(id, "model.rds"))
    write_rds(dataset, dataset_file(id, "dataset.rds"))
  })
}




