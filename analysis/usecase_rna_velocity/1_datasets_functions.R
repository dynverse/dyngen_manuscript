create_dataset <- function(id = "test", seed = 10, backbone = "linear_simple") {
  requireNamespace("dyngen")

  assertthat::assert_that(is.function(dataset_file))

  file_dataset <- dataset_file(id, "dataset.rds")

  reread(file_dataset, function() {
    set.seed(seed)
    wsr_multiplier <- 5
    backbone <- dyngen::list_backbones()[[backbone]]()
    # backbone$expression_patterns$time <- backbone$expression_patterns$time * wsr_multiplier
    model <-
      initialise_model(
        num_tfs = 20,
        num_targets = 50,
        num_hks = 15,
        backbone = backbone,
        verbose = TRUE,
        num_cells = 1000,
        download_cache_dir = "~/.cache/dyngen",
        kinetics_params = kinetics_default(
          # sample_wsr = function(n) rep(1, n)#rnorm_bounded(n, 2, 1, min = 1)
        ),
        gold_standard_params = gold_standard_default(
          # tau = .001 * wsr_multiplier, census_interval = 0.01 * wsr_multiplier
        ),
        simulation_params = simulation_default(
          # ssa_algorithm = GillespieSSA2::ssa_etl(tau = 0.001 * wsr_multiplier),
          # census_interval = .01 * wsr_multiplier,
          burn_time = simtime_from_backbone(backbone, burn = TRUE),
          total_time = simtime_from_backbone(backbone),
          perform_dimred = F,
          # experiment_params = simulation_type_wild_type(num_simulations = 10),
          store_reaction_propensities = TRUE
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
