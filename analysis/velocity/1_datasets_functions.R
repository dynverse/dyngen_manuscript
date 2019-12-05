create_dataset <- function(id = "test", seed = 10) {
  requireNamespace("dyngen")

  assertthat::assert_that(is.function(dataset_file))

  file_dataset <- dataset_file(id, "dataset.rds")

  reread(file_dataset, function() {
    set.seed(seed)
    wsr_multiplier <- 2
    backbone <- backbone_linear_simple()
    backbone$expression_patterns$time <- backbone$expression_patterns$time * wsr_multiplier
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
          sample_wsr = function(n) rep(1, n)#rnorm_bounded(n, 2, 1, min = 1)
        ),
        gold_standard_params = gold_standard_default(
          tau = .001 * wsr_multiplier, census_interval = 0.01 * wsr_multiplier
        ),
        simulation_params = simulation_default(
          ssa_algorithm = GillespieSSA2::ssa_etl(tau = 0.001 * wsr_multiplier),
          census_interval = .01 * wsr_multiplier,
          burn_time = simtime_from_backbone(backbone, burn = TRUE),
          total_time = simtime_from_backbone(backbone),
          perform_dimred = F,
          experiment_params = simulation_type_wild_type(num_simulations = 10),
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


load_dataset <- function(id) {
  dataset <- read_rds(dataset_file(id, "dataset.rds"))
  dataset
}


load_model <- function(id) {
  model <- read_rds(dataset_file(id, "model.rds"))
  model
}





extract_groundtruth_velocity <- function(model) {
  if(is.null(model$simulations$reaction_propensities)) {
    stop("model has to contain reaction_propensities for estimation of velocity ground truth")
  }

  # extract for each spliced mRNA its relevant reactions
  feature_ids <- model$feature_info$feature_id
  relevant_molecules <- paste0("x_", feature_ids)

  relevant_reactions <- map2_dfr(seq_along(model$simulation_system$reactions), model$simulation_system$reactions, function(reaction_ix, reaction) {
    if(names(reaction$effect)[1] %in% relevant_molecules) {
      tibble(reaction_ix = reaction_ix, feature_id = names(reaction$effect)[1], effect = reaction$effect[1], name = reaction$name)
    }
  })

  propensities <- model$simulations$reaction_propensities[model$experiment$cell_info$step_ix, relevant_reactions$reaction_ix]
  rownames(propensities) <- model$experiment$cell_info$cell_id
  colnames(propensities) <- relevant_reactions$name

  formulas <- relevant_reactions %>%
    group_by(feature_id) %>%
    mutate(effect = ifelse(effect == 1, "production", "degradation")) %>%
    pivot_wider("feature_id", names_from = "effect", values_from = "name")

  propensity_ratios <- formulas %>%
    mutate(data = map2(production, degradation, function(production, degradation) {
      propensities[,production] / propensities[,degradation]
    })
    ) %>%
    pull(data) %>%
    do.call(rbind, .)
  propensity_ratios[is.na(propensity_ratios)] <- runif(sum(is.na(propensity_ratios)), 1-1e-10, 1+1e-10)

  rownames(propensity_ratios) <- str_replace_all(formulas$feature_id, "x_", "")
  t(propensity_ratios)
}

dataset_file <- dynamic_file(derived_file("datasets"))

load_dataset <- function(id) {
  dataset <- read_rds(dataset_file(id, "dataset.rds"))
  dataset
}
