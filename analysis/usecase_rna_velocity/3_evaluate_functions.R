pairwise_correlations <- function(x, y) {
  assertthat::assert_that(all(colnames(x) == colnames(y)))
  assertthat::assert_that(all(nrow(x) == nrow(y)))
  map_dbl(seq_len(nrow(x)), function(i) cor(x[i, ], y[i, ], method = "spearman"))
  # map_dbl(seq_len(nrow(x)), function(i) mean(sign(x[i, ]) == sign(y[i, ])))
  # map_dbl(seq_len(nrow(x)), function(i) sum(x[i, ] * y[i, ])/(sqrt(sum(x[i, ]^2)) * sqrt(sum(y[i, ]^2))))
}



run_metric <- function(
  dataset_id = "initial",
  method_id = "scvelo",
  params = list(),
  params_id = "default",
  metric_id = "pairwise_correlations"
) {
  dataset <- load_dataset_velocity(dataset_id, method_id, params_id)

  experiment_file <- dynamic_file(metric_file(dataset_id, method_id, params_id))

  reread(
    experiment_file("scores.rds"),
    function() {
      if(metric_id == "pairwise_correlations") {

      }

      write_rds(scores, experiment_file("scores.rds"))
    }
  )
}


metric_file <- function(dataset_id, method_id, params_id, metric_id) {
  dynamic_file(derived_file("metrics", dataset_id, method_id, params_id, metric_id))
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
