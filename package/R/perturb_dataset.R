
#' Generate difference in abundance
#' @export
generate_diff_abundance <- function(model, froms, tos, ratios, normalise_gs = T){
  m <- model
  if(normalise_gs){
    m <- normalise_goldstandard(m)
  }
  for(i in seq(from = 1, to = length(froms))){
    print(froms[i])
    m <- adapt_gold_standard_network(m, froms[i], tos[i], ratios[i])
  }

  m <- m %>% dyngen:::generate_cells()

}

#' @export
normalise_goldstandard <- function(model){
  model$gold_standard$network <- model$gold_standard$network %>% mutate_at("length", ~sapply(model$gold_standard$network$length, function(z) z / sum(model$gold_standard$network$length)))
  model
}

#' @export
adapt_gold_standard_network <- function(model, from, to, ratio){
  netw <- model$gold_standard$network
  id1 <- which((to == netw$to & from == netw$from))
  id2 <- which((to == netw$to & from == netw$from))
  sq <- seq(from = id1, to = id2)

  old <- netw$length[id1:id2]
  newl <- sapply(seq(from = 1, to = length(netw$length)), function(i) ifelse(is.element(i, sq), netw$length[[i]] * ratio, netw$length[[i]]))

  model$gold_standard$network$length <- newl
  model
}

#' @export
combine_models <- function(model1, model2){
  model12 <- model1

  m1gs <- model1$gold_standard
  m2gs <- model2$gold_standard
  # No need to normalize to 100% here --> if stuff breaks, this may be why though.
  model12$gold_standard <- list(
    meta = bind_rows(
      m1gs$meta %>% mutate_at(c("from", "to", "from_", "to_"), ~paste0("left_", .)),
      m2gs$meta %>% mutate_at(c("from", "to", "from_", "to_"), ~paste0("right_", .))
    ),
    counts = rbind(
      m1gs$counts,
      m2gs$counts
    ),
    network = bind_rows(
      m1gs$network %>% mutate_at(c("from", "to"), ~paste0("left_", .)),
      m2gs$network %>% mutate_at(c("from", "to"), ~paste0("right_", .))

    )
  )
  m1sim <- model1$simulations
  m2sim <- model2$simulations
  num_m1sim <- max(m1sim$meta$simulation_i)
  model12$simulations <- list(
    meta = bind_rows(
      m1sim$meta %>% mutate_at(c("from", "to"), ~paste0("left_", .)),
      m2sim$meta %>% mutate_at(c("from", "to"), ~paste0("right_", .)) %>% mutate(simulation_i = simulation_i + num_m1sim)
    ),
    counts = rbind(
      m1sim$counts,
      m2sim$counts
    ),
    regulation = rbind(
      m1sim$regulation,
      m2sim$regulation
    ),
    reaction_firings = rbind(
      m1sim$reaction_firings,
      m2sim$reaction_firings
    ),
    reaction_propensities = rbind(
      m1sim$reaction_propensities,
      m2sim$reaction_propensities
    )
  )
  m1simsys <- model1$simulation_system
  m2simsys <- model2$simulation_system
  model12$simulation_system <- list(
    reactions = rbind(
      m1simsys$reactions,
      m2simsys$reactions
    ),
    molecule_ids = rbind(
      m1simsys$molecule_ids,
      m2simsys$molecule_ids
    ),
    initial_state = rbind(
      m1simsys$initial_state,
      m2simsys$initial_state
    ),
    parameters = rbind(
      m1simsys$parameters,
      m2simsys$parameters
    ),
    burn_variables = rbind(
      m1simsys$burn_variables,
      m2simsys$burn_variables
    )

  )

  model12 <- model12 %>%
    dyngen:::calculate_dimred()  %>%
    dyngen:::generate_experiment()
  model12
}
