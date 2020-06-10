
#' Generate difference in abundance
#'
#' Multiplies the abundance of cells on particular edges by a given ratio
#'
#' @param model A dyngen model for which to alter the abundance in generated cells
#' @param froms The start milestone of each edge to alter
#' @param tos The end milestone of each edge to alter
#' @param ratios The ratio by which to alter the abundances
#'
#' @export
generate_diff_abundance <- function(model, froms, tos, ratios) {
  tib <- tibble(
    from = froms,
    to = tos,
    ratio = ratios
  )

  model$gold_standard$network <-
    model$gold_standard$network %>%
    full_join(tib, by = c("from", "to")) %>%
    mutate(
      ratio = ratio %|% 1,
      length = length * ratio
    ) %>%
    select(-ratio)


  model
}

#' Combine two dyngen models
#'
#' @param model1 the first model
#' @param model2 the second model
#'
#' @export
combine_models <- function(model1, model2){
  model12 <- model1

  # combine gold standard
  m1gs <- model1$gold_standard
  m2gs <- model2$gold_standard
  model12$gold_standard <- list(
    mod_changes = bind_rows(
      m1gs$mod_changes %>% mutate_at(c("from", "to", "from_", "to_"), ~paste0("left_", .)),
      m2gs$mod_changes %>% mutate_at(c("from", "to", "from_", "to_"), ~paste0("right_", .))
    ),
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

  # combine simulations
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
    )
    # could also propensities, rna velocity, etc
  )

  # recalculate the dimred
  model12 <- model12 %>%
    dyngen:::calculate_dimred()

  model12
}
