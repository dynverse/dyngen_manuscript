run_velocity <- function(dataset, method_id, param_id, params = list(), store = NULL) {
  id <- paste0(method_id, "-", param_id)
  if (!is.null(store) && store$exists(id)) {
    store$get(id)
  } else {
    if(method_id == "scvelo") {
      dataset <- add_velocity(dataset, mode = params$mode)

      # get differences
      differences <- dataset$velocity$layers[["velocity"]]
      dimnames(differences) <- dimnames(dataset$expression)
      differences[is.na(differences)] <- 0

      result <- list(
        expression_projected = dataset$expression_projected,
        differences = differences,
        raw = dataset$velocity
      )
    } else if (method_id == "current") {
      result <- list(
        expression_projected = dataset$expression,
        differences = dataset$expression,
        raw = NULL
      )
    } else if (method_id == "future") {
      expression_projected <- find_future_expression(dataset, step = params$step)

      result <- lst(
        expression_projected,
        differences = expression_projected - dataset$expression
      )
    } else if (method_id == "unspliced") {
      expression_projected <- dataset$expression_unspliced

      result <- lst(
        expression_projected,
        differences = expression_projected - dataset$expression
      )
    } else {
      stop("unknown method")
    }

    if(!is.null(store)) {
      store$set(id, result)
    }

    result
  }
}





find_future_expression <- function(dataset, step) {
  # add waypoints
  dataset <- dataset %>% add_waypoints(n_waypoints = 500)
  dataset$waypoints$geodesic_distances_reverse <- dynwrap::calculate_geodesic_distances(
    dataset,
    waypoint_milestone_percentages = dataset$waypoints$milestone_percentages,
    directed = "reverse"
  )


  # find closest waypoints
  find_future_waypoint <- function(dataset, step, cell_id) {
    names(which.min(abs(dataset$waypoints$geodesic_distances_reverse[dataset$waypoints$waypoints$waypoint_id, cell_id] - step)))
  }

  cell_ids <- dataset$cell_ids
  closest_waypoints <- map_chr(cell_ids, find_future_waypoint, dataset = dataset, step = step)

  if (any(is.na(closest_waypoints))) warning("SOme closest future waypoints are NA")

  # calculate expression at each waypoint
  bw <- 0.15
  weights <- dnorm(dataset$waypoints$geodesic_distances, sd = bw)
  weights <- weights/rowSums(weights)
  expression_waypoints <- weights %*% dataset$expression

  # return expression at requested waypoints
  expression_future <- expression_waypoints[closest_waypoints, ]
  rownames(expression_future) <- cell_ids

  expression_future
}



#
#
#
# find_future_expression <- function(model, step) {
#   closest_ix <- model$experiment$cell_info %>% pmap_int(function(simulation_i, sim_time, ...) {
#     model$simulations$meta %>%
#       mutate(ix = row_number()) %>%
#       filter(simulation_i == !!simulation_i) %>%
#       arrange(abs(sim_time - !!sim_time - step)) %>%
#       slice(10) %>%
#       pull(ix)
#   })
#
#   counts_future <- model$simulations$counts[closest_ix, model$feature_info$x]
#   colnames(counts_future) <- model$feature_info$feature_id
#   rownames(counts_future) <- model$experiment$cell_info$cell_id
#   expression_future <- as(log2(counts_future + 1), "dgCMatrix")
#   expression_future <- expression_future[, colnames(dataset$expression_projected)]
#
#   expression_future
# }
#
