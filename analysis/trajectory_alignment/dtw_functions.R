library(tidyverse)
library(dyngen)
library(dynplot)
library(gtools)
library(dtw)

calculate_correct_pseudotime <- function(dataset){
  length_pieces <- dataset$milestone_network$length/sum(dataset$milestone_network$length)
  nr_piece <- 0
  start <- "sA"
  traj <- numeric()
  names_traj <- list()
  # TODO convert to for for each piece
  while(start %in% dataset$milestone_network[["from"]]){
    idx <- which(start == dataset$milestone_network[["from"]])[[1]]
    end <- dataset$milestone_network[["to"]][idx]
    to_add <- sum(length_pieces[0:nr_piece])

    perc <-(dataset[["progressions"]] %>% filter(from == start & to == end))$percentage * length_pieces[nr_piece + 1] + to_add
    nams <-(dataset[["progressions"]] %>% filter(from == start & to == end))$cell_id

    traj <- c(traj, perc)
    names_traj <- c(names_traj, nams)

    start <- end
    nr_piece <- nr_piece + 1

  }

  names(traj) <- names_traj
  traj <- traj[mixedorder(names(traj))]
  traj
}

get_geodesic_distances_from_progressions <- function(dataset, waypoint_progressions){
  waypoint_milestone_percentages <- dynwrap::convert_progressions_to_milestone_percentages(cell_id, dataset$milestone_ids, dataset$milestone_network, waypoint_progressions)
  waypoint_milestone_percentages <- dplyr::rename(waypoint_milestone_percentages, waypoint_id = cell_id)

  geodesic_distances <- dynwrap::calculate_geodesic_distances(dataset, waypoint_milestone_percentages = waypoint_milestone_percentages, directed=FALSE)
  geodesic_distances
}

get_waypoint_progression <- function(dataset, nr){
  cell_id <- character()
  froms <- character()
  tos <- character()
  percentage <- numeric()

  length_pieces <- dataset$milestone_network$length/sum(dataset$milestone_network$length)
  idx <- 1
  for(piece in length_pieces){
    amount <- round(nr * piece)
    wps <- modify(list(1:amount), function(x){x/amount})[[1]]
    ids <- modify(list(1:amount), function(x){sprintf("w%d_%03d", idx, x)})[[1]]
    from <- rep(dataset$milestone_network$from[idx], amount)
    to <- rep(dataset$milestone_network$to[idx], amount)

    cell_id <- c(cell_id, ids)
    froms <- c(froms, from)
    tos <- c(tos, to)
    percentage <- c(percentage, wps)

    idx <- idx + 1

  }
  waypoint_progressions <- tibble(cell_id=cell_id, from=froms, to=tos, percentage=percentage)
  waypoint_progressions
}

interpolate_expression <- function(dataset, wps, gds){
  weighted_expr <- do.call('cbind', lapply(seq_along(wps), function(idx){
    dist <- gds[idx,]
    weighted <- exp(-(dist^2)/(0.1)^2)
    weighted <- weighted/sum(weighted)
    weighted_exp <- t(as.matrix(dataset$counts)) %*% weighted
    weighted_exp
  }))
  weighted_expr
}
