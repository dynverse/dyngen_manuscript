#' Get gene expression of inferred waypoints
#'
#' @param dataset A dyngen dataset
#' @param amount The amount of pseudocells wanted
#' @param ws The windowsize of the Gaussian kernel
#'
#' @export
get_waypoint_expression <- function(dataset, amount, ws = 0.05){
  wpp1 <- get_waypoint_progression(dataset, amount)
  gd1 <- get_geodesic_distances_from_progressions(dataset, wpp1)
  wexpr1 <- interpolate_expression(dataset, wpp1$percentage, gd1, ws)
  pt <- wpp1$cumulative_percentages
  names(pt) <- wpp1$cell_id

  return(list(pseudotime = pt, expression = t(wexpr1)))
}

# get_waypoint_pseudotime <- function(dataset, amount){
#   wpp1 <- get_waypoint_progression(dataset, amount)
#   # mst_perc <- dynwrap::convert_progressions_to_milestone_percentages(cell_id, dataset$milestone_ids, dataset$milestone_network, wpp1)
#
#   ids <- c()
#   percentages <- c()
#   milestone <- "sA"
#   prev_val <- 0
#   toadd <- 0
#   for(i in seq(1, dim(wpp1)[1])){
#     print(i)
#     row <- mst_perc[i,]
#     if(row$from != milestone){
#       toadd <- toadd+1
#       milestone <- row$from
#     }
#
#     percentages <- c(percentages, row$percentage + toadd)
#     ids <- c(ids, row$cell_id)
#   }
#   names(percentages) <- ids
#   percentages
# }

get_waypoint_progression <- function(dataset, nr){
  cell_id <- character()
  froms <- character()
  tos <- character()
  percentage <- numeric()
  cumulative_percentages <- numeric()

  length_pieces <- dataset$milestone_network$length/sum(dataset$milestone_network$length)
  idx <- 1
  toadd <- 0
  for(piece in length_pieces){
    amount <- round(nr * piece)
    wps <- lapply(list(1:amount), function(x){x/amount})[[1]]
    cumulative_perc <- lapply(list(1:amount), function(x){(x/nr) + toadd})[[1]]
    ids <- lapply(list(1:amount), function(x){sprintf("w%d_%03d", idx, x)})[[1]]
    from <- rep(dataset$milestone_network$from[idx], amount)
    to <- rep(dataset$milestone_network$to[idx], amount)

    cell_id <- c(cell_id, ids)
    froms <- c(froms, from)
    tos <- c(tos, to)
    percentage <- c(percentage, wps)
    cumulative_percentages <- c(cumulative_percentages, cumulative_perc)

    idx <- idx + 1
    toadd <- toadd + piece

  }
  waypoint_progressions <- tibble(cell_id=cell_id, from=froms, to=tos, percentage=percentage, cumulative_percentages=cumulative_percentages)
  waypoint_progressions
}

get_geodesic_distances_from_progressions <- function(dataset, waypoint_progressions){
  waypoint_milestone_percentages <- dynwrap::convert_progressions_to_milestone_percentages(cell_id, dataset$milestone_ids, dataset$milestone_network, waypoint_progressions)
  waypoint_milestone_percentages <- dplyr::rename(waypoint_milestone_percentages, waypoint_id = cell_id)

  geodesic_distances <- dynwrap::calculate_geodesic_distances(dataset, waypoint_milestone_percentages = waypoint_milestone_percentages, directed=FALSE)
  geodesic_distances
}

interpolate_expression <- function(dataset, wps, gds, ws = 0.05){
  dc <- t(as.matrix(dataset$counts))
  weighted_expr <- do.call('cbind', lapply(seq_along(wps), function(idx){
    dist <- gds[idx,]
    weighted <- exp(-(dist^2)/(ws)^2)
    weighted <- weighted/sum(weighted)
    weighted_exp <- dc %*% weighted
    weighted_exp
  }))
  weighted_expr
}
