#' Get cell expression, but sorted according to pseudotime
#'
#' @param dataset A dyngen dataset
#'
#' @export
get_cell_expression <- function(dataset, milestone_network, start){
  pseudotime <- sort(calculate_correct_pseudotime(dataset, start))
  nor2 <- names(pseudotime)
  volg2 <- sapply(nor2, function(i) strtoi(substr(i, 5, nchar(i))))
  cnt2 <- as.matrix(dataset$counts)
  cnts2 <- cnt2[volg2,]
  return(list(pseudotime = pseudotime, expression = cnts2))
}

#   function(dataset){
#   pseudotime <- sort(calculate_correct_pseudotime(dataset))
#   nor2 <- names(pseudotime)
#   volg2 <- sapply(nor2, function(i) strtoi(substr(i, 5, nchar(i))))
#   cnt2 <- as.matrix(dataset$counts)
#   cnts2 <- cnt2[volg2,]
#
#   return(list(pseudotime = pseudotime, expression = cnts2))
# }

#' Get corrected pseudotime
#'
#' @param dataset A dyngen dataset
#' @param start The start milestone
#'
#' @importFrom gtools mixedorder
#'
#' @export
calculate_correct_pseudotime <- function(dataset, start, normalized = TRUE){
  if (normalized) {
    dataset$milestone_network <-
      dataset$milestone_network %>%
      mutate(length = length / sum(length))
  }
  pct <- dynwrap::calculate_geodesic_distances(
    dataset,
    waypoint_milestone_percentages = tibble(
      waypoint_id = "wp",
      milestone_id = start,
      percentage = 1
    )
  )[1,]
}

#' Get uncorrected pseudotime
#'
#' @param dataset A dyngen dataset
#'
#' @importFrom gtools mixedorder
#'
#' @export
calculate_pseudotime <- function(dataset, milestone_network, start){
  length_pieces <- milestone_network$length/sum(milestone_network$length)
  nr_piece <- 0
  traj <- numeric()
  names_traj <- list()
  # TODO convert to for for each piece
  while(start %in% milestone_network[["from"]]){
    idx <- which(start == milestone_network[["from"]])[[1]]
    end <- milestone_network[["to"]][idx]
    to_add <- sum(length_pieces[0:nr_piece])

    perc <-(dataset[["progressions"]] %>% filter(from == start & to == end))$percentage + to_add
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

