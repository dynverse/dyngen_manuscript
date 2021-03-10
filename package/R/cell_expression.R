#' Get cell expression, but sorted according to pseudotime
#'
#' @param dataset A dyngen dataset
#'
#' @export
get_cell_expression <- function(dataset) {
  pseudotime <- sort(compute_pseudotime_from_root(dataset))
  list(
    pseudotime = pseudotime,
    counts = as.matrix(dataset$counts)[names(pseudotime),],
    expression = as.matrix(dataset$expression)[names(pseudotime),]
  )
}

#' Get corrected pseudotime
#'
#' @param dataset A dyngen dataset
#' @param root The root milestone to compute the pseudotime from.
#'
#' @importFrom gtools mixedorder
#'
#' @export
compute_pseudotime_from_root <- function(
  dataset,
  root = unique(setdiff(dataset$milestone_network$from, dataset$milestone_network$to)),
  normalized = TRUE
) {
  root_pct <- tibble(
    waypoint_id = "wp",
    milestone_id = root,
    percentage = 1
  )

  pseudotime <- dynwrap::calculate_geodesic_distances(
    dataset,
    waypoint_milestone_percentages = root_pct
  )[1,]

  if (normalized) {
    pseudotime <- pseudotime / sum(dataset$milestone_network$length)
  }

  pseudotime
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

