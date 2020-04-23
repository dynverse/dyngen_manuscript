#' Get cell expression, but sorted according to pseudotime
#'
#' @param dataset A dyngen dataset
#'
#' @export
get_cell_expression <- function(dataset){
  pseudotime <- sort(calculate_correct_pseudotime(dataset))
  nor2 <- names(pseudotime)
  volg2 <- sapply(nor2, function(i) strtoi(substr(i, 5, nchar(i))))
  cnt2 <- as.matrix(dataset$counts)
  cnts2 <- cnt2[volg2,]

  return(list(pseudotime = pseudotime, expression = cnts2))
}

#' Get corrected pseudotime
#'
#' @param dataset A dyngen dataset
#'
#' @importFrom gtools mixedorder
#'
#' @export
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
