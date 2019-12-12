get_wexpr <- function(dataset, amount){
  wpp1 <- get_waypoint_progression(dataset, amount)
  gd1 <- get_geodesic_distances_from_progressions(dataset, wpp1)
  wexpr1 <- interpolate_expression(dataset, wpp1$percentage, gd1)
  t(wexpr1)
}

get_cell_expression <- function(dataset){
  nor2 <- names(sort(calculate_correct_pseudotime(dataset)))
  volg2 <- sapply(nor2, function(i) strtoi(substr(i, 5, 7)))
  cnt2 <- as.matrix(dataset$counts)
  cnts2 <- cnt2[volg2,]
  cnts2
}
