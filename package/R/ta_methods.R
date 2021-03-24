#' @importFrom cellAlign interWeights globalAlign
#' @importFrom dtw dtw symmetric2
#'
#' @export
ta_methods <- list(
  "cellAlign" = function(dataset1, dataset2) {
    pt_expr1 <- get_cell_expression(dataset1)
    pt_expr2 <- get_cell_expression(dataset2)

    interp1 <- cellAlign::interWeights(
      expDataBatch = t(pt_expr1$expression),
      trajCond = pt_expr1$pseudotime,
      winSz = 0.1,
      numPts = 200
    )

    interp2 <- cellAlign::interWeights(
      expDataBatch = t(pt_expr2$expression),
      trajCond = pt_expr2$pseudotime,
      winSz = 0.1,
      numPts = 200
    )

    alignment <- cellAlign::globalAlign(
      interp1$interpolatedVals,
      interp2$interpolatedVals,
      scores = list(query = interp1$traj, ref = interp2$traj),
      sigCalc = FALSE
    )

    align <- alignment$align[[1]]
    align$pt1_aligned <- interp1$traj[align$index1]
    align$pt2_aligned <- interp2$traj[align$index2]
    align$costMatrix <- alignment$costMatrix

    align
  },
  "DTW" = function(dataset1, dataset2) {
    res1 <- get_cell_expression(dataset1)
    res2 <- get_cell_expression(dataset2)

    dtw_alignment <- dtw::dtw(
      res2$expression,
      res1$expression,
      step.pattern = dtw::symmetric2,
      keep.internals = TRUE
    )

    align <- dtw_alignment[c("index1", "index2", "index1s", "index2s", "stepsTaken", "costMatrix")]
    align$pt1_aligned <- res1$pseudotime[align$index2]
    align$pt2_aligned <- res2$pseudotime[align$index1]

    align
  },
  "DTW+smoothing" = function(dataset1, dataset2) {
    res1 <- get_waypoint_expression(dataset1, 100, ws = 0.125)
    res2 <- get_waypoint_expression(dataset2, 100, ws = 0.125)

    dtw_alignment <- dtw::dtw(res2$expression, res1$expression, step.pattern = dtw::symmetric2, keep.internals = TRUE)
    # dtwPlotDensity(dtw_alignment)

    align <- dtw_alignment[c("index1", "index2", "index1s", "index2s", "stepsTaken", "costMatrix")]
    align$pt1_aligned <- res1$pseudotime[align$index2]
    align$pt2_aligned <- res2$pseudotime[align$index1]

    align
  }
)
