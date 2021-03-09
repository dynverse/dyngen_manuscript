#' @importFrom cellAlign interWeights globalAlign
#' @importFrom dtw dtw symmetric2
#'
#' @export
ta_methods <- list(
  "cellAlign" = function(dataset1, dataset2) {

    pt_expr1 <- get_cell_expression(dataset1, dataset1$milestone_network, "sA")
    pt_expr2 <- get_cell_expression(dataset2, dataset2$milestone_network, "sA")

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

    pt1_aligned <- interp1$traj[alignment$align[[1]]$index1]
    pt2_aligned <- interp2$traj[alignment$align[[1]]$index2]

    lst(pt1_aligned, pt2_aligned)
  },
  "DTW" = function(dataset1, dataset2) {
    res1 <- get_cell_expression(dataset1, dataset1$milestone_network, "sA")
    res2 <- get_cell_expression(dataset2, dataset2$milestone_network, "sA")

    dtw_alignment <- dtw::dtw(res2$expression, res1$expression, step.pattern = dtw::symmetric2, keep.internals = TRUE)
    # dtwPlotDensity(dtw_alignment)

    pt1_aligned <- res1$pseudotime[dtw_alignment$index2]
    pt2_aligned <- res2$pseudotime[dtw_alignment$index1]

    lst(pt1_aligned, pt2_aligned)
  },
  "DTW+smoothing" = function(dataset1, dataset2) {
    res1 <- get_waypoint_expression(dataset1, 100, ws = 0.125)
    res2 <- get_waypoint_expression(dataset2, 100, ws = 0.125)

    dtw_alignment <- dtw::dtw(res2$expression, res1$expression, step.pattern = dtw::symmetric2, keep.internals = TRUE)
    # dtwPlotDensity(dtw_alignment)

    pt1_aligned <- res1$pseudotime[dtw_alignment$index2]
    pt2_aligned <- res2$pseudotime[dtw_alignment$index1]

    lst(pt1_aligned, pt2_aligned)
  }
)
