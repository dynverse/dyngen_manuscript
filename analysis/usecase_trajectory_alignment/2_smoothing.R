library(dyngen.manuscript)
library(tidyverse)
library(dtw)
library(ggbeeswarm)
library(viridis)

exp <- start_analysis("usecase_trajectory_alignment")

# For each pair, use 2 different kind of processing
design_smoothing <- exp$result("design_smooting.rds") %cache% {
  crossing(
    read_rds(exp$result("design_grouped.rds")),
    method = forcats::fct_inorder(c("DTW", "DTW+smoothing", "cellAlign"))
  )
}


#' design_smoothing %>% mutate(rn = row_number()) %>% dynutils::extract_row_to_list(1) %>% list2env(.GlobalEnv)
result_smoothing <- exp$result("results.rds") %cache% pmap_dfr(
  design_smoothing %>% mutate(rn = row_number()),
  function(rn, group, base1, base2, id, alpha, method,...) {

    datasets <- read_rds(exp$dataset_file(group))

    dataset1 <- datasets$dataset1 # read_rds(exp$dataset_file(id1))
    dataseti <- datasets$dataseti # read_rds(exp$dataset_file(id2))

    scores <- map2_dbl(id, seq_along(id), function(id_val, index){
      dataset2 <- dataseti[[index]]

      if(method == "cellAlign"){
          pt_expr1 <- get_cell_expression(dataset1, dataset1$milestone_network, "sA")
          interp1 <- cellAlign::interWeights(expDataBatch = t(pt_expr1$expression), trajCond = pt_expr1$pseudotime,
                                             winSz = 0.1, numPts = 200)
          # interscaled1 <- cellAlign::scaleInterpolate(interp1)

          pt_expr2 <- get_cell_expression(dataset2, dataset2$milestone_network, "sA")
          interp2 <- cellAlign::interWeights(expDataBatch = t(pt_expr2$expression), trajCond = pt_expr2$pseudotime,
                                             winSz = 0.1, numPts = 200)
          # interscaled2 <- cellAlign::scaleInterpolate(interp2)


          alignment <- cellAlign::globalAlign(interp1$interpolatedVals, interp2$interpolatedVals,
                                  scores = list(query = interp1$traj,
                                                ref = interp2$traj),
                                  sigCalc = F)

          pt1_aligned <- interp1$traj[alignment$align[[1]]$index1]
          pt2_aligned <- interp2$traj[alignment$align[[1]]$index2]

      } else {
        if (method == "DTW+smoothing") {
          res1 <- get_waypoint_expression(dataset1, 100, ws = 0.125)
          res2 <- get_waypoint_expression(dataset2, 100, ws = 0.125)
        } else{
          res1 <- get_cell_expression(dataset1, dataset1$milestone_network, "sA")
          res2 <- get_cell_expression(dataset2, dataset2$milestone_network, "sA")
        }

        pt1 <- res1$pseudotime
        pt2 <- res2$pseudotime

        expr1 <- res1$expression
        expr2 <- res2$expression

        dtw_alignment <- dtw(expr2, expr1, step.pattern = symmetric2, keep.internals = TRUE)
        # dtwPlotDensity(dtw_alignment)

        pt1_aligned <- pt1[dtw_alignment$index2]
        pt2_aligned <- pt2[dtw_alignment$index1]

      }

      score <- mean(abs(pt1_aligned - pt2_aligned))

      print(score)

      score
    })

    print(rn)
    design_smoothing[rn, ] %>% mutate(scores = paste(scores,collapse=" "))
  }
)
