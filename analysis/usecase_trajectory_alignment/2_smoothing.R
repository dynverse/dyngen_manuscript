library(dyngen.manuscript)
library(tidyverse)
library(dtw)
library(ggbeeswarm)
library(viridis)

exp <- start_analysis("usecase_trajectory_alignment")

# For each pair, use 3 different kind of processing
design_smoothing <- exp$result("design_methods.rds") %cache% {
  crossing(
    read_rds(exp$result("design_datasets.rds")),
    method = forcats::fct_inorder(c("DTW", "DTW+smoothing"))
  )
}

#' design_smoothing %>% mutate(rn = row_number()) %>% dynutils::extract_row_to_list(1) %>% list2env(.GlobalEnv)
result_smoothing <- exp$result("results.rds") %cache% pmap_dfr(
  design_smoothing %>% mutate(rn = row_number()),
  function(base_id1, base_id2, id1, id2, method, rn, noise, ...) {

    dataset1 <- read_rds(exp$dataset_file(id1))
    dataset2 <- read_rds(exp$dataset_file(id2))

    if (method == "DTW+smoothing") {
      res1 <- get_waypoint_expression(dataset1, 100)
      res2 <- get_waypoint_expression(dataset2, 100)
    } else {
      res1 <- get_cell_expression(dataset1, dataset1$milestone_network, "sA")
      res2 <- get_cell_expression(dataset2, dataset2$milestone_network, "sA")
    }

    pt1 <- res1$pseudotime
    pt2 <- res2$pseudotime

    expr1 <- res1$expression
    expr2 <- res2$expression

    dtw_alignment <- dtw(expr2, expr1, step.pattern = symmetric2, keep.internals = TRUE)

    pt1_aligned <- pt1[dtw_alignment$index2]
    pt2_aligned <- pt2[dtw_alignment$index1]

    score <- mean(abs(pt1_aligned - pt2_aligned))

    design_smoothing[rn, ] %>% mutate(score)
  }
)
