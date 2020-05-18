library(dyngen.manuscript)
library(tidyverse)
library(dtw)
library(ggbeeswarm)
library(viridis)

exp <- start_analysis("usecase_trajectory_alignment")

# For each pair, use 3 different kind of processing
design_smoothing <- exp$result("design_smoothing.rds") %cache% {
  crossing(
    read_rds(exp$result("design_datasets.rds")),
    smooth = c("smoothed", "subsampled", "original cells")
  )
}

#' design_smoothing %>% mutate(rn = row_number()) %>% dynutils::extract_row_to_list(1) %>% list2env(.GlobalEnv)
result_smoothing <- exp$result("result_smoothing.rds") %cache% pmap_dfr(
  design_smoothing %>% mutate(rn = row_number()),
  function(base_id1, base_id2, id1, id2, smooth, rn, noise, ...) {

    dataset1 <- read_rds(exp$dataset_file(id1))
    dataset2 <- read_rds(exp$dataset_file(id2))

    if (smooth == "smoothed") {
      res1 <- get_waypoint_expression(dataset1, 100)
      res2 <- get_waypoint_expression(dataset2, 100)
    } else {
      res1 <- get_cell_expression(dataset1)
      res2 <- get_cell_expression(dataset2)
    }

    pt1 <- res1$pseudotime
    pt2 <- res2$pseudotime

    expr1 <- res1$expression
    expr2 <- res2$expression

    if (smooth == "subsampled") {
      smp1 <- seq(from = 1, to = nrow(expr1), by = 10)
      smp2 <- seq(from = 1, to = nrow(expr2), by = 10)
      pt1 <- pt1[smp1]
      pt2 <- pt2[smp2]
      expr1 <- expr1[names(pt1),]
      expr2 <- expr2[names(pt2),]
    }

    dtw_alignment <- dtw(expr2, expr1, step.pattern = symmetric2, keep.internals = TRUE)

    pt1_aligned <- pt1[dtw_alignment$index2]
    pt2_aligned <- pt2[dtw_alignment$index1]

    score <- sum(abs(pt1_aligned - pt2_aligned))

    design_smoothing[rn, ] %>% mutate(score)
  }
)
