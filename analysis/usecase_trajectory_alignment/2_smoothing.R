library(dyngen.manuscript)
library(tidyverse)
library(dtw)
library(ggbeeswarm)
library(viridis)

exp <- start_analysis("usecase_trajectory_alignment")

# For each pair, use 3 different kind of processing: smoothing, subsampling, original cells
design_smoothing <- exp$result("design_smoothing.rds") %cache% {
  read_rds(exp$result("design_datasets.rds")) %>%
    select("base_id1","base_id2", "id1", "id2", "noise") %>%
    mutate("smooth" = c(rep("smoothed", 25), rep("subsampled", 25), rep("original cells", 25))) %>%
    expand(nesting(base_id1, base_id2, id1, id2, noise), smooth)
}


alignment_results <- pmap(design_smoothing %>% mutate(rn = row_number()),
      function(base_id1, base_id2, id1, id2, smooth, rn, noise) {

        dataset1 <- read_rds(exp$dataset_file(id1))
        dataset2 <- read_rds(exp$dataset_file(id2))

        if(smooth == "smoothed"){
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

        if(smooth == "subsampled"){
          smp1 <- seq(from = 1, to = 1000, by = 10)
          smp2 <- seq(from = 1, to = 1000, by = 10)
          pt1_smp <- pt1[smp1]
          pt2_smp <- pt2[smp2]
          expr1_smp <- expr1[names(pt1_smp),]
          expr2_smp <- expr2[names(pt2_smp),]

          dtw_alignment <- dtw(expr2_smp, expr1_smp, step.pattern=symmetric2, keep.internals=T)
          dtwPlotAlignment(dtw_alignment)

          pt1_aligned_smp <- pt1_smp[dtw_alignment$index2]
          pt2_aligned_smp <- pt2_smp[dtw_alignment$index1]

          res <- sum(abs(pt1_aligned_smp - pt2_aligned_smp))
          cat(id1, "&", id2, "&", smooth, "=", res, "\n", sep=" ")

          return(res)

        } else {

          dtw_alignment <- dtw(expr2, expr1, step.pattern=symmetric2, keep.internals=T)
          dtwPlotAlignment(dtw_alignment)

          pt1_aligned <- pt1[dtw_alignment$index2]
          pt2_aligned <- pt2[dtw_alignment$index1]

          res <- sum(abs(pt1_aligned - pt2_aligned))
          cat(id1, "&", id2, "&", smooth, "=", res, "\n", sep=" ")
          return(res)
        }
      })

result_smoothing <- exp$result("result_smoothing.rds") %cache% {
  design_smoothing %>% mutate(result = as.numeric(alignment_results))
  }

