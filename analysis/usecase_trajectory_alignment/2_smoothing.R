library(dyngen.manuscript)
library(tidyverse)

exp <- start_analysis("usecase_trajectory_alignment")

# cross datasets and methods
design_smoothing <-
  crossing(
    read_rds(exp$result("design_datasets.rds")) %>% unnest(c(id, alpha)),
    method = factor(names(ta_methods), levels = c("DTW", "DTW+smoothing", "cellAlign"))
  )

#' design_smoothing %>% mutate(rn = row_number()) %>% dynutils::extract_row_to_list(1) %>% list2env(.GlobalEnv)
result_smoothing <- exp$result("results.rds") %cache% pmap_dfr(
  design_smoothing %>% mutate(rn = row_number()),
  function(rn, group, id, alpha, method, ...) {
    datasets <- read_rds(exp$dataset_file(id))

    dataset1 <- datasets$left
    dataset2 <- datasets$right

    method_fun <- ta_methods[[method]]
    out <- method_fun(dataset1, dataset2)
    pt1_aligned <- out$pt1_aligned
    pt2_aligned <- out$pt2_aligned

    distance <- mean(abs(pt1_aligned - pt2_aligned))

    design_smoothing[rn, ] %>% mutate(distance)
  }
)
