library(dyngen.manuscript)
library(tidyverse)
library(ggbeeswarm)
library(viridis)

exp <- start_analysis("usecase_trajectory_alignment")

datasets <- read_rds(exp$dataset_file("test"))

# cross datasets and methods
design_smoothing <- exp$result("design_smooting.rds") %cache% {
  crossing(
    read_rds(exp$result("design_grouped.rds")),
    method = forcats::fct_inorder(names(ta_methods))
  )
}

#' design_smoothing %>% mutate(rn = row_number()) %>% dynutils::extract_row_to_list(1) %>% list2env(.GlobalEnv)
result_smoothing <- exp$result("results.rds") %cache% pmap_dfr(
  design_smoothing %>% mutate(rn = row_number()),
  function(rn, group, base1, base2, id, alpha, method,...) {

    datasets <- read_rds(exp$dataset_file(group))

    dataset1 <- datasets$dataset1 # read_rds(exp$dataset_file(id1))
    dataset2 <- datasets$dataset2 # read_rds(exp$dataset_file(id2))

    scores <- map2_dbl(id, seq_along(id), function(id_val, index){
      dataset2 <- dataseti[[index]]

      method_fun <- ta_methods[[method]]
      out <- method_fun(dataset1, dataset2)
      pt1_aligned <- out$pt1_aligned
      pt2_aligned <- out$pt2_aligned

      score <- mean(abs(pt1_aligned - pt2_aligned))

      print(score)

      score
    })

    print(rn)
    design_smoothing[rn, ] %>% mutate(scores = paste(scores,collapse=" "))
  }
)
