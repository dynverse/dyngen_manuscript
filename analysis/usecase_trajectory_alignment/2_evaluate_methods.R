library(dyngen.manuscript)
library(tidyverse)

exp <- start_analysis("usecase_trajectory_alignment")

# cross datasets and methods
design_smoothing <-
  crossing(
    read_rds(exp$result("design_datasets.rds")),
    method = c("DTW", "cellAlign")
  )

#' design_smoothing %>% mutate(rn = row_number()) %>% dynutils::extract_row_to_list(1) %>% list2env(.GlobalEnv)
results <- exp$result("results.rds") %cache% pmap_dfr(
  design_smoothing %>% mutate(rn = row_number()),
  function(rn, group, id, method, ...) {
    file <- exp$dataset_file(id)
    if (file.exists(file)) {
      dataset <- read_rds(file)

      method_fun <- ta_methods[[method]]

      dataset1 <- dataset %>% filter_cells(model == "left")
      dataset2 <- dataset %>% filter_cells(model == "right")

      out <- method_fun(dataset1, dataset2)
      pt1 <- out$pt1_aligned
      pt2 <- out$pt2_aligned

      distance <- mean(abs(pt1 - pt2))

      abwap <- ta_abwap(pt1, pt2)

      design_smoothing[rn, ] %>% mutate(
        mean_distance = distance,
        abwap
      )
    }
  }
)

ggstatsplot::ggwithinstats(
  results,
  x = method,
  y = abwap,
  type = "np",
  pairwise.display = "all"
)
