library(dyngen.manuscript)
library(tidyverse)

exp <- start_analysis("usecase_trajectory_alignment")

# cross datasets and methods
design_smoothing <-
  crossing(
    read_rds(exp$result("design_datasets.rds")),
    method = forcats::fct_inorder(c("DTW", "cellAlign"))
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
      pt1_aligned <- out$pt1_aligned
      pt2_aligned <- out$pt2_aligned

      distance <- mean(abs(pt1_aligned - pt2_aligned))

      aupt <- pracma::trapz(
        pt1_aligned + pt2_aligned,
        abs(pt1_aligned - pt2_aligned)
      )

      score <- 1 - aupt

      design_smoothing[rn, ] %>% mutate(
        distance,
        aupt,
        score
      )
    }
  }
)

ggstatsplot::ggwithinstats(
  results,
  x = method,
  y = score,
  type = "np",
  pairwise.display = "all"
)

z <- ggstatsplot::ggwithinstats(
  results,
  x = method,
  y = distance,
  type = "np",
  pairwise.display = "all"
)
z$labels$subtitle
