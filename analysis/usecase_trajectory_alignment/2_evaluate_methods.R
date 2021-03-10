library(dyngen.manuscript)
library(tidyverse)

exp <- start_analysis("usecase_trajectory_alignment")

# cross datasets and methods
design_smoothing <-
  crossing(
    read_rds(exp$result("design_datasets.rds")),
    method = forcats::fct_inorder(c("DTW", "cellAlign"))
  )

filter_cells <- function(dataset, ...) {
  df <-
    dataset$cell_info %>%
    left_join(dataset$progressions, by = "cell_id")

  df <- df %>% filter(...)

  cids <- df$cell_id

  dataset$cell_ids <- cids
  dataset$cell_info <- dataset$cell_info %>% slice(match(cids, .data$cell_id))
  dataset$counts_spliced <- dataset$counts_spliced[cids, , drop = FALSE]
  dataset$counts_protein <- dataset$counts_protein[cids, , drop = FALSE]
  dataset$counts_unspliced <- dataset$counts_unspliced[cids, , drop = FALSE]
  dataset$counts <- dataset$counts[cids, , drop = FALSE]
  dataset$expression <- dataset$expression[cids, , drop = FALSE]
  dataset$milestone_percentages <- dataset$milestone_percentages %>% filter(.data$cell_id %in% cids)
  dataset$progressions <- dataset$progressions %>% slice(match(cids, .data$cell_id))
  dataset$milestone_network <- dataset$milestone_network %>%
    inner_join(dataset$progressions %>% select(from, to) %>% unique, by = c("from", "to"))

  dataset
}

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

      design_smoothing[rn, ] %>% mutate(distance)
    }
  }
)

ggstatsplot::ggwithinstats(
  results %>% filter(method != "DTW+smoothing"),
  x = method,
  y = distance,
  type = "np",
  pairwise.display = "all"
)
