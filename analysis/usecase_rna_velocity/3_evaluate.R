library(tidyverse)
library(dyngen.manuscript)

exp <- start_analysis("usecase_rna_velocity")

design_datasets <- read_rds(exp$result("design_datasets.rds"))
design_velocity <- read_rds(exp$result("design_velocity.rds"))

#' @examples
#' design_velocity %>% dynutils::extract_row_to_list(46) %>% list2env(.GlobalEnv)

# Calculate scores -----
scores <- exp$result("scores_individual.rds") %cache% {
  pmap_dfr(
    design_velocity,
    function(method_id, params_id, dataset_id, ...) {
      dataset <- read_rds(exp$dataset_file(dataset_id))
      groundtruth_velocity <- dataset$propensity_ratios

      velocity <- read_rds(exp$velocity_file(dataset_id, method_id, params_id))

      velocity_differences <- velocity$expression_future - dataset$expression
      velocity_differences[velocity_differences == 0] <- NA
      velocity_differences[is.na(velocity_differences)] <- runif(sum(is.na(velocity_differences)), -1e-10, 1e-10)

      paired_simil(
        velocity_differences,
        groundtruth_velocity,
        method = "spearman",
        margin = 2
      ) %>%
        enframe("feature_id", "score") %>%
        mutate(method_id, params_id, dataset_id)
    }
  ) %>% left_join(design_datasets, c("dataset_id" = "id"))
}

mean_scores <- exp$result("scores_aggregated.rds") %cache% {
  scores %>%
    group_by(dataset_id, method_id, params_id) %>%
    summarise(score = mean(score)) %>%
    left_join(design_datasets, c("dataset_id" = "id"))
}

