library(tidyverse)
library(dyngen.manuscript)

exp <- start_analysis("usecase_rna_velocity_net")

design_datasets <- read_rds(exp$result("design_datasets.rds"))
design_velocity <- read_rds(exp$result("design_velocity.rds"))

#' @examples
#' design_velocity %>% dynutils::extract_row_to_list(16) %>% list2env(.GlobalEnv)
# file.remove(exp$result("scores_individual.rds"))
# file.remove(exp$result("scores_aggregated.rds"))

# Calculate scores -----
scores <- exp$result("scores_individual.rds") %cache% {
  pmap_dfr(
    design_velocity,
    function(method_id, params_id, dataset_id, ...) {
      dataset <- read_rds(exp$dataset_file(dataset_id))
      groundtruth_velocity <- dataset$log_propensity_ratios
      groundtruth_velocity[groundtruth_velocity == 0] <- runif(sum(groundtruth_velocity == 0), -1e-10, 1e-10)
      groundtruth_velocity1 <- dataset$net_propensity
      groundtruth_velocity1[groundtruth_velocity1 == 0] <- runif(sum(groundtruth_velocity1 == 0), -1e-10, 1e-10)

      velocity_file <- exp$velocity_file(dataset_id, method_id, params_id)
      if (!file.exists(velocity_file)) {
        return(NULL)
      }
      velocity <- read_rds(velocity_file)

      velocity_differences <- velocity$expression_future - dataset$expression
      velocity_differences[velocity_differences == 0] <- NA
      velocity_differences[is.na(velocity_differences)] <- runif(sum(is.na(velocity_differences)), -1e-10, 1e-10)

      # paired_simil(
      #   velocity_differences,
      #   groundtruth_velocity,
      #   method = "spearman",
      #   margin = 2
      # ) %>%
      #   enframe("feature_id", "score") %>%
      #   mutate(method_id, params_id, dataset_id)

      score <- paired_simil(
        velocity_differences,
        groundtruth_velocity,
        method = "spearman",
        margin = 2
      ) %>%
        enframe("feature_id", "score")
      score_new <- paired_simil(
        velocity_differences,
        groundtruth_velocity1,
        method = "spearman",
        margin = 2
      ) %>%
        enframe("feature_id", "score_new")
      inner_join(score, score_new, by = c("feature_id")) %>%
        mutate(method_id, params_id, dataset_id)
    }
  ) %>% left_join(design_datasets, c("dataset_id" = "id"))
}

mean_scores <- exp$result("scores_aggregated.rds") %cache% {
  scores %>%
    group_by(dataset_id, method_id, params_id) %>%
    summarise(score = mean(score), score_new = mean(score_new)) %>%
    left_join(design_datasets, c("dataset_id" = "id")) %>%
    ungroup()
}

ggplot(mean_scores) + geom_point(aes(score, score_new, colour = params_id))

