library(tidyverse)
library(dyngen.manuscript)

exp <- start_analysis("usecase_rna_velocity")

design_datasets <- read_rds(exp$result("design_datasets.rds"))
design_velocity <- read_rds(exp$result("design_velocity.rds"))

#' @examples
#' design_velocity %>% dynutils::extract_row_to_list(16) %>% list2env(.GlobalEnv)
# file.remove(exp$result("scores.rds"))

# Calculate scores -----
scores <- exp$result("scores.rds") %cache% {
  outs <- pmap(
    design_velocity %>% mutate(i = row_number()),
    function(method_id, params_id, dataset_id, i, ...) {
      cat(i, "/", nrow(design_velocity), ": ", dataset_id, "\n", sep = "")
      if (!file.exists(exp$dataset_file(dataset_id))) return(NULL)
      dataset <- read_rds(exp$dataset_file(dataset_id))
      groundtruth_velocity <- dataset$rna_velocity

      velocity_file <- exp$velocity_file(dataset_id, method_id, params_id)
      if (!file.exists(velocity_file)) return(NULL)
      velocity <- read_rds(velocity_file)
      predicted_velocity <- velocity$velocity_vector

      per_cell <- tibble(
        method_id, params_id, dataset_id,
        feature_id = dataset$feature_ids,
        cor = paired_simil(predicted_velocity, groundtruth_velocity, method = "spearman", margin = 2)[feature_id],
      ) %>%
        left_join(design_datasets, c("dataset_id" = "id"))

      summ <- tibble(
        method_id, params_id, dataset_id,
        cor = cor(as.vector(predicted_velocity), as.vector(groundtruth_velocity)),
        mean_corr = mean(per_cell$cor),
      ) %>% left_join(design_datasets, c("dataset_id" = "id"))

      lst(per_cell, summ)
    }
  )
  list(
    per_cell = map_dfr(outs, "per_cell"),
    summ = map_dfr(outs, "summ")
  )
}
