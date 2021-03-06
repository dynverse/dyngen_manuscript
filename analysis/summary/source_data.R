library(tidyverse)
library(dyngen.manuscript)

exp <- start_analysis("summary")

out <- list()

# Trajectory alignment
exp_a <- start_analysis("usecase_trajectory_alignment")
scores_a <- read_rds(exp_a$result("results.rds"))
stats_a <- read_rds(exp_a$result("statistics.rds"))

out$trajectory_alignment_summary <- scores_a %>% select(-mean_distance) %>% rename(dataset_id = id, method_id = method)
out$trajectory_alignment_statistics <- stats_a %>% select(-label)

# RNA velocity
exp_b <- start_analysis("usecase_rna_velocity")
scores_b <- read_rds(exp_b$result("scores.rds"))
stats_b <- read_rds(exp_b$result("statistics.rds"))

design_names <-
  tribble(
    ~method_id, ~params_id, ~method_label,
    "velocyto", "constant_velocity", "velocyto",
    "scvelo", "stochastic", "scvelo stochastic",
    "scvelo", "dynamical", "scvelo dynamical"
  )

out$RNA_velocity_perfeature <- scores_b$per_cell %>% inner_join(design_names, by = c("method_id", "params_id")) %>% select(-params_id, -method_id) %>% rename(method_id = method_label) %>% select(method_id, everything())
out$RNA_velocity_summary <- scores_b$summ %>% inner_join(design_names, by = c("method_id", "params_id")) %>% select(-params_id, -method_id, -mean_corr) %>% rename(method_id = method_label) %>% select(method_id, everything())
out$RNA_velocity_statistics <- stats_b %>% select(-label)

# Network inference
exp_c <- start_analysis("usecase_network_inference")
scores_c <- read_rds(exp_c$result("scores.rds"))
stats_c <- read_rds(exp_c$result("statistics.rds"))

out$CSNI_percell <- scores_c$aucs %>% filter(method == "casewise_casewise") %>% select(-method, -method_label, cni_method_id) %>% rename(method_id = cni_method_name)
out$CSNI_summary <- scores_c$summ %>% filter(method == "casewise_casewise") %>% select(-method, -method_label, cni_method_id) %>% rename(method_id = cni_method_name)
out$CSNI_statistics <- stats_c %>% select(-label)

# write output
openxlsx::write.xlsx(out, exp$result("Source Data.xlsx"))
