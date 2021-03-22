library(tidyverse)
library(dyngen.manuscript)

exp <- start_analysis("summary")

out <- list()

# Trajectory alignment
exp_a <- start_analysis("usecase_trajectory_alignment")
scores_a <- read_rds(exp_a$result("results.rds"))

out$trajectory_alignment_summary <- scores_a

# RNA velocity
exp_b <- start_analysis("usecase_rna_velocity")
scores_b <- read_rds(exp_b$result("scores.rds"))

out$RNA_velocity_percell <- scores_b$per_cell
out$RNA_velocity_summary <- scores_b$summ

# Network inference
exp_c <- start_analysis("usecase_network_inference")

scores_c <- read_rds(exp_c$result("scores.rds"))

out$CSNI_percell <- scores_c$aucs %>% filter(method == "casewise_casewise") %>% select(-method, -method_label)
out$CSNI_summary <- scores_c$summ %>% filter(method == "casewise_casewise") %>% select(-method, -method_label)

# write output
openxlsx::write.xlsx(out, exp$result("scores.xlsx"))
