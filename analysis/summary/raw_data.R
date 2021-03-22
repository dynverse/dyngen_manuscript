library(tidyverse)
library(dyngen.manuscript)

exp <- start_analysis("summary")

out <- list()

# Network inference
exp_c <- start_analysis("usecase_network_inference")

exp_c_scores <- read_rds(exp_c$result("scores.rds"))

out$CSNI_raw <- exp_c_scores$aucs %>% filter(method == "casewise_casewise") %>% select(-method, -method_label)
out$CSNI_summary <- exp_c_scores$summ %>% filter(method == "casewise_casewise") %>% select(-method, -method_label)

# write output
openxlsx::write.xlsx(out, exp$result("scores.xlsx"))
