library(tidyverse)
library(dynwrap)
library(dyneval)
library(dyngenanalysis)

dataset <- read_rds("derived_files/network_inference/small_bifurcating/dataset.rds")

expression <- dataset$expression
regulators <- dataset$regulators
targets <- dataset$targets
parameters <- list(num_int_per_cell = 10000L)


dataset$prior_information$regulators <- dataset$regulators
dataset$prior_information$targets <- dataset$targets
#model <- infer_trajectory(dataset, cni_deltacor(), verbose = TRUE)
model <- infer_trajectory(dataset, cni_pyscenic_sgbm(), verbose = TRUE)
eval <- calculate_metrics(dataset, model, metrics = list(auc = cni_auc))
eval
ggplot(eval$evals[[1]]) + geom_point(aes(auroc, aupr, colour = method))
