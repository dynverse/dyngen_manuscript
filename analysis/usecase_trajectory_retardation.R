library(tidyverse)
library(dyngen)
library(dynplot)
library(gtools)
library(dtw)

setwd("/home/louisedc/Work/dyngen_analysis/analysis")
source("dtw_functions.R")

# # Function that will change the speed of the targeted gene.
# # TODO: what is wpr & xdr again?
# change_speed <- function(model, target, rate) {
#   param_id <- which(target == model$feature_info$feature_id)[[1]]
#
#   model$feature_info$wpr[param_id] <- model$feature_info$wpr[param_id] * rate
#   model$feature_info$xdr[param_id] <- model$feature_info$xdr[param_id] * rate
#
#   param_wpr <- paste0("wpr_", target)
#   param_xdr <- paste0("xdr_", target)
#
#   model$simulation_system$parameters[param_wpr] <- model$simulation_system$parameters[param_wpr] * rate
#   model$simulation_system$parameters[param_xdr] <- model$simulation_system$parameters[param_xdr] * rate
#
#   model
# }
#
# generate_retardation <- function(model, id){
#   model1 <- model
#   model1 <- change_speed(model1, "A4_TF1", 0.5)
#   model1 <- generate_gold_standard(model1)
#   model1 <- generate_cells(model1)
#   model1 <- generate_experiment(model1)
#   dataset2 <- wrap_dataset(model1)
#
#   # save model & dataset
#   # saveRDS(dataset2, paste0("data/retardation_ds_", id))
#   # saveRDS(model1, paste0("analysis/data/retardation_model_", id))
#
#   # for(i in 1:9){
#   #   generate_replication(model1, i)
#   # }
#
#   dataset2
# }
#
# generate_replication <- function(model, id){
#   model <- generate_cells(model)
#   model <- generate_experiment(model)
#   saveRDS(model, paste0("analysis/data/retardation_replication_", id))
# }
#
# # Construct a synthetic model
# # A very simple linear model
# backbone <- bblego(
#   bblego_start("A", type = "simple", num_modules = 4),
#   bblego_linear("A", "B", type = "simple", num_modules = 5),
#   bblego_linear("B", "C", type = "simple", num_modules = 6),
#   bblego_end("C", type = "simple", num_modules = 5)
# )
#
# # Simpler model to check the outline of the backbone
# model <- initialise_model(
#     backbone = backbone,
#     num_cells = 1000,
#     num_tfs = 30,
#     num_targets = 500,
#     num_hks = 500,
#     simulation_params = simulation_default(num_simulations=5)
#   )
#
# model <- generate_tf_network(model)
# model <- generate_feature_network(model)
# model <- generate_kinetics(model)
# model <- generate_gold_standard(model)
# model <- generate_cells(model)
# model <- generate_experiment(model)
# dataset1 <- wrap_dataset(model)
#
# model1 <- model
# # TODO: check if "A4_TF1" is in model$feature_network$from, else choose another tf?
# model1 <- change_speed(model1, "A5_TF1", 0.25)
# model1 <- generate_gold_standard(model1)
# model1 <- generate_cells(model1)
# model1 <- generate_experiment(model1)
# dataset2 <- wrap_dataset(model1)

# saveRDS(model, paste0("analysis/data/orig_model_", id))

# dataset2 <- generate_retardation(model, 1)

# read in orig model, retardation & 9 replicates?

# orig: A5_TF1 delayed
# 2: B6_TF1 op 0
# ctrl: nothing happened

am1 <- get_alignment("data/orig_model", "data/retardation_model_1")
am2 <- get_alignment("data/2_orig_model", "data/2_retardation_model_1")
am3 <- get_alignment("data/3_orig_model", "data/3_retardation_model_1")
am4 <- get_alignment("data/4_orig_model", "data/4_retardation_model_1")
am5 <- get_alignment("data/5_orig_model", "data/5_retardation_model_1")
am6 <- get_alignment("data/6_orig_model", "data/6_retardation_model_1")
am7 <- get_alignment("data/7_orig_model", "data/7_retardation_model_1")
am7 <- get_alignment("data/7_orig_model", "data/7_retardation_model_1")
am8 <- get_alignment("data/8_orig_model", "data/8_retardation_model_1")

am7b <- get_alignment2("data/7_orig_model", "data/7_retardation_model_1")


dtwPlotDensity(am7b)
abline(coef = c(0,1))


get_alignment <- function(file1, file2){
  ref_model <- wrap_dataset(readRDS(file1))
  ret_model <- wrap_dataset(readRDS(file2))

  traj1 <- calculate_correct_pseudotime(ref_model)
  traj2 <- calculate_correct_pseudotime(ret_model)

  wpp1 <- get_waypoint_progression(ref_model, 100)
  wpp2 <- get_waypoint_progression(ret_model, 100)

  gd1 <- get_geodesic_distances_from_progressions(ref_model, wpp1)
  gd2 <- get_geodesic_distances_from_progressions(ret_model, wpp2)

  wexpr1 <- interpolate_expression(ref_model, wpp1$percentage, gd1)
  wexpr2 <- interpolate_expression(ret_model, wpp2$percentage, gd2)

  pr1 <- to_python_readable(ref_model, gd1, wpp1, traj1, wpp1)
  pr2 <- to_python_readable(ret_model, gd2, wpp2, traj2, wpp2)

  saveRDS(pr1, "../../trajectory_alignment/data/python_readable_big9.rds")
  saveRDS(pr2, "../../trajectory_alignment/data/python_readable_big10.rds")


  am1 <- calculate_alignment(wexpr1, wexpr2)
  am1
}

get_alignment2 <- function(file1, file2){
  ref_model <- wrap_dataset(readRDS(file1))
  ret_model <- wrap_dataset(readRDS(file2))

  traj1 <- calculate_correct_pseudotime(ref_model)
  traj2 <- calculate_correct_pseudotime(ret_model)

  wpp1 <- get_waypoint_progression(ref_model, 100)
  wpp2 <- get_waypoint_progression(ret_model, 100)

  gd1 <- get_geodesic_distances_from_progressions(ref_model, wpp1)
  gd2 <- get_geodesic_distances_from_progressions(ret_model, wpp2)

  wexpr1 <- interpolate_expression(ref_model, wpp1$percentage, gd1)
  wexpr2 <- interpolate_expression(ret_model, wpp2$percentage, gd2)

  wexpr1 <- wexpr1[1:30,]
  wexpr2 <- wexpr2[1:30,]

  pr1 <- to_python_readable(ref_model, gd1, wpp1, traj1, wpp1)
  pr2 <- to_python_readable(ret_model, gd2, wpp2, traj2, wpp2)

  saveRDS(pr1, "../../trajectory_alignment/data/python_readable_big11.rds")
  saveRDS(pr2, "../../trajectory_alignment/data/python_readable_big12.rds")


  am1 <- calculate_alignment(wexpr1, wexpr2)
  am1
}

calculate_alignment <- function(wexpr1, wexpr2){
  d1 <- dist(as.matrix(t(wexpr1)), as.matrix(t(wexpr2)), dist.method="Euclidean")
  d1 <- (d1 - min(d1)) / (max(d1) - min(d1))
  d1 <- d1/2
  am <- dtw(d1, step.pattern = symmetric2, keep.internals = TRUE, open.end = TRUE)
  am
}



d2 <- dist(as.matrix(t(wexpr1)), as.matrix(t(wexpr2)), dist.method="Euclidean")
d2 <- (d2 - min(d2)) / (max(d2) - min(d2))
d2 <- d2/2

am2 <- dtw(d2, step.pattern = symmetric2, keep.internals = TRUE)
dtwPlotDensity(am2)
abline(coef = c(0,1))

d3 <- dist(as.matrix(t(wexpr1)), as.matrix(t(wexpr2)), dist.method="Euclidean")
d3 <- (d3 - min(d2)) / (max(d2) - min(d2))
d3 <- d2/2

am2 <- dtw(d2, step.pattern = symmetric2, keep.internals = TRUE)
dtwPlotDensity(am2)
abline(coef = c(0,1))

rep_models <- list()
wexprs <- list()
dists <- list()
cum_dist <- d1
for(i in 1:9){
  rep_model <- wrap_dataset(readRDS(paste0("data/2_retardation_replication_", i)))
  rep_models[[i]] <- rep_model
  wpp <- get_waypoint_progression(rep_model, 100)
  gd <- get_geodesic_distances_from_progressions(rep_model, wpp)
  wexpr <- interpolate_expression(rep_model, wpp$percentage, gd)
  wexprs[[i]] <- wexpr
  d <- dist(as.matrix(t(wexpr1)), as.matrix(t(wexpr)), dist.method="Euclidean")
  d <- (d - min(d)) / (max(d) - min(d))
  dists[[i]] <- d
  cum_dist <- cum_dist + d
}

d <- cum_dist / 10


am <- dtw(d, step.pattern = symmetric2, keep.internals = TRUE)
dtwPlotDensity(am)
abline(coef = c(0,1))

traj1 <- calculate_correct_pseudotime(dataset1)
traj2 <- calculate_correct_pseudotime(dataset2)

wpp1 <- get_waypoint_progression(dataset1, 100)
wpp2 <- get_waypoint_progression(dataset2, 100)
gd1 <- get_geodesic_distances_from_progressions(dataset1, wpp1)
gd2 <- get_geodesic_distances_from_progressions(dataset2, wpp2)

wexpr1 <- interpolate_expression(dataset1, wpp1$percentage, gd1)
wexpr2 <- interpolate_expression(dataset2, wpp2$percentage, gd2)

# calculate own costs & normalize
d <- dist(as.matrix(t(wexpr1)), as.matrix(t(wexpr2)), dist.method="Euclidean")
d <- (d - min(d)) / (max(d) - min(d))
am <- dtw(d, step.pattern = symmetric2, keep.internals = TRUE)
dtwPlotDensity(am)
abline(coef = c(0,1))

to_python_readable <- function(dataset, gd, wp, traj, wp_whole){
  pr <- list(
    milestone_network = dataset$milestone_network,
    counts = as.matrix(dataset$counts),
    gd = gd,
    cell_ids = rownames(dataset$counts),
    genes = colnames(dataset$counts),
    wp = wp,
    progressions = dataset$progressions,
    traj_percentages = traj,
    wp_whole_traj = wp_whole
  )

  pr
}

pr1 <- to_python_readable(dataset1, gd1, wpp1, traj1, wpp1)
pr2 <- to_python_readable(dataset2, gd2, wpp2, traj2, wpp2)

saveRDS(pr1, "../trajectory_alignment/data/python_readable_big1.rds")
saveRDS(pr2, "../trajectory_alignment/data/python_readable_big2.rds")


library(ggplot2)
library(reshape2)
library(pheatmap)
library(cellAlign)

exp1 <- t(as.matrix(dataset1$counts))
exp2 <- t(as.matrix(dataset2$counts))

numPts <- 100
inter1 <- cellAlign::interWeights(expDataBatch = exp1, trajCond = traj1, winSz = 0.1, numPts = numPts)
inter2 <- cellAlign::interWeights(expDataBatch = exp2, trajCond = traj2, winSz = 0.1, numPts = numPts)

alignment <- globalAlign(inter1$interpolatedVals, inter2$interpolatedVals,
                         scores = list(query = interScaled1$traj, ref = interScaled2$traj),
                         sigCalc = F, numPerm = 20)
plotAlign(alignment)
abline(coef = c(0,1))
alignment2 <- localAlign(inter1$interpolatedVals,inter2$interpolatedVals,threshPercent = 0.25)
plotAlign(alignment2)
