library(tidyverse)
library(dyngen)
library(dynplot)

set.seed(1)

out_dir <- "fig/usecase/trajectory_alignment/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Merge two models, useful if we want to plot the models together
merge_models <- function(model, model1, model2){
  model12 <- model
  model12$feature_info <- model$feature_info %>%
    mutate(
      w = paste0("w_", feature_id),
      x = paste0("x_", feature_id),
      y = paste0("y_", feature_id)
    )
  m1gs <- model1$gold_standard
  m2gs <- model2$gold_standard
  model12$gold_standard <- list(
    meta = bind_rows(
      m1gs$meta %>% mutate_at(c("from", "to", "from_", "to_"), ~paste0("left_", .)),
      m2gs$meta %>% mutate_at(c("from", "to", "from_", "to_"), ~paste0("right_", .))
    ),
    counts = rbind(
      m1gs$counts,
      m2gs$counts
    ),
    network = bind_rows(
      m1gs$network %>% mutate_at(c("from", "to"), ~paste0("left_", .)),
      m2gs$network %>% mutate_at(c("from", "to"), ~paste0("right_", .))
    )
  )
  m1sim <- model1$simulations
  m2sim <- model2$simulations
  num_m1sim <- max(m1sim$meta$simulation_i)
  model12$simulations <- list(
    meta = bind_rows(
      m1sim$meta %>% mutate_at(c("from", "to"), ~paste0("left_", .)),
      m2sim$meta %>% mutate_at(c("from", "to"), ~paste0("right_", .)) %>% mutate(simulation_i = simulation_i + num_m1sim)
    ),
    counts = rbind(
      m1sim$counts,
      m2sim$counts
    ),
    regulation = rbind(
      m1sim$regulation,
      m2sim$regulation
    ),
    reaction_firings = rbind(
      m1sim$reaction_firings,
      m2sim$reaction_firings
    ),
    reaction_propensities = rbind(
      m1sim$reaction_propensities,
      m2sim$reaction_propensities
    )
  )
  model12
}


# Generate base model -----------------------------------------------------
backbone <- bblego(
  bblego_start("A", type = "simple", num_modules = 4),
  bblego_linear("A", "B", type = "doublerep1", num_modules = 6),
  bblego_linear("B", "C", type = "doublerep2", num_modules = 6),
  bblego_end("C", type = "simple", num_modules = 4)
)
model <-
  initialise_model(
    num_tfs = 100,
    num_targets = 200,
    num_hks = 200,
    num_cells = 1000,
    backbone = backbone,
    verbose = TRUE,
    download_cache_dir = "~/.cache/dyngen",
    num_cores = 8,
    distance_metric = "pearson",
    tf_network_params = tf_network_default(min_tfs_per_module = 2, sample_num_regulators = function() 1),
    simulation_params = simulation_default(census_interval = .01)
  ) %>%
  generate_tf_network() %>%
  generate_feature_network()


# Show GRN ----------------------------------------------------------------
plot_backbone_statenet(model)
plot_backbone_modulenet(model)
plot_feature_network(model, show_targets = FALSE)
plot_feature_network(model)
plot_feature_network(model, show_hks = TRUE)


# Generate 2 variations with slightly randomised kinetics -----------------
model1 <- model
model1$feature_info <- model1$feature_info %>%
  mutate_at(c("ba", "ind"), ~ . * rnorm(length(.), mean = 1, sd = .01)) %>%
  mutate_at(c("ba", "ind"), ~ pmin(., 1))
model1$feature_network <- model1$feature_network %>%
  mutate_at(c("strength", "cooperativity"), ~ . * rnorm(length(.), mean = 1, sd = .01))
model1 <- model1 %>%
  generate_kinetics() %>%
  generate_gold_standard() %>%
  generate_cells()

model1_ <- generate_experiment(model1)
dataset1 <- wrap_dataset(model1_)

model2 <- model
# by_backbone <- model2$feature_network$from %>% grepl("_TF", .)
# fn <- model2$feature_network
# model2$feature_network$from[by_backbone] <-
#   model2$feature_network[by_backbone, ] %>%
#   mutate(from_modmod = gsub("_TF.*", "", from)) %>%
#   group_by(from_modmod) %>%
#   mutate(
#     from = sample(from)
#   ) %>%
#   ungroup() %>%
#   select(-from_modmod) %>%
#   pull(from)
model2$feature_info <- model2$feature_info %>%
  mutate_at(c("ba", "ind"), ~ . * rnorm(length(.), mean = 1, sd = .01)) %>%
  mutate_at(c("ba", "ind"), ~ pmin(., 1))
model2$feature_network <- model2$feature_network %>%
  mutate_at(c("strength", "cooperativity"), ~ . * rnorm(length(.), mean = 1, sd = .01))
model2 <- model2 %>%
  generate_kinetics() %>%
  generate_gold_standard() %>%
  generate_cells()
model2_ <- generate_experiment(model2)
dataset2 <- wrap_dataset(model2_)

model12 <- merge_models(model, model1, model2)
model12 <- model12 %>%
  dyngen:::calculate_dimred(dimred_premrna = FALSE) %>%
  generate_experiment()

plot_gold_simulations(model12)

# View result -------------------------------------------------------------
dataset <- wrap_dataset(model12)

library(dynplot)
g <- plot_dimred(dataset)
ggsave(paste0(out_dir, "traj_dimred.pdf"), g, width = 8, height = 6)
g <- plot_graph(dataset)
ggsave(paste0(out_dir, "traj_graph.pdf"), g, width = 8, height = 6)
g <- plot_heatmap(dataset, features_oi = 40)
ggsave(paste0(out_dir, "traj_heatmap.pdf"), g, width = 8, height = 6)

# Apply dynamic time warping --------------------------------------------------

library(dtw)
setwd("/home/louisedc/Work/dyngen_analysis/analysis")
source("dtw_functions.R")

traj_order <- calculate_correct_pseudotime(dataset1)

wpp1 <- get_waypoint_progression(dataset1, 200)
wpp2 <- get_waypoint_progression(dataset2, 200)

gd1 <- get_geodesic_distances_from_progressions(dataset1, wpp1)
gd2 <- get_geodesic_distances_from_progressions(dataset2, wpp2)

wexpr1 <- interpolate_expression(dataset1, wpp1$percentage, gd1)
wexpr2 <- interpolate_expression(dataset2, wpp2$percentage, gd2)

wexpr1 <- wexpr1[1:40,]
wexpr2 <- wexpr2[1:40,]


d1 <- dist(as.matrix(t(wexpr1)), as.matrix(t(wexpr2)), dist.method="Euclidean")
# d1 <- (d1 - min(d1)) / (max(d1) - min(d1))
am <- dtw(d1, step.pattern = symmetric2, keep.internals = TRUE, open.end = FALSE)

# Plot the dynamic time warping -----------------------------------------------

dtwPlotDensity(am)
abline(coef = c(0,1))

# Is dtw robust? --------------------------------------------------------------
## Make initial, unperturbed model --------------------------------------------

backbone <- bblego(
  bblego_start("A", type = "simple", num_modules = 4),
  bblego_linear("A", "B", type = "doublerep1", num_modules = 6),
  bblego_linear("B", "C", type = "doublerep2", num_modules = 6),
  bblego_end("C", type = "simple", num_modules = 4)
)
model_init <-
  initialise_model(
    num_tfs = 100,
    num_targets = 200,
    num_hks = 200,
    num_cells = 1000,
    backbone = backbone,
    verbose = TRUE,
    download_cache_dir = "~/.cache/dyngen",
    num_cores = 8,
    distance_metric = "pearson",
    tf_network_params = tf_network_default(min_tfs_per_module = 2, sample_num_regulators = function() 1),
    simulation_params = simulation_default(num_simulations = 10, census_interval = .01)
  ) %>%
  generate_tf_network() %>%
  generate_feature_network()

model_unperturbed <- generate_kinetics(model_init) %>%
  generate_gold_standard() %>%
  generate_cells()

model_unperturbed_ <- generate_experiment(model_unperturbed)
dataset_unperturbed <- wrap_dataset(model_unperturbed_)

## Perform generate_cells & generate_experiment again -------------------------
## This will give us different cells compared to the previous model, which
## contain different noise & will lie on different parts of the trajectory

model_replicate1 <- generate_cells(model_unperturbed)
model_replicate1_ <- generate_experiment(model_replicate1)
dataset_replicate1 <- wrap_dataset(model_replicate1_)

merged_model <- merge_models(model_init, model_unperturbed, model_replicate1)

merged_model <- merged_model %>%
  dyngen:::calculate_dimred(dimred_premrna = FALSE) %>%
  generate_experiment()

plot_gold_simulations(merged_model)
dataset_merged <- wrap_dataset(merged_model)

library(dynplot)
plot_dimred(dataset_merged)

traj_order <- calculate_correct_pseudotime(dataset_unperturbed)

# Apply dtw to these

wpp1 <- get_waypoint_progression(dataset_unperturbed, 1000)
wpp2 <- get_waypoint_progression(dataset_replicate1, 1000)

gd1 <- get_geodesic_distances_from_progressions(dataset_unperturbed, wpp1)
gd2 <- get_geodesic_distances_from_progressions(dataset_replicate1, wpp2)

wexpr1 <- interpolate_expression(dataset_unperturbed, wpp1$percentage, gd1)
wexpr2 <- interpolate_expression(dataset_replicate1, wpp2$percentage, gd2)

# wexpr1 <- wexpr1[1:40,]
# wexpr2 <- wexpr2[1:40,]


d1 <- dist(as.matrix(t(wexpr1)), as.matrix(t(wexpr2)), dist.method="Euclidean")
# d1 <- (d1 - min(d1)) / (max(d1) - min(d1))
am1 <- dtw(d1, step.pattern = symmetric2, keep.internals = TRUE, open.end = FALSE)

# Plot the dynamic time warping -----------------------------------------------

dtwPlotDensity(am1)
abline(coef = c(0,1))

# Compare different step patterns ---------------------------------------------

am2 <- dtw(d1, step.pattern = typeIds, keep.internals = TRUE, open.end = FALSE)
am3 <- dtw(d1, step.pattern = symmetric2, keep.internals = TRUE, open.end = TRUE)
am4 <- dtw(d1, step.pattern = asymmetric, keep.internals = TRUE, open.end = FALSE)

dtwPlotDensity(am2)
abline(coef = c(0,1))

# Manually perturb genes

dataset_perturbed1 <- dataset_unperturbed
dataset_perturbed1$counts[sample(length(dataset_perturbed1$counts), 500000 * 0.2, replace=FALSE)] <- 0

wpp1 <- get_waypoint_progression(dataset_unperturbed, 1000)
wpp2 <- get_waypoint_progression(dataset_perturbed1, 1000)

gd1 <- get_geodesic_distances_from_progressions(dataset_unperturbed, wpp1)
gd2 <- get_geodesic_distances_from_progressions(dataset_perturbed1, wpp2)

wexpr1 <- interpolate_expression(dataset_unperturbed, wpp1$percentage, gd1)
wexpr2 <- interpolate_expression(dataset_perturbed1, wpp2$percentage, gd2)

# wexpr1 <- wexpr1[1:40,]
# wexpr2 <- wexpr2[1:40,]


d2 <- dist(as.matrix(t(wexpr1)), as.matrix(t(wexpr2)), dist.method="Euclidean")
# d1 <- (d1 - min(d1)) / (max(d1) - min(d1))
am2 <- dtw(d2, step.pattern = symmetric2, keep.internals = TRUE, open.end = TRUE)

# Plot the dynamic time warping -----------------------------------------------

dtwPlotDensity(am2)
abline(coef = c(0,1))

# index1 -> x index
# index2 -> y index
xymap <- rep(list(), length = 100)
for(i in 1:length(am2$index1)){
  x_index = am2$index1[[i]]
  y_index = am2$index2[[i]]
  xymap[[x_index]] <- c(xymap[[x_index]], y_index)
}

y_indices <- rep(0, 100)
for(i in 1:100){
  poss_indices <- xymap[[i]]
  y_indices[i] <- max(poss_indices)
}

eucl_dist <- function(x1, x2, y1, y2){
  sqrt((x1 - x2)^2 + (y1 - y2)^2)
}

area <- MESS::auc(sequence(100), y_indices - sequence(100), absolutearea = TRUE)

dist <- sum(eucl_dist(sequence(100), sequence(100), y_indices, sequence(100)))

d3 <- dist(as.matrix(dataset_unperturbed$counts), as.matrix(dataset_replicate1$counts), dist.method="Euclidean")
am3 <- dtw(d3, step.pattern = symmetric2, keep.internals = TRUE, open.end = TRUE)
dtwPlotDensity(am3)
abline(coef = c(0,1))

# Topology change -------------------------------------------------------------

change_speed <- function(model, target, rate) {
  param_id <- which(target == model$feature_info$feature_id)[[1]]

  model$feature_info$wpr[param_id] <- model$feature_info$wpr[param_id] * rate
  model$feature_info$xdr[param_id] <- model$feature_info$xdr[param_id] * rate

  param_wpr <- paste0("wpr_", target)
  param_xdr <- paste0("xdr_", target)

  model$simulation_system$parameters[param_wpr] <- model$simulation_system$parameters[param_wpr] * rate
  model$simulation_system$parameters[param_xdr] <- model$simulation_system$parameters[param_xdr] * rate

  model
}

## Generate base model -----------------------------------------------------
backbone <- bblego(
  bblego_start("A", type = "simple", num_modules = 4),
  bblego_linear("A", "B", type = "simple", num_modules = 6),
  bblego_linear("B", "C", type = "simple", num_modules = 6),
  bblego_end("C", type = "simple", num_modules = 4)
)
model1 <-
  initialise_model(
    num_tfs = 20,
    num_targets = 200,
    num_hks = 200,
    num_cells = 1000,
    backbone = backbone,
    verbose = TRUE,
    download_cache_dir = "~/.cache/dyngen",
    num_cores = 8,
    distance_metric = "pearson",
    tf_network_params = tf_network_default(min_tfs_per_module = 1),
    simulation_params = simulation_default(num_simulations = 32, census_interval = .01),
  ) %>%
  generate_tf_network() %>%
  generate_feature_network() %>%
  generate_kinetics() %>%
  generate_gold_standard() %>%
  generate_cells()

# # Show GRN ----------------------------------------------------------------
# plot_backbone_statenet(model0)
# plot_backbone_modulenet(model0)
# plot_feature_network(model0, show_targets = FALSE)
# plot_feature_network(model0)
# plot_feature_network(model0, show_hks = TRUE)


# Generate 2 variations with slightly randomised kinetics -----------------

model1_ <- generate_experiment(model1)
dataset1 <- wrap_dataset(model1_)

model2 <- change_speed(model1, "C1_TF1", 0)
model2 <- generate_gold_standard(model2)
model2 <- generate_cells(model2)

model2_ <- generate_experiment(model2)
dataset2 <- wrap_dataset(model2_)

model12 <- merge_models(model1, model1, model2)
model12 <- model12 %>%
  dyngen:::calculate_dimred(dimred_premrna = FALSE) %>%
  generate_experiment()

plot_gold_simulations(model12)

# View result -------------------------------------------------------------
dataset <- wrap_dataset(model12)
plot_dimred(dataset)
plot_graph(dataset)
setwd("/home/louisedc/Work/dyngen_analysis/analysis")
source("dtw_functions.R")
nr_wps <- 1000

wpp1 <- get_waypoint_progression(dataset1, nr_wps)
wpp2 <- get_waypoint_progression(dataset2, nr_wps)

gd1 <- get_geodesic_distances_from_progressions(dataset1, wpp1)
gd2 <- get_geodesic_distances_from_progressions(dataset2, wpp2)

wexpr1 <- interpolate_expression(dataset1, wpp1$percentage, gd1)
wexpr2 <- interpolate_expression(dataset2, wpp2$percentage, gd2)

# wexpr1 <- wexpr1[1:100,]
# wexpr2 <- wexpr2[1:100,]

d2 <- dist(as.matrix(t(wexpr2)), as.matrix(t(wexpr1)), dist.method="pearson")
am2 <- dtw(as.matrix(t(wexpr2)), as.matrix(t(wexpr1)), step.pattern = symmetric2, keep.internals = TRUE, open.end = TRUE)

nor1 <- names(sort(calculate_correct_pseudotime(dataset1)))
volg1 <- sapply(nor1, function(i) strtoi(substr(i, 5, 7)))
cnt1 <- as.matrix(dataset1$counts)
cnts1 <- cnt1[volg1,]

nor2 <- names(sort(calculate_correct_pseudotime(dataset2)))
volg2 <- sapply(nor2, function(i) strtoi(substr(i, 5, 7)))
cnt2 <- as.matrix(dataset2$counts)
cnts2 <- cnt2[volg2,]

# cnts1 <- cnts1[,1:100]
# cnts2 <- cnts2[,1:100]

d3 <- dist(cnts2, cnts1, dist.method="pearson")
am3 <- dtw(d3, step.pattern = symmetric2, keep.internals = TRUE, open.end = TRUE)

dtwPlotDensity(am2)
abline(coef= c(0, sum(dataset1$milestone_network$length[1:2]) / sum(dataset1$milestone_network$length)))
dtwPlotDensity(am3)
abline(coef = c(0, sum(dataset1$milestone_network$length[1:2]) / sum(dataset1$milestone_network$length)))




library(MESS)

get_y_indices <- function(nr_points, index1, index2){
  xymap <- rep(list(), length = nr_points)
  for(i in 1:length(index1)){
    x_index = index1[[i]]
    y_index = index2[[i]]
    xymap[[x_index]] <- c(xymap[[x_index]], y_index)
  }

  sapply(sequence(nr_points), function(i) max(xymap[[i]]))
}

eucl_dist <- function(x1, x2, y1, y2){
  sqrt((x1 - x2)^2 + (y1 - y2)^2)
}

nr_x <- length(unique(am2$index1))

y_indices_smooth <- get_y_indices(length(unique(am2$index1)), am2$index1, am2$index2)
y_indices_no_smooth <- get_y_indices(length(unique(am2$index1)), am3$index1, am3$index2)

area_smooth <- MESS::auc(sequence(nr_x), y_indices_smooth - sequence(nr_x), absolutearea = TRUE)
dist_smooth <- sum(eucl_dist(sequence(nr_x), sequence(nr_x), y_indices_smooth, sequence(nr_x)))

area_no_smooth <- MESS::auc(sequence(nr_x), y_indices_no_smooth - sequence(nr_x), absolutearea = TRUE)
dist_no_smooth <- sum(eucl_dist(sequence(nr_x), sequence(nr_x), y_indices_no_smooth, sequence(nr_x)))

res <- data.frame("Smooth" = c("yes", "no"), "Area" = c(area_smooth, area_no_smooth), "Eucl" = c(dist_smooth, dist_no_smooth))

print(res)

