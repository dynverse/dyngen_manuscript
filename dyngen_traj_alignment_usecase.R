library(tidyverse)
library(dyngen)
library(dynplot)
library(dtw)

set.seed(276)

setwd("/home/louisedc/Work/dyngen_analysis/analysis")
source("dtw_functions.R")

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

get_wexpr <- function(dataset, amount){
  wpp1 <- get_waypoint_progression(dataset, amount)
  gd1 <- get_geodesic_distances_from_progressions(dataset, wpp1)
  wexpr1 <- interpolate_expression(dataset, wpp1$percentage, gd1)
  wexpr1
}

surpress_module <- function(model, module){
  fts_to_surpress <- which(module == model$feature_network$to_module)
  for(index in fts_to_surpress){
    ft <- model$feature_network$from[[index]]
    model <- change_speed(model, ft, 0)
  }
  model <- model %>% generate_gold_standard() %>% generate_cells() %>% generate_experiment()
  model
}

backbone <- bblego(
  bblego_start("A", type = "simple", num_modules = 4),
  bblego_linear("A", "B", type = "simple", num_modules = 6),
  bblego_linear("B", "C", type = "simple", num_modules = 6),
  bblego_end("C", type = "simple", num_modules = 6)
)

runmodel <- function(){
  model <- initialise_model(
    num_tfs = 50,
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
    generate_feature_network() %>%
    generate_kinetics() %>%
    generate_gold_standard() %>%
    generate_cells() %>%
    generate_experiment()
  model
}

dists_100 <- list()
dists_200 <- list()
dists_1000 <- list()
dists_cells <- list()

for(exp_i in seq(1:2)){

  # lm1 <- runmodel()
  # lm2 <- lm1 %>% generate_cells() %>% generate_experiment()

  # Generate first two linear datasets
  # lm1 <- initialise_model(
  #   num_tfs = 50,
  #   num_targets = 200,
  #   num_hks = 200,
  #   num_cells = 1000,
  #   backbone = backbone,
  #   verbose = TRUE,
  #   num_cores = 8,
  #   distance_metric = "pearson",
  #   tf_network_params = tf_network_default(min_tfs_per_module = 2, sample_num_regulators = function() 1),
  #   simulation_params = simulation_default(census_interval = .01)
  # ) %>%
  #   generate_tf_network() %>%
  #   generate_feature_network()
  #
  # lm2 <- lm1

  lm1 <- runmodel()
  lm2 <- lm1 %>% generate_cells() %>% generate_experiment()

  # lm2$feature_info <- lm1$feature_info %>%
  #   mutate_at(c("basal", "independence"), ~ . * rnorm(length(.), mean = 1, sd = .01)) %>%
  #   mutate_at(c("basal", "independence"), ~ pmin(., 1))
  # lm2$feature_network <- lm1$feature_network %>%
  #   mutate_at(c("strength", "cooperativity"), ~ . * rnorm(length(.), mean = 1, sd = .01))
  # lm2 <- lm2 %>%
  #   generate_kinetics() %>%
  #   generate_gold_standard() %>%
  #   generate_cells() %>%
  #   generate_experiment()
  #
  # lm1 <- lm1 %>% generate_kinetics() %>%
  #   generate_gold_standard() %>%
  #   generate_cells() %>%
  #   generate_experiment

  # lm2 <- lm1 %>% generate_cells() %>% generate_experiment()
  dlm1 <- wrap_dataset(lm1)
  dlm2 <- wrap_dataset(lm2)

  # Generate datasets with C cut off
  rmC1 <- surpress_module(lm1, "C1")
  rmC2 <- surpress_module(lm2, "C1")

  drmCD1 <- wrap_dataset(rmC1)
  drmCD2 <- wrap_dataset(rmC2)

  names <- list("lm1", "lm2", "rmC1", "rmC2")
  datasets <- list(dlm1, dlm2, drmCD1, drmCD2)
  expressions100 <- list()
  expressions200 <- list()
  expressions1000 <- list()
  for(j in seq(1:length(datasets))){
    dataset <- datasets[[j]]
    expressions100[[length(expressions100) + 1]] <- get_wexpr(dataset, 100)
    expressions200[[length(expressions200) + 1]] <- get_wexpr(dataset, 200)
    expressions1000[[length(expressions1000) + 1]] <- get_wexpr(dataset, 1000)

    for(i in c(0.1, 0.2)){
      ds <- dataset
      ds$counts[sample(length(ds$counts), length(ds$counts) * i, replace=FALSE)] <- 0

      expressions100[[length(expressions100) + 1]] <- get_wexpr(ds, 100)
      expressions200[[length(expressions200) + 1]] <- get_wexpr(ds, 200)
      expressions1000[[length(expressions1000) + 1]] <- get_wexpr(ds, 1000)

      datasets <- append(datasets, list(ds))
      name <- paste0(names[[j]], i)
      names <- c(names, name)

      print(name)
    }
  }

  lm3 <- runmodel()
  dlm3 <- wrap_dataset(lm3)
  expr3_100 <- get_wexpr(dlm3, 100)
  expr3_200 <- get_wexpr(dlm3, 200)
  expr3_1000 <- get_wexpr(dlm3, 1000)
  datasets <- append(datasets, list(dlm3))
  names <- c(names, "lm3")

  expressions100[[length(expressions100) + 1]] <- expr3_100
  expressions200[[length(expressions200) + 1]] <- expr3_200
  expressions1000[[length(expressions1000) + 1]] <- expr3_1000

  dist_dtw_100 <- matrix(NA, nrow=13, ncol=13)
  dist_dtw_200 <- matrix(NA, nrow=13, ncol=13)
  dist_dtw_1000 <- matrix(NA, nrow=13, ncol=13)
  dist_cells <- matrix(NA, nrow=13, ncol=13)

  for(i in seq(1:length(datasets))){
    expr1_100 <- expressions100[[i]]
    expr1_200 <- expressions200[[i]]
    expr1_1000 <- expressions1000[[i]]
    dataset1 <- datasets[[i]]

    for(j in seq(1:length(datasets))){
      expr2_100 <- expressions100[[j]]
      expr2_200 <- expressions200[[j]]
      expr2_1000 <- expressions1000[[j]]

      if(i > j){

        expr1 <- expr1_100[1:50,]
        expr2 <- expr2_100[1:50,]

        am <- dtw(dist(as.matrix(t(expr1)), as.matrix(t(expr2)), dist.method="correlation"), step.pattern = symmetric2, keep.internals = TRUE, open.end = FALSE)
        dist_dtw_100[i, j] <- am$normalizedDistance

        expr1 <- expr1_200[1:50,]
        expr2 <- expr2_200[1:50,]

        am <- dtw(dist(as.matrix(t(expr1)), as.matrix(t(expr2)), dist.method="correlation"), step.pattern = symmetric2, keep.internals = TRUE, open.end = FALSE)
        dist_dtw_200[i, j] <- am$normalizedDistance

        expr1 <- expr1_1000[1:50,]
        expr2 <- expr2_1000[1:50,]

        am <- dtw(dist(as.matrix(t(expr1)), as.matrix(t(expr2)), dist.method="correlation"), step.pattern = symmetric2, keep.internals = TRUE, open.end = FALSE)
        dist_dtw_1000[i, j] <- am$normalizedDistance


        dataset2 <- datasets[[j]]
        nor1 <- names(sort(calculate_correct_pseudotime(dataset1)))
        volg1 <- sapply(nor1, function(i) strtoi(substr(i, 5, 7)))
        cnt1 <- as.matrix(dataset1$counts)
        cnts1 <- cnt1[volg1,]

        nor2 <- names(sort(calculate_correct_pseudotime(dataset2)))
        volg2 <- sapply(nor2, function(i) strtoi(substr(i, 5, 7)))
        cnt2 <- as.matrix(dataset2$counts)
        cnts2 <- cnt2[volg2,]

        cnts1 <- cnts1[1:50,]
        cnts2 <- cnts2[1:50,]

        am3 <- dtw(dist(cnts1, cnts2, dist.method="correlation"), step.pattern = symmetric2, keep.internals = TRUE, open.end = FALSE)

        dist_cells[i, j] <- am3$normalizedDistance


        print(j)

      }
    }
  }


  colnames(dist_dtw_100) <- names
  rownames(dist_dtw_100) <- names
  colnames(dist_dtw_200) <- names
  rownames(dist_dtw_200) <- names
  colnames(dist_dtw_1000) <- names
  rownames(dist_dtw_1000) <- names
  colnames(dist_cells) <- names
  rownames(dist_cells) <- names

  dists_100[[length(dists_100) + 1]] <- dist_dtw_100
  dists_200[[length(dists_200) + 1]] <- dist_dtw_200
  dists_1000[[length(dists_1000) + 1]] <- dist_dtw_1000
  dists_cells[[length(dists_cells) + 1]] <- dist_cells

  # saveRDS(dist_dtw_100, paste0("2dist_dtw_100_", exp_i))
  # saveRDS(dist_dtw_200, paste0("2dist_dtw_200_", exp_i))
  # saveRDS(dist_dtw_1000, paste0("2dist_dtw_1000_", exp_i))
  # saveRDS(dist_cells, paste0("2dist_cells_", exp_i))

  print(paste0("DONE WITH ROUND ", exp_i))

}

saveRDS(dists_100, "2dists_100")
saveRDS(dists_200, "2dists_200")
saveRDS(dists_1000, "2dists_1000")
saveRDS(dists_cells, "2_dists_cells")

d_1_100 <- readRDS("dists_100")
d_1_200 <- readRDS("dists_200")
d_1_1000 <- readRDS("dists_1000")
d_1_cell <- lapply(seq(1:10), function(i) readRDS(paste0("dist_cells_", i)))

avg_100 <- Reduce('+', d_1_100) / length(d_1_100)
avg_200 <- Reduce('+', d_1_200) / length(d_1_200)
avg_1000 <- Reduce('+', d_1_1000) / length(d_1_1000)
avg_cells <- Reduce('+', d_1_cell) / length(d_1_cell)


expr1 <- expressions100[[1]][1:50,]
expr2 <- expressions100[[2]][1:50,]

am <- dtw(dist(as.matrix(t(expr1)), as.matrix(t(expr2)), dist.method="correlation"), step.pattern = symmetric2, keep.internals = TRUE, open.end = FALSE)
dtwPlotTwoWay(am, expr1, expr2)

get_ratios1 <- function(data){
  verh <- list(0)
  verh[[1]] <- data[[3]] / data[[4]]
  verh[[2]] <- data[[5]] / data[[7]]
  verh[[3]] <- data[[6]] / data[[8]]
  verh[[4]] <- data[[9]] / data[[11]]
  verh[[5]] <- data[[10]] / data[[12]]
  verh
}

get_ratios2 <- function(data){
  verh <- list(0)
  verh[[1]] <- data[[4]] / data[[3]]
  verh[[2]] <- data[[7]] / data[[5]]
  verh[[3]] <- data[[8]] / data[[6]]
  verh[[4]] <- data[[11]] / data[[9]]
  verh[[5]] <- data[[12]] / data[[10]]
  verh
}


plot_data <- function(data, model){

}

lm1s_100 <- lapply(d_1_100, function(data) get_ratios2(data[,1]))
lm2s_100 <- lapply(d_1_100, function(data) get_ratios(data[,2]))
lm1s_200 <- lapply(d_1_200, function(data) get_ratios2(data[,1]))
lm2s_200 <- lapply(d_1_200, function(data) get_ratios(data[,2]))
lm1s_1000 <- lapply(d_1_1000, function(data) get_ratios2(data[,1]))
lm2s_1000 <- lapply(d_1_1000, function(data) get_ratios(data[,2]))
lm1s_cells <- lapply(d_1_cell, function(data) get_ratios2(data[,1]))
lm2s_cells <- lapply(d_1_cell, function(data) get_ratios(data[,2]))

# df <- data.frame(c(unlist(lm1s_100), unlist(lm2s_100), unlist(lm1s_1000), unlist(lm2s_1000), unlist(lm1s_cells), unlist(lm2s_cells)))
# df$model <- c(rep("100 wp, traj 1", 50), rep("100 wp, traj 2", 50), rep("1000 wp,  traj 1", 50), rep("1000 wp, traj 2", 50), rep("all cells,  traj 1", 50), rep("all cells,  traj 2", 50))
# df$color <- rep(c("remove milestone C", "10% noise", "20% noise", "remove milestone C + 10% noise", "remove milestone C + 20% noise"), 60)
# df <- df %>% rename(distance.ratio = c.unlist.lm1s_100...unlist.lm2s_100...unlist.lm1s_1000...unlist.lm2s_1000...)
# ggplot2::qplot(model, distance.ratio, data=df, aes(color = color), geom='beeswarm') + scale_y_log10() + theme_minimal()

df <- data.frame(c(unlist(lm1s_100), unlist(lm2s_100), unlist(lm1s_cells), unlist(lm2s_cells)))
df$Model <- c(rep("100 wp, traj 1", 50), rep("100 wp, traj 2", 50), rep("all cells,  traj 1", 50), rep("all cells,  traj 2", 50))
df$Noise <- rep(c("remove milestone C", "10% noise", "20% noise", "remove milestone C + 10% noise", "remove milestone C + 20% noise"), 40)
df <- df %>% rename(distance.ratio = c.unlist.lm1s_100...unlist.lm2s_100...unlist.lm1s_cells...unlist.lm2s_cells..)
ggplot2::qplot(Model, distance.ratio, data=df, aes(color = Noise), geom='beeswarm') + scale_y_log10() + theme_bw() + ylab("Discrimination power (higher is better)") + xlab("")

plot(dtw(dist(cnts1, cnts2, dist.method="correlation"), step.pattern = symmetric2, keep.internals = TRUE, open.end = FALSE), type="twoway")
