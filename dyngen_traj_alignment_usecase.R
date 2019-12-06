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

for(exp_i in seq(1:10)){

  # Generate first two linear datasets
  lm1 <- runmodel()
  lm2 <- lm1 %>% generate_cells() %>% generate_experiment()
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

        expr1 <- expr1_100[,1:50]
        expr2 <- expr2_100[,1:50]

        am <- dtw(dist(as.matrix(t(expr1)), as.matrix(t(expr2)), dist.method="correlation"), step.pattern = symmetric2, keep.internals = TRUE, open.end = FALSE)
        dist_dtw_100[i, j] <- am$normalizedDistance

        expr1 <- expr1_200[,1:50]
        expr2 <- expr2_200[,1:50]

        am <- dtw(dist(as.matrix(t(expr1)), as.matrix(t(expr2)), dist.method="correlation"), step.pattern = symmetric2, keep.internals = TRUE, open.end = FALSE)
        dist_dtw_200[i, j] <- am$normalizedDistance

        expr1 <- expr1_1000[,1:50]
        expr2 <- expr2_1000[,1:50]

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

        cnts1 <- cnts1[,1:50]
        cnts2 <- cnts2[,1:50]

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
  dists_cells[[length(dists_cells) + 1]] <- dist_cellss

  print(paste0("DONE WITH ROUND ", exp_i))

}




