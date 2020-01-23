library(tidyverse)
library(dyngen)
library(dynplot)
library(dtw)
library(ggbeeswarm)
library(reshape2)

setwd("/home/louisedc/Work/dyngen_analysis/analysis/trajectory_alignment")
output_dir <- paste0(getwd(), "/data/")
source("dtw_functions.R")
source("1_runmodels.R")
source("2_get_expression.R")
source("3_plot.R")

set.seed(123)

## GENERAL
lm1 <- runmodel()
dlm2 <- lm1 %>% generate_cells() %>% generate_experiment() %>% wrap_dataset()
dlm1 <- wrap_dataset(lm1)

# VALIDATE DISCRIMINATION POWER METRIC
## 1. INDIVIDUAL GENES?

## 2. DIMINISH WINDOW SIZE
get_dataset <- function(new_dataset = T, i = 1){
  if(new_dataset){
    lm1 <- runmodel()
    dlm2 <- lm1 %>% generate_cells() %>% generate_experiment() %>% wrap_dataset()
    dlm1 <- wrap_dataset(lm1)
    
    saveRDS(dlm1, paste0(output_dir, paste0("dlm1_", i)))
    saveRDS(dlm2, paste0(output_dir, paste0("dlm2_", i)))
  } else{
    dlm1 <- readRDS(paste0(output_dir, paste0("dlm1_", i)))
    dlm2 <- readRDS(paste0(output_dir, paste0("dlm2_", i)))
  }
  return(list(dlm1=dlm1, dlm2=dlm2))
}

get_normalized_dist <- function(expr1, expr2){
  dm <- proxy::dist(expr1[,1:50], expr2[,1:50], method="euclidean")
  am <- dtw(dm, step.pattern=symmetric2, keep.internals=T)
  
  return(list(dtw = am$normalizedDistance))
  # return(list(dtw = am$normalizedDistance, eucl = sum(diag(dm)) / length(diag(dm))))
}

get_disc_distances <- function(models, noisy_models){
  mapply(function(lms, noisys){
    distances <- list()
    for(i in 1:length(lms)){
      distances[[i]] <- list()
      for(j in 1:length(noisys)){
        m1 <- lms[[i]]
        m2 <- noisys[[j]]
        
        distances[[i]][[j]] <- get_normalized_dist(m1, m2)
      }
    }
    distances
  }, models, noisy_models, SIMPLIFY = F)
}

get_ratio_lm1_lm2 <- function(data, lm1_on_lm2, metric="dtw"){
  if(lm1_on_lm2){
    return(sapply(seq(1:15), function(i) data[[i]][[metric]] / data[[i+15]][[metric]]))
  } else {
    return(sapply(seq(1:15), function(i) data[[i+15]][[metric]] / data[[i]][[metric]]))
  }
}

get_ratios <- function(dists){
  ratio_dtw1 <- lapply(dists, function(data) get_ratio_lm1_lm2(data[[1]], F))
  ratio_dtw2 <- lapply(dists, function(data) get_ratio_lm1_lm2(data[[2]], T))
  # ratio_euc1 <- lapply(dists, function(data) get_ratio_lm1_lm2(data[[1]], F, metric="eucl"))
  # ratio_euc2 <- lapply(dists, function(data) get_ratio_lm1_lm2(data[[2]], T, metric="eucl"))
  return(list(dtw1 = ratio_dtw1, dtw2 = ratio_dtw2))
}

cl_expr_dss <- lapply(dss, lapply, get_cell_expression)
cl_expr_noisy <- lapply(noisy_data, lapply, lapply, get_cell_expression) %>% lapply(unlist, recursive=F)
cl_disc_power <- get_disc_distances(cl_expr_dss, cl_expr_noisy)
cl_ratio_dp <- get_ratios(cl_disc_power) %>% melt() %>% rename(experiment = L2) %>% 
  add_column(metric = c(rep("dtw", 75), rep("euc", 75))) %>%
  add_column(noise = rep(seq(from = 0.05, to = 0.75, by = 0.05), 10)) %>% add_column(wp = F) %>% add_column(window_size = 1)

generate_experiment <- function(ws, nr_wp){
  # Generate two datasets & the noisy versions of those two.
  dss <- lapply(seq(from=1, to=5), get_dataset, new_dataset = F)
  noisy_data <- lapply(dss, lapply, get_all_noisy_data, start=0.05, stop=0.75, step=0.05)
  
  print("DATA")
  
  # Get expressions of waypoints & actual cells of the two datasets. Then, get them for all the noisy
  # versions as well.
  wp_expr_dss <- lapply(dss, lapply, get_wexpr, amount = nr_wp)
  
  wp_expr_noisy <- lapply(noisy_data, lapply, lapply, get_wexpr, amount = nr_wp, ws = ws) %>% lapply(unlist, recursive=F)
  
  print("EXPRESSION")
  
  # Calculate DTW between datasets. Calculate distance between model 1 and all noisy models,
  # and distance between model 2 and all noisy models.
  wp_disc_power <- get_disc_distances(wp_expr_dss, wp_expr_noisy)
  
  print("DISC DISTANCES")
  
  wp_ratio_dp <- get_ratios(wp_disc_power) %>% melt() %>% rename(experiment = L2) %>% 
    add_column(metric = c(rep("dtw", 75), rep("euc", 75))) %>%
    add_column(noise = rep(seq(from = 0.05, to = 0.75, by = 0.05), 10)) %>% add_column(wp = T)

  ratio_dp <- wp_ratio_dp %>% add_column(window_size = ws) %>% add_column(waypoints = nr_wp)
  
  print("DONE")
  
  return(ratio_dp)
  
}

exps <- lapply(c(0.1, 0.05, 0.025, 0.01, 0.005), generate_experiment, nr_wp = 100)
res <- exps %>% reduce(full_join) %>% full_join(cl_ratio_dp)
res$window_size <- as.factor(res$window_size)
saveRDS(exps, paste0(output_dir, "results_ws.RDS"))

plot_graph(res, aes(x = noise, y = value, group = experiment, color = window_size))
ggsave("plot_windowsize.png", path = output_dir)

exps2 <- lapply(c(100, 250, 500, 750, 1000), generate_experiment, ws = 0.05)
res2 <- exps2 %>% reduce(full_join) %>% full_join(cl_ratio_dp)
res2$waypoints <- as.factor(res2$waypoints)
saveRDS(exps2, paste0(output_dir, "results_wp.RDS"))


plot_graph(res2, aes(x = noise, y = value, group = experiment, color = waypoints))
ggsave("plot_waypoints.png", path = output_dir)


exps3 <- lapply(c(0.1, 0.05, 0.025, 0.01, 0.005), function(i) lapply(c(100, 250, 500, 750, 1000), function(j) generate_experiment(i, j)))
exps3_100 <- lapply(c(0.1, 0.05, 0.025, 0.01, 0.005), generate_experiment, nr_wp = 100)
exps3_250 <- lapply(c(0.1, 0.05, 0.025, 0.01, 0.005), generate_experiment, nr_wp = 250)
exps3_500 <- lapply(c(0.1, 0.05, 0.025, 0.01, 0.005), generate_experiment, nr_wp = 500)
exps3_750 <- lapply(c(0.1, 0.05, 0.025, 0.01, 0.005), generate_experiment, nr_wp = 750)
exps3_1000 <- lapply(c(0.1, 0.05, 0.025, 0.01, 0.005), generate_experiment, nr_wp = 1000)

saveRDS(exps3_100, "exps3_100.RDS")
saveRDS(exps3_250, "exps3_250.RDS")
saveRDS(exps3_500, "exps3_500.RDS")

res3_100 <- exps3_100 %>% reduce(full_join) %>% full_join(cl_ratio_dp)
res3_100$window_size <- as.factor(res3_100$window_size)
plot_graph(res3_100, aes(x = noise, y = value, group = experiment, color = window_size))

res3_250 <- exps3_250 %>% reduce(full_join) %>% full_join(cl_ratio_dp)
res3_250$window_size <- as.factor(res3_250$window_size)
plot_graph(res3_250, aes(x = noise, y = value, group = experiment, color = window_size))

res3_500 <- exps3_500 %>% reduce(full_join) %>% full_join(cl_ratio_dp)
res3_500$window_size <- as.factor(res3_500$window_size)
plot_graph(res3_500, aes(x = noise, y = value, group = experiment, color = window_size))

## PARALLEL TRY
library(parallel)
cl <- makeCluster(5, type="FORK")
res_par_750 <- parLapply(cl, c(0.1, 0.05, 0.025, 0.01, 0.005), generate_experiment, nr_wp = 750)
res_par_1000 <- parLapply(cl, c(0.1, 0.05, 0.025, 0.01, 0.005), generate_experiment, nr_wp = 1000)
stopCluster(cl)

saveRDS(res_par_750, "exps3_750.RDS")
saveRDS(res_par_1000, "exps3_1000.RDS")


saveRDS(exps3, paste0(output_dir, "results_wp_ws.RDS"))
# TODO: exps3 plotten & impact van beide zien
res3 <- exps3 %>% reduce(full_join)


# wp_dp_df <- melt(wp_disc_power) %>% rename(experiment = L1, model = L2, noise = L3, metric = L4) %>% add_column(wp = T)
# cl_dp_df <- melt(cl_disc_power) %>% rename(experiment = L1, model = L2, noise = L3, metric = L4) %>% add_column(wp = F)
# dp_df <- full_join(wp_dp_df, cl_dp_df)

## 3. INCREASE WAYPOINTS

