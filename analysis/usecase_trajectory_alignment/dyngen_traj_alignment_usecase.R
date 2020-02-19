library(tidyverse)
library(dyngen)
library(dynplot)
library(dtw)
library(ggbeeswarm)

set.seed(42)

setwd("/home/louisedc/Work/dyngen_analysis/analysis/trajectory_alignment")
output_dir <- paste0(getwd(), "/data/")
source("dtw_functions.R")
source("1_runmodels.R")
source("2_get_expression.R")
source("3_plot.R")

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
  
  return(list(dtw = am$normalizedDistance, eucl = sum(diag(dm)) / length(diag(dm))))
}

get_disc_distances <- function(model, noisy_models){
  mapply(function(lms, noisys){
    mapply(function(lm, noisy){
      lapply(noisy, get_normalized_dist, expr2=lm)
    }, lms, noisys, SIMPLIFY = F)
  }, model, noisy_models, SIMPLIFY = F)
}

dss <- lapply(seq(from=1, to=5), get_dataset, new_dataset = F)
noisy_data <- lapply(dss, lapply, get_all_noisy_data, start=0.05, stop=0.75, step=0.05)

wp_expr_dss <- lapply(dss, lapply, get_wexpr, amount = 100)
cl_expr_dss <- lapply(dss, lapply, get_cell_expression)

wp_expr_noisy <- lapply(noisy_data, lapply, lapply, get_wexpr, amount = 100)
cl_expr_noisy <- lapply(noisy_data, lapply, lapply, get_cell_expression)

wp_disc_power <- get_disc_distances(wp_expr_dss, wp_expr_noisy)
cl_disc_power <- get_disc_distances(cl_expr_dss, cl_expr_noisy)

wp_noisy_metric2 <- lapply(wp_expr_noisy, function(m) mapply(get_normalized_dist, m$dlm1, m$dlm2, SIMPLIFY = F))
wp_dss_metric2 <- lapply(wp_expr_dss, function(r) get_normalized_dist(r$dlm1, r$dlm2))
cl_noisy_metric2 <- lapply(cl_expr_noisy, function(m) mapply(get_normalized_dist, m$dlm1, m$dlm2, SIMPLIFY = F))
cl_dss_metric2 <- lapply(cl_expr_dss, function(r) get_normalized_dist(r$dlm1, r$dlm2))


