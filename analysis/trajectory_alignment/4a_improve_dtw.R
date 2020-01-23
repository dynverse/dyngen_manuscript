library(tidyverse)
library(dyngen)
library(dynplot)
library(dtw)
library(ggbeeswarm)

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
gen <- 6

## BASELINE DTW ALIGNMENT
expression1 <- get_wexpr(dlm1, 500, 0.1)
expression2 <- get_wexpr(dlm2, 500, 0.1)
e1 <- expression1[,gen]
e2 <- expression2[,gen]
a <- dtw(e1, e2, step.pattern = symmetric2, keep.internals = TRUE, open.end = FALSE)
dtwPlotTwoWay(a)

## HOW TO IMPROVE DTW RESULTS:
## 1. SCALE GENE EXPRESSION

scale_expression <- function(expr){
  scaled_expression <- do.call('cbind', lapply(1:ncol(expr), function(ind){
    minv <- min(expr[,ind])
    maxv <- max(expr[,ind])
    scaled <- (expr[,ind] - minv) / (maxv - minv)
  }))
}

scaled_expression1 <- scale_expression(expression1)
scaled_expression2 <- scale_expression(expression2)
e3 <- scaled_expression1[,gen]
e4 <- scaled_expression2[,gen]
a_scale <- dtw(e3, e4, step.pattern = symmetric2, keep.internals = TRUE, open.end = FALSE)
dtwPlotTwoWay(a_scale)

## 2. PERTURB & AVERAGE

rnorm_bounded <- function(n, mean = 0, sd = 1, min = -Inf, max = Inf) {
  unif_min <- pnorm(min, mean = mean, sd = sd)
  unif_max <- pnorm(max, mean = mean, sd = sd)
  quan <- runif(n, unif_min, unif_max)
  qnorm(quan, mean = mean, sd = sd)
}


noisy_model <- function(model, mean, min, sd){
  m <- model
  m$counts <- m$counts * rnorm_bounded(length(m$counts), mean = mean, min = min, sd = sd)
  m
}

aggregate_models <- function(models){
  Reduce("+", models) / length(models)
}

gen_avg_noisy_models <- function(model, it){
  ms <- lapply(seq(from = 1, to = it), function(i){
    m1 <- noisy_model(model, 1, 0, 0.1)
    get_wexpr(m1, 500, 0.1)
  })
  aggregate_models(ms)
}

mm1 <- gen_avg_noisy_models(dlm1, 5)
mm2 <- gen_avg_noisy_models(dlm2, 5)

e1 <- mm1[,gen]
e2 <- mm2[,gen]
a_pert <- dtw(e1, e2, step.pattern = symmetric1, keep.internals = TRUE, open.end = FALSE)
dtwPlotTwoWay(a_pert)


