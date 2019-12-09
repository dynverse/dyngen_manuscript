library(tidyverse)
library(dyngen)
library(dynplot)
library(dtw)

set.seed(42)

setwd("/home/louisedc/Work/dyngen_analysis/analysis")
source("dtw_functions.R")
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

lm1 <- runmodel()
lm2 <- lm1 %>% generate_cells() %>% generate_experiment()
dlm1 <- wrap_dataset(lm1)
dlm2 <- wrap_dataset(lm2)

surpress_module <- function(model, module){
  fts_to_surpress <- which(module == model$feature_network$to_module)
  for(index in fts_to_surpress){
    ft <- model$feature_network$from[[index]]
    model <- change_speed(model, ft, 0)
  }
  model <- model %>% generate_gold_standard() %>% generate_cells() %>% generate_experiment()
  model
}

fts_to_surpress <- which("C1" == lm1$feature_network$to_module)
for(index in fts_to_surpress){
  ft <- lm1
}

rmC1 <- change_speed(lm1, "B6_TF2", 0) %>% generate_gold_standard() %>% generate_cells() %>% generate_experiment()
rmC2 <- change_speed(lm2, "B6_TF2", 0) %>% generate_gold_standard() %>% generate_cells() %>% generate_experiment()

drmCD1 <- wrap_dataset(rmCD1)
drmCD2 <- wrap_dataset(rmCD2)

# dlm1_50 <- dlm1
# dlm1_50$counts[sample(length(dlm1_50$counts), length(dlm1_50$counts) * 0.5, replace=FALSE)] <- 0


#--------
# names <- list("lm1", "lm2", "rmD1", "rmD2", "rmCD1", "rmCD2")
# datasets <- list(dlm1, dlm2, drmD1, drmD2, drmCD1, drmCD2)
names <- list("lm1", "lm2", "rmCD1", "rmCD2")
datasets <- list(dlm1, dlm2, drmCD1, drmCD2)
expressions200 <- list()
expressions500 <- list()
expressions1000 <- list()
expressions5000 <- list()
for(j in seq(1:length(datasets))){
  dataset <- datasets[[j]]
  expressions200[[length(expressions200) + 1]] <- get_wexpr(dataset, 200)
  # expressions500[[length(expressions) + 1]] <- get_wexpr(dataset, 500)
  expressions1000[[length(expressions1000) + 1]] <- get_wexpr(dataset, 1000)
  # expressions5000[[length(expressions) + 1]] <- get_wexpr(dataset, 5000)

  for(i in c(0.1, 0.2)){
    ds <- dataset
    ds$counts[sample(length(ds$counts), length(ds$counts) * i, replace=FALSE)] <- 0

    expressions200[[length(expressions200) + 1]] <- get_wexpr(ds, 200)
    # expressions500[[length(expressions) + 1]] <- get_wexpr(ds, 500)
    expressions1000[[length(expressions1000) + 1]] <- get_wexpr(ds, 1000)
    # expressions5000[[length(expressions) + 1]] <- get_wexpr(ds, 5000)

    datasets <- append(datasets, list(ds))
    name <- paste0(names[[j]], i)
    names <- c(names, name)

    print(name)
  }
}

dist_dtw <- matrix(0, nrow=12, ncol=12)
dist_eucl <- matrix(0, nrow=12, ncol=12)
# vector("list", 18)
# for(i in seq(1:18)){
#   # dist_dtw[[i]] <- vector("list", 18)
# }
# dist_eucl <- vector("list", 18)
# for(i in seq(1:18)){
#   dist_eucl[[i]] <- vector("list", 18)
# }

dist_dtw_200 <- matrix(0, nrow=12, ncol=12)
dist_euc_200 <- matrix(0, nrow=12, ncol=12)
dist_dtw_1000 <- matrix(0, nrow=12, ncol=12)
dist_euc_1000 <- matrix(0, nrow=12, ncol=12)
dist_dtw_5000 <- matrix(0, nrow=12, ncol=12)
dist_euc_5000 <- matrix(0, nrow=12, ncol=12)

for(i in seq(1:length(datasets))){
  expr1_200 <- expressions200[[i]]
  # expr1_500 <- expressions500[[i]]
  expr1_1000 <- expressions1000[[i]]
  # expr1_5000 <- expressions5000[[i]]

  for(j in seq(1:length(datasets))){
    expr2_200 <- expressions200[[j]]
    # expr2_500 <- expressions500[[j]]
    expr2_1000 <- expressions1000[[j]]
    # expr2_5000 <- expressions5000[[j]]

    if(i > j){

      expr1 <- expr1_200[,1:50]
      expr2 <- expr2_200[,1:50]

      d1 <- dist(as.matrix(t(expr1)), as.matrix(t(expr2)), dist.method="correlation")
      am <- dtw(d1, step.pattern = symmetric2, keep.internals = TRUE, open.end = FALSE)
      dist_dtw_200[i, j] <- am$normalizedDistance

      euc <- sum(diag(d1))
      dist_euc_200[i, j] <- euc[[1]]

      expr1 <- expr1_1000[,1:50]
      expr2 <- expr2_1000[,1:50]

      d1 <- dist(as.matrix(t(expr1)), as.matrix(t(expr2)), dist.method="correlation")
      am <- dtw(d1, step.pattern = symmetric2, keep.internals = TRUE, open.end = FALSE)
      dist_dtw_1000[i, j] <- am$normalizedDistance

      euc <- sum(diag(d1))
      dist_euc_1000[i, j] <- euc[[1]]


      # expr1 <- expr1_5000[,1:50]
      # expr2 <- expr2_5000[,1:50]
      #
      # d1 <- dist(as.matrix(t(expr1)), as.matrix(t(expr2)), dist.method="correlation")
      # am <- dtw(d1, step.pattern = symmetric2, keep.internals = TRUE, open.end = FALSE)
      # dist_dtw_5000[i, j] <- am$normalizedDistance
      #
      # euc <- sum(diag(d1))
      # dist_eucl_5000[i, j] <- euc[[1]]

      print(j)

    }
  }
}

dist_nosmooth <- matrix(NA, nrow=12, ncol=12)
dist_nosmooth_eucl <- matrix(NA, nrow=12, ncol=12)

for(i in seq(1:length(datasets))){
  dataset1 <- datasets[[i]]
  for(j in seq(1:length(datasets))){
    if(i > j){
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

      d3 <- dist(cnts1, cnts2, dist.method="correlation")
      am3 <- dtw(d3, step.pattern = symmetric2, keep.internals = TRUE, open.end = FALSE)

      dist_nosmooth[i, j] <- am3$normalizedDistance

      euc <- sum(diag(d3))
      dist_nosmooth_eucl[i, j] <- euc[[1]]

    }
  }
}


colnames(dist_euc_200) <- names
rownames(dist_euc_200) <- names
colnames(dist_euc_1000) <- names
rownames(dist_euc_1000) <- names
colnames(dist_dtw_200) <- names
rownames(dist_dtw_200) <- names
colnames(dist_dtw_1000) <- names
rownames(dist_dtw_1000) <- names
colnames(dist_nosmooth) <- names
rownames(dist_nosmooth) <- names
colnames(dist_nosmooth_eucl) <- names
rownames(dist_nosmooth_eucl) <- names

dist_dtw_format <- dist_dtw %>%
  as_tibble() %>%
  rowid_to_column(var="X") %>%
  gather(key="Y", value="Z", -1)

dist_eucl_format <- dist_eucl %>%
  as_tibble() %>%
  rowid_to_column(var="X") %>%
  gather(key="Y", value="Z", -1)

ggplot(dist_dtw_format, aes(X, Y, fill=Z)) + geom_tile()
ggplot(dist_eucl_format, aes(X, Y, fill=Z)) + geom_tile()

d1 <- dist(as.matrix(t(expressions1000[[1]])), as.matrix(t(expressions1000[[4]])), dist.method="Euclidean")
am <- dtw(d1, step.pattern = symmetric2, keep.internals = TRUE, open.end = FALSE)
dtwPlotDensity(am)
#--------


dlm1_expr <- get_wexpr(dlm1, 1000)
dlm1_50_expr <- get_wexpr(dlm1_50, 1000)
dlm2_expr <- get_wexpr(dlm2, 1000)

d_1_50 <- dist(as.matrix(t(dlm1_expr)), as.matrix(t(dlm1_50_expr)))
a_1_50 <- dtw(d_1_50, step.pattern = symmetric2, keep.internals = TRUE, open.end = FALSE)

dtwPlotDensity(a_1_50)
abline(coef = c(0,1))

d_1_2 <- dist(as.matrix(t(dlm1_expr)), as.matrix(t(dlm2_expr)))
a_1_2 <- dtw(d_1_2, step.pattern = symmetric2, keep.internals = TRUE, open.end = FALSE)

dtwPlotDensity(a_1_2)
abline(coef = c(0,1))

drmD1_expr <- get_wexpr(drmD1, 1000)

d_1_rmD <- dist(as.matrix(t(dlm1_expr)), as.matrix(t(drmD1_expr)))
a_1_rmD <- dtw(d_1_rmD, step.pattern = symmetric2, keep.internals = TRUE, open.end = FALSE)

dtwPlotDensity(a_1_rmD)
abline(coef = c(0,1))

d_1_50_rmD <- dist(as.matrix(t(dlm1_50_expr)), as.matrix(t(drmD1_expr)))
a_1_50_rmD <- dtw(d_1_50_rmD, step.pattern = symmetric2, keep.internals = TRUE, open.end = FALSE)

dtwPlotDensity(a_1_50_rmD)
abline(coef = c(0,1))

thing <- dist(as.matrix(t(dlm1_50_expr)), as.matrix(t(drmD1_expr)))

wpp1 <- get_waypoint_progression(dataset1, 1000)
wpp11 <- get_waypoint_progression(dataset11, 1000)
wpp2 <- get_waypoint_progression(dataset2, 1000)
wpp3 <- get_waypoint_progression(dataset3, 1000)

gd1 <- get_geodesic_distances_from_progressions(dataset1, wpp1)
gd11 <- get_geodesic_distances_from_progressions(dataset11, wpp11)
gd2 <- get_geodesic_distances_from_progressions(dataset2, wpp2)
gd3 <- get_geodesic_distances_from_progressions(dataset3, wpp3)

wexpr1 <- interpolate_expression(dlm1, wpp1$percentage, get_geodesic_distances_from_progressions(dlm1))
wexpr11 <- interpolate_expression(dataset11, wpp11$percentage, gd11)
wexpr2 <- interpolate_expression(dataset2, wpp2$percentage, gd2)
wexpr3 <- interpolate_expression(dataset3, wpp3$percentage, gd3)

# wexpr1 <- wexpr1[1:50,]
# wexpr2 <- wexpr2[1:50,]


d1 <- dist(as.matrix(t(wexpr1)), as.matrix(t(wexpr11)), dist.method="Euclidean")
# d1 <- (d1 - min(d1)) / (max(d1) - min(d1))
am <- dtw(d1, step.pattern = symmetric2, keep.internals = TRUE, open.end = FALSE)
dtwPlotDensity(a)
abline(coef = c(0,1))

#------------------------------------------------------------------------------

model2 <- change_speed(model1, "C1_TF1", 0)
model2 <- generate_gold_standard(model2)
model2 <- generate_cells(model2)

model2_ <- generate_experiment(model2)
dataset2 <- wrap_dataset(model2_)


model1 <- runmodel()
model2 <- runmodel()
model3 <- runmodel()
model4 <- runmodel()
model5 <- runmodel()

dataset1 <- wrap_dataset(model1)
dataset2 <- wrap_dataset(model2)
dataset3 <- wrap_dataset(model3)
dataset4 <- wrap_dataset(model4)
dataset5 <- wrap_dataset(model5)

model11 <- model1 %>% generate_cells() %>% generate_experiment()
dataset11 <- wrap_dataset(model11)

model1 <-
  initialise_model(
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
  generate_cells()

library(dtw)
setwd("/home/louisedc/Work/dyngen_analysis/analysis")
source("dtw_functions.R")

wpp1 <- get_waypoint_progression(dataset1, 1000)
wpp11 <- get_waypoint_progression(dataset11, 1000)
wpp2 <- get_waypoint_progression(dataset2, 1000)
wpp3 <- get_waypoint_progression(dataset3, 1000)

gd1 <- get_geodesic_distances_from_progressions(dataset1, wpp1)
gd11 <- get_geodesic_distances_from_progressions(dataset11, wpp11)
gd2 <- get_geodesic_distances_from_progressions(dataset2, wpp2)
gd3 <- get_geodesic_distances_from_progressions(dataset3, wpp3)

wexpr1 <- interpolate_expression(dataset1, wpp1$percentage, gd1)
wexpr11 <- interpolate_expression(dataset11, wpp11$percentage, gd11)
wexpr2 <- interpolate_expression(dataset2, wpp2$percentage, gd2)
wexpr3 <- interpolate_expression(dataset3, wpp3$percentage, gd3)

# wexpr1 <- wexpr1[1:50,]
# wexpr2 <- wexpr2[1:50,]


d1 <- dist(as.matrix(t(wexpr1)), as.matrix(t(wexpr11)), dist.method="Euclidean")
# d1 <- (d1 - min(d1)) / (max(d1) - min(d1))
am <- dtw(d1, step.pattern = symmetric2, keep.internals = TRUE, open.end = FALSE)

d2 <- dist(as.matrix(t(wexpr2)), as.matrix(t(wexpr3)), dist.method="Euclidean")
am <- dtw(d2, step.pattern = symmetric2, keep.internals = TRUE, open.end = FALSE)

# Plot the dynamic time warping -----------------------------------------------

dtwPlotDensity(am)
abline(coef = c(0,1))

