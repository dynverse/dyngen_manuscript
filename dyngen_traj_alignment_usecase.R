library(tidyverse)
library(dyngen)
library(dynplot)
library(dtw)

set.seed(276)
# set.seed(15)

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
    simulation_params = simulation_default(census_interval = .01, experiment_params = bind_rows(simulation_type_wild_type(num_simulations = 32), simulation_type_knockdown(num_simulations = 0)))
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
dists_cells <- list()

for(exp_i in seq(1:10)){
#
  lm1 <- runmodel()
  lm2 <- lm1 %>% generate_cells() %>% generate_experiment()

  # dlm1 <- wrap_dataset(lm1)
  # dlm2 <- wrap_dataset(lm2)
  # dlm1 <- readRDS(paste0("2shuffle_lm1", exp_i))
  # dlm2 <- readRDS(paste0("2shuffle_lm2", exp_i))

  # Generate datasets with C cut off
  # rmC1 <- surpress_module(lm1, "C1")
  # rmC2 <- surpress_module(lm2, "C1")
  #
  # drmCD1 <- wrap_dataset(rmC1)
  # drmCD2 <- wrap_dataset(rmC2)

  names <- list("lm1", "lm2") #, "rmC1", "rmC2")
  datasets <- list(dlm1, dlm2) #, drmCD1, drmCD2)
  expressions100 <- list()
  expressions200 <- list()
  expressions1000 <- list()
  expressionscells <- list()

  new_datasets <- list()
  new_names <- list()

  for(j in seq(1:length(names))){
    dataset <- datasets[[j]]
    expressions100[[length(expressions100) + 1]] <- get_wexpr(dataset, 100)

    nor2 <- names(sort(calculate_correct_pseudotime(dataset)))
    volg2 <- sapply(nor2, function(i) strtoi(substr(i, 5, 7)))
    cnt2 <- as.matrix(dataset$counts)
    cnts2 <- cnt2[volg2,]

    expressionscells[[length(expressionscells) + 1]] <- cnts2

    name <- names[[j]]
    new_datasets <- append(new_datasets, list(ds))
    new_names <- append(new_names, name)
    saveRDS(dataset, paste0(paste0("shuffle_", name), exp_i))

    # for(i in c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)){
    for(i in c(0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2)){
      ds <- dataset

      # ds$counts[sample(length(ds$counts), length(ds$counts) * i, replace=FALSE)] <- 0
      smp <- sample(length(ds$counts), size = length(ds$counts) * i, replace=FALSE)
      res <- sample(ds$counts[smp])
      ds$counts[smp] <- res

      expressions100[[length(expressions100) + 1]] <- get_wexpr(ds, 100)

      nor2 <- names(sort(calculate_correct_pseudotime(ds)))
      volg2 <- sapply(nor2, function(i) strtoi(substr(i, 5, 7)))
      cnt2 <- as.matrix(ds$counts)
      cnts2 <- cnt2[volg2,]
      expressionscells[[length(expressionscells) + 1]] <- cnts2

      new_datasets <- append(new_datasets, list(ds))
      name <- paste0(names[[j]], i)
      new_names <- c(new_names, name)

      print(name)

      saveRDS(ds, paste0(paste0("shuffle_", name), exp_i))
    }
  }

  dist_dtw_100 <- matrix(NA, nrow=22, ncol=12)
  dist_cells <- matrix(NA, nrow=22, ncol=12)

  # i = lm1 of lm2
  for(i in seq(1:length(new_datasets))){
    expr1_100 <- expressions100[[i]]
    expr1_cells <- expressionscells[[i]]

    for(j in c(1, 12)){
      # if(i > j){
        expr2_100 <- expressions100[[j]]
        expr2_cells <- expressionscells[[j]]

        expr1 <- expr1_100[1:50,]
        expr2 <- expr2_100[1:50,]

        am <- dtw(dist(as.matrix(t(expr1)), as.matrix(t(expr2)), dist.method="correlation"), step.pattern = symmetric2, keep.internals = TRUE, open.end = FALSE)
        dist_dtw_100[i, j] <- am$normalizedDistance

        # somehow, this is saved differently, so transposing isn't needed
        cnts1 <- expr1_cells[,1:50]
        cnts2 <- expr2_cells[,1:50]

        am3 <- dtw(dist(cnts1, cnts2, dist.method="correlation"), step.pattern = symmetric2, keep.internals = TRUE, open.end = FALSE)
        dist_cells[i, j] <- am3$normalizedDistance

        print(i)

      # }
    }
  }


  colnames(dist_dtw_100) <- new_names[1:12]
  rownames(dist_dtw_100) <- new_names
  colnames(dist_cells) <- new_names[1:12]
  rownames(dist_cells) <- new_names

  dists_100[[length(dists_100) + 1]] <- dist_dtw_100
  dists_cells[[length(dists_cells) + 1]] <- dist_cells

  saveRDS(dist_dtw_100, paste0("shuffle_dist_dtw_100_", exp_i))
  saveRDS(dist_cells, paste0("shuffle_dist_cells_", exp_i))

  print(paste0("DONE WITH ROUND ", exp_i))

}

saveRDS(dists_100, "shuffle_dists_100")
# saveRDS(dists_200, "dists_200")
# saveRDS(dists_1000, "dists_1000")
saveRDS(dists_cells, "shuffle_dists_cells")

d_1_100 <- readRDS("shuffle_dists_100")
d_1_200 <- readRDS("dists_200")
d_1_1000 <- readRDS("dists_1000")
d_1_cell <- readRDS("shuffle_dists_cells") #

dists_100 <- lapply(seq(1:10), function(i) readRDS(paste0("5shuffle_dist_dtw_100_", i)))
dists_cells <- lapply(seq(1:10), function(i) readRDS(paste0("5shuffle_dist_cells_", i)))

avg_100 <- Reduce('+', d_1_100) / length(d_1_100)
avg_200 <- Reduce('+', d_1_200) / length(d_1_200)
avg_1000 <- Reduce('+', d_1_1000) / length(d_1_1000)
avg_cells <- Reduce('+', d_1_cell) / length(d_1_cell)


expr1 <- expressions100[[1]][1:50,]
expr2 <- expressions100[[2]][1:50,]

am <- dtw(dist(as.matrix(t(expr1)), as.matrix(t(expr2)), dist.method="correlation"), step.pattern = symmetric2, keep.internals = TRUE, open.end = FALSE)
dtwPlotTwoWay(am, expr1, expr2)

get_ratio_lm1_lm2 <- function(data){
  verh <- list(0)
  for(i in seq(1:15)){
    verh[[i]] <- data[[i+1]] / data[[22+i]]
  }
  verh
}

get_ratio_lm2_lm1 <- function(data){
  verh <- list(0)
  for(i in seq(1:15)){
    verh[[i]] <- data[[i+22]] / data[[1+i]]
  }

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

df <- data.frame(c(unlist(lm1s_100), unlist(lm2s_100), unlist(lm1s_200), unlist(lm2s_200), unlist(lm1s_cells), unlist(lm2s_cells)))
df$model <- c(rep("100 wp, traj 1", 50), rep("100 wp, traj 2", 50), rep("200 wp,  traj 1", 50), rep("200 wp, traj 2", 50), rep("all cells,  traj 1", 50), rep("all cells,  traj 2", 50))
df$Noise <- rep(c("remove milestone C", "10% noise", "20% noise", "remove milestone C + 10% noise", "remove milestone C + 20% noise"), 60)
df <- df %>% rename(distance.ratio = c.unlist.lm1s_100...unlist.lm2s_100...unlist.lm1s_200...unlist.lm2s_200...)
ggplot2::qplot(model, distance.ratio, data=df, aes(color = Noise), geom='beeswarm') + scale_y_log10() + theme_bw() + ylab("Discrimination power (higher is better)") + xlab("")
library(ggbeeswarm)
df <- data.frame(c(unlist(lm1s_100), unlist(lm2s_100), unlist(lm1s_cells), unlist(lm2s_cells)))
df$Model <- c(rep("100 wp, traj 1", 50), rep("100 wp, traj 2", 50), rep("all cells,  traj 1", 50), rep("all cells,  traj 2", 50))
df$Noise <- rep(c("remove milestone C", "10% noise", "20% noise", "remove milestone C + 10% noise", "remove milestone C + 20% noise"), 40)
df <- df %>% rename(distance.ratio = c.unlist.lm1s_100...unlist.lm2s_100...unlist.lm1s_cells...unlist.lm2s_cells..)
ggplot2::qplot(Model, distance.ratio, data=df, aes(color = Noise), geom='beeswarm') + scale_y_log10() + theme_bw() + ylab("Discrimination power (higher is better)") + xlab("")

plot(dtw(dist(cnts1, cnts2, dist.method="correlation"), step.pattern = symmetric2, keep.internals = TRUE, open.end = FALSE), type="twoway")

lm1s_100 <- lapply(d_1_100, function(data) get_ratios(data[,1]))
lm2s_100 <- lapply(d_1_100, function(data) get_ratios2(data[,2]))

# long scatter plot noise levels
df <- data.frame(c(unlist(lm1s_100), unlist(lm2s_100), unlist(lm1s_cells), unlist(lm2s_cells)))
df$m <- c(rep("smoothed", 200), rep("unsmoothed", 200))
df$x <- rep(c("0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1"), 40)
df$s <- c(rep("traj1_sm", 100), rep("traj2_sm", 100), rep("traj1_nosm", 100), rep("traj2_nosm", 100))
# df$Model <- c(rep("100 wp, traj 1", 50), rep("100 wp, traj 2", 50), rep("all cells,  traj 1", 50), rep("all cells,  traj 2", 50))
df$Noise <- rep(c("remove milestone C", "10% noise", "20% noise", "remove milestone C + 10% noise", "remove milestone C + 20% noise"), 40)
df <- df %>% rename(distance.ratio = c.unlist.lm1s_100...unlist.lm2s_100...unlist.lm1s_cells...unlist.lm2s_cells..)
ggplot2::qplot(x, distance.ratio, data=df, geom='beeswarm', aes(color=s)) + scale_y_log10() + theme_bw() + ylab("Discrimination power (higher is better)") + xlab("")


lm1a <- get_ratios2(avg_100[, 1])
lm2a <- get_ratios(avg_100[, 2])

lm1acells <- get_ratios2(avg_cells[, 1])
lm2acells <- get_ratios(avg_cells[, 2])
# line plot
df <- data.frame(c(unlist(lm1a), unlist(lm2a), unlist(lm1acells), unlist(lm2acells)))
df$m <- c(rep("smoothed", 20), rep("unsmoothed", 20))
df$x <- rep(c("0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1"), 4)
df$s <- c(rep("traj1_sm", 10), rep("traj2_sm", 10), rep("traj1_nosm", 10), rep("traj2_nosm", 10))
df$t <- c(rep("traj1_sm", 20), rep("traj1_nosm", 20))
# df$Model <- c(rep("100 wp, traj 1", 50), rep("100 wp, traj 2", 50), rep("all cells,  traj 1", 50), rep("all cells,  traj 2", 50))
# df$Noise <- rep(c("remove milestone C", "10% noise", "20% noise", "remove milestone C + 10% noise", "remove milestone C + 20% noise"), 40)
df <- df %>% rename(distance.ratio = c.unlist.lm1a...unlist.lm2a...unlist.lm1acells...unlist.lm1acells..)
ggplot2::qplot(x, distance.ratio, data=df, geom='beeswarm', aes(color=t)) + scale_y_log10() + theme_bw() + ylab("Discrimination power (higher is better)") + xlab("")

lm1a <- get_ratios2(avg_100[, 1])
lm2a <- get_ratios(avg_100[, 2])

lm1acells <- get_ratios2(avg_cells[, 1])
lm2acells <- get_ratios(avg_cells[, 2])

lm1a <- get_ratios2(avg_100[, 1])
lm2a <- get_ratios(avg_100[, 2])

dists_100 <- readRDS("2shuffle_dists_100")
dists_cells <- readRDS("2shuffle_dists_cells")

lm1s_100 <- lapply(dists_100, function(data) get_ratios2(data[,1]))
lm2s_100 <- lapply(dists_100, function(data) get_ratios(data[,2]))

lm1acells <- lapply(dists_cells, function(data) get_ratios2(data[,1]))
lm2acells <- lapply(dists_cells, function(data) get_ratios(data[,2]))

# line plot
df <- data.frame(c(unlist(lm1s_100), unlist(lm2s_100), unlist(lm1acells), unlist(lm2acells)))
df$m <- c(rep("smoothed", 200), rep("unsmoothed", 200))
df$x <- rep(c("0.02", "0.04", "0.06", "0.08", "0.1", "0.12", "0.14", "0.16", "0.18", "0.2"), 40)
df$s <- c(rep("traj1_sm", 100), rep("traj2_sm", 100), rep("traj1_nosm", 100), rep("traj2_nosm", 100))
df$t <- c(rep("traj1_sm", 200), rep("traj1_nosm", 200))
# df$Model <- c(rep("100 wp, traj 1", 50), rep("100 wp, traj 2", 50), rep("all cells,  traj 1", 50), rep("all cells,  traj 2", 50))
# df$Noise <- rep(c("remove milestone C", "10% noise", "20% noise", "remove milestone C + 10% noise", "remove milestone C + 20% noise"), 40)
df <- df %>% rename(distance.ratio = c.unlist.lm1s_100...unlist.lm2s_100...unlist.lm1acells...unlist.lm2acells..)
ggplot2::qplot(x, distance.ratio, data=df, geom='beeswarm', aes(color=t)) + scale_y_log10() + theme_bw() + ylab("Discrimination power (higher is better)") + xlab("")


library(ggbeeswarm)
dists_100 <- readRDS("shuffle_005_1_dists_100")
dists_cells <- readRDS("shuffle_005_1_dists_cells")

lm1s_100 <- lapply(dists_100, function(data) get_ratio_lm2_lm1(data[,1]))
lm2s_100 <- lapply(dists_100, function(data) get_ratio_lm1_lm2(data[,22]))

lm1acells <- lapply(dists_cells, function(data) get_ratio_lm2_lm1(data[,1]))
lm2acells <- lapply(dists_cells, function(data) get_ratio_lm1_lm2(data[,22]))
# line plot
df <- data.frame(c(unlist(lm1s_100), unlist(lm2s_100), unlist(lm1acells), unlist(lm2acells)))
df$m <- c(rep("smoothed", 300), rep("unsmoothed", 300))
# df$x <- rep(c("0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9"), 40)
df$x <- rep(sapply(seq(0.05,0.75,0.05), toString), 40)
df$s <- c(rep("traj1_sm", 150), rep("traj2_sm", 150), rep("traj1_nosm", 150), rep("traj2_nosm", 150))
# df$Model <- c(rep("100 wp, traj 1", 50), rep("100 wp, traj 2", 50), rep("all cells,  traj 1", 50), rep("all cells,  traj 2", 50))
# df$Noise <- rep(c("remove milestone C", "10% noise", "20% noise", "remove milestone C + 10% noise", "remove milestone C + 20% noise"), 40)
df <- df %>% rename(distance.ratio = c.unlist.lm1s_100...unlist.lm2s_100...unlist.lm1acells...unlist.lm2acells..)

ggplot(data = df, aes(x = x, y = distance.ratio, color = m)) +
  geom_beeswarm(dodge.width=0.2) +
  theme_bw() +
  ylab("Discrimination power (higher is better)") +
  xlab("Amount of noise added") +
  scale_y_log10(breaks = scales::pretty_breaks(n = 10)) +
  theme(legend.position = "bottom", text = element_text(size=20)) +
  scale_color_discrete(name = "Smoothing method", labels = c("100 smoothed pseudocells", "original cells"))
  labs(colour = "Smoothing method")
# ggplot2::qplot(x, distance.ratio, data=df, geom='beeswarm', aes(color=m)) +
#   theme_bw() +
#   ylab("Discrimination power (higher is better)") +
#   xlab("Amount of noise added") +
#   scale_y_log10(breaks = scales::pretty_breaks(n = 10)) +
#   theme(legend.position = "bottom")



avg_100 <- Reduce('+', dists_100) / length(dists_100)
avg_cells <- Reduce('+', dists_cells) / length(dists_cells)
ra100_1 <- get_ratio_lm2_lm1(avg_100[,1])
ra100_2 <- get_ratio_lm1_lm2(avg_100[,2])

rac_1 <- get_ratio_lm2_lm1(avg_cells[,1])
rac_2 <- get_ratio_lm1_lm2(avg_cells[,2])

df <- data.frame(c(unlist(ra100_1), unlist(ra100_2), unlist(rac_1), unlist(rac_2)))
df$m <- c(rep("smoothed", 18), rep("unsmoothed", 18))
# df$x <- rep(c("0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9"), 4)
df$x <- rep(c("0.02", "0.04", "0.06", "0.08", "0.1", "0.12", "0.14", "0.16", "0.18"), 4)
df$s <- c(rep("traj1_sm", 9), rep("traj2_sm", 9), rep("traj1_nosm", 9), rep("traj2_nosm", 9))
# df$t <- c(rep("traj1_sm", 200), rep("traj1_nosm", 200))
# df$Model <- c(rep("100 wp, traj 1", 50), rep("100 wp, traj 2", 50), rep("all cells,  traj 1", 50), rep("all cells,  traj 2", 50))
# df$Noise <- rep(c("remove milestone C", "10% noise", "20% noise", "remove milestone C + 10% noise", "remove milestone C + 20% noise"), 40)
df <- df %>% rename(distance.ratio = c.unlist.ra100_1...unlist.ra100_2...unlist.rac_1...unlist.rac_2..)
ggplot2::qplot(x, distance.ratio, data=df, geom='beeswarm', aes(color=s)) + scale_y_log10() + theme_bw() + ylab("Discrimination power (higher is better)") + xlab("")


ra100_1 <- get_ratio_lm2_lm1(dist_dtw_100[,1])
ra100_2 <- get_ratio_lm1_lm2(dist_dtw_100[,12])

rac_1 <- get_ratio_lm2_lm1(dist_cells[,1])
rac_2 <- get_ratio_lm1_lm2(dist_cells[,12])

e1 <- expressions100[[1]][15,]
e2 <- expressions100[[4]][15,]
a <- dtw(e1, e2, step.pattern = symmetric2, keep.internals = TRUE, open.end = FALSE)
dtwPlotTwoWay(a, e1, e2, match.indices = 100, offset = 5, xlab = "Cell", ylab = "Gene expression")
legend("bottomright",c("Trajectory 1","Trajectory 2 (right axis)           "), col=1:2, lt=c(1,2), bty="n")

d1 <- data.frame(e1) %>% tibble::rownames_to_column("X")
d1$X <- as.numeric(as.character(d1$X))
d2 <- data.frame(e2) %>% tibble::rownames_to_column("X")
d2$X <- as.numeric(as.character(d2$X))
p1 <- ggplot(data = d1, aes(y=e1, x=X, group = 1)) + geom_line() + theme_void()
p2 <- ggplot(data = d2, aes(y=e2, x=X, group = 1)) + geom_line() + theme_void()

# Generate twowayplot
X1_edges <- a$index1
X2_edges <- a$index2
dtwoway_edge1 <- data.frame(X1_edges)
dtwoway_edge2 <- data.frame(X2_edges)

total_edge1 <- merge(d1, dtwoway_edge1, by.x = "X", by.y = "X1_edges")
total_edge2 <- merge(d2, dtwoway_edge2, by.x = "X", by.y = "X2_edges")


edges <- total_edge1
edges$X2 <- total_edge2$X
edges$Y1 <- total_edge1$e1
edges$Y2 <- total_edge2$e2 + 5
edges$X1 <- edges$X


dtwoway <- data.frame(e1, e2) %>% tibble::rownames_to_column("X")
dtwoway$X <-as.numeric(as.character(dtwoway$X))
dtwoway$e2 <- dtwoway$e2 + 5

p_lines <- ggplot(data=dtwoway) +
  scale_color_brewer(palette = "Dark2") +
  geom_line(aes(x=X, y=e1, group = 1, color="a"), size=2) +
  geom_line(aes(x=X, y=e2, group= 1, color = "b"), size=2) +
  geom_segment(data = edges, aes(x = X1, y = Y1, xend = X2, yend = Y2), colour="grey", lty=1, size=0.7) +
  theme_void() +
  theme(legend.position = "none")


dtwoway_1 <- dtwoway
dtwoway_1$e2 <- dtwoway_1$e2 + 20
p1 <- ggplot(data=dtwoway_1) +
  geom_line(aes(x=X, y=e1, group = 1), colour = "#1B9E77", size=2) +
  geom_line(aes(x=X, y=e2, group = 1), colour = "#D95F02", size=2) +
  theme_void() +
  labs(col="") +
  theme(legend.position = "bottom")

p2 <- ggplot(data=dtwoway_1) +
  geom_line(aes(x=X, y=e2, group = 1), colour = "#D95F02", size=2) +
  theme_void() +
  labs(col="") +
  theme(legend.position = "bottom")
