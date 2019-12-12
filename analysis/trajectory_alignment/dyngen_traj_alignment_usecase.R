library(tidyverse)
library(dyngen)
library(dynplot)
library(dtw)
library(ggbeeswarm)

set.seed(276)

output_dir <- paste0(getwd(), "/data/")
source("dtw_functions.R")
source("1_runmodels.R")
source("2_get_expression.R")

noise_experiment <- function(i){
  # Generate two datasets containing a trajectory
  lm1 <- runmodel()
  dlm2 <- lm1 %>% generate_cells() %>% generate_experiment() %>% wrap_dataset()
  dlm1 <- wrap_dataset(lm1)

  saveRDS(dlm1, paste0(output_dir, paste0("dlm1_", i)))
  saveRDS(dlm2, paste0(output_dir, paste0("dlm2_", i)))

  # Generate noisy datasets for each base dataset
  # noise -> [dlm1, dlm2, dlm1_0.05, dlm1_0.1, dlm1_0.15, ...]
  noisy <- lapply(list(dlm1, dlm2), get_all_noisy, start=0.05, stop=1, step=0.05)
  noisy <- c(list(dlm1), list(dlm2), unlist(noisy, recursive = F))

  # Get expression of pseudocells or for actual cells
  expr100 <- lapply(noisy, get_wexpr, amount = 100)
  exprcells <- lapply(noisy, get_cell_expression)

  # Get distances between expr_base & all noisy ones
  get_normalized_dists <- function(expr_base, expr_list){
    sapply(expr_list, get_normalized_dist, expr1=expr_base)
  }

  get_normalized_dist <- function(expr1, expr2){
    am <- dtw(dist(expr1[,1:50], expr2[,1:50], dist.method="correlation"), step.pattern=symmetric2, keep.internals=T)
    am$normalizedDistance
  }

  dist100 <- lapply(list(expr100[[1]], expr100[[2]]), get_normalized_dists, expr_list=expr100[3:length(expr100)])
  distcells <- lapply(list(exprcells[[1]], exprcells[[2]]), get_normalized_dists, expr_list=exprcells[3:length(exprcells)])

  return(list(wp=dist100, orig=distcells, dlm1=dlm1, dlm2=dlm2))
}

nr_it <- 10
res <- lapply(seq(1:nr_it), noise_experiment)
x <-list()
x$dist100 <- lapply(res, "[[", 1)
x$distcells <- lapply(res, "[[", 2)

get_ratio_lm1_lm2 <- function(data, lm1_on_lm2){
  verh <- list()
  if(lm1_on_lm2){
    verh <- sapply(seq(1:15), function(i) data[[i]] / data[[20+i]])
  } else {
    verh <- sapply(seq(1:15), function(i) data[[i+20]] / data[[i]])
  }

  verh
}

get_ratios <- function(dists){
  ratio100_1 <- lapply(dists, function(data) get_ratio_lm1_lm2(data[[1]], F))
  ratio100_2 <- lapply(dists, function(data) get_ratio_lm1_lm2(data[[2]], T))
  return(list(ratio100_1, ratio100_2))
}

# Get ratio data in a dataframe
ratios <- c(unlist(get_ratios(x$dist100)), unlist(get_ratios(x$distcells)))
df <- data.frame(ratios)
df$m <- c(rep("smoothed", 15 * 2 * nr_it), rep("unsmoothed", 15 * 2 * nr_it))
df$x <- rep(sapply(seq(0.05, 0.75, 0.05), toString), 4 * nr_it)
df$s <- c(rep("traj1_sm", 15 * nr_it), rep("traj2_sm", 15 * nr_it), rep("traj1_nosm", 15 * nr_it), rep("traj2_nosm", 15 * nr_it))

# Plot ratios (discrimination power) over noise
ggplot(data = df, aes(x = x, y = ratios, color = m)) +
  geom_beeswarm(dodge.width=0.2) +
  theme_bw() +
  ylab("Discrimination power (higher is better)") +
  xlab("Amount of noise added") +
  scale_y_log10(breaks = scales::pretty_breaks(n = 10)) +
  theme(legend.position = "bottom", text = element_text(size=20)) +
  scale_color_brewer(palette = "Dark2", name = "Smoothing method", labels = c("100 smoothed pseudocells", "original cells"))

# Plot lines & edges
e1 <- get_wexpr(res[[1]]$dlm1, 100)[,15]
e2 <- get_wexpr(res[[1]]$dlm2, 100)[,15]
a <- dtw(e1, e2, step.pattern = symmetric2, keep.internals = TRUE, open.end = FALSE)

d1 <- data.frame(e1) %>% tibble::rownames_to_column("X")
d1$X <- as.numeric(as.character(d1$X))
d2 <- data.frame(e2) %>% tibble::rownames_to_column("X")
d2$X <- as.numeric(as.character(d2$X))

X1_edges <- a$index1
X2_edges <- a$index2
total_edge1 <- merge(d1, data.frame(X1_edges), by.x = "X", by.y = "X1_edges")
total_edge2 <- merge(d2, data.frame(X2_edges), by.x = "X", by.y = "X2_edges")

edges <- total_edge1
edges$X2 <- total_edge2$X
edges$Y1 <- total_edge1$e1
edges$Y2 <- total_edge2$e2
edges$X1 <- edges$X
edges$X <- NULL
edges$e1 <- NULL

p_lines <- ggplot(data=edges) +
  geom_line(aes(x=X1, y=Y1, group = 1), colour = "#1B9E77", size=2) +
  geom_line(aes(x=X2, y=Y2+5, group = 1), colour = "#D95F02", size=2) +
  geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2+5), colour="grey", lty=1, size=0.7) +
  theme_void() +
  theme(legend.position = "none")

p1 <- ggplot(data=edges) +
  geom_line(aes(x=X1, y=Y1, group = 1), colour = "#1B9E77", size=2) +
  geom_line(aes(x=X2, y=Y2 + 15, group = 1), colour = "#D95F02", size=2) +
  theme_void() +
  labs(col="")
