get_ratio_lm1_lm2 <- function(data, lm1_on_lm2){
  if(lm1_on_lm2){
    return(sapply(seq(1:15), function(i) data[[i]] / data[[15+i]]))
  } else {
    return(sapply(seq(1:15), function(i) data[[i+15]] / data[[i]]))
  }
}

get_ratios <- function(dists){
  ratio100_1 <- lapply(dists, function(data) get_ratio_lm1_lm2(data[[1]], F))
  ratio100_2 <- lapply(dists, function(data) get_ratio_lm1_lm2(data[[2]], T))
  return(list(ratio100_1, ratio100_2))
}

get_2nd_ratio <- function(data){
  return(sapply(seq(from = 2, to = 16), function(i) data[[i]] / data[[1]]))
}

plot_graph <- function(df, aest){
  ggplot(data = df, aest) +
    ggbeeswarm::geom_beeswarm() + 
    theme_bw() + 
    xlab("Amount of noise added") +
    scale_y_log10(breaks = scales::pretty_breaks(n = 10)) +
    theme(legend.position = "bottom", text = element_text(size=20))
}

ratios <- c(unlist(lapply(x$distnoisy, get_2nd_ratio)), 
            unlist(lapply(x$distnoisyorig, get_2nd_ratio)), 
            unlist(lapply(x$distnoisyeucl, get_2nd_ratio)), 
            unlist(lapply(x$distnoisyorigeucl, get_2nd_ratio)))
df <- data.frame(ratios)
df$x <- rep(sapply(seq(0.05, 0.75, 0.05), toString), 4 * nr_it)
df$exp <- c(unlist(lapply(seq(1,nr_it), function(i) rep(paste0("sm", i), 15))), 
            unlist(lapply(seq(1,nr_it), function(i) rep(paste0("og", i), 15))),
            unlist(lapply(seq(1,nr_it), function(i) rep(paste0("sm_eucl", i), 15))), 
            unlist(lapply(seq(1,nr_it), function(i) rep(paste0("og_eucl", i), 15))))
df$s <- c(rep("smooth", 15 * nr_it), rep("orig", 15 * nr_it), rep("smooth_eucl", 15 * nr_it), rep("orig_eucl", 15 * nr_it))

ggplot(data = df, aes(x = x, y = ratios, group = exp, color=s)) +
  geom_line() +
  theme_bw() +
  ylab("d(t1n, t2n)/d(t1, t2)") +
  xlab("Amount of noise added") +
  scale_y_log10(breaks = scales::pretty_breaks(n = 10)) +
  theme(legend.position = "bottom", text = element_text(size=20))


# Get ratio data in a dataframe
ratios <- c(unlist(get_ratios(x$dist100)), unlist(get_ratios(x$distorig)), unlist(get_ratios(x$dist100eucl)), unlist(get_ratios(x$distorigeucl)))
df <- data.frame(ratios)
df$color <- c(rep("smoothed", 15 * 2 * nr_it), rep("unsmoothed", 15 * 2 * nr_it), rep("smoothed_eucl", 15 * 2 * nr_it), rep("unsmoothed_eucl", 15 * 2 * nr_it))
df$x <- rep(sapply(seq(0.05, 0.75, 0.05), toString), 8 * nr_it)
df$s <- c(rep("traj1_sm", 15 * nr_it), rep("traj2_sm", 15 * nr_it), rep("traj1_nosm", 15 * nr_it), rep("traj2_nosm", 15 * nr_it), rep("traj1_sm", 15 * nr_it), rep("traj2_sm", 15 * nr_it), rep("traj1_nosm", 15 * nr_it), rep("traj2_nosm", 15 * nr_it))

# Plot ratios (discrimination power) over noise
ggplot(data = df, aes(x = x, y = ratios, color = color)) +
  geom_beeswarm() +
  theme_bw() +
  ylab("Discrimination power (higher is better)") +
  xlab("Amount of noise added") +
  scale_y_log10(breaks = scales::pretty_breaks(n = 13)) +
  theme(legend.position = "bottom", text = element_text(size=20))
# scale_color_brewer(palette = "Dark2", name = "Smoothing method", labels = c("100 smoothed pseudocells 1", "100 smoothed pseudocells 2", "original cells 1", "original cells 2"))



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