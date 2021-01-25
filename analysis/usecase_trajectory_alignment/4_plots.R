library(tidyverse)
library(dyngen.manuscript)
library(patchwork)
library(tidyr)
library(reshape2)
library(dyngen)
library(dynplot)
library(gtools)
library(dtw)

exp <- start_analysis("usecase_trajectory_alignment")

# PART 3: Scores ----------------------------------------------------------
results <- read_rds(exp$result("results.rds"))  %>%
  separate_rows(id, alpha, scores, sep = " ", convert = T) %>%
  mutate(noise = factor(alpha))

g <-
  results %>%
  group_by(method, noise) %>%
  summarise_at(vars(scores), list(
    min = ~quantile(., .05),
    lower = ~quantile(., .25),
    mean = ~mean(.),
    upper = ~quantile(., .75),
    max = ~quantile(., .95)
  )) %>%
  ggplot() +
  geom_boxplot(
    aes(noise, ymin = min, lower = lower, middle = mean, upper = upper, max = max, fill = method),
    stat = "identity", width = 0.5, size = 0.45
  ) +
  theme_bw() +
  theme_common() +
  labs(x = "Amount of added noise", y = "Distance (lower is better)", fill = "Processing method") +
  scale_fill_brewer(palette = "Set2")
g

# PART 2: dtw heatmaps ----------------------------------------------------
d1 <- readRDS(exp$dataset_file("linear1_1_0.5"))
d2 <- readRDS(exp$dataset_file("linear1_2_0.5"))

res1 <- get_cell_expression(d1, d1$milestone_network, "sA")
res2 <- get_cell_expression(d2, d2$milestone_network, "sA")

res1_sm <- get_waypoint_expression(d1, 100)
res2_sm <- get_waypoint_expression(d2, 100)

pt1 <- res1$pseudotime
pt2 <- res2$pseudotime
expr1 <- res1$expression
expr2 <- res2$expression

alignment_original <- dtw(expr2, expr1, step.pattern=symmetric2, keep.internals=T)
ao <- plot_density(alignment_original, title = "DTW")

smp1 <- seq(from = 1, to = 1000, by = 10) #sample(1000, size = 100, replace = FALSE)
pt1_smp <- pt1[seq(from = 1, to = 1000, by = 10)]
pt2_smp <- pt2[seq(from = 1, to = 1000, by = 10)]
expr1_smp <- expr1[names(pt1_smp),]
expr2_smp <- expr2[names(pt2_smp),]

alignment_smooth <- dtw(res2_sm$expression, res1_sm$expression, step.pattern=symmetric2, keep.internals=T)
a_sm <- plot_density(alignment_smooth, title = "DTW+smoothing")


# PART 1: Explanation plot ------------------------------------------------
part1plots <- read_rds(exp$result("usecase_separateplots.rds"))

all_plots <- wrap_plots(
  wrap_plots(
    part1plots$t2_pt + labs(tag = "A", colour = "Healthy\nprogression   "),
    part1plots$t1_pt + labs(tag = "B", colour = "Diseased\nprogression   "),
    part1plots$prediction + theme(legend.position = "bottom") + labs(tag = "C")
  ),
  wrap_plots(
    part1plots$plot_dens + labs(tag = "D") + theme(legend.position = "none"),
    ao + theme(legend.position = "none") + labs(tag = "E"),
    a_sm + theme(legend.position = "none") + labs(tag = "F"),
    ggpubr::as_ggplot(ggpubr::get_legend(ao + theme(legend.position = "right") + labs(fill = "Alignment\ndistance"))),
    widths = c(1, 1, 1, .5)
  ),
  g + labs(tag = "G", y = "Alignment distance"),
  heights = c(.5, 1.4, 1.5)
)

ggsave(exp$result("usecase.pdf"), all_plots, height = 9, width = 10, useDingbats = FALSE)
ggsave(exp$result("usecase.png"), all_plots, height = 9, width = 10)
