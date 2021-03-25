library(tidyverse)
library(dyngen.manuscript)
library(patchwork)
library(tidyr)
library(reshape2)
library(dyngen)
library(dynplot)
library(gtools)
library(dtw)
library(ggstatsplot)

exp <- start_analysis("usecase_trajectory_alignment")

# PART 3: Scores ----------------------------------------------------------
results <- read_rds(exp$result("results.rds")) %>%
  mutate(
    method = factor(as.character(method), c("DTW", "cellAlign"))
  )

g <-
  results %>%
  group_by(method) %>%
  summarise_at(vars(distance), list(
    min = ~quantile(., .05),
    lower = ~quantile(., .25),
    mean = ~median(.),
    upper = ~quantile(., .75),
    max = ~quantile(., .95)
  )) %>%
  ggplot() +
  geom_boxplot(
    aes(1, ymin = min, lower = lower, middle = mean, upper = upper, max = max, fill = method),
    stat = "identity", width = 0.5, size = 0.45
  ) +
  theme_bw() +
  theme_common() +
  labs(x = "Amount of added noise", y = "Distance (lower is better)", fill = "Processing method") +
  scale_fill_brewer(palette = "Set2")
g

# normalized minimum global distance computed
res_plot <- ggstatsplot::ggwithinstats(
  results,
  x = method,
  y = distance,
  type = "np",
  pairwise.comparisons = TRUE,
  pairwise.display = TRUE
)
res_plot

# results %>%
#   ggplot() +
#   geom_path(
#     aes(noise, score, colour = method, group = paste0(method, "_", seed))
#   ) +
#   theme_bw() +
#   theme_common() +
#   labs(x = "Amount of added noise", y = "Distance (lower is better)", fill = "Processing method") +
#   scale_colour_brewer(palette = "Set2")
#
# results %>% spread(method, score) %>%
#   ggplot() + geom_point(aes(DTW, `DTW+smoothing`)) + geom_abline(intercept = 0, slope = 1)
# results %>% spread(method, score) %>%
#   ggplot() + geom_point(aes(`DTW+smoothing`, `cellAlign`)) + geom_abline(intercept = 0, slope = 1)

# PART 2: dtw heatmaps
# read dataset
design_dataset <- read_rds(exp$result("design_datasets.rds"))
file <- exp$dataset_file(design_dataset$id[[1]])
dataset <- read_rds(file) %>% dynwrap::add_dimred(dimred = dyndimred::dimred_landmark_mds)

# split up in two
dataset1 <- dataset %>% filter_cells(model == "left", keep_dimred = TRUE)
dataset2 <- dataset %>% filter_cells(model == "right", keep_dimred = TRUE)

cal_alignment <- ta_methods$cellAlign(dataset1, dataset2)
dtw_alignment <- ta_methods$DTW(dataset1, dataset2)

cal_dens <- plot_density(cal_alignment)
cal_dens
dtw_dens <- plot_density(dtw_alignment)
dtw_dens

# PART 2: dtw heatmaps ----------------------------------------------------
# d1 <- readRDS(exp$dataset_file("linear1_1_0.5"))
# d2 <- readRDS(exp$dataset_file("linear1_2_0.5"))
#
# res1 <- get_cell_expression(d1)
# res2 <- get_cell_expressio#n(d2)
#
# res1_sm <- get_waypoint_expression(d1, 100)
# res2_sm <- get_waypoint_expression(d2, 100)
#
# pt1 <- res1$pseudotime
# pt2 <- res2$pseudotime
# expr1 <- res1$expression
# expr2 <- res2$expression
#
# alignment_original <- dtw::dtw(expr2, expr1, step.pattern = dtw::symmetric2, keep.internals = TRUE)
# ao <- plot_density(alignment_original, title = "DTW")
#
# smp1 <- seq(from = 1, to = 1000, by = 10) #sample(1000, size = 100, replace = FALSE)
# pt1_smp <- pt1[seq(from = 1, to = 1000, by = 10)]
# pt2_smp <- pt2[seq(from = 1, to = 1000, by = 10)]
# expr1_smp <- expr1[names(pt1_smp),]
# expr2_smp <- expr2[names(pt2_smp),]
#
# alignment_smooth <- dtw::dtw(res2_sm$expression, res1_sm$expression, step.pattern = dtw::symmetric2, keep.internals = TRUE)
# a_sm <- plot_density(alignment_smooth, title = "DTW+smoothing")
t1 <- (plot_spacer() + ground_truth + plot_spacer() + prediction + plot_spacer() + plot_layout(nrow = 1, widths = c(1, 2, 1, 2, 1))) / (dtw_dens + cal_dens + results)
t1
t2 <- (ground_truth + prediction) / (dtw_dens + results)
t2
ggsave(exp$result("usecase1.pdf"), t1, height = 7, width = 9)

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
