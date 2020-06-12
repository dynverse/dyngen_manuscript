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
results <- read_rds(exp$result("results.rds")) %>%
  mutate(noise = factor(noise))

g <-
  results %>%
  group_by(method, noise) %>%
  summarise_at(vars(score), list(
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
ao <- plot_density(alignment_original) + labs(title = "Original cells", subtitle = NULL)

smp1 <- seq(from = 1, to = 1000, by = 10) #sample(1000, size = 100, replace = FALSE)
pt1_smp <- pt1[seq(from = 1, to = 1000, by = 10)]
pt2_smp <- pt2[seq(from = 1, to = 1000, by = 10)]
expr1_smp <- expr1[names(pt1_smp),]
expr2_smp <- expr2[names(pt2_smp),]

alignment_subsample <- dtw(expr2_smp, expr1_smp, step.pattern=symmetric2, keep.internals=T)
asubs <- plot_density(alignment_subsample) + labs(title = "Subsampled", subtitle = NULL)

alignment_smooth <- dtw(res2_sm$expression, res1_sm$expression, step.pattern=symmetric2, keep.internals=T)
a_sm <- plot_density(alignment_smooth) + labs(title = "Smoothed", subtitle = NULL)


# PART 1: Explanation plot ------------------------------------------------

part1 <- readRDS(exp$result("explanation_flat.rds"))
part2 <- (ao + plot_spacer() + a_sm) + plot_layout(widths = c(0.35, 0.3, 0.35))
all_plots <- part1 / part2 / g
all_plots <- all_plots + plot_annotation(tag_levels = c('A')) + plot_layout(heights = c(0.5, 1.5, 2))

ggsave(exp$result("usecase.pdf"), all_plots, height = 11, width = 11, useDingbats = FALSE)
ggsave(exp$result("usecase.png"), all_plots, height = 11, width = 11)
