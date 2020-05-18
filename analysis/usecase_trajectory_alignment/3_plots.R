library(tidyverse)
library(dyngen.manuscript)
library(patchwork)

exp <- start_analysis("usecase_trajectory_alignment")

result_smoothing <- read_rds(exp$result("result_smoothing.rds")) %>%
  mutate(
    noise = ordered(as.factor(noise)),
    smooth = factor(smooth, levels = c("smoothed", "original cells", "subsampled")),
    score_scaled = ifelse(smooth == "original cells", score / 10, score)
  )

g <-
  ggplot(data = result_smoothing, aes(noise, score_scaled, fill = smooth)) +
  geom_boxplot(width = 0.5, size=0.45,outlier.size=0.5) +
  theme_bw() +
  theme_common() +
  labs(x = "Amount of added noise", y = "Distance (lower is better)", fill = "Processing method")

# Plot the three different alignments of a single pair of trajectories
plot_density <- function(a1, title = "Distance matrix + warping path", subtitle = ""){
  x <- as_tibble(a1$index1)
  x$Y <- a1$index2

  gga_1 <- reshape2::melt(a1$costMatrix)
  p_heat1 <- ggplot(gga_1, aes(Var1, Var2, fill= value)) +
    geom_raster(show.legend = FALSE) +
    scale_fill_distiller(palette = "RdYlGn") +
    scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
    geom_line(data=x, aes(x=value, y=Y)) +
    theme_classic() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    labs(title = title, fill = "Accumulated\ndistance", x = "Sample 2", y = "Sample 1")

  p_heat1
}

d1 <- readRDS(exp$dataset_file("linear1_1_0.5"))
d2 <- readRDS(exp$dataset_file("linear1_2_0.5"))

res1 <- get_cell_expression(d1)
res2 <- get_cell_expression(d2)

res1_sm <- get_waypoint_expression(d1, 100)
res2_sm <- get_waypoint_expression(d2, 100)

pt1 <- res1$pseudotime
pt2 <- res2$pseudotime
expr1 <- res1$expression
expr2 <- res2$expression

alignment_original <- dtw(expr2, expr1, step.pattern=symmetric2, keep.internals=T)
ao <- plot_density(alignment_original, title = "Original cells")

smp1 <- seq(from = 1, to = 1000, by = 10) #sample(1000, size = 100, replace = FALSE)
pt1_smp <- pt1[seq(from = 1, to = 1000, by = 10)]
pt2_smp <- pt2[seq(from = 1, to = 1000, by = 10)]
expr1_smp <- expr1[names(pt1_smp),]
expr2_smp <- expr2[names(pt2_smp),]

alignment_subsample <- dtw(expr2_smp, expr1_smp, step.pattern=symmetric2, keep.internals=T)
asubs <- plot_density(alignment_subsample, title = "Subsampled")

alignment_smooth <- dtw(res2_sm$expression, res1_sm$expression, step.pattern=symmetric2, keep.internals=T)
a_sm <- plot_density(alignment_smooth, title = "Smoothed")

# Combine all plots together
part1 <- cowplot::ggdraw() +
  cowplot::draw_image(exp$result("dyngen_ta_plot_flat.svg"))

all_plots <- part1 / (ao + asubs + a_sm) / g
all_plots <- all_plots + plot_annotation(tag_levels = c('A')) + plot_layout(heights = c(1, 1.1, 1.5))

ggsave(exp$result("usecase.pdf"), all_plots, height = 11, width = 11, useDingbats = FALSE)
ggsave(exp$result("usecase.png"), all_plots, height = 11, width = 11)
