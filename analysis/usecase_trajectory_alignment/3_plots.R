library(tidyverse)
library(dyngen.manuscript)
library(patchwork)
library(tidyr)
library(reshape2)
library(dyngen)
library(dynplot)
library(gtools)

exp <- start_analysis("usecase_trajectory_alignment")

result_smoothing <- read_rds(exp$result("result_smoothing.rds")) %>%
  mutate(
    noise = ordered(as.factor(noise)),
    smooth = factor(smooth, levels = c("smoothed", "original cells", "subsampled")),
    score_scaled = ifelse(smooth == "original cells", score / 10, score)
  )

g <- ggplot(data = result_smoothing, aes(noise, score_scaled, fill = smooth)) +
  geom_boxplot(width = 0.5, size=0.45,outlier.size=0.5) +
  theme_bw() +
  theme_common() +
  labs(x = "Amount of added noise", y = "Distance (lower is better)", fill = "Processing method")

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

# Plot explanation plot
get_explanation_plot <- function(save = F){
  set.seed(42)
  exp <- start_analysis("usecase_trajectory_alignment")

  backbone <- bblego(
    bblego_start("A", type = "simple", num_modules = 4),
    bblego_linear("A", "B", type = "simple", num_modules = 6),
    bblego_linear("B", "C", type = "simple", num_modules = 6),
    bblego_linear("C", "D", type = "simple", num_modules = 6),
    bblego_end("D", type = "simple", num_modules = 4)
  )

  runmodel <- function(backbone, num_tfs = 50){
    model <- initialise_model(
      num_tfs = num_tfs,
      num_targets = 250,
      num_hks = 200,
      num_cells = 1000,
      backbone = backbone,
      verbose = TRUE,
      num_cores = 8,
    ) %>%
      generate_tf_network() %>%
      generate_feature_network() %>%
      generate_kinetics() %>%
      generate_gold_standard() %>%
      generate_cells()
  }
  model1 <- runmodel(backbone) %>% normalise_goldstandard() %>% generate_cells()
  model2D <- generate_diff_abundance(model1, c("sC", "sD"), c("sD", "sEndD"), c(0, 0))
  model2D <- generate_cells(model2D)
  modelD <- combine_models(model1, model2D)
  datasetD <- wrap_dataset(modelD)
  plot_dimred(datasetD)

  milestone_colors <- c("#ff681c", "#ff681c", "#ff681c", "#ff681c", "#3abbba","#3abbba", "#3abbba", "#3abbba", "#ff681c", "#3abbba")
  together <- plot_trajectory_in_color(datasetD, milestone_colors, c("sample 1", "sample 2"), c("#3abbba", "#ff681c"), plot_sd = 0.11)

  dataset1 <- model1 %>% generate_experiment() %>% wrap_dataset()
  dataset2D <- model2D %>% generate_experiment() %>% wrap_dataset()

  t1_pt <- plot_pseudotime(dataset1, high = "#ff681c", low = "#66290B", trajectory_projection_sd = 0.11, plot_trajectory = T)
  t2_pt <- plot_pseudotime(dataset2D, trajectory_projection_sd = 0.11, plot_trajectory = T)

  library(dtw)
  ds <- datasetD

  expr2 <- get_cell_expression(ds, ds$milestone_network[1:4,], "left_sA")$expression
  expr1 <- get_cell_expression(ds, ds$milestone_network[5:8,], "right_sA")$expression
  align_normal <- dtw(expr2, expr1, keep.internals = T, open.end = F)
  plota <- dtwPlotAlignment(align_normal)
  plotd <- dtwPlotDensity(align_normal)

  plot_dens <- plot_density(align_normal, title = "", show_legend = T)
  plot_dens

  # ggsave(plot_dens, filename="ts_dtw.png", bg="transparent", width = 4, height = 4)

  dataset <- ds
  a <- align_normal
  dr2 <- dynwrap::get_dimred(dataset, dyndimred::dimred_landmark_mds)
  dr2 <- data.frame(dr2)

  e2 <- calculate_correct_pseudotime(dataset, dataset$milestone_network[1:4,], "left_sA")
  e1 <- calculate_correct_pseudotime(dataset, dataset$milestone_network[5:8,], "right_sA")

  traj1 <- dr2[rownames(expr1),]
  traj2 <- dr2[rownames(expr2),]
  traj1$color <- 1
  traj2$color <- 2
  traj1$color2 <- sort(e1)
  traj2$color2 <- sort(e2)

  traj2$comp_1 <- traj2$comp_1 + 1.1
  comb_traj <- rbind(traj1, traj2)
  comb_traj$ID <- seq.int(nrow(comb_traj))

  comb_traj1 <- comb_traj[1:nrow(expr1),] #%>% sample_frac(size = 0.7)
  comb_traj2 <- comb_traj[nrow(expr1):5000,] #%>% sample_frac(size = 0.7)
  comb_traj_less <- rbind(comb_traj1, comb_traj2)

  ggplot() + geom_point(data = comb_traj_less, mapping = aes(x = color2, y = comp_1, colour = as.factor(color)), size = 2) +
    scale_colour_manual(values = c("#3abbba", "#ff681c")) +
    theme_void()

  traj1_warped <- traj1[a$index2,]
  traj2_warped <- traj2[a$index1,]

  segm_traj <- traj1_warped
  segm_traj$comp_3 <- traj2_warped$comp_1
  segm_traj$comp_4 <- traj2_warped$comp_2
  segm_traj$color22 <- traj2_warped$color2
  segm_traj$color <- as.factor(segm_traj$color)

  comp1t <- comb_traj_less$comp_1
  comp2t <- comb_traj_less$comp_2
  segm_traj_test <- segm_traj %>% filter(comp_1 %in% comp1t & comp_2 %in% comp2t & comp_3 %in% comp1t & comp_4 %in% comp2t)

  cellmappings <- ggplot() + geom_point(data = comb_traj_less, mapping = aes(x = color2, y = comp_1, colour = as.factor(color)), size = 3.5, show.legend = F) +
    scale_colour_manual(values = c("#ff681c", "#3abbba"), label = "") +
    geom_segment(data = segm_traj_test, mapping = aes(x = color2, y = comp_1, xend = color22, yend = comp_3), alpha = 0.25) +
    geom_segment(aes(x = 0, y = -0.5, xend = 1, yend = -0.5), arrow = arrow(length = unit(0.25, "cm"))) +
    annotate(geom = "text", color = "black", x = 0.5, y = -0.75, label = "Pseudotime", size = 4) +
    theme_void()
  cellmappings


  part1 <- t2_pt | t1_pt |plot_spacer()| cellmappings | plot_dens

  if(save){
    saveRDS(part1, file = exp$result("explanation_flat.rds"))
    ggsave(part1, filename = paste0(exp$result(), "/explanation_flat.png"), bg = 'transparent', width = 20, height = 4)
  }

  part1 <- part1 + plot_layout(widths = c(1,1,.2,1,1))
  part1
}

part1 <- get_explanation_plot(save = T) #readRDS(exp$result("explanation_flat.rds"))

all_plots <- part1 / (ao + asubs + a_sm) / g
all_plots <- all_plots + plot_annotation(tag_levels = c('A')) + plot_layout(heights = c(1, 1.1, 1.5))

ggsave(exp$result("usecase.pdf"), all_plots, height = 11, width = 11, useDingbats = FALSE)
ggsave(exp$result("usecase.png"), all_plots, height = 11, width = 11)
