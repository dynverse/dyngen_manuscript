library(tidyverse)
library(dyngen.manuscript)
library(patchwork)
library(tidyr)
library(reshape2)
library(dyngen)
library(dynplot)
library(gtools)
library(dtw)

set.seed(1)
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
    simulation_params = simulation_default(census_interval = 10),
    verbose = TRUE,
    num_cores = 8
  ) %>%
    generate_tf_network() %>%
    generate_feature_network() %>%
    generate_kinetics() %>%
    generate_gold_standard()
  }
model1 <- runmodel(backbone) %>% normalise_goldstandard() %>% generate_cells()
model2D <- generate_diff_abundance(model1, c("sC", "sD"), c("sD", "sEndD"), c(0, 0))
model2D <- generate_cells(model2D)
modelD <- combine_models(model1, model2D)
datasetD <- wrap_dataset(modelD)
plot_dimred(datasetD)

dataset1 <- model1 %>% generate_experiment() %>% wrap_dataset()
dataset2D <- model2D %>% generate_experiment() %>% wrap_dataset()

## Save all 3 datasets
write_rds(datasetD, paste0(exp$dataset_folder("expplot"), "datasetD.rds"), compress = "gz")
write_rds(dataset1, paste0(exp$dataset_folder("expplot"), "dataset1.rds"), compress = "gz")
write_rds(dataset2D, paste0(exp$dataset_folder("expplot"), "dataset2D.rds"), compress = "gz")

milestone_colors <- c("#ff681c", "#ff681c", "#ff681c", "#ff681c", "#3abbba","#3abbba", "#3abbba", "#3abbba", "#ff681c", "#3abbba")
together <- plot_trajectory_in_color(datasetD, milestone_colors, c("sample 1", "sample 2"), c("#3abbba", "#ff681c"), plot_sd = 0.11)

t1_pt <- plot_pseudotime(dataset1, palette = "Oranges", trajectory_projection_sd = 0.11, plot_trajectory = T)
t2_pt <- plot_pseudotime(dataset2D, palette = "Blues", trajectory_projection_sd = 0.11, plot_trajectory = T)

expr1 <- get_cell_expression(dataset1, dataset1$milestone_network, "sA")$expression
expr2 <- get_cell_expression(dataset2D, dataset2D$milestone_network, "sA")$expression
a <- dtw(expr2, expr1, keep.internals = T, open.end = F)
dtwPlotAlignment(a)

plot_dens <- plot_density(a, show_legend = TRUE) + labs(title = NULL, subtitle = NULL)
plot_dens


e1 <- calculate_correct_pseudotime(dataset1, dataset1$milestone_network, "sA", normalized = F)
e2 <- calculate_correct_pseudotime(dataset2D, dataset2D$milestone_network, "sA", normalized = F)

# warped segmenten
traj1 <- data.frame(x = sort(e1))
traj2 <- data.frame(x = sort(e2))
traj1$color <- 1
traj2$color <- 2
traj1$y <- 0
traj2$y <- 1

traj1_warped <- traj1[a$index2,]
traj2_warped <- traj2[a$index1,]

combinatie <- traj1_warped
combinatie$x2 <- traj2_warped$x
combinatie$y2 <- traj2_warped$y
leftover_combinatie <- combinatie[seq(1, nrow(combinatie), 20), ]

alles <- data.frame(x1 = c(leftover_combinatie$x, leftover_combinatie$x2))
alles$y <- c(leftover_combinatie$y, leftover_combinatie$y2)
alles$color <- c(rep(1, nrow(leftover_combinatie)), rep(2, nrow(leftover_combinatie)))

cell_mappings <- ggplot() +
  geom_segment(data = leftover_combinatie, mapping = aes(x = x, y = y, xend = x2, yend = y2), alpha = 1) +
  geom_point(data = alles, mapping = aes(x = x1, y = y, colour = as.factor(color)), size = 3.5, show.legend = F) +
  scale_colour_manual(values = c("#fd8d3c", "#6baed6"), label = "") +
  scale_x_continuous(breaks = c(0, 1), labels = c("Start", "End")) +
  theme_void() +
  theme(
    axis.title = element_text(),
    axis.title.y = element_blank(),
    axis.line = element_line(),
    axis.line.y = element_blank(),
    axis.ticks = element_line(),
    axis.ticks.y = element_blank(),
    axis.text = element_text(),
    axis.text.y = element_blank()
  ) +
  labs(x = "Simulation time")

cell_mappings

write_rds(lst(leftover_combinatie, alles, combinatie), exp$result("cell_mappings_data.rds"), compress = "gz")


# For each 50th cell -> 26
# For each 30th cell -> 41
# For each 20th cell -> 63
demarcation_line <- 58

straight_part <- leftover_combinatie[1:demarcation_line,]
straight_part$x2 <- straight_part$x
partial_part <- leftover_combinatie[demarcation_line+1:nrow(leftover_combinatie),]
partial_part$x2 <- partial_part$x
partial_part$y2 <- 0.5
ground_truth <- rbind(straight_part, partial_part)

ground_truth_mappings <- ggplot() +
  geom_segment(data = ground_truth, mapping = aes(x = x, y = y, xend = x2, yend = y2), alpha = 1) +
  geom_point(data = alles, mapping = aes(x = x1, y = y, colour = as.factor(color)), size = 3.5, show.legend = F) +
  scale_colour_manual(values = c("#fd8d3c", "#6baed6"), label = "") +
  scale_x_continuous(breaks = c(0, 1), labels = c("Start", "End")) +
  theme_void() +
  theme(
    axis.title = element_text(),
    axis.title.y = element_blank(),
    axis.line = element_line(),
    axis.line.y = element_blank(),
    axis.ticks = element_line(),
    axis.ticks.y = element_blank(),
    axis.text = element_text(),
    axis.text.y = element_blank()
  ) +
  labs(x = "Simulation time")

ground_truth_mappings

write_rds(lst(ground_truth, alles, combinatie), exp$result("ground_truth_mappings_data.rds"), compress = "gz")
#
#
# # segm_traj <- traj1_warped
# # segm_traj$comp_3 <- traj2_warped$comp_1
# # segm_traj$comp_4 <- traj2_warped$comp_2
# # segm_traj$x_pt2 <- traj2_warped$x
# # segm_traj$x_pt1 <-segm_traj$x
#
# comb_traj1 <- comb_traj1[seq(1, nrow(comb_traj1), 20), ]
# comb_traj_part <- rbind(comb_traj1, comb_traj2)
#
#
# comp1t <- comb_traj_part$comp_1
# comp2t <- comb_traj_part$comp_2
# segm_traj_test <- segm_traj %>% filter(comp_1 %in% comp1t & comp_2 %in% comp2t & comp_3 %in% comp1t & comp_4 %in% comp2t)
# segm_traj_test$fake_val <- segm_traj_test$comp_1 - 0.5
# segm_traj_test$fake_x <- segm_traj_test$x_pt1
#
# st2 <- rbind(segm_traj_test, segm_traj_test)
# st2$fake_val
#
#
#
# cell_mappings <- ggplot() +
#   geom_segment(data = segments, mapping = aes(x = x, y = y1, xend = x2, yend = y2), alpha = 1) +
#   geom_point(data = comb_traj_less, mapping = aes(x = x, y = comp_1, colour = as.factor(color)), size = 3.5, show.legend = F) +
#   scale_colour_manual(values = c("#fd8d3c", "#6baed6"), label = "") +
#   scale_x_continuous(breaks = c(0, 1), labels = c("Start", "End")) +
#   theme_void() +
#   theme(
#     axis.title = element_text(),
#     axis.title.y = element_blank(),
#     axis.line = element_line(),
#     axis.line.y = element_blank(),
#     axis.ticks = element_line(),
#     axis.ticks.y = element_blank(),
#     axis.text = element_text(),
#     axis.text.y = element_blank()
#   ) +
#   labs(x = "Simulation time")
#
# cell_mappings
#
#
#
# ############
# #
# # ds <- datasetD
# #
# # expr2 <- get_cell_expression(ds, ds$milestone_network[1:4,], "left_sA")$expression
# # expr1 <- get_cell_expression(ds, ds$milestone_network[5:8,], "right_sA")$expression
# # align_normal <- dtw(expr2, expr1, keep.internals = T, open.end = F)
# # plota <- dtwPlotAlignment(align_normal)
# # plotd <- dtwPlotDensity(align_normal)
# #
# # plot_dens <- plot_density(align_normal, show_legend = TRUE) + labs(title = NULL, subtitle = NULL)
# # plot_dens
# #
# # # ggsave(plot_dens, filename="ts_dtw.png", bg="transparent", width = 4, height = 4)
# #
# # dataset <- ds
# # a <- align_normal
# # dr2 <- dynwrap::get_dimred(dataset, dyndimred::dimred_landmark_mds)
# # dr2 <- data.frame(dr2)
# #
# # e2 <- calculate_correct_pseudotime(dataset, dataset$milestone_network[1:4,], "left_sA")
# # e1 <- calculate_correct_pseudotime(dataset, dataset$milestone_network[5:8,], "right_sA")
# #
# # traj1 <- dr2[rownames(expr1),]
# # traj2 <- dr2[rownames(expr2),]
# # traj1$color <- 1
# # traj2$color <- 2
# # traj1$color2 <- sort(e1)
# # traj2$color2 <- sort(e2)
# #
# # traj2$comp_1 <- traj2$comp_1 + 1.1
# # comb_traj <- rbind(traj1, traj2)
# # comb_traj$ID <- seq.int(nrow(comb_traj))
# #
# # comb_traj1 <- comb_traj[1:nrow(expr1),] #%>% sample_frac(size = 0.7)
# # comb_traj2 <- comb_traj[nrow(expr1):1000,] #%>% sample_frac(size = 0.7)
# # comb_traj_less <- rbind(comb_traj1, comb_traj2)
# #
# # ggplot() + geom_point(data = comb_traj_less, mapping = aes(x = color2, y = comp_1, colour = as.factor(color)), size = 2) +
# #   scale_colour_manual(values = c("#3abbba", "#ff681c")) +
# #   theme_void()
# #
# # traj1_warped <- traj1[a$index2,]
# # traj2_warped <- traj2[a$index1,]
# #
# # segm_traj <- traj1_warped
# # segm_traj$comp_3 <- traj2_warped$comp_1
# # segm_traj$comp_4 <- traj2_warped$comp_2
# # segm_traj$color22 <- traj2_warped$color2
# # segm_traj$color <- as.factor(segm_traj$color)
# #
# # comp1t <- comb_traj_less$comp_1
# # comp2t <- comb_traj_less$comp_2
# # segm_traj_test <- segm_traj %>% filter(comp_1 %in% comp1t & comp_2 %in% comp2t & comp_3 %in% comp1t & comp_4 %in% comp2t)
# #
# # cellmappings <- ggplot() +
# #   geom_segment(data = segm_traj_test, mapping = aes(x = color2, y = comp_1, xend = color22, yend = comp_3), alpha = 0.25) +
# #   geom_point(data = comb_traj_less, mapping = aes(x = color2, y = comp_1, colour = as.factor(color)), size = 3.5, show.legend = F) +
# #   scale_colour_manual(values = c("#fd8d3c", "#6baed6"), label = "") +
# #   scale_x_continuous(breaks = c(0, 1), labels = c("Start", "End")) +
# #   theme_void() +
# #   theme(
# #     axis.title = element_text(),
# #     axis.title.y = element_blank(),
# #     axis.line = element_line(),
# #     axis.line.y = element_blank(),
# #     axis.ticks = element_line(),
# #     axis.ticks.y = element_blank(),
# #     axis.text = element_text(),
# #     axis.text.y = element_blank()
# #   ) +
# #   labs(x = "Simulation time")
# ###############
#
# write_rds(lst(segm_traj_test, comb_traj_less), exp$result("explanation_plot_data.rds"), compress = "gz")
# cellmappings


part1 <-
  patchwork::wrap_plots(
    t2_pt,
    t1_pt,
    plot_spacer(),
    cell_mappings,
    plot_dens,
    widths = c(1, 1, .2, 1, 1)
  )

saveRDS(part1, file = exp$result("explanation_flat.rds"))
ggsave(part1, filename = exp$result("explanation_flat.png"), bg = 'transparent', width = 20, height = 4)
