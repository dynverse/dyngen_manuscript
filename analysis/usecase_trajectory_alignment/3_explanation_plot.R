library(tidyverse)
library(dyngen.manuscript)
library(patchwork)
library(tidyr)
library(reshape2)
library(dyngen)
library(dynplot)
library(gtools)
library(dtw)

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
    simulation_params = simulation_default(census_interval = 10),
    verbose = TRUE,
    num_cores = 8
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

t1_pt <- plot_pseudotime(dataset1, palette = "Oranges", trajectory_projection_sd = 0.11, plot_trajectory = T)
t2_pt <- plot_pseudotime(dataset2D, palette = "Blues", trajectory_projection_sd = 0.11, plot_trajectory = T)

ds <- datasetD

expr2 <- get_cell_expression(ds, ds$milestone_network[1:4,], "left_sA")$expression
expr1 <- get_cell_expression(ds, ds$milestone_network[5:8,], "right_sA")$expression
align_normal <- dtw(expr2, expr1, keep.internals = T, open.end = F)
plota <- dtwPlotAlignment(align_normal)
plotd <- dtwPlotDensity(align_normal)

plot_dens <- plot_density(align_normal, show_legend = TRUE) + labs(title = NULL, subtitle = NULL)
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

cellmappings <- ggplot() +
  geom_segment(data = segm_traj_test, mapping = aes(x = color2, y = comp_1, xend = color22, yend = comp_3), alpha = 0.25) +
  geom_point(data = comb_traj_less, mapping = aes(x = color2, y = comp_1, colour = as.factor(color)), size = 3.5, show.legend = F) +
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
  labs(x = "Pseudotime")
write_rds(lst(segm_traj_test, comb_traj_less), exp$result("explanation_plot_data.rds"), compress = "gz")
cellmappings


part1 <-
  patchwork::wrap_plots(
    t2_pt,
    t1_pt,
    plot_spacer(),
    cellmappings,
    plot_dens,
    widths = c(1, 1, .2, 1, 1)
  )

saveRDS(part1, file = exp$result("explanation_flat.rds"))
ggsave(part1, filename = exp$result("explanation_flat.png"), bg = 'transparent', width = 20, height = 4)
