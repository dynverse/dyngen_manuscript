library(tidyverse)
library(dyngen.manuscript)
library(dynplot2)
library(rgl)

exp <- start_analysis("usecase_rna_velocity")

design_velocity <- read_rds(exp$result("design_velocity.rds"))

design_names <-
  tribble(
    ~method_id, ~params_id, ~method_label,
    "velocyto", "constant_velocity", "velocyto",
    "scvelo", "stochastic", "scvelo stochastic",
    "scvelo", "dynamical", "scvelo dynamical"
  ) %>%
  mutate(method_label = forcats::fct_inorder(method_label))

# PART A: Illustration of velocity ----------------------------------------
dataset_id <- "bifurcating_seed3"

dataset <- read_rds(exp$dataset_file(dataset_id))
model <- read_rds(exp$model_file(dataset_id))
dyngen::plot_backbone_modulenet(model)
# # choosing a good rotation from 3D to 2D
# mil_cols <- dynplot2::define_milestone_colors(milestone_colors = NULL, milestone_ids = dataset$milestone_ids)
# mil_pct_mat <- dataset$milestone_percentages %>% reshape2::acast(cell_id~milestone_id, value.var = "percentage", fill = 0) %>% .[rownames(dataset$dimred), dataset$milestone_ids]
# mil_rgb <- mil_pct_mat %*% mil_cols
# cols <- rgb(mil_rgb[,1], mil_rgb[,2], mil_rgb[,3], maxColorValue = 256)
# rgl::plot3d(dataset$dimred, col = cols, size = 10)
# userMatrix <- rgl::par3d()$userMatrix[1:3,1:3]
# userMatrix %>% as.vector %>% paste0(collapse = ", ")

## apply previous rotation
if (dataset_id == "bifurcating_seed3") {
  userVec <- c(-0.0365856438875198, -0.96261727809906, -0.268382370471954, -0.998166739940643, 0.0222430042922497, 0.0562904588878155, -0.0482163727283478, 0.269949942827225, -0.961666285991669)
  feature_oi <- "B5_TF1"
}
userMatrix <- matrix(userVec, nrow = 3)

dimred <- dataset$dimred %*% t(userMatrix)
colnames(dimred) <- colnames(dataset$dimred)
dataset <- dataset %>% dynwrap::add_dimred(dimred = dimred)

# Plot 1, trajectory
plot_trajectory <- {dynplot_dimred(dataset) +
    geom_cell_point(aes(color = milestone_percentages), size = 1) +
    scale_milestones_colour() +
    geom_trajectory_segments(size = 1, color = "#333333") +
    geom_milestone_label(aes(label = gsub("^s", "", label)), color = "black", fill = "#EEEEEE") +
    theme_common(legend.position = "none") +
    ggtitle("Trajectory")} %>%
  reduce_size()

# Plot 2, expression of a gene that goes up and down
expression_plot <- {dynplot_dimred(dataset) +
    geom_cell_point(aes(color = select_feature_expression(feature_oi, .data)), size = 1) +
    geom_trajectory_segments(size = 1, color = "#333333") +
    scale_expression_color(breaks = c(0, 1), labels = c("min", "max")) +
    theme_common() +
    ggtitle(paste0("Expression of gene ", feature_oi))}# %>%
# reduce_size()

# Plot 3, ground truth velocity
transform_groundtruth_velocity <- function(x) {
  scales::squish(x, c(-1, 1), only.finite = FALSE)
}
RdGyBu <- RColorBrewer::brewer.pal(9, "RdBu")
RdGyBu <- c(RdGyBu[1:3], "lightgray", RdGyBu[7:9])

gs_plot <- {dynplot_dimred(dataset) +
    geom_cell_point(aes(color = transform_groundtruth_velocity(dataset$rna_velocity[,feature_oi])), size = 1) +
    scale_colour_gradientn(colours = RdGyBu, limits = c(-1, 1), name = "Velocity", guide = dynplot2:::common_colorbar_legend) +
    ggtitle(paste0("Ground truth velocity of ", feature_oi)) +
    theme_common()} %>%
  reduce_size()


# Combine plots
plot_part_A <- patchwork::wrap_plots(plot_trajectory, expression_plot, gs_plot, nrow = 1, guides = "collect")


# PART B: RNA velocity estimates of different methods ---------------------
design_velocity_oi <- design_velocity %>%
  filter(dataset_id == !!dataset_id) %>%
  inner_join(design_names, by = c("method_id", "params_id")) %>%
  arrange(method_label)

#' design_velocity_oi %>% dynutils::extract_row_to_list(1) %>% list2env(.GlobalEnv)
plot_part_B <- pmap(design_velocity_oi, function(dataset_id, method_id, params_id, method_label, ...) {
  velocity <- read_rds(exp$velocity_file(dataset_id, method_id, params_id))
  velo_vec <- velocity$velocity_vector[,feature_oi]
  maxv <- quantile(abs(velo_vec), .8)
  g <- dynplot_dimred(dataset) +
    geom_cell_point(aes(color = transform_groundtruth_velocity(velo_vec / maxv)), size = 1) +
    scale_colour_gradientn(colours = RdGyBu, limits = c(-1, 1), name = "Velocity", guide = "none") +
    ggtitle(method_label) +
    theme_common()

  rm(velocity)

  reduce_size(g)
}) %>% patchwork::wrap_plots(nrow = 1)


# PART C: Embedded RNA velocity estimates of different methods ------------
plot_part_C <- pmap(design_velocity_oi, function(dataset_id, method_id, params_id, params, method_label, ...) {
  velocity_file <- exp$velocity_file(dataset_id, method_id, params_id)
  velocity <- read_rds(velocity_file)
  if (method_id == "scvelo") {
    velocity$scvelo <- reticulate::py_load_object(paste0(dirname(velocity_file), "/scvelo.pkl"))
  }

  dataset2 <- dataset %>%
    scvelo::add_velocity(velocity = velocity)
  dataset2 <- dataset2 %>%
    scvelo::add_dimred_future()

  g <- {dynplot_dimred(dataset2) +
      geom_cell_point(aes(color = milestone_percentages), size = 1) +
      scale_milestones_colour() +
      geom_velocity_stream(
        size = .8,
        color = "#333333",
        stat = stat_velocity_stream(grid_bandwidth = 1, filter = rlang::quo(mass > max(mass) * 0.07)),
        arrow = arrow(length = unit(0.3, "cm"), type = "closed")
      ) +
      ggtitle(method_label) +
      theme_common() +
      theme(legend.position = "none")} %>%
    reduce_size()

  rm(velocity, dataset2)
  g
}) %>%
  patchwork::wrap_plots(nrow = 1)

# PART D: Average scores --------------------------------------------------
scores <- read_rds(exp$result("scores.rds"))

metric_labels <- c(cor = "Velocity correlation", mean_cosine = "Velocity arrow cosine")

# create boxplots
summ <- scores$summ %>%
  inner_join(design_names, by = c("method_id", "params_id"))
results <- summ %>%
  gather(metric, score, cor, mean_cosine) %>%
  mutate(
    metric_label = factor(metric_labels[metric], levels = metric_labels)
  )

# preview
ggstatsplot::grouped_ggwithinstats(
  results,
  method_label,
  score,
  grouping.var = metric_label,
  type = "nonparametric"
)

# compute statistics
quants <- results %>%
  group_by(method_label, metric) %>%
  summarise_at(vars(score), list(
    min = ~quantile(., 0),
    lower = ~quantile(., .25),
    median = ~median(.),
    upper = ~quantile(., .75),
    max = ~quantile(., 1)
  ))

pairwise <-
  results %>%
  group_by(metric) %>%
  do({
    pairwiseComparisons::pairwise_comparisons(
      data = .,
      x = method_label,
      y = score,
      type = "nonparametric",
      paired = TRUE
    )
  }) %>%
  ungroup() %>%
  mutate(label = gsub("\\[Holm-corrected\\]", "", label))

# create metric plots
g_metrics <- map2(names(metric_labels), metric_labels, function(metric, metric_name) {
  data <- results %>% filter(metric == !!metric)

  df_pairwise <- pairwise %>%
    filter(metric == !!metric) %>%
    mutate(groups = pmap(.l = list(group1, group2), .f = c)) %>%
    arrange(group1, group2)
  df_pairwise$stat_y_pos <- ggstatsplot:::ggsignif_xy(data$method_label, data$score)
  df_pairwise <- df_pairwise %>% filter(p.value < .05)

  g <- ggplot(data, aes(method_label, score)) +
    # geom_violin() +
    geom_violin(colour = NA, aes(fill = method_label), alpha = .3) +
    # geom_violin(aes(colour = cni_method_name, fill = cni_method_name), alpha = .2) +
    geom_boxplot(
      aes(x = method_label, y = NULL, ymin = min, lower = lower, middle = median, upper = upper, ymax = max),
      quants %>% filter(metric == !!metric),
      stat = "identity", width = 0.35, size = 0.45, fill = NA
    ) +
    geom_point(aes(colour = method_label), size = 1) +
    # expand_limits(y = c(auroc = .95, aupr = 0.26)[[metric]]) +
    theme_classic() +
    theme_common(legend.position = "none") +
    theme(axis.text.x = element_text(angle = 20, hjust = 1)) +
    labs(
      x = NULL,
      y = metric_name,
      colour = "Method"
    ) +
    scale_colour_brewer(palette = "Set1") +
    scale_fill_brewer(palette = "Set1")
  if (nrow(df_pairwise) > 0) {
    g <- g + ggsignif::geom_signif(
      comparisons = df_pairwise$groups,
      map_signif_level = TRUE,
      y_position = df_pairwise$stat_y_pos,
      annotations = df_pairwise$label,
      test = NULL,
      parse = TRUE,
      textsize = 2.5
    )
  }
  g
})
names(g_metrics) <- names(metric_labels)
g_metrics$cor
g_metrics$mean_cosine <- g_metrics$mean_cosine + scale_y_continuous(breaks = c(0, .25, .5, .75, 1))


# create pairwise plots
g_pairwise <-
  ggplot(summ) +
  geom_point(aes(cor, mean_cosine, colour = method_label), size = 1) +
  theme_classic() +
  theme_common(legend.position = "right") +
  labs(x = metric_labels[["cor"]], y = metric_labels[["mean_cosine"]], colour = "Method") +
  scale_color_brewer(palette = "Set1")


plot_part_D <- patchwork::wrap_plots(
  g_metrics$cor,
  g_metrics$mean_cosine,
  g_pairwise + theme(axis.title.x = element_text(margin = margin(t = -30, unit = "pt"))), # perform some manual alignment
  nrow = 1,
  guides = "collect"
)



# COMBINE ALL PARTS -------------------------------------------------------
plot_part_A[[1]] <- plot_part_A[[1]] + labs(tag = "A")
plot_part_B[[1]] <- plot_part_B[[1]] + labs(tag = "B")
plot_part_C[[1]] <- plot_part_C[[1]] + labs(tag = "C")
plot_part_D[[1]] <- plot_part_D[[1]] + labs(tag = "D")

g <- patchwork::wrap_plots(
  plot_part_A,
  plot_part_B,
  (plot_part_C & theme(plot.title = element_blank(), plot.subtitle = element_blank())),
  plot_part_D,
  ncol = 1,
  heights = c(1, 1, 1, 1)
)
ggsave(exp$result("supp_fig.pdf"), g, height = 12, width = 12, useDingbats = FALSE)
# convert for previewing in gdoc
pdftools::pdf_convert(
  pdf = exp$result("supp_fig.pdf"),
  filenames = exp$result("supp_fig.png"),
  dpi = 120
)

# ggsave(exp$result("usecase.png"), g, height = 12, width = 12)

write_rds(list(groundtruth = g[[1]][[1]], prediction = g[[3]][[1]]), exp$result("usecase_separateplots.rds"), compress = "gz")



