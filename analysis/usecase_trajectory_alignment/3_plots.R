library(tidyverse)
library(dyngen.manuscript)
library(patchwork)
library(dynplot2)

exp <- start_analysis("usecase_trajectory_alignment")


# PART 0: Load example dataset -------------------------------------------
# read dataset
design_dataset <- read_rds(exp$result("design_datasets.rds"))
file <- exp$dataset_file(design_dataset$id[[1]])
dataset <- read_rds(file) %>% dynwrap::add_dimred(dimred = dyndimred::dimred_landmark_mds)

# split up in two
dataset1 <- dataset %>% filter_cells(model == "left", keep_dimred = TRUE)
dataset2 <- dataset %>% filter_cells(model == "right", keep_dimred = TRUE)

# add healthy/diseased grouping
dataset$cell_info <- dataset$cell_info %>% mutate(
  group = c("left" = "Healthy", "right" = "Diseased")[model]
)
group_palette <- c(Diseased = "#fd8d3c", Healthy = "#6baed6")

# run methods
dtw_out <- ta_methods$DTW(dataset1, dataset2)


# PART 1: create explanation plots ----------------------------------------
# construct gold standard alignment
step_progressions <- crossing(
  dataset$milestone_network %>% select(from, to),
  percentage = seq(0, 1, by = .2)
) %>%
  mutate(cell_id = paste0("step", row_number()))
step_waypoints <- list(
  milestone_percentages = dynwrap::convert_progressions_to_milestone_percentages(
    cell_id = step_progressions$waypoint_id,
    milestone_ids = dataset$milestone_ids,
    milestone_network = dataset$milestone_network,
    progressions = step_progressions
  ) %>% rename(waypoint_id = cell_id)
)
step_proj <- dynwrap::project_waypoints(trajectory = dataset, space = dataset$dimred, waypoints = step_waypoints)

ds <- data.frame(step_progressions, step_proj[step_progressions$cell_id, ]) %>%
  as_tibble() %>%
  mutate(from_ = gsub("^(left|right)_", "", from), to_ = gsub("^(left|right)_", "", to))

# remove large matrices
dataset$counts <-
  dataset$counts_spliced <-
  dataset$counts_unspliced <-
  dataset$counts_protein <-
  dataset$expression <-
  NULL

# create plot
ground_truth <- dynplot_dimred(dataset) +
  geom_path(aes(comp_1, comp_2, group = paste(from_, to_, percentage, sep = "~")), linetype = "dashed", colour = "darkgray", ds) +
  geom_cell_point(colour = "black", size = 3) +
  geom_cell_point(aes(colour = group), size = 2.5) +
  geom_trajectory_segments(size = 2, arrow_size = .2) +
  geom_milestone_point(size = 4, colour = "black") +
  labs(colour = "Sample") +
  scale_colour_manual(values = group_palette) +
  theme_common()
ground_truth


ground_truth_small <- dynplot_dimred(dataset) +
  geom_path(aes(comp_1, comp_2, group = paste(from_, to_, percentage, sep = "~")), linetype = "dashed", colour = "darkgray", ds) +
  geom_cell_point(colour = "black", size = 2.5) +
  geom_cell_point(aes(colour = group), size = 2) +
  geom_trajectory_segments(size = 1, arrow_size = .2, grid::arrow(type = "closed", length = unit(0.06, "inches"))) +
  geom_milestone_point(size = 2, colour = "black") +
  labs(colour = "Sample") +
  scale_colour_manual(values = group_palette) +
  theme_common()

# run method
res1 <- get_cell_expression(dataset1)
res2 <- get_cell_expression(dataset2)

df <- tibble(
  i = dtw_out$index1,
  j = dtw_out$index2,
  ci = names(res1$pseudotime)[i],
  cj = names(res2$pseudotime)[j],
  pi = dtw_out$pt1_aligned,
  pj = dtw_out$pt2_aligned,
  rn = seq_along(i)
) %>%
  filter(rn %% round(n()/15) == 1)
df2 <- df %>%
  select(pos = rn, left = ci, right = cj) %>%
  gather(group, cell_id, left, right) %>%
  mutate(pred_id = paste0("pred_", cell_id))

pred_waypoints <- list(
  milestone_percentages =
    dataset$milestone_percentages %>%
    filter(cell_id %in% df2$cell_id) %>%
    mutate(cell_id = paste0("pred_", cell_id)) %>%
    rename(waypoint_id = cell_id)
)
pred_proj <- dynwrap::project_waypoints(trajectory = dataset, space = dataset$dimred, waypoints = pred_waypoints)

pred_df <- data.frame(df2, pred_proj[df2$pred_id, ]) %>%
  as_tibble()

prediction <- dynplot_dimred(dataset) +
  geom_path(aes(comp_1, comp_2, group = pos), linetype = "dashed", colour = "darkgray", pred_df) +
  geom_cell_point(colour = "black", size = 3) +
  geom_cell_point(aes(colour = group), size = 2.5) +
  geom_trajectory_segments(size = 2, arrow_size = .2) +
  geom_milestone_point(size = 4, colour = "black") +
  labs(colour = "Sample") +
  scale_colour_manual(values = group_palette) +
  theme_common()
prediction

prediction_small <- dynplot_dimred(dataset) +
  geom_path(aes(comp_1, comp_2, group = pos), linetype = "dashed", colour = "darkgray", pred_df) +
  geom_cell_point(colour = "black", size = 2.5) +
  geom_cell_point(aes(colour = group), size = 2) +
  geom_trajectory_segments(size = 1, arrow_size = .2, grid::arrow(type = "closed", length = unit(0.06, "inches"))) +
  geom_milestone_point(size = 2, colour = "black") +
  labs(colour = "Sample") +
  scale_colour_manual(values = group_palette) +
  theme_common()



# save plots for figure 1 and 2
plots <- lst(
  ground_truth,
  ground_truth_small,
  prediction,
  prediction_small
)

write_rds(plots, exp$result("usecase_separateplots.rds"), compress = "xz")


# PART 2: heatmaps --------------------------------------------------------

heatmap_paths <- dtw_out[c("index1", "index2")] %>% as_tibble()

heatmap_raster <-
  reshape2::melt(dtw_out$costMatrix, varnames = c("index1", "index2"), value.name = "distance") %>%
  as_tibble() %>%
  mutate(
    pt1 = res1$pseudotime[index1],
    pt2 = res2$pseudotime[index2],
  )

dtw_dens <-
  ggplot(mapping = aes(index1, index2)) +
  geom_raster(aes(fill = distance), heatmap_raster, alpha = .9) +
  geom_path(data = heatmap_paths) +
  scale_fill_distiller(palette = "RdYlGn", breaks = range(heatmap_raster$distance), labels = c("min", "max")) +
  scale_x_continuous(expand = c(0, 0), breaks = quantile(heatmap_paths$index1, c(0, .5, 1)), labels = c("0.0", "0.5", "1.0")) +
  scale_y_continuous(expand = c(0, 0), breaks = quantile(heatmap_paths$index2, c(0, .5, 1)), labels = c("0.0", "0.5", "1.0")) +
  theme_classic() +
  theme(
    axis.line = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    legend.position = "right"
  ) +
  labs(x = "Dataset 2 index", y = "Dataset 1 index", fill = "Accumulated\ndistance")

dtw_dens

# PART 3: ABWAP explanation --------------------------------------------
paths <- bind_rows(
  tibble(pt1_aligned = dtw_out$pt1_aligned, pt2_aligned = dtw_out$pt2_aligned, method = "Prediction"),
  tibble(pt1_aligned = c(0, 1, 1), pt2_aligned = c(0, 0, 1), method = "Worst"),
  tibble(pt1_aligned = c(0, 1), pt2_aligned = c(0, 1), method = "Best")
) %>%
  mutate(method = factor(method, levels = c("Prediction", "Best", "Worst")))

# safety check
paths %>% group_by(method) %>% summarise(abwap = ta_abwap(pt1_aligned, pt2_aligned))

fill <- bind_rows(
  tibble(pt1_aligned = dtw_out$pt1_aligned, pt2_aligned = dtw_out$pt2_aligned),
  tibble(pt1_aligned = c(1, 1, 0), pt2_aligned = c(1, 0, 0))
) %>% mutate(method = "Prediction")

text <- tribble(
  ~x, ~y, ~label, ~vjust, ~method,
  1.3, .025, "Best possible\nalignment", 0, "Best",
  1.4, 1, "Worst possible\nalignment", 1, "Worst",
  1.4, .75, paste0("Area Between Worst\nAnd Prediction =\nABWAP = ", round(ta_abwap(dtw_out$pt1_aligned, dtw_out$pt2_aligned), 2)), .5, "Prediction"
)
arrows <- tribble(
  ~x0, ~y0, ~x1, ~y1, ~method,
  1.25, .075, 1.1, 0.02, "Best",
  1.35, .95, 1.1, 0.98, "Worst",
  1.35, .75, 1.1, .5, "Prediction"
)

g_abwap <-
  ggplot() +
  geom_polygon(aes(pt1_aligned + pt2_aligned, abs(pt1_aligned - pt2_aligned)), fill, colour = NA, fill = "darkgray", alpha = .2) +
  geom_path(aes(pt1_aligned + pt2_aligned, abs(pt1_aligned - pt2_aligned), group = method), paths %>% filter(method != "Prediction"), linetype = "dashed") +
  geom_path(aes(pt1_aligned + pt2_aligned, abs(pt1_aligned - pt2_aligned), group = method), paths %>% filter(method == "Prediction"), linetype = "solid") +
  geom_segment(aes(x = x0, y = y0, xend = x1, yend = y1), arrows, arrow = grid::arrow(type = "closed", length = unit(0.06, "inches"))) +
  geom_text(aes(x, y, label = label, vjust = vjust), text, hjust = 0, size = 3) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "pt1 + pt2", y = "|pt1 - pt2|", colour = "Method", fill = "Method")
g_abwap

# PART 4: scores ----------------------------------------------------------
results <- read_rds(exp$result("results.rds")) %>%
  gather(metric, score, abwap) %>%
  mutate(
    method = factor(as.character(method), c("DTW", "cellAlign"))
  )

# preview
ggstatsplot::grouped_ggwithinstats(
  results,
  method,
  score,
  grouping.var = metric,
  type = "nonparametric"
)

# compute statistics
quants <- results %>%
  group_by(method, metric) %>%
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
      x = method,
      y = score,
      type = "nonparametric",
      paired = TRUE
    )
  }) %>%
  ungroup() %>%
  mutate(label = gsub("\\[Holm-corrected\\]", "", label))

write_rds(pairwise, exp$result("statistics.rds"), compress = "gz")

metric_labels <- c(abwap = "ABWAP score")

# create metric plots
g_metrics <- map2(names(metric_labels), metric_labels, function(metric, metric_name) {
  data <- results %>% filter(metric == !!metric)

  df_pairwise <- pairwise %>%
    filter(metric == !!metric) %>%
    mutate(groups = pmap(.l = list(group1, group2), .f = c)) %>%
    arrange(group1, group2)
  df_pairwise$stat_y_pos <- ggstatsplot:::ggsignif_xy(data$method, data$score)
  df_pairwise <- df_pairwise %>% filter(p.value < .05)

  g <- ggplot(data, aes(method, score)) +
    # geom_violin() +
    geom_violin(colour = NA, aes(fill = method), alpha = .3) +
    # geom_violin(aes(colour = cni_method_name, fill = cni_method_name), alpha = .2) +
    geom_boxplot(
      aes(x = method, y = NULL, ymin = min, lower = lower, middle = median, upper = upper, ymax = max),
      quants %>% filter(metric == !!metric),
      stat = "identity", width = 0.35, size = 0.45, fill = NA
    ) +
    geom_point(aes(colour = method), size = 1) +
    expand_limits(y = list(abwap = c(0,1))[[metric]]) +
    theme_classic() +
    theme_common(legend.position = "none") +
    labs(
      x = NULL,
      y = metric_name,
      colour = "Method"
    ) +
    scale_colour_brewer(palette = "Set2") +
    scale_fill_brewer(palette = "Set2")
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



# FINAL: Combine plots ----------------------------------------------------
g <- wrap_plots(
  wrap_plots(
    ground_truth_small + coord_flip() + theme(legend.position = "right") + labs(tag = "A"),
    prediction_small + coord_flip() + theme(legend.position = "right") + labs(tag = "B"),
    dtw_dens + theme(legend.position = "right") + labs(tag = "C"),
    ncol = 1,
    heights = c(1, 1, 1.5)
  ),
  wrap_plots(
    g_abwap + labs(tag = "D"),
    g_metrics$abwap + expand_limits(y = 1.05) + scale_y_continuous(breaks = c(0, .25, .5, .75, 1)) + labs(tag = "E"),
    ncol = 1
  ),
  nrow = 1,
  widths = c(1, 1.5)
)
ggsave(exp$result("supp_fig.pdf"), g, height = 7, width = 9, useDingbats = FALSE)


pdftools::pdf_convert(
  pdf = exp$result("supp_fig.pdf"),
  filenames = exp$result("supp_fig.png"),
  dpi = 120
)
