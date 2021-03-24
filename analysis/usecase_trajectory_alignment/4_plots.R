library(tidyverse)
library(dyngen.manuscript)
library(patchwork)
library(dynplot)

exp <- start_analysis("usecase_trajectory_alignment")

# PART 1: Explanation plot ------------------------------------------------
part1plots <- read_rds(exp$result("usecase_separateplots.rds"))

# PART 2: heatmaps --------------------------------------------------------

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


# PART 3: scores ----------------------------------------------------------
results <- read_rds(exp$result("results.rds")) %>%
  mutate(
    metric = "score",
    method = factor(as.character(method), c("DTW", "cellAlign"))
  )

# preview
ggstatsplot::ggwithinstats(
  results,
  method,
  score,
  type = "nonparametric"
)

# compute statistics
quants <- results %>%
  group_by(method) %>%
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

metric_labels <- c(score = "Score")

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
g_metrics$score
# g_metrics$cor
# g_metrics$mean_cosine <- g_metrics$mean_cosine + scale_y_continuous(breaks = c(0, .25, .5, .75, 1)) + expand_limits(y = 1.28)



# FINAL: Combine plots ----------------------------------------------------
g <- wrap_plots(
  wrap_plots(
    plot_spacer(),
    part1plots$ground_truth + theme(legend.position = "right") + labs(tag = "A"),
    plot_spacer(),
    part1plots$prediction + theme(legend.position = "right") + labs(tag = "B"),
    plot_spacer(),
    nrow = 1,
    widths = c(1, 2, 1, 2, 1),
    guides = "collect"
  ),
  wrap_plots(
    dtw_dens + theme(legend.position = "right", panel.background = element_blank()) + expand_limits(x = 1000 * 1.01) + labs(tag = "C", title = NULL),
    cal_dens + theme(legend.position = "right", panel.background = element_blank()) + expand_limits(x = 200 * 1.01) + labs(tag = "D", title = NULL),
    g_metrics$score + expand_limits(y = 1.05) + scale_y_continuous(breaks = c(0, .25, .5, .75, 1)) + labs(tag = "E"),
    guides = "collect"
  ),
  ncol = 1
)

ggsave(exp$result("supp_fig.pdf"), g, height = 6, width = 10, useDingbats = FALSE)
pdftools::pdf_convert(
  pdf = exp$result("supp_fig.pdf"),
  filenames = exp$result("supp_fig.png"),
  dpi = 120
)
