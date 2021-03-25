library(tidyverse)
library(dyngen.manuscript)
library(patchwork)
library(dynplot)

exp <- start_analysis("usecase_trajectory_alignment")

# PART 1: Explanation plot ------------------------------------------------
part1plots <- read_rds(exp$result("usecase_separateplots.rds"))

# PART 2: heatmaps --------------------------------------------------------

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
dataPastel2 <- dataset %>% filter_cells(model == "left", keep_dimred = TRUE)
dataset2 <- dataset %>% filter_cells(model == "right", keep_dimred = TRUE)

cal_alignment <- ta_methods$cellAlign(dataPastel2, dataset2)
dtw_alignment <- ta_methods$DTW(dataPastel2, dataset2)

cal_dens <- plot_density(cal_alignment, title = NULL)
cal_dens
dtw_dens <- plot_density(dtw_alignment, title = NULL)
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


# PART 2bis: naupa explanation --------------------------------------------

# run method
paths <- bind_rows(
  tibble(pt1_aligned = cal_alignment$pt1_aligned, pt2_aligned = cal_alignment$pt2_aligned, method = "cellAlign"),
  tibble(pt1_aligned = dtw_alignment$pt1_aligned, pt2_aligned = dtw_alignment$pt2_aligned, method = "DTW"),
  tibble(pt1_aligned = c(0, 1, 1), pt2_aligned = c(0, 0, 1), method = "Worst possible"),
  tibble(pt1_aligned = c(0, 1), pt2_aligned = c(0, 1), method = "Best possible")
) %>%
  mutate(method = factor(method, levels = c("cellAlign", "DTW", "Best possible", "Worst possible")))
g_naupa <-
  ggplot(paths) +
  geom_path(aes(pt1_aligned + pt2_aligned, abs(pt1_aligned - pt2_aligned), colour = method)) +
  scale_colour_brewer(palette = "Set2") +
  theme_classic() +
  coord_equal() +
  labs(x = "pt1 + pt2", y = "|pt1 - pt2|", colour = "Method")




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

metric_labels <- c(score = "NAUPA score")

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

# alternative layout, to be discussed
galt <- wrap_plots(
  wrap_plots(
    part1plots$ground_truth + coord_flip() + theme(legend.position = "right") + labs(tag = "A"),
    part1plots$prediction + coord_flip() + theme(legend.position = "right") + labs(tag = "B"),
    nrow = 1,
    guides = "collect"
  ),
  wrap_plots(
    g_naupa + labs(tag = "C"),
    g_metrics$score + coord_flip() + expand_limits(y = 1.05) + scale_y_continuous(breaks = c(0, .25, .5, .75, 1)) + labs(tag = "D"),
    guides = "collect"
  ),
  heights = c(1, 1.3),
  ncol = 1
)
ggsave(exp$result("supp_fig.pdf"), galt, height = 6, width = 10, useDingbats = FALSE)

pdftools::pdf_convert(
  pdf = exp$result("supp_fig.pdf"),
  filenames = exp$result("supp_fig.png"),
  dpi = 120
)
