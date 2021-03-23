library(tidyverse)
library(dyngen.manuscript)

exp <- start_analysis("usecase_network_inference")

# load example plots and data from previous script
g_example <- read_rds(exp$result("casewise_grn.rds"))
list2env(read_rds(exp$result("scores.rds")), .GlobalEnv)


# create boxplots
results <- summ %>%
  filter(method == "casewise_casewise") %>%
  gather(metric, score, auroc, aupr) %>%
  mutate(metric = factor(metric, levels = c("auroc", "aupr")))

# preview
ggstatsplot::grouped_ggwithinstats(
  results,
  cni_method_name,
  score,
  grouping.var = metric,
  type = "np"
)

# compute statistics
quants <- results %>%
  group_by(cni_method_id, cni_method_name, metric) %>%
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
      x = cni_method_name,
      y = score,
      type = "np",
      paired = TRUE
    )
  }) %>%
  ungroup() %>%
  mutate(label = gsub("\\[Holm-corrected\\]", "", label))

# create metric plots
g_metrics <- map2(c("auroc", "aupr"), c("mean AUROC", "mean AUPR"), function(metric, metric_name) {
  data <- results %>% filter(metric == !!metric)

  df_pairwise <- pairwise %>%
    filter(metric == !!metric) %>%
    mutate(groups = pmap(.l = list(group1, group2), .f = c)) %>%
    filter(p.value < 0.05) %>%
    arrange(group1, group2)
  stat_y_pos <- ggstatsplot:::ggsignif_xy(data$cni_method_name, data$score)

  ggplot(data, aes(cni_method_name, score)) +
    # geom_violin() +
    geom_violin(colour = NA, aes(fill = cni_method_name), alpha = .3) +
    # geom_violin(aes(colour = cni_method_name, fill = cni_method_name), alpha = .2) +
    geom_boxplot(
      aes(x = cni_method_name, y = NULL, ymin = min, lower = lower, middle = median, upper = upper, ymax = max),
      quants %>% filter(metric == !!metric),
      stat = "identity", width = 0.35, size = 0.45, fill = NA
    ) +
    geom_point(aes(colour = cni_method_name), size = 1) +
    ggsignif::geom_signif(
      comparisons = df_pairwise$groups,
      map_signif_level = TRUE,
      y_position = stat_y_pos,
      annotations = df_pairwise$label,
      test = NULL,
      parse = TRUE,
      textsize = 2.5
    ) +
    expand_limits(y = c(auroc = .95, aupr = 0.26)[[metric]]) +
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
})
names(g_metrics) <- c("auroc", "aupr")
g_metrics$aupr


# create pairwise plots
g_pairwise <-
  ggplot(summ %>% filter(method == "casewise_casewise")) +
  geom_point(aes(auroc, aupr, colour = cni_method_name), size = 1) +
  theme_classic() +
  theme_common(legend.position = "right") +
  labs(x = "mean AUROC", y = "mean AUPR", colour = "Method") +
  scale_color_brewer(palette = "Set1")



grnh <- 3.5
bch <- 3.5
g <- patchwork::wrap_plots(
  g_example + labs(tag = "A"),
  patchwork::wrap_plots(
    g_metrics$auroc + labs(tag = "B"),
    g_metrics$aupr,
    g_pairwise + theme(axis.title.x = element_text(margin = margin(t = -35, unit = "pt"))), # perform some manual alignment
    nrow = 1,
    guides = "collect"
  ),
  heights = c(grnh, bch),
  ncol = 1
) & theme(plot.tag.position = c(0,1))

# g
ggsave(exp$result("supp_fig.pdf"), g, width = 10, height = (grnh + bch) * .9)

