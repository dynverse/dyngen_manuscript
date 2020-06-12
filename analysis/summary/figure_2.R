library(tidyverse)
library(dyngen)
library(dyngen.manuscript)
library(patchwork)

set.seed(1)

exp <- start_analysis("summary")

#################
# TRAJECTORY ALIGNMENT
#################
exp_a <- start_analysis("usecase_trajectory_alignment")
plots_a <- read_rds(exp_a$result("usecase_separateplots.rds"))

ga1 <- plots_a$ground_truth +
  coord_cartesian() +
  theme_classic() +
  theme(
    legend.position = "none",
    text = element_text(family = "Helvetica"),
    axis.title = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank()
  ) +
  labs(title = NULL, subtitle = NULL, tag = NULL)

ga2 <- plots_a$prediction +
  coord_cartesian() +
  theme_classic() +
  theme(
    legend.position = "none",
    text = element_text(family = "Helvetica"),
    axis.title = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank()
  ) +
  labs(title = NULL, subtitle = NULL, tag = NULL)

scores_a <-
  read_rds(exp_a$result("results.rds")) %>%
  mutate(noise = factor(noise))

plota_data <-
  scores_a %>%
  group_by(method, noise) %>%
  summarise_at(vars(score), list(
    min = ~quantile(., .05, na.rm = TRUE),
    lower = ~quantile(., .25, na.rm = TRUE),
    mean = ~mean(., na.rm = TRUE),
    upper = ~quantile(., .75, na.rm = TRUE),
    max = ~quantile(., .95, na.rm = TRUE)
  )) %>%
  ungroup()

ga3 <-
  ggplot(plota_data) +
  geom_boxplot(
    aes(noise, ymin = min, lower = lower, middle = mean, upper = upper, max = max, fill = method, group = paste0(noise, method)),
    stat = "identity", width = 0.5, size = 0.45
  ) +
  theme_classic() +
  theme_common() +
  labs(x = "Noise levels", y = "Distance (lower is better)", colour = "Method") +
  scale_fill_brewer(palette = "Dark2") +
  theme(legend.position = "right")

#################
# RNA VELOCITY
#################
exp_b <- start_analysis("usecase_rna_velocity")
plots_b <- read_rds(exp_b$result("usecase_separateplots.rds"))
gb1 <- plots_b$groundtruth +
  coord_cartesian() +
  theme_classic() +
  theme(
    legend.position = "none",
    text = element_text(family = "Helvetica"),
    axis.title = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank()
  ) +
  labs(title = NULL, subtitle = NULL, tag = NULL)

gb2 <- plots_b$prediction +
  coord_cartesian() +
  theme_classic() +
  theme(
    legend.position = "none",
    text = element_text(family = "Helvetica"),
    axis.title = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank()
  ) +
  labs(title = NULL, subtitle = NULL, tag = NULL)

scores_b <- read_rds(exp_b$result("scores.rds"))

metric_labels_b <- c(cor = "Velocity correlation", mean_cosine = "Velocity arrow cosine")
plotb_data <-
  scores_b$summ %>%
  gather(metric, score, cor, mean_cosine) %>%
  mutate(
    metric_label = factor(metric_labels_b[metric], levels = metric_labels_b),
    group = paste0(method_id, " ", params_id)
  ) %>%
  group_by(group, method_id, params_id, metric_label) %>%
  summarise_at(vars(score), list(
    min = ~quantile(., .05, na.rm = TRUE),
    lower = ~quantile(., .25, na.rm = TRUE),
    mean = ~mean(., na.rm = TRUE),
    upper = ~quantile(., .75, na.rm = TRUE),
    max = ~quantile(., .95, na.rm = TRUE)
  ))

limits_tib <- tibble(
  metric_label = forcats::fct_inorder(metric_labels_b),
  x = 1,
  y0 = c(0, .75),
  y1 = c(1, 1)
) %>% gather(l, y, y0, y1)

gb3 <-
  patchwork::wrap_plots(
    list = map(metric_labels_b, function(met_lab) {
      g <- ggplot(plotb_data %>% filter(metric_label == met_lab)) +
        geom_boxplot(
          aes(forcats::fct_rev(group), ymin = min, lower = lower, middle = mean, upper = upper, max = max, fill = group),
          stat = "identity", width = 0.5, size = 0.45
        ) +
        geom_point(aes(x, y), limits_tib %>% filter(metric_label == met_lab), alpha = 0) +
        theme_classic() +
        theme_common() +
        coord_flip() +
        labs(x = "", y = met_lab, colour = "Method") +
        scale_fill_brewer(palette = "Dark2") +
        theme(legend.position = "none")
      if (met_lab != metric_labels_b[[1]]) {
        g <- g + theme(
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()
        )
      }
      g
    }),
    ncol = 2
  )

#################
# CSNI
#################
exp_c <- start_analysis("usecase_network_inference")
plots_c <- read_rds(exp_c$result("usecase_separateplots.rds"))
gc1 <- plots_c$groundtruth +
  coord_cartesian() +
  theme_classic() +
  labs(colour = "Regulation strength\nof cell 1") +
  theme(
    legend.position = "bottom",
    text = element_text(family = "Helvetica"),
    axis.title = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank()
  ) +
  labs(title = NULL, subtitle = NULL, tag = NULL)

gc2 <- plots_c$prediction +
  coord_cartesian() +
  theme_classic() +
  labs(colour = "Regulation strength\nof cell 1") +
  theme(
    legend.position = "bottom",
    text = element_text(family = "Helvetica"),
    axis.title = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank()
  ) +
  labs(title = NULL, subtitle = NULL, tag = NULL)

scores_c <- read_rds(exp_c$result("scores.rds"))

metric_labels_c <- c(auroc = "mean AUROC", aupr = "mean AUPR")
plotc_data <-
  scores_c$summ %>%
  filter(method == "casewise_casewise") %>%
  gather(metric, score, auroc, aupr) %>%
  mutate(
    metric_label = factor(metric_labels_c[metric], levels = metric_labels_c)
  ) %>%
  group_by(cni_method_name, metric, metric_label) %>%
  summarise_at(vars(score), list(
    min = ~quantile(., .05, na.rm = TRUE),
    lower = ~quantile(., .25, na.rm = TRUE),
    mean = ~mean(., na.rm = TRUE),
    upper = ~quantile(., .75, na.rm = TRUE),
    max = ~quantile(., .95, na.rm = TRUE)
  ))

gc3 <-
  patchwork::wrap_plots(
    list = map(metric_labels_c, function(met_lab) {
      g <- ggplot(plotc_data %>% filter(metric_label == met_lab)) +
        geom_boxplot(
          aes(forcats::fct_rev(cni_method_name), ymin = min, lower = lower, middle = mean, upper = upper, max = max, fill = cni_method_name),
          stat = "identity", width = 0.5, size = 0.45
        ) +
        theme_classic() +
        theme_common() +
        coord_flip() +
        labs(x = "", y = met_lab, colour = "Method") +
        scale_fill_brewer(palette = "Dark2") +
        theme(legend.position = "none")
      if (met_lab != metric_labels_c[[1]]) {
        g <- g + theme(
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()
        )
      }
      g
    }),
    nrow = 1
  )

#################
# COMBINE
#################
g <- wrap_plots(
  plot_spacer(),
  grid::textGrob("Ground-truth"),
  grid::textGrob("Prediction"),
  grid::textGrob("Evaluation"),
  grid::textGrob("Trajectory alignment", rot = 90), ga1, ga2, ga3,
  grid::textGrob("RNA velocity", rot = 90), gb1, gb2, gb3,
  grid::textGrob("Cell-specific\nnetwork inference", rot = 90), gc1, gc2, gc3,
  ncol = 4,
  heights = c(.25, 1, 1, 1.2),
  widths = c(.25, 1, 1, 2),
  byrow = TRUE
)

g

ggsave(exp$result("figure_2.pdf"), g, width = 12, height = 9, useDingbats = FALSE)
ggsave(exp$result("figure_2.png"), g, width = 12, height = 9)



