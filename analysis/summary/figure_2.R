library(tidyverse)
# library(dyngen)
library(dyngen.manuscript)
library(patchwork)

set.seed(1)

exp <- start_analysis("summary")

theme_title <- function() {
  theme(
    plot.title = element_text(hjust = .5),
    axis.title = element_blank(),
    axis.text = element_blank(),
    line = element_blank(),
    rect = element_blank()
  )
}
plot_title <- function(value, angle = 0, size = 10) {
  ggplot() + geom_text(aes(1, 1, label = value), angle = angle, size = size) + theme_title()
  # ggplot() + geom_point(aes(1,1), alpha = 0) + labs(title = title, subtitle = subtitle, tag = tag) + theme_title()
}


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
    legend.position = "right",
    text = element_text(family = "Helvetica"),
    axis.title = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank()
  ) +
  labs(title = NULL, subtitle = NULL, tag = NULL)

scores_a <-
  read_rds(exp_a$result("results.rds")) %>%
  gather(metric, score, abwap)

plota_data <-
  scores_a %>%
  group_by(method) %>%
  summarise_at(vars(score), list(
    min = ~quantile(., .05, na.rm = TRUE),
    lower = ~quantile(., .25, na.rm = TRUE),
    mean = ~median(., na.rm = TRUE),
    upper = ~quantile(., .75, na.rm = TRUE),
    max = ~quantile(., .95, na.rm = TRUE)
  )) %>%
  ungroup()

ga3 <-
  ggplot(plota_data) +
  geom_boxplot(
    aes(method, ymin = min, lower = lower, middle = mean, upper = upper, max = max, fill = method),
    stat = "identity", width = 0.5, size = 0.45
  ) +
  coord_flip() +
  theme_classic() +
  theme_common() +
  labs(x = NULL, y = "ABWAP score") +
  scale_fill_brewer(palette = "Pastel2") +
  theme(legend.position = "none") +
  expand_limits(y = c(0, 1))


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

scores_b <- read_rds(exp_b$result("scores.rds"))$summ %>%
  inner_join(
    tribble(
      ~method_id, ~params_id, ~method_label,
      "velocyto", "constant_velocity", "velocyto",
      "scvelo", "stochastic", "scvelo\nstochastic",
      "scvelo", "dynamical", "scvelo\ndynamical"
    ) %>%
    mutate(method_label = forcats::fct_inorder(method_label)),
    by = c("method_id", "params_id")
  )


metric_labels_b <- c(cor = "Velocity correlation", mean_cosine = "Velocity arrow cosine")

plotb_data <-
  scores_b %>%
  gather(metric, score, cor, mean_cosine) %>%
  mutate(
    metric_label = factor(metric_labels_b[metric], levels = metric_labels_b)
  ) %>%
  group_by(method_label, metric_label) %>%
  summarise_at(vars(score), list(
    min = ~quantile(., .05, na.rm = TRUE),
    lower = ~quantile(., .25, na.rm = TRUE),
    mean = ~median(., na.rm = TRUE),
    upper = ~quantile(., .75, na.rm = TRUE),
    max = ~quantile(., .95, na.rm = TRUE)
  ))

limits_tib_b <- tibble(
  metric_label = forcats::fct_inorder(metric_labels_b),
  x = 1,
  y0 = c(0, 0),
  y1 = c(.4, 1)
) %>% gather(l, y, y0, y1)

gb3 <-
  patchwork::wrap_plots(
    list = pmap(
      crossing(metric = forcats::fct_inorder(metric_labels_b)),
      function(metric) {
      g <- ggplot(plotb_data %>% filter(metric_label == metric)) +
        geom_boxplot(
          aes(forcats::fct_rev(method_label), ymin = min, lower = lower, middle = mean, upper = upper, max = max, fill = method_label),
          stat = "identity", width = 0.5, size = 0.45
        ) +
        geom_blank(aes(x, y), limits_tib_b %>% filter(metric_label == metric)) +
        theme_classic() +
        theme_common() +
        coord_flip() +
        labs(x = NULL, y = metric, colour = "Method") +
        scale_fill_brewer(palette = "Pastel2") +
        theme(legend.position = "none")
      if (metric != metric_labels_b[[1]]) {
        g <- g +
          theme(
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank()
          ) +
          labs(x = NULL)
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
  # labs(colour = "Regulatory strength\nin cell 1") +
  theme(
    legend.position = "none",
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
  labs(colour = "Regulatory\nstrength in cell 1") +
  theme(
    legend.position = "right",
    text = element_text(family = "Helvetica"),
    axis.title = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank()
  ) +
  labs(title = NULL, subtitle = NULL, tag = NULL)

scores_c <- read_rds(exp_c$result("scores.rds"))

metric_labels_c <- c(auroc = "mean AUROC", aupr = "mean AUPR")
method_labels_c <- c(lionesspearson = "LIONESS\n+ Pearson", pyscenic = "pySCENIC", ssn = "SSN*")
plotc_data <-
  scores_c$summ %>%
  mutate(cni_method_name = method_labels_c[cni_method_id]) %>%
  filter(method == "casewise_casewise") %>%
  gather(metric, score, auroc, aupr) %>%
  mutate(
    metric_label = factor(metric_labels_c[metric], levels = metric_labels_c)
  ) %>%
  group_by(cni_method_name, metric, metric_label) %>%
  summarise_at(vars(score), list(
    min = ~quantile(., .05, na.rm = TRUE),
    lower = ~quantile(., .25, na.rm = TRUE),
    mean = ~median(., na.rm = TRUE),
    upper = ~quantile(., .75, na.rm = TRUE),
    max = ~quantile(., .95, na.rm = TRUE)
  ))


limits_tib_c <- tibble(
  metric_label = forcats::fct_inorder(metric_labels_c),
  x = 1,
  y0 = c(.45, 0),
  y1 = c(.75, .12)
) %>% gather(l, y, y0, y1)
gc3 <-
  patchwork::wrap_plots(
    list = map(metric_labels_c, function(met_lab) {
      g <- ggplot(plotc_data %>% filter(metric_label == met_lab)) +
        geom_boxplot(
          aes(forcats::fct_rev(cni_method_name), ymin = min, lower = lower, middle = mean, upper = upper, max = max, fill = cni_method_name),
          stat = "identity", width = 0.5, size = 0.45
        ) +
        geom_blank(aes(x, y), limits_tib_c %>% filter(metric_label == met_lab)) +
        theme_classic() +
        theme_common() +
        coord_flip() +
        labs(x = "", y = met_lab, colour = "Method") +
        scale_fill_brewer(palette = "Pastel2") +
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
gc3[[1]] <- gc3[[1]] + scale_y_continuous(breaks = c(.45, .55, .65, .75))
gc3[[2]] <- gc3[[2]] + scale_y_continuous(breaks = c(0, .04, .08, .12))

#################
# COMBINE
#################
g <- wrap_plots(
  plot_title("", angle = 90, size = 5) + labs(title = "A") + theme(plot.title = element_text(hjust = 0)),
  ga1 + labs(title = "Ground-truth") + theme(plot.title = element_text(hjust = .5)),
  ga2 + labs(title = "Prediction") + theme(plot.title = element_text(hjust = .5)),
  ga3 + labs(title = "Evaluation") + theme(plot.title = element_text(hjust = .5)),
  plot_title("", angle = 90, size = 5) + labs(title = "B") + theme(plot.title = element_text(hjust = 0)),
  gb1, gb2, gb3,
  plot_title("", angle = 90, size = 5) + labs(title = "C") + theme(plot.title = element_text(hjust = 0)),
  gc1, gc2, gc3,
  ncol = 4,
  heights = c(1, 1, 1.2),
  widths = c(.1, 1, 1, 1.5),
  byrow = TRUE
)

# g


ggsave(exp$result("figure_2.pdf"), g, width = 12, height = 9, useDingbats = FALSE)
ggsave(exp$result("figure_2.png"), g, width = 12, height = 9)



