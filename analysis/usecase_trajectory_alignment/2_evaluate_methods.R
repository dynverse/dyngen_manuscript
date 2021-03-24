library(dyngen.manuscript)
library(tidyverse)

exp <- start_analysis("usecase_trajectory_alignment")

# cross datasets and methods
design_smoothing <-
  crossing(
    read_rds(exp$result("design_datasets.rds")),
    method = forcats::fct_inorder(c("DTW", "cellAlign"))
  )

#' design_smoothing %>% mutate(rn = row_number()) %>% dynutils::extract_row_to_list(1) %>% list2env(.GlobalEnv)
results <- exp$result("results.rds") %cache% pmap_dfr(
  design_smoothing %>% mutate(rn = row_number()),
  function(rn, group, id, method, ...) {
    file <- exp$dataset_file(id)
    if (file.exists(file)) {
      dataset <- read_rds(file)

      method_fun <- ta_methods[[method]]

      dataset1 <- dataset %>% filter_cells(model == "left")
      dataset2 <- dataset %>% filter_cells(model == "right")

      out <- method_fun(dataset1, dataset2)
      pt1_aligned <- out$pt1_aligned
      pt2_aligned <- out$pt2_aligned

      distance <- mean(abs(pt1_aligned - pt2_aligned))

      aupt <- pracma::trapz(
        pt1_aligned + pt2_aligned,
        abs(pt1_aligned - pt2_aligned)
      )

      score <- 1 - aupt

      design_smoothing[rn, ] %>% mutate(
        distance,
        aupt,
        score
      )
    }
  }
)

ggstatsplot::ggwithinstats(
  results,
  x = method,
  y = score,
  type = "np",
  pairwise.display = "all"
)

z <- ggstatsplot::ggwithinstats(
  results,
  x = method,
  y = distance,
  type = "np",
  pairwise.display = "all"
)
z$labels$subtitle

# PART 3: Scores ----------------------------------------------------------
results <- read_rds(exp$result("results.rds")) %>%
  mutate(
    method = factor(as.character(method), c("DTW", "cellAlign"))
  )

quants <- results %>%
  group_by(method) %>%
  summarise_at(vars(distance), list(
    min = ~quantile(., 0),
    lower = ~quantile(., .25),
    mean = ~median(.),
    upper = ~quantile(., .75),
    max = ~quantile(., 1)
  ))

g <- ggplot() +
  geom_violin(aes(method, distance), results) +
  geom_path(aes(method, distance, group = id), results, linetype = "dashed", colour = "gray") +
  geom_boxplot(
    aes(method, ymin = min, lower = lower, middle = mean, upper = upper, max = max),
    quants,
    stat = "identity", width = 0.35, size = 0.45, fill = NA
  ) +
  geom_point(aes(method, distance, colour = method), results) +
  theme_classic() +
  theme_common() +
  labs(
    x = NULL, y = "Distance (lower is better)",
    colour = "Method", subtitle = z$labels$subtitle) +
  scale_colour_brewer(palette = "Set2")

ggsave(exp$result("score_summary.pdf"), g, width = 6, height = 6)

write_rds(g, exp$result("score_summary.rds"), compress = "gz")
