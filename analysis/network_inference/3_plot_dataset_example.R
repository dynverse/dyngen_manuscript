library(tidyverse)
library(dyngen)
library(dyngenanalysis)
library(dynutils)
library(tidygraph)
library(ggraph)

dataset_file <- paste0("derived_files/network_inference/datasets/small_disconnected/dataset.rds")
dataset <- read_rds(dataset_file)

set.seed(1)

feature_info <- dataset$feature_info
feature_network <-
  dataset$regulatory_network %>%
  arrange(as.character(regulator) == as.character(target))

# add extra edges invisible between regulators from the same module
feature_network_tmp <-
  bind_rows(
    feature_network,
    feature_info %>%
      filter(is_tf) %>%
      select(module_id, feature_id) %>%
      group_by(module_id) %>%
      do({
        crossing(regulator = .$feature_id, target = .$feature_id) %>%
          mutate(effect = -2)
      }) %>%
      ungroup() %>%
      filter(regulator < target)
  )

gr <- tbl_graph(nodes = feature_info, edges = feature_network_tmp)
layout <- igraph::layout_with_fr(gr) %>%
  dynutils::scale_minmax() %>%
  magrittr::set_rownames(feature_info$feature_id) %>%
  magrittr::set_colnames(c("x", "y")) %>%
  as.data.frame()


cell_wps <- do.call(dynwrap::select_waypoint_cells, c(dataset[c("milestone_ids", "milestone_network", "milestone_percentages", "progressions", "divergence_regions")], list(num_cells_selected = 1)))
cells <- sample(cell_wps, 5) %>% {.[order(match(., dataset$cell_ids))]}

node_df <- bind_cols(feature_info %>% select(-w:-y), layout)

cap <- .015
capend <- .02
edge_df <-
  bind_rows(
    feature_network %>% mutate(group = "Global GRN"),
    dataset$regulatory_network_sc %>%
      filter(cell_id %in% cells) %>%
      arrange(as.character(regulator) == as.character(target)) %>%
      left_join(feature_network %>% select(regulator, target, effect), by = c("regulator", "target")) %>%
      mutate(effect = strength, group = paste0("GRN of cell ", match(cell_id, cells)))
  ) %>%
  mutate(
    regulator = as.character(regulator),
    target = as.character(target),
    x_ = layout[regulator, "x"],
    y_ = layout[regulator, "y"],
    xend_ = layout[target, "x"],
    yend_ = layout[target, "y"],
    len = sqrt( (x_ - xend_)^2 + (y_ - yend_)^2 ),
    al = (len - cap) / len,
    alend = (len - capend) / len,
    x = al * x_ + (1 - al) * xend_,
    y = al * y_ + (1 - al) * yend_,
    xend = alend * xend_ + (1 - alend) * x_,
    yend = alend * yend_ + (1 - alend) * y_
  )

arrow_up <- grid::arrow(type = "closed", angle = 30, length = grid::unit(1.3, "mm"))
arrow_down <- grid::arrow(type = "closed", angle = 89, length = grid::unit(1.3, "mm"))
g <- ggplot() +
  # geom_segment(aes(x = x, y = y, xend = xend, yend = yend), edge_df %>% filter(is.na(cell_id)) %>% select(-group), colour = "lightgray") +
  geom_point(aes(x, y), colour = "gray", node_df) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend, colour = effect), edge_df %>% filter(effect >= 0), arrow = arrow_up) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend, colour = effect), edge_df %>% filter(effect < 0), arrow = arrow_down) +
  theme_graph(base_family = "Helvetica") +
  coord_equal() +
  labs(colour = "Regulatory\nactivity") +
  scale_colour_gradientn(colours = c(rev(RColorBrewer::brewer.pal(5, "Reds")), RColorBrewer::brewer.pal(5, "Greens"))) +
  facet_wrap(~group) +
  theme(
    plot.margin = margin(0, 0, 0, 0, "cm")
  )

ggsave("fig/network_inference/casewise_grn.pdf", g, width = 8, height = 5)
write_rds(g, "fig/network_inference/casewise_grn.rds")

