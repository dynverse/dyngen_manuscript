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
feature_network <-
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

gr <- tbl_graph(nodes = feature_info, edges = feature_network)
layout <- igraph::layout_with_fr(gr) %>%
  dynutils::scale_minmax() %>%
  magrittr::set_rownames(feature_info$feature_id) %>%
  magrittr::set_colnames(c("x", "y")) %>%
  as.data.frame()
gr <- gr %>% activate(edges) %>% filter(is.na(effect) | effect != -2)
feature_network <- feature_network %>% filter(effect != -2)

cap <- circle(1.3, "mm")
str <- .1
arrow_up <- grid::arrow(type = "closed", angle = 30, length = grid::unit(1.3, "mm"))
arrow_down <- grid::arrow(type = "closed", angle = 89, length = grid::unit(1.3, "mm"))

g1 <- ggraph(gr, layout = "manual", x = layout$x, y = layout$y) +
  geom_edge_loop(aes(strength = str, filter = !is.na(effect) & effect >= 0 & from == to, colour = sign(effect)), arrow = arrow_up, start_cap = cap, end_cap = cap) +
  geom_edge_loop(aes(strength = str, filter = !is.na(effect) & effect < 0 & from == to, colour = sign(effect)), arrow = arrow_down, start_cap = cap, end_cap = cap) +
  geom_edge_loop(aes(strength = str, filter = is.na(effect), colour = sign(effect))) +
  geom_edge_fan(aes(filter = !is.na(effect) & effect >= 0 & from != to, colour = sign(effect)), arrow = arrow_up, start_cap = cap, end_cap = cap) +
  geom_edge_fan(aes(filter = !is.na(effect) & effect < 0, colour = sign(effect)), arrow = arrow_down, start_cap = cap, end_cap = cap) +
  geom_edge_fan(aes(filter = is.na(effect))) +
  geom_node_point(size = 1, colour = "gray") +
  theme_graph(base_family = "Helvetica") +
  coord_equal() +
  labs(size = "is TF", color = "Module group", edge_color = "Regulatory\nactivity") +
  scale_edge_colour_gradientn(colours = c(rev(RColorBrewer::brewer.pal(5, "Reds")), RColorBrewer::brewer.pal(5, "Greens")))

cell_wps <- do.call(dynwrap::select_waypoint_cells, c(dataset[c("milestone_ids", "milestone_network", "milestone_percentages", "progressions", "divergence_regions")], list(num_cells_selected = 1)))
cells <- sample(cell_wps, 6) %>% {.[order(match(., dataset$cell_ids))]}

g2plots <- map(seq_along(cells), function(i) {
  cell_id <- cells[[i]]
  feature_info <- dataset$feature_info
  feature_network_sc <-
    dataset$regulatory_network_sc %>%
    filter(cell_id == !!cell_id) %>%
    arrange(as.character(regulator) == as.character(target)) %>%
    left_join(feature_network %>% select(regulator, target, effect), by = c("regulator", "target")) %>%
    mutate(effect = strength)
  fn <- bind_rows(
    feature_network_sc,
    anti_join(feature_network, feature_network_sc, by = c("regulator", "target")) %>% mutate(effect = NA)
  ) %>%
    select(regulator, target, everything())
  grsc <- tbl_graph(nodes = feature_info, edges = fn)


  # cap <- circle(1, "mm")
  # str <- .05
  # arrow_up <- grid::arrow(type = "closed", angle = 30, length = grid::unit(1, "mm"))
  # arrow_down <- grid::arrow(type = "closed", angle = 89, length = grid::unit(1, "mm"))

  ggraph(grsc, layout = "manual", x = layout$x, y = layout$y) +
    geom_edge_loop(aes(strength = str, filter = is.na(effect)), colour = "#EEEEEE") +
    geom_edge_fan(aes(filter = is.na(effect)), colour = "#EEEEEE") +
    geom_node_point(size = 1, colour = "gray") +
    geom_edge_loop(aes(strength = str, filter = !is.na(effect) & effect >= 0 & from == to, colour = effect), arrow = arrow_up, start_cap = cap, end_cap = cap) +
    geom_edge_loop(aes(strength = str, filter = !is.na(effect) & effect < 0 & from == to, colour = effect), arrow = arrow_down, start_cap = cap, end_cap = cap) +
    geom_edge_fan(aes(filter = !is.na(effect) & effect >= 0 & from != to, colour = effect), arrow = arrow_up, start_cap = cap, end_cap = cap) +
    geom_edge_fan(aes(filter = !is.na(effect) & effect < 0, colour = effect), arrow = arrow_down, start_cap = cap, end_cap = cap) +
    theme_graph(base_family = "Helvetica") +
    coord_equal() +
    labs(size = "is TF", color = "Module group", title = paste0("GRN of cell ", i)) +
    scale_edge_colour_gradientn(colours = c(rev(RColorBrewer::brewer.pal(5, "Reds")), RColorBrewer::brewer.pal(5, "Greens"))) +
    theme(legend.position = "none")
})
#
# g <- patchwork::wrap_plots(g1 + labs(title = "Global GRN"), patchwork::wrap_plots(g2plots, nrow = 2), nrow = 1, widths = c(2, 3.5)) &
#   theme(
#     plot.margin = margin(0, 0, 0, 0, "cm"),
#     plot.title = element_text(size = 10),
#     legend.text = element_text(size = 6),
#     legend.title = element_text(size = 8)
#   )
# ggsave("fig/network_inference/casewise_grn.pdf", width = 11, height = 4)



g <- patchwork::wrap_plots(
  c(
    list(g1 + labs(title = "Global GRN") + theme(legend.position = "none")),
    g2plots[1:2],
    list(ggpubr::as_ggplot(ggpubr::get_legend(g1 + theme(legend.text = element_text(size = 6), legend.title = element_text(size = 8))))),
    g2plots[3:5],
    list(patchwork::plot_spacer())
  ),
  nrow = 2,
  widths = c(3, 3, 3, .75)
) &
  theme(
    plot.margin = margin(0, 0, 0, 0, "cm"),
    plot.title = element_text(size = 10)
  )
ggsave("fig/network_inference/casewise_grn.pdf", width = 8, height = 6)

