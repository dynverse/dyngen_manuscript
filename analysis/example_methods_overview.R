library(tidyverse)
library(dyngen)

set.seed(4)

dir.create("fig/example_methods", showWarnings = FALSE, recursive = TRUE)

backbone <- backbone(
  module_info = tribble(
    ~module_id, ~ba, ~burn,
    "A", 1, TRUE,
    "B", 0, FALSE,
    "C", 0, FALSE,
    "D", 0, FALSE,
  ),
  module_network = tribble(
    ~from, ~to, ~effect, ~strength, ~cooperativity,
    "A", "B", 1, 1, 2,
    "B", "C", 1, 1, 2,
    "B", "D", 1, 1, 2,
    "C", "D", -1, 10, 2,
    "D", "C", -1, 10, 2
  ),
  expression_patterns = tribble(
    ~from, ~to, ~module_progression, ~start, ~burn, ~time,
    "Burn0", "S1", "+A", TRUE, TRUE, .5,
    "S1", "S2", "+B", FALSE, FALSE, .5,
    "S2", "S3", "+C", FALSE, FALSE, 1,
    "S2", "S4", "+D", FALSE, FALSE, 1,
  )
)
backbone$module_info$color <- c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3")
model <-
  initialise_model(
    num_tfs = 8,
    num_targets = 10,
    num_hks = 0,
    num_cells = 100,
    backbone = backbone,
    verbose = TRUE,
    download_cache_dir = "~/.cache/dyngen",
    num_cores = 8,
    distance_metric = "euclidean",
    tf_network_params = tf_network_default(min_tfs_per_module = 2, sample_num_regulators = function() 1),
    simulation_params = simulation_default(ssa_algorithm = GillespieSSA2::ssa_exact(), num_simulations = 8, census_interval = .01)
  ) %>%
  generate_tf_network() %>%
  generate_feature_network() %>%
  generate_kinetics()

## grn
plot_backbone_statenet(model) %>% ggsave(filename = "fig/example_methods/statenet.pdf", width = 6, height = 6)
plot_backbone_modulenet(model) %>% ggsave(filename = "fig/example_methods/modulenet.pdf", width = 8, height = 6)

plot_feature_network(model, show_targets = FALSE)

plot_feature_network(model)
plot_feature_network(model, show_hks = TRUE)
plot_feature_network(model)
plot_feature_network(model, show_hks = TRUE) %>% ggsave(filename = "fig/example_methods/featnet.pdf", width = 8, height = 6)
plot_feature_network(model, show_targets = FALSE) %>% ggsave(filename = "fig/example_methods/featnet_onlytfs.pdf", width = 8, height = 6)


## simulate
model2 <- model %>%
  generate_gold_standard() %>%
  generate_cells() %>%
  generate_experiment()

model2 <- dyngen:::calculate_dimred(model2, dimred_premrna = FALSE)
write_rds(model2, "fig/example_methods/model.rds", compress = "gz")

g <- plot_gold_simulations(model2) + scale_colour_brewer(palette = "Dark2") + theme_classic() + theme(legend.position = "none") + labs(x = "Comp1", y = "Comp2") +
  scale_x_continuous(breaks = 1000) +
  scale_y_continuous(breaks = 1000) +
  coord_equal()
ggsave("fig/example_methods/gold_simulations.pdf", g, width = 2, height = 2)

g <- plot_gold_expression(model2, what = "x", label_changing = FALSE) + facet_wrap(~edge, nrow = 1) +
  theme_classic() + theme(legend.position = "none") + labs(x = "Simulation time", y = "mRNA expression") +
  scale_x_continuous(breaks = c(0, 1))
maxy <- max(g$data$value)
maxy <- ceiling(maxy / 2.5) * 2.5
g <- g + scale_y_continuous(breaks = c(0, maxy), limits = c(0, maxy))
ggsave("fig/example_methods/gold_mrna.pdf", g, width = 6, height = 2)

plot_gold_expression(model2, what = "x")

g <- plot_gold_expression(model2, label_changing = FALSE) # premrna, mrna, and protein
ggsave("fig/example_methods/gold_expression.pdf", g, width = 8, height = 6)

g <- plot_simulations(model2) + theme_classic() + theme(axis.ticks = element_blank(), axis.line = element_blank(), axis.text = element_blank(), axis.title = element_blank()) + labs(colour = "Simulation\ntime")
ggsave("fig/example_methods/simulations.pdf", g, width = 8, height = 6)

g <- plot_gold_mappings(model2, do_facet = FALSE) + scale_colour_brewer(palette = "Dark2")
ggsave("fig/example_methods/simulations_mapping.pdf", g, width = 8, height = 6)

g <- plot_simulation_expression(model2, c(7), what = "x") +
  facet_wrap(~simulation_i, ncol = 1) +
  theme_classic() + theme(legend.position = "none") +
  labs(x = "Simulation time", y = "mRNA expression") +
  scale_x_continuous(breaks = c(0, 10), limits = c(0, 10))
maxy <- max(g$data$value)
maxy <- ceiling(maxy / 5) * 5
g <- g + scale_y_continuous(breaks = c(0, maxy), limits = c(0, maxy))
ggsave("fig/example_methods/simulation_expression.pdf", g, width = 6, height = 2)




g <- plot_simulation_expression(model2, 1:6, what = "x") +
  facet_wrap(~simulation_i, ncol = 3) +
  theme_classic() + theme(legend.position = "none", strip.background = element_blank(), strip.text = element_blank()) +
  labs(x = NULL, y = NULL) +
  scale_x_continuous(breaks = 1000, limits = c(0, 10)) +
  scale_y_continuous(breaks = 1000)
ggsave("fig/example_methods/simulation_expression_many.pdf", g, width = 5, height = 2)


dataset <- wrap_dataset(model2)

library(dynplot)
g <- plot_dimred(dataset)
ggsave("fig/example_methods/traj_dimred.pdf", g, width = 8, height = 6)
g <- plot_graph(dataset)
ggsave("fig/example_methods/traj_graph.pdf", g, width = 8, height = 6)
g <- plot_heatmap(dataset, features_oi = 40)
ggsave("fig/example_methods/traj_heatmap.pdf", g, width = 8, height = 6)


gold_counts <- model2$simulations$counts[model2$experiment$cell_info$step_ix, paste0("x_", colnames(dataset$expression))]
colnames(gold_counts) <- colnames(dataset$expression)

ph <- pheatmap::pheatmap(
  t(as.matrix(gold_counts)),
  filename = "fig/example_methods/heatmap_orig.pdf",
  width = 10,
  height = 8,
  border_color = NA,
  show_colnames = FALSE,
  show_rownames = FALSE,
  treeheight_col = 0,
  treeheight_row = 0
)

pheatmap::pheatmap(
  t(as.matrix(dataset$counts)),
  filename = "fig/example_methods/heatmap_sampled.pdf",
  # cluster_rows = ph$tree_row,
  # cluster_cols = ph$tree_col,
  width = 10,
  height = 8,
  border_color = NA,
  show_colnames = FALSE,
  show_rownames = FALSE,
  treeheight_col = 0,
  treeheight_row = 0
)


library(dyno)
pred <- infer_trajectory(dataset, ti_slingshot())
g <- plot_dimred(pred)
ggsave("fig/example_methods/slingshot.pdf", g, width = 8, height = 6)


regs <- unique(model$feature_network$from)
tars <- model$feature_info$feature_id
g3 <- GENIE3::run_genie3(
  data = as.matrix(dataset$expression),
  regulators = regs,
  targets = tars,
  parallel_type = 8
)

library(tidygraph)
library(ggraph)
gr <- tbl_graph(edges = g3 %>% head(200))
layout <- igraph::layout_with_fr(gr) %>%
  dynutils::scale_minmax() %>%
  magrittr::set_rownames(names(igraph::V(gr))) %>%
  magrittr::set_colnames(c("x", "y")) %>%
  as.data.frame()

cap <- circle(3, "mm")
g <- ggraph(gr, layout = "manual", x = layout$x, y = layout$y) +
  geom_edge_fan(arrow = grid::arrow(type = "closed", angle = 30, length = grid::unit(3, "mm")), start_cap = cap, end_cap = cap) +
  geom_node_circle(aes(r = .02), fill = "white") +
  theme_graph(base_family = "Helvetica") +
  coord_equal()

ggsave("fig/example_methods/genie3.pdf", g, width = 8, height = 6)

dimred <- dyndimred::dimred_umap(dataset$expression, n_neighbors = 70)
g <- ggplot(as.data.frame(dimred)) + geom_point(aes(comp_1, comp_2)) + coord_equal() + theme_classic() +
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  labs(x = "Comp 1", y = "Comp 2")
ggsave("fig/example_methods/dimred.pdf", g, width = 8, height = 6)





set.seed(1)

model2 <-
  initialise_model(
    num_tfs = 10,
    num_targets = 10,
    num_hks = 10,
    num_cells = 100,
    backbone = backbone,
    verbose = TRUE,
    download_cache_dir = "~/.cache/dyngen",
    num_cores = 8,
    distance_metric = "euclidean",
    tf_network_params = tf_network_default(min_tfs_per_module = 2, sample_num_regulators = function() 1),
    simulation_params = simulation_default(ssa_algorithm = GillespieSSA2::ssa_exact(), num_simulations = 8, census_interval = .01)
  ) %>%
  generate_tf_network() %>%
  generate_feature_network() %>%
  generate_kinetics()




## custom plots for grn
feature_info <-
  model2$feature_info %>%
  mutate(
    color_by = case_when(
      !is.na(module_id) ~ module_id,
      is_tf ~ "TF",
      !is_hk ~ "Target",
      is_hk ~ "HK"
    )
  )

color_legend <- c(
  "TF" = "black",
  "Target" = "darkgray",
  "HK" = "lightgray",
  model2$backbone$module_info %>% select(module_id, color) %>% deframe()
)

# filter feature network
feature_network <-
  model2$feature_network %>%
  arrange(from == to) %>%
  left_join(feature_info %>% select(to = feature_id, to_tf = is_tf, to_hk = is_hk), by = "to")

# # add extra edges invisible between regulators from the same module
# feature_network <-
#   bind_rows(
#     feature_network,
#     feature_info %>%
#       filter(is_tf) %>%
#       select(module_id, feature_id) %>%
#       group_by(module_id) %>%
#       do({
#         crossing(from = .$feature_id, to = .$feature_id) %>%
#           mutate(effect = -2)
#       }) %>%
#       ungroup() %>%
#       filter(from < to)
#   )

library(tidygraph)
library(ggraph)
gr <- tbl_graph(nodes = feature_info, edges = feature_network)

# layout <- igraph::layout_with_fr(gr) %>%
#   dynutils::scale_minmax() %>%
#   magrittr::set_rownames(feature_info$feature_id) %>%
#   magrittr::set_colnames(c("x", "y")) %>%
#   as.data.frame()

ggplot(layout %>% rownames_to_column("id")) +
  geom_text(aes(x, y, label = id))

layout <- tribble(
  ~id, ~x, ~y,
  "A_TF2", 107, 150,
  "A_TF1", 115, 150,
  "B_TF1", 107, 130,
  "B_TF2", 115, 130,
  "C_TF3", 91, 111,
  "C_TF2", 96, 105,
  "C_TF1", 101, 100,
  "D_TF1", 121, 100,
  "D_TF2", 126, 105,
  "D_TF3", 131, 110,
  "Target2", 97, 142,
  "Target3", 125, 142,
  "Target4", 97, 128,
  "Target8", 89, 91,
  "Target7", 96, 86,
  "Target9", 105, 85,
  "Target1", 136, 96,
  "Target10", 144, 88,
  "Target5", 144, 100,
  "Target6", 147, 110,
  "HK2", 166, 136,
  "HK1", 159, 149,
  "HK3", 170, 149,
  "HK4", 179, 141,
  "HK5", 180, 132,
  "HK6", 174, 124,
  "HK7", 165, 122,
  "HK8", 155, 126,
  "HK9", 149, 134,
  "HK10", 152, 143
) %>%
  slice(order(match(id, feature_info$feature_id))) %>%
  column_to_rownames("id") %>%
  dynutils::scale_minmax()

gr <- gr %>% activate(edges) %>% filter(is.na(effect) | effect != -2)

cap <- circle(2.5, "mm")
str <- .2

arrow_up <- grid::arrow(type = "closed", angle = 30, length = grid::unit(3, "mm"))
arrow_down <- grid::arrow(type = "closed", angle = 89, length = grid::unit(3, "mm"))

g1 <-
  ggraph(gr, layout = "manual", x = layout$x, y = layout$y) +
  geom_node_point(aes(colour = color_by, size = as.character(is_tf)), alpha = 0) +
  geom_node_point(aes(colour = color_by, size = as.character(is_tf), filter = is_tf)) +
  theme_graph(base_family = "Helvetica") +
  scale_colour_manual(values = color_legend) +
  scale_size_manual(values = c("TRUE" = 5, "FALSE" = 3)) +
  coord_equal() +
  theme(legend.position = "none")

g2 <-
  ggraph(gr, layout = "manual", x = layout$x, y = layout$y) +
  geom_node_point(aes(colour = color_by, size = as.character(is_tf)), alpha = 0) +
  geom_edge_fan(aes(filter = !is.na(effect) & effect >= 0  & to_tf), arrow = arrow_up, start_cap = cap, end_cap = cap) +
  geom_edge_fan(aes(filter = !is.na(effect) & effect < 0 & to_tf), arrow = arrow_down, start_cap = cap, end_cap = cap) +
  geom_node_point(aes(colour = color_by, size = as.character(is_tf), filter = is_tf)) +
  theme_graph(base_family = "Helvetica") +
  scale_colour_manual(values = color_legend) +
  scale_size_manual(values = c("TRUE" = 5, "FALSE" = 3)) +
  coord_equal() +
  theme(legend.position = "none")

g3 <-
  ggraph(gr, layout = "manual", x = layout$x, y = layout$y) +
  geom_node_point(aes(colour = color_by, size = as.character(is_tf)), alpha = 0) +
  geom_edge_fan(aes(filter = !is.na(effect) & effect >= 0  & !to_hk), arrow = arrow_up, start_cap = cap, end_cap = cap) +
  geom_edge_fan(aes(filter = !is.na(effect) & effect < 0 & !to_hk), arrow = arrow_down, start_cap = cap, end_cap = cap) +
  geom_node_point(aes(colour = color_by, size = as.character(is_tf), filter = !is_hk)) +
  theme_graph(base_family = "Helvetica") +
  scale_colour_manual(values = color_legend) +
  scale_size_manual(values = c("TRUE" = 5, "FALSE" = 3)) +
  coord_equal() +
  theme(legend.position = "none")

g4 <-
  ggraph(gr, layout = "manual", x = layout$x, y = layout$y) +
  geom_node_point(aes(colour = color_by, size = as.character(is_tf)), alpha = 0) +
  geom_edge_fan(aes(filter = !is.na(effect) & effect >= 0), arrow = arrow_up, start_cap = cap, end_cap = cap) +
  geom_edge_fan(aes(filter = !is.na(effect) & effect < 0), arrow = arrow_down, start_cap = cap, end_cap = cap) +
  geom_node_point(aes(colour = color_by, size = as.character(is_tf))) +
  theme_graph(base_family = "Helvetica") +
  scale_colour_manual(values = color_legend) +
  scale_size_manual(values = c("TRUE" = 5, "FALSE" = 3)) +
  coord_equal() +
  theme(legend.position = "none")

g <- patchwork::wrap_plots(
  g1, g2, g3, g4, nrow = 1
)
ggsave("fig/example_methods/gen_feature_network.pdf", g, width = 15, height = 3)


