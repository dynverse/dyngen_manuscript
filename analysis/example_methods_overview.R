library(tidyverse)
library(dyngen)

set.seed(4)

dir.create("fig/example_methods", showWarnings = FALSE, recursive = TRUE)

backbone <- backbone(
  module_info = tribble(
    ~module_id, ~a0, ~burn,
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


