library(tidyverse)
library(dyngen)

set.seed(10)

backbone <- bblego(
  bblego_start("A", type = "simple", num_modules = 1),
  bblego_linear("A", "B", type = "simple", num_modules = 3),
  bblego_branching("B", c("C", "D")),
  bblego_end("C", type = "simple", num_modules = 2),
  bblego_end("D", type = "simple", num_modules = 2)
)
model <-
  initialise_model(
    num_tfs = 12,
    num_targets = 30,
    num_hks = 15,
    num_cells = 150,
    backbone = backbone,
    verbose = TRUE,
    download_cache_dir = "~/.cache/dyngen",
    num_cores = 8,
    distance_metric = "euclidean",
    simulation_params = simulation_default(ssa_algorithm = GillespieSSA2::ssa_exact(), num_simulations = 8, census_interval = .01)
  ) %>%
  generate_tf_network() %>%
  generate_feature_network() %>%
  generate_kinetics()

## grn
plot_backbone_statenet(model) %>% ggsave(filename = "dyngen/statenet.pdf", width = 6, height = 6)
plot_backbone_modulenet(model) %>% ggsave(filename = "dyngen/modulenet.pdf", width = 8, height = 6)

plot_feature_network(model, show_targets = FALSE)

plot_feature_network(model)
plot_feature_network(model, show_hks = TRUE)
plot_feature_network(model)
plot_feature_network(model, show_hks = TRUE) %>% ggsave(filename = "dyngen/featnet.pdf", width = 8, height = 6)
plot_feature_network(model, show_targets = FALSE) %>% ggsave(filename = "dyngen/featnet_onlytfs.pdf", width = 8, height = 6)

## simulate
model2 <- model %>%
  generate_gold_standard() %>%
  generate_cells() %>%
  generate_experiment()

# recalculate dimred
space <- prcomp(model2$simulations$counts %>% as.matrix, rank. = 2)$x
colnames(space) <- c("comp_1", "comp_2")
model2$simulations$dimred <- space

write_rds(model2, "dyngen/model.rds", compress = "gzip")

g <- plot_gold_simulations(model2) + scale_colour_brewer(palette = "Dark2")
ggsave("dyngen/gold_simulations.pdf", g, width = 8, height = 6)

g <- plot_gold_expression(model2, what = "x") # mrna
ggsave("dyngen/gold_mrna.pdf", g, width = 8, height = 6)
g <- plot_gold_expression(model2, label_changing = FALSE) # premrna, mrna, and protein
ggsave("dyngen/gold_expression.pdf", g, width = 8, height = 6)

g <- plot_simulations(model2) + theme_classic() + theme(axis.ticks = element_blank(), axis.line = element_blank(), axis.text = element_blank(), axis.title = element_blank()) + labs(colour = "Simulation\ntime")
ggsave("dyngen/simulations.pdf", g, width = 8, height = 6)

g <- plot_gold_mappings(model2, do_facet = FALSE) + scale_colour_brewer(palette = "Dark2")
ggsave("dyngen/simulations_mapping.pdf", g, width = 8, height = 6)

g <- plot_simulation_expression(model2, 1:4, what = "x")
ggsave("dyngen/simulation_expression.pdf", g, width = 8, height = 6)


dataset <- wrap_dataset(model2)

library(dynplot)
g <- plot_dimred(dataset)
ggsave("dyngen/traj_dimred.pdf", g, width = 8, height = 6)
g <- plot_graph(dataset)
ggsave("dyngen/traj_graph.pdf", g, width = 8, height = 6)
g <- plot_heatmap(dataset, features_oi = 40)
ggsave("dyngen/traj_heatmap.pdf", g, width = 8, height = 6)

pheatmap::pheatmap(
  t(as.matrix(dataset$expression[1:150,])),
  filename = "dyngen/heatmap.pdf",
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
ggsave("dyngen/slingshot.pdf", g, width = 8, height = 6)


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

ggsave("dyngen/genie3.pdf", g, width = 8, height = 6)

dimred <- dyndimred::dimred_umap(dataset$expression, n_neighbors = 70)
g <- ggplot(as.data.frame(dimred)) + geom_point(aes(comp_1, comp_2)) + coord_equal() + theme_classic() +
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  labs(x = "Comp 1", y = "Comp 2")
ggsave("dyngen/dimred.pdf", g, width = 8, height = 6)

