library(tidyverse)
library(dyngen)

set.seed(3)

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
    "Burn0", "S1", "+A", TRUE, TRUE, 2,
    "S1", "S2", "+B", FALSE, FALSE, 2,
    "S2", "S3", "+C", FALSE, FALSE, 2,
    "S2", "S4", "+D", FALSE, FALSE, 2,
  )
)
model <-
  initialise_model(
    num_tfs = nrow(backbone$module_info),
    num_targets = 10,
    num_hks = 10,
    num_cells = 0,
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
plot_backbone_statenet(model) #%>% ggsave(filename = "dyngen/statenet.pdf", width = 6, height = 6)
plot_backbone_modulenet(model) #%>% ggsave(filename = "dyngen/modulenet.pdf", width = 8, height = 6)

plot_feature_network(model, show_targets = FALSE)

plot_feature_network(model)
plot_feature_network(model, show_hks = TRUE)
plot_feature_network(model)
plot_feature_network(model, show_hks = TRUE) #%>% ggsave(filename = "dyngen/featnet.pdf", width = 8, height = 6)
plot_feature_network(model, show_targets = FALSE) #%>% ggsave(filename = "dyngen/featnet_onlytfs.pdf", width = 8, height = 6)

## simulate
model2 <- model %>%
  generate_gold_standard() %>%
  generate_cells() %>%
  generate_experiment()

model2 <- dyngen:::calculate_dimred(model2, dimred_premrna = FALSE)
# write_rds(model2, "dyngen/model.rds", compress = "gzip")

g <- plot_gold_simulations(model2) + scale_colour_brewer(palette = "Dark2")
g # ggsave("dyngen/gold_simulations.pdf", g, width = 8, height = 6)

g <- plot_gold_expression(model2, what = "x") # mrna
g # ggsave("dyngen/gold_mrna.pdf", g, width = 8, height = 6)
g <- plot_gold_expression(model2, label_changing = FALSE) # premrna, mrna, and protein
g # ggsave("dyngen/gold_expression.pdf", g, width = 8, height = 6)

g <- plot_simulations(model2) + theme_classic() + theme(axis.ticks = element_blank(), axis.line = element_blank(), axis.text = element_blank(), axis.title = element_blank()) + labs(colour = "Simulation\ntime")
g # ggsave("dyngen/simulations.pdf", g, width = 8, height = 6)

g <- plot_gold_mappings(model2, do_facet = FALSE) + scale_colour_brewer(palette = "Dark2")
g # ggsave("dyngen/simulations_mapping.pdf", g, width = 8, height = 6)

g <- plot_simulation_expression(model2, 1:4, what = "x")
g # ggsave("dyngen/simulation_expression.pdf", g, width = 8, height = 6)


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

