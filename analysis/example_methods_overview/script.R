library(tidyverse)
library(dyngen)
library(dyngen.manuscript)

exp <- start_analysis("example_methods")

# construct backbone
backbone <- backbone(
  module_info = tribble(
    ~module_id, ~basal, ~burn, ~independence,
    "A", 1, TRUE, 0,
    "B", 0, FALSE, 0,
    "C", 0, FALSE, 0,
    "D", 0, FALSE, 0,
  ),
  module_network = tribble(
    ~from, ~to, ~effect, ~strength, ~hill,
    "A", "B", 1L, 1, 2,
    "B", "C", 1L, 1, 2,
    "B", "D", 1L, 1, 2,
    "C", "D", -1L, 50, 2,
    "D", "C", -1L, 50, 2
  ),
  expression_patterns = tribble(
    ~from, ~to, ~module_progression, ~start, ~burn, ~time,
    "Burn0", "S1", "+A", TRUE, TRUE, 20,
    "S1", "S2", "+B", FALSE, FALSE, 10,
    "S2", "S3", "+C", FALSE, FALSE, 30,
    "S2", "S4", "+D", FALSE, FALSE, 30,
  )
)
backbone$module_info$color <- c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3")

model <- exp$result("model.rds") %cache% {
  set.seed(10)
  initialise_model(
    num_tfs = 8,
    num_targets = 16,
    num_hks = 0,
    num_cells = 100,
    backbone = backbone,
    verbose = TRUE,
    distance_metric = "euclidean",
    tf_network_params = tf_network_default(min_tfs_per_module = 2, sample_num_regulators = function() 1),
    simulation_params = simulation_default(
      total_time = 120,
      ssa_algorithm = GillespieSSA2::ssa_etl(tau = .1),
      experiment_params = simulation_type_wild_type(num_simulations = 16),
      census_interval = 1
    )
  ) %>%
  generate_tf_network() %>%
  generate_feature_network() %>%
  generate_kinetics() %>%
  generate_gold_standard() %>%
  generate_cells() %>%
  generate_experiment()
}

## show grn
plot_backbone_statenet(model) %>% ggsave(filename = exp$result("statenet.pdf"), width = 6, height = 6)
plot_backbone_modulenet(model) %>% ggsave(filename = exp$result("modulenet.pdf"), width = 8, height = 6)

plot_feature_network(model, show_targets = FALSE)

plot_feature_network(model)
plot_feature_network(model, show_hks = TRUE)
plot_feature_network(model)
plot_feature_network(model, show_hks = TRUE) %>% ggsave(filename = exp$result("featnet.pdf"), width = 8, height = 6)
plot_feature_network(model, show_targets = FALSE) %>% ggsave(filename = exp$result("featnet_onlytfs.pdf"), width = 8, height = 6)


## show simulations
g <- plot_gold_simulations(model) +
  scale_colour_brewer(palette = "Dark2") +
  theme_classic() + theme(legend.position = "none") + labs(x = "Comp1", y = "Comp2") +
  scale_x_continuous(breaks = 1000) +
  scale_y_continuous(breaks = 1000) +
  coord_equal()
ggsave(exp$result("gold_simulations.pdf"), g, width = 2, height = 2)

g <- plot_gold_expression(model, what = "mol_mrna", label_changing = FALSE) + facet_wrap(~edge, nrow = 1) +
  theme_classic() + theme(legend.position = "none") + labs(x = "Simulation time", y = "mRNA expression") +
  scale_x_continuous(breaks = c(0, 1))
maxy <- max(g$data$value)
maxy <- ceiling(maxy / 2.5) * 2.5
g <- g + scale_y_continuous(breaks = c(0, maxy), limits = c(0, maxy))
ggsave(exp$result("gold_mrna.pdf"), g, width = 6, height = 2)

plot_gold_expression(model, what = "mol_mrna")

g <- plot_gold_expression(model, what = "mol_mrna", label_changing = FALSE) # premrna, mrna, and protein
ggsave(exp$result("gold_expression.pdf"), g, width = 8, height = 6)

g <- plot_simulations(model) + theme_classic() + theme(axis.ticks = element_blank(), axis.line = element_blank(), axis.text = element_blank(), axis.title = element_blank()) + labs(colour = "Simulation\ntime")
ggsave(exp$result("simulations.pdf"), g, width = 8, height = 6)

g <- plot_gold_mappings(model, do_facet = FALSE) + scale_colour_brewer(palette = "Dark2")
ggsave(exp$result("simulations_mapping.pdf"), g, width = 8, height = 6)

g <- plot_simulation_expression(model, c(7), what = "mol_mrna") +
  facet_wrap(~simulation_i, ncol = 1) +
  theme_classic() + theme(legend.position = "none") +
  labs(x = "Simulation time", y = "mRNA expression")
ggsave(exp$result("simulation_expression.pdf"), g, width = 6, height = 2)




g <- plot_simulation_expression(model, 1:6, what = "mol_mrna") +
  facet_wrap(~simulation_i, ncol = 3) +
  theme_classic() + theme(legend.position = "none", strip.background = element_blank(), strip.text = element_blank()) +
  labs(x = NULL, y = NULL) +
  scale_x_continuous(breaks = 1000) +
  scale_y_continuous(breaks = 1000)
ggsave(exp$result("simulation_expression_many.pdf"), g, width = 5, height = 2)


# convert to dyno object for more plotting functionality
dataset <- as_dyno(model)

library(dynplot)
g <- plot_dimred(dataset)
ggsave(exp$result("traj_dimred.pdf"), g, width = 8, height = 6)
g <- plot_graph(dataset)
ggsave(exp$result("traj_graph.pdf"), g, width = 8, height = 6)
g <- plot_heatmap(dataset, features_oi = 40)
ggsave(exp$result("traj_heatmap.pdf"), g, width = 8, height = 6)


gold_counts <- model$simulations$counts[model$experiment$cell_info$step_ix, model$feature_info$mol_mrna]
colnames(gold_counts) <- colnames(dataset$expression)

ph <- pheatmap::pheatmap(
  t(as.matrix(gold_counts)),
  filename = exp$result("heatmap_orig.pdf"),
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
  filename = exp$result("heatmap_sampled.pdf"),
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


regs <- unique(model$feature_network$from)
tars <- model$feature_info$feature_id
g3 <- GENIE3bis::run_genie3(
  data = as.matrix(dataset$expression),
  regulators = regs,
  targets = tars,
  parallel_type = 8
)

library(tidygraph)
library(ggraph)
gr <- tbl_graph(edges = g3 %>% head(50))
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

ggsave(exp$result("genie3.pdf"), g, width = 8, height = 6)

dimred <- dyndimred::dimred_knn_fr(dataset$expression)
g <- ggplot(as.data.frame(dimred)) + geom_point(aes(comp_1, comp_2)) + coord_equal() + theme_classic() +
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  labs(x = "Comp 1", y = "Comp 2")
ggsave(exp$result("dimred.pdf"), g, width = 8, height = 6)





set.seed(1)

model <-
  initialise_model(
    num_tfs = 10,
    num_targets = 10,
    num_hks = 10,
    num_cells = 100,
    backbone = backbone,
    verbose = TRUE,
    distance_metric = "euclidean",
    tf_network_params = tf_network_default(min_tfs_per_module = 2, sample_num_regulators = function() 1),
    simulation_params = simulation_default(
      ssa_algorithm = GillespieSSA2::ssa_etl(tau = 300 / 3600),
      experiment_params = simulation_type_wild_type(num_simulations = 8),
      census_interval = .01
    )
  ) %>%
  generate_tf_network() %>%
  generate_feature_network() %>%
  generate_kinetics()




## custom plots for grn
feature_info <-
  model$feature_info %>%
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
  model$backbone$module_info %>% select(module_id, color) %>% deframe()
)

# filter feature network
feature_network <-
  model$feature_network %>%
  arrange(from == to) %>%
  left_join(feature_info %>% select(to = feature_id, to_tf = is_tf, to_hk = is_hk), by = "to")

library(tidygraph)
library(ggraph)
gr <- tbl_graph(nodes = feature_info, edges = feature_network)

# first run, open cytoscape and position nodes manually:
# write_csv(feature_network, "feature_network.csv")
# write_csv(feature_info, "feature_info.csv")

layout <- tribble(
  ~id, ~x, ~y,
  "A_TF1", 204.43, 133.4,
  "A_TF2", 237.84, 133.40,
  "B_TF1", 253.84, 184.20,
  "B_TF2", 221.13, 184.20,
  "B_TF3", 188.42, 184.20,
  "C_TF1", 183.17, 250.28,
  "C_TF2", 197.30, 273.30,
  "D_TF1", 259.09, 234.98,
  "D_TF2", 249.11, 253.86,
  "D_TF3", 239.13, 272.74,
  "Target1", 170.42, 202.20,
  "Target2", 257.13, 290.74,
  "Target3", 275.84, 198.20,
  "Target4", 267.84, 206.20,
  "Target5", 176.43, 133.40,
  "Target6", 180.43, 143.41,
  "Target7", 186.43, 151.41,
  "Target8", 203.13, 202.20,
  "Target9", 165.17, 268.28,
  "Target10", 267.11, 271.86,
  "HK1", 355.05, 212.89,
  "HK2", 388.47, 196.22,
  "HK3", 375.17, 231.19,
  "HK4", 352.91, 183.09,
  "HK5", 365.97, 165.91,
  "HK6", 401.77, 231.19,
  "HK7", 388.47, 157.47,
  "HK8", 421.88, 212.89,
  "HK9", 410.97, 165.91,
  "HK10", 424.03, 183.09
  # "A_TF1", 115, 150,
  # "A_TF2", 107, 150,
  # "B_TF1", 103, 130,
  # "B_TF2", 111, 130,
  # "B_TF3", 119, 130,
  # "C_TF1", 96, 105,
  # "C_TF2", 101, 100,
  # "D_TF1", 121, 100,
  # "D_TF2", 126, 105,
  # "D_TF3", 131, 110,
  # "Target2", 97, 142,
  # "Target3", 125, 142,
  # "Target4", 97, 128,
  # "Target8", 89, 91,
  # "Target7", 96, 86,
  # "Target9", 105, 85,
  # "Target1", 136, 96,
  # "Target10", 144, 88,
  # "Target5", 144, 100,
  # "Target6", 147, 110,
  # "HK2", 166, 136,
  # "HK1", 159, 149,
  # "HK3", 170, 149,
  # "HK4", 179, 141,
  # "HK5", 180, 132,
  # "HK6", 174, 124,
  # "HK7", 165, 122,
  # "HK8", 155, 126,
  # "HK9", 149, 134,
  # "HK10", 152, 143
) %>%
  mutate(y = -y) %>%
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
ggsave(exp$result("gen_feature_network.pdf"), g, width = 15, height = 3)


