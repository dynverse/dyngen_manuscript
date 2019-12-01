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
    simulation_params = simulation_default(
      kinetics_noise_function = kinetics_noise_none(),
      ssa_algorithm = GillespieSSA2::ssa_exact(),
      experiment_params = simulation_type_wild_type(num_simulations = 8),
      census_interval = .01,
      store_grn = TRUE
    )
  ) %>%
  generate_tf_network() %>%
  generate_feature_network() %>%
  generate_kinetics()

## grn
plot_backbone_statenet(model) %>% ggsave(filename = "fig/small_bifurc/statenet.pdf", width = 6, height = 6)
plot_backbone_modulenet(model) %>% ggsave(filename = "fig/small_bifurc/modulenet.pdf", width = 8, height = 6)

plot_feature_network(model, show_targets = FALSE)

plot_feature_network(model)
plot_feature_network(model, show_hks = TRUE)
plot_feature_network(model)
plot_feature_network(model, show_hks = TRUE) %>% ggsave(filename = "fig/small_bifurc/featnet.pdf", width = 8, height = 6)
plot_feature_network(model, show_targets = FALSE) %>% ggsave(filename = "fig/small_bifurc/featnet_onlytfs.pdf", width = 8, height = 6)

## simulate
model2 <- model %>%
  generate_gold_standard()

g <- plot_gold_expression(model2, what = "x") # mrna
ggsave("fig/small_bifurc/gold_mrna.pdf", g, width = 8, height = 6)
g <- plot_gold_expression(model2, label_changing = FALSE) # premrna, mrna, and protein
ggsave("fig/small_bifurc/gold_expression.pdf", g, width = 8, height = 6)

model2 <- model2 %>%
  generate_cells() %>%
  generate_experiment()

# recalculate dimred
co1 <- model2$gold_standard$counts %>% as.matrix
co2 <- model2$simulations$counts[,colnames(co1)] %>% as.matrix
space <- prcomp(rbind(co1, co2), rank. = 2)$x
colnames(space) <- c("comp_1", "comp_2")
model2$gold_standard$dimred <- space[seq_len(nrow(co1)),]
model2$simulations$dimred <- space[-seq_len(nrow(co1)),]


g <- plot_gold_simulations(model2) + scale_colour_brewer(palette = "Dark2")
ggsave("fig/small_bifurc/gold_simulations.pdf", g, width = 8, height = 6)


write_rds(model2, "fig/small_bifurc/model.rds", compress = "gz")


g <- plot_simulations(model2) + theme_classic() + theme(axis.ticks = element_blank(), axis.line = element_blank(), axis.text = element_blank(), axis.title = element_blank()) + labs(colour = "Simulation\ntime")
ggsave("fig/small_bifurc/simulations.pdf", g, width = 8, height = 6)

g <- plot_gold_mappings(model2, do_facet = FALSE) + scale_colour_brewer(palette = "Dark2")
ggsave("fig/small_bifurc/simulations_mapping.pdf", g, width = 8, height = 6)

g <- plot_simulation_expression(model2, 1:4, what = "x")
ggsave("fig/small_bifurc/simulation_expression.pdf", g, width = 8, height = 6)


dataset <- wrap_dataset(model2, store_grn = TRUE)

library(dynplot)
g <- plot_dimred(dataset)
ggsave("fig/small_bifurc/traj_dimred.pdf", g, width = 8, height = 6)
g <- plot_graph(dataset)
ggsave("fig/small_bifurc/traj_graph.pdf", g, width = 8, height = 6)
g <- plot_heatmap(dataset, features_oi = 40)
ggsave("fig/small_bifurc/traj_heatmap.pdf", g, width = 8, height = 6)

pheatmap::pheatmap(
  t(as.matrix(dataset$expression[1:150,])),
  filename = "fig/small_bifurc/heatmap.pdf",
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
ggsave("fig/small_bifurc/slingshot.pdf", g, width = 8, height = 6)

expr <- dataset$expression[, Matrix::colMeans(dataset$expression > 0) > 0, drop = FALSE]
regs <- unique(model$feature_network$from) %>% intersect(colnames(expr))
tars <- model$feature_info$feature_id %>% intersect(colnames(expr))
g3 <- GENIE4::run_genie3(
  data = as.matrix(dataset$expression),
  regulators = regs,
  targets = tars,
  parallel_type = 8,
  rf_package = "ranger"
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

ggsave("fig/small_bifurc/genie3.pdf", g, width = 8, height = 6)

dimred <- dyndimred::dimred_umap(dataset$expression, n_neighbors = 70)
g <- ggplot(as.data.frame(dimred)) + geom_point(aes(comp_1, comp_2)) + coord_equal() + theme_classic() +
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  labs(x = "Comp 1", y = "Comp 2")
ggsave("fig/small_bifurc/dimred.pdf", g, width = 8, height = 6)

# single cell grn
feature_info <- dataset$feature_info %>%
  mutate(gene_group = factor(ifelse(!is.na(module_id), gsub("[0-9]*", "", module_id), gsub("[0-9]*", "", feature_id)))) %>%
  select(feature_id, gene_group)
interaction_info <-
  dataset$feature_network %>%
  transmute(interaction = paste0(from, "->", to), from, to) %>%
  left_join(feature_info %>% rename(to = feature_id, `Target module` = gene_group), by = "to") %>%
  left_join(feature_info %>% rename(from = feature_id, `Regulator module` = gene_group), by = "from") %>%
  select(-from, -to) %>%
  as.data.frame() %>%
  column_to_rownames("interaction")
cell_info <-
  data.frame(`Cell group` = dynwrap::group_onto_nearest_milestones(dataset), check.names = FALSE)
gene_colours <- model$backbone$module_info %>% transmute(feature_group = gsub("[0-9]*", "", module_id), color) %>% group_by(feature_group) %>% slice(n()) %>% deframe() %>% c(HK = "lightgray", Target = "darkgray")
ann_colours <- list(
  `Cell group` = dynplot:::add_milestone_coloring(tibble(id = levels(unique(cell_info$`Cell group`)))) %>% deframe(),
  `Regulator module` = gene_colours,
  `Target module` = gene_colours
)
# color <- c("white", grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "Blues")[-1:-3])(99))
color <- c(
  rev(grDevices::colorRampPalette(RColorBrewer::brewer.pal(7, "Reds"))(49)),
  "#BBBBBB",
  grDevices::colorRampPalette(RColorBrewer::brewer.pal(7, "Blues"))(49)
)
breaks <- seq(-1, 1, length.out = length(color) + 1)

cell_df <- tibble(
  id = dataset$cell_ids,
  pseudotime = dynwrap::calculate_pseudotime(dataset),
  group = dynwrap::group_onto_nearest_milestones(dataset),
  group_ord = match(group, c("sA", "sB", "sC", "sEndC", "sD", "sEndD"))
) %>%
  arrange(group_ord, pseudotime)
pheatmap::pheatmap(
  t(as.matrix(dataset$feature_network_sc[cell_df$id,Matrix::colMeans(dataset$feature_network_sc != 0) > 0])),
  cluster_cols = FALSE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  annotation_col = cell_info,
  annotation_row = interaction_info,
  annotation_colors = ann_colours,
  color = color,
  breaks = breaks,
  border = NA,
  filename = "fig/small_bifurc/scgrn.pdf",
  width = 8,
  height = 6,
  angle_col = 315,
  treeheight_col = 0,
  treeheight_row = 0,
  gaps_col = which(diff(cell_df$group_ord) != 0),
  clustering_distance_rows = "correlation"
)

pheatmap::pheatmap(
  t(as.matrix(dataset$feature_network_sc[,Matrix::colMeans(dataset$feature_network_sc != 0) > 0])),
  show_rownames = FALSE,
  show_colnames = FALSE,
  color = color,
  breaks = breaks,
  border = NA,
  filename = "fig/small_bifurc/scgrn_noann.pdf",
  width = 8,
  height = 6,
  angle_col = 315,
  # treeheight_col = 0,
  # treeheight_row = 0,
  gaps_col = which(diff(cell_df$group_ord) != 0),
  clustering_distance_rows = "correlation"
)

