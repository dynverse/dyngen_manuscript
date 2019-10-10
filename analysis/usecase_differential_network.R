library(tidyverse)
library(dyngen)

set.seed(2)

out_dir <- "fig/usecase/differential_network/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Generate base model -----------------------------------------------------
backbone <- bblego(
  bblego_start("A", type = "simple", num_modules = 4),
  bblego_linear("A", "B", type = "doublerep1", num_modules = 6),
  bblego_linear("B", "C", type = "doublerep2", num_modules = 6),
  bblego_end("C", type = "simple", num_modules = 4)
)
model <-
  initialise_model(
    num_tfs = 30,
    num_targets = 50,
    num_hks = 50,
    num_cells = 1000,
    backbone = backbone,
    verbose = TRUE,
    download_cache_dir = "~/.cache/dyngen",
    num_cores = 8,
    distance_metric = "pearson",
    tf_network_params = tf_network_default(min_tfs_per_module = 1, sample_num_regulators = function() 1),
    simulation_params = simulation_default(num_simulations = 10, census_interval = .01)
  ) %>%
  generate_tf_network() %>%
  generate_feature_network() %>%
  generate_kinetics()


# Show GRN ----------------------------------------------------------------
plot_backbone_statenet(model)
plot_backbone_modulenet(model)
plot_feature_network(model, show_targets = FALSE)
plot_feature_network(model)
plot_feature_network(model, show_hks = TRUE)


# Generate 2 variations with slightly randomised kinetics -----------------
model1 <- model %>%
  generate_gold_standard() %>%
  generate_cells()

model2 <- model
# flip effect
ix <- sample.int(nrow(model$feature_network), 5)
model2$feature_network$effect[ix] <- 1 - model2$feature_network$effect[ix]

model2$simulation_system$reactions <- dyngen:::.kinetics_generate_formulae(model2)
generate_kinetics
model2$simulation_system$reactions
model2$simulation_system$initial_state
model2 <- model2 %>%
  generate_gold_standard() %>%
  generate_cells()


# Aggregate back into model -----------------------------------------------
model12 <- model
model12$feature_info <- model$feature_info %>%
  mutate(
    w = paste0("w_", feature_id),
    x = paste0("x_", feature_id),
    y = paste0("y_", feature_id)
  )
m1gs <- model1$gold_standard
m2gs <- model2$gold_standard
model12$gold_standard <- list(
  meta = bind_rows(
    m1gs$meta %>% mutate_at(c("from", "to", "from_", "to_"), ~paste0("left_", .)),
    m2gs$meta %>% mutate_at(c("from", "to", "from_", "to_"), ~paste0("right_", .))
  ),
  counts = rbind(
    m1gs$counts,
    m2gs$counts
  ),
  network = bind_rows(
    m1gs$network %>% mutate_at(c("from", "to"), ~paste0("left_", .)),
    m2gs$network %>% mutate_at(c("from", "to"), ~paste0("right_", .))
  )
)
m1sim <- model1$simulations
m2sim <- model2$simulations
num_m1sim <- max(m1sim$meta$simulation_i)
model12$simulations <- list(
  meta = bind_rows(
    m1sim$meta %>% mutate_at(c("from", "to"), ~paste0("left_", .)),
    m2sim$meta %>% mutate_at(c("from", "to"), ~paste0("right_", .)) %>% mutate(simulation_i = simulation_i + num_m1sim)
  ),
  counts = rbind(
    m1sim$counts,
    m2sim$counts
  ),
  regulation = rbind(
    m1sim$regulation,
    m2sim$regulation
  ),
  reaction_firings = rbind(
    m1sim$reaction_firings,
    m2sim$reaction_firings
  ),
  reaction_propensities = rbind(
    m1sim$reaction_propensities,
    m2sim$reaction_propensities
  )
)
model12 <- model12 %>%
  dyngen:::calculate_dimred(dimred_premrna = FALSE) %>%
  generate_experiment()

plot_gold_simulations(model12)

# View result -------------------------------------------------------------
dataset <- wrap_dataset(model12)

library(dynplot)
g <- plot_dimred(dataset)
ggsave(paste0(out_dir, "traj_dimred.pdf"), g, width = 8, height = 6)
g <- plot_graph(dataset)
ggsave(paste0(out_dir, "traj_graph.pdf"), g, width = 8, height = 6)
g <- plot_heatmap(dataset, features_oi = 40)
ggsave(paste0(out_dir, "traj_heatmap.pdf"), g, width = 8, height = 6)




show_tfs <- show_targets <- show_hks <- TRUE

# get feature info
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

# remove unwanted features and convert NA into "NA"
feature_info <-
  feature_info %>%
  filter(
    show_tfs & is_tf |
      show_targets & !is_tf & !is_hk |
      show_hks & is_hk
  ) %>%
  mutate(
    color_by = ifelse(is.na(color_by), "NA", color_by)
  )

# filter feature network
fn1 <- model1$feature_network
fn2 <- model2$feature_network
feature_network <-
  bind_rows(
    inner_join(fn1, fn2, by = colnames(fn1)) %>% mutate(type = "common"),
    anti_join(fn1, fn2, by = colnames(fn1)) %>% mutate(type = "normal"),
    anti_join(fn2, fn1, by = colnames(fn1)) %>% mutate(type = "perturbed")
  )

# add extra edges invisible between regulators from the same module
feature_network <-
  bind_rows(
    feature_network,
    feature_info %>%
      filter(is_tf) %>%
      select(module_id, feature_id) %>%
      group_by(module_id) %>%
      do({
        crossing(from = .$feature_id, to = .$feature_id) %>%
          mutate(effect = -2)
      }) %>%
      ungroup() %>%
      filter(from < to)
  )

library(tidygraph)
library(ggraph)

gr <- tbl_graph(nodes = feature_info, edges = feature_network)
layout <- igraph::layout_with_fr(gr) %>%
  dynutils::scale_minmax() %>%
  magrittr::set_rownames(feature_info$feature_id) %>%
  magrittr::set_colnames(c("x", "y")) %>%
  as.data.frame()
gr <- gr %>% activate(edges) %>% filter(is.na(effect) | effect != -2)

cap <- circle(2.5, "mm")
str <- .2

arrow_up <- grid::arrow(type = "closed", angle = 30, length = grid::unit(3, "mm"))
arrow_down <- grid::arrow(type = "closed", angle = 89, length = grid::unit(3, "mm"))

g <-
  ggraph(gr, layout = "manual", x = layout$x, y = layout$y) +
  geom_edge_fan(aes(filter = !is.na(effect) & effect >= 0 & from != to, colour = type, width = type), arrow = arrow_up, start_cap = cap, end_cap = cap) +
  geom_edge_fan(aes(filter = !is.na(effect) & effect < 0, colour = type, width = type), arrow = arrow_down, start_cap = cap, end_cap = cap) +
  geom_node_point(aes()) +
  theme_graph(base_family = "Helvetica") +
  scale_edge_color_manual(values = c(common = "darkgray", normal = "darkgreen", perturbed = "darkred")) +
  scale_edge_width_manual(values = c(common = .5, normal = 2, perturbed = 2)) +
  coord_equal() +
  labs(edge_color = "Interaction type", edge_width = "Interaction type")

ggsave(paste0(out_dir, "feature_network.pdf"), g, width = 8, height = 6)
