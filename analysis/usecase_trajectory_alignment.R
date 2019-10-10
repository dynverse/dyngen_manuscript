library(tidyverse)
library(dyngen)

set.seed(1)

out_dir <- "fig/usecase/trajectory_alignment/"
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
    num_tfs = 100,
    num_targets = 200,
    num_hks = 200,
    num_cells = 1000,
    backbone = backbone,
    verbose = TRUE,
    download_cache_dir = "~/.cache/dyngen",
    num_cores = 8,
    distance_metric = "pearson",
    tf_network_params = tf_network_default(min_tfs_per_module = 2, sample_num_regulators = function() 1),
    simulation_params = simulation_default(num_simulations = 10, census_interval = .01)
  ) %>%
  generate_tf_network() %>%
  generate_feature_network()


# Show GRN ----------------------------------------------------------------
plot_backbone_statenet(model)
plot_backbone_modulenet(model)
plot_feature_network(model, show_targets = FALSE)
plot_feature_network(model)
plot_feature_network(model, show_hks = TRUE)


# Generate 2 variations with slightly randomised kinetics -----------------
model1 <- model
model1$feature_info <- model1$feature_info %>%
  mutate_at(c("ba", "ind"), ~ . * rnorm(length(.), mean = 1, sd = .01)) %>%
  mutate_at(c("ba", "ind"), ~ pmin(., 1))
model1$feature_network <- model1$feature_network %>%
  mutate_at(c("strength", "cooperativity"), ~ . * rnorm(length(.), mean = 1, sd = .01))
model1 <- model1 %>%
  generate_kinetics() %>%
  generate_gold_standard() %>%
  generate_cells()

model2 <- model
# by_backbone <- model2$feature_network$from %>% grepl("_TF", .)
# fn <- model2$feature_network
# model2$feature_network$from[by_backbone] <-
#   model2$feature_network[by_backbone, ] %>%
#   mutate(from_modmod = gsub("_TF.*", "", from)) %>%
#   group_by(from_modmod) %>%
#   mutate(
#     from = sample(from)
#   ) %>%
#   ungroup() %>%
#   select(-from_modmod) %>%
#   pull(from)
model2$feature_info <- model2$feature_info %>%
  mutate_at(c("ba", "ind"), ~ . * rnorm(length(.), mean = 1, sd = .01)) %>%
  mutate_at(c("ba", "ind"), ~ pmin(., 1))
model2$feature_network <- model2$feature_network %>%
  mutate_at(c("strength", "cooperativity"), ~ . * rnorm(length(.), mean = 1, sd = .01))
model2 <- model2 %>%
  generate_kinetics() %>%
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


