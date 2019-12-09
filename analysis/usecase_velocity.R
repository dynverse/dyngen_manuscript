library(tidyverse)
devtools::load_all('~/projects/dynverse/dynwrap/')
library(dyngen)
devtools::load_all('~/projects/dynverse/dynplot2/')
devtools::load_all('~/projects/dynverse/libraries/scvelo/')

set.seed(10)
wsr_multiplier <- 1
backbone <- backbone_linear_simple()
backbone$expression_patterns$time <- backbone$expression_patterns$time * wsr_multiplier
model <-
  initialise_model(
    num_tfs = 12,
    num_targets = 30,
    num_hks = 15,
    backbone = backbone,
    verbose = TRUE,
    num_cells = 1000,
    download_cache_dir = "~/.cache/dyngen",
    # kinetics_params = kinetics_default(sample_wsr = function(n) rnorm(n, 2, 1) %>% pmax(0.5)),
    gold_standard_params = gold_standard_default(
      # tau = .001 * wsr_multiplier, census_interval = 0.01 * wsr_multiplier
    ),
    simulation_params = simulation_default(
      # ssa_algorithm = GillespieSSA2::ssa_etl(tau = 0.001 * wsr_multiplier),
      census_interval = .01 * wsr_multiplier,
      burn_time = simtime_from_backbone(backbone, burn = TRUE),
      total_time = simtime_from_backbone(backbone),
      perform_dimred = F,
      experiment_params = simulation_type_wild_type(num_simulations = 1)
    )
  )
model <- generate_tf_network(model)
model <- generate_feature_network(model)
model <- generate_kinetics(model)
model <- generate_gold_standard(model)
model <- generate_cells(model)
model <- generate_experiment(model)
dataset <- wrap_dataset(model)



# plot_gold_simulations(model)
# plot_simulations(model)

plot_backbone_modulenet(model)

dataset <- dataset %>% add_velocity(mode = "dynamical")

# Plot the spliced vs unspliced changes
# feature_id <- model$feature_info$feature_id[[2]]
feature_id <- names(sort(apply(dataset$expression, 2, sd), decreasing = T)) %>% .[11]
meta <- model$simulations$meta %>%
  mutate(step_ix = row_number()) %>%
  filter(simulation_i == 1) %>%
  identity()
counts <- model$simulations$counts[meta$step_ix, ]
plotdata <- bind_rows(
  tibble(
    expression =  counts[,paste0("w_", feature_id)],
    step_ix = meta$step_ix,
    simulaton_i = meta$simulation_i,
    sim_time = meta$sim_time,
    molecule = "unspliced"
  ),
  tibble(
    expression =  counts[,paste0("x_", feature_id)],
    step_ix = meta$step_ix,
    simulaton_i = meta$simulation_i,
    sim_time = meta$sim_time,
    molecule = "spliced"
  )
)
plotdata %>%
  pivot_wider(names_from = "molecule", values_from = "expression") %>%
  ggplot(aes(spliced, unspliced)) + geom_path(aes(color = sim_time))

# plot dataset
plotdata <- left_join(
  dataset$dimred %>% as.data.frame() %>% rownames_to_column("cell_id"),
  tibble(
    cell_id = dataset$cell_ids,
    spliced = dataset$expression[, feature_id],
    unspliced = dataset$expression_unspliced[, feature_id]
  )
) %>%
  left_join(
    dataset %>% group_onto_trajectory_edges() %>% enframe("cell_id", "group")
  )

plotdata %>%
  ggplot(aes(spliced, unspliced)) +
  geom_point(aes(color = group)) +
  coord_equal()


# Basic plotting --------------------------------------------------------------------------------------------------
datasets_after_velocity <- map(

)

dataset <- dataset %>% add_velocity(mode = "dynamical")
dataset <- dataset %>% add_dimred(dyndimred::dimred_landmark_mds, pair_with_velocity = F)
dataset <- dataset %>% add_dimred(dimred = dataset$dimred, dimred_projected = embed_velocity(dataset$velocity$scvelo, dataset$dimred))

# Plot the dimensionality reduction
feature_id <- model$feature_info$feature_id[[10]]
dynplot_dimred(dataset) +
  geom_cell_point(color = "grey") +
  # geom_cell_point(aes(color = select_feature_expression(feature_id, d = .data))) + scale_expression_fillcolour() +
  geom_cell_point(aes(color = select_feature_velocity(feature_id, d = .data))) + scale_velocity_fillcolor() +
  geom_velocity_arrow(color = "black", size = 1, stat = stat_velocity_grid(grid_bandwidth=1.5)) +
  # geom_trajectory_segments() +
  NULL

dynplot_dimred(dataset) +
  geom_cell_point(color = "grey") +
  # geom_cell_point(aes(color = select_feature_expression(feature_id, d = .data))) + scale_expression_fillcolour() +
  geom_cell_point(aes(color = select_feature_expression(feature_id, d = .data))) + scale_expression_fillcolour() +
  geom_velocity_arrow(color = "black", size = 1, stat = stat_velocity_grid(grid_bandwidth=1.5)) +
  geom_trajectory_segments() +
  NULL



plot_backbone_modulenet(model)



# Evaluation --------------------------------------------------------------
# get for each cell its state at time + step
get_scvelo_difference <- function(dataset, ...) {
  velocity <- get_velocity(dataset$expression, dataset$expression_unspliced, ...)
  differences <- velocity$layers[["velocity"]]
  dimnames(differences) <- dimnames(dataset$expression)
  differences[is.na(differences)] <- 0
  differences
}

# expression_future <- find_future_expression(model, 0.5)
expression_current <- dataset$expression
# expression_current <- find_future_expression(dataset, 0)
expression_future <- find_future_expression(dataset, 0.2)
difference_future <- expression_future - expression_current

velocity_store <- storr::storr(storr::driver_rds("output/velocities"))

designs <- list()
designs$current <- list(
  method_id = "current",
  param_id = "default",
  params = list()
)
designs$unspliced <- list(
  method_id = "unspliced",
  param_id = "default",
  params = list()
)
designs$future <- list(
  method_id = "future",
  param_id = "default",
  params = list(step = 0.5)
)
designs$scvelo_deterministic <- list(
  method_id = "scvelo",
  param_id = "deterministic",
  params = list(mode = "deterministic")
)
designs$scvelo_dynamical <- list(
  method_id = "scvelo",
  param_id = "dynamical",
  params = list(mode = "dynamical")
)
designs$scvelo_stochastic <- list(
  method_id = "scvelo",
  param_id = "stochastic",
  params = list(mode = "stochastic")
)

design <- dynutils::list_as_tibble(designs) %>%
  mutate(id = paste0(method_id, "-", param_id))
models <- design

models <- models %>%
  mutate(model = pmap(lst(method_id, param_id, params), run_velocity, dataset = dataset, store = velocity_store))

# velocity_store$clear()
# velocity_store$del(design %>% filter(method_id == "future") %>% pull(id))
# velocity_store$del(design %>% filter(method_id == "scvelo") %>% pull(id))

models$differences <- map(models$model, "differences")

## Score based on correlation between scaled expression of estimate and future expression

score <- function(x, y) {
  assertthat::assert_that(all(colnames(x) == colnames(y)))
  assertthat::assert_that(all(nrow(x) == nrow(y)))
  # map_dbl(seq_len(nrow(x)), function(i) cor(x[i, ], y[i, ], method = "spearman"))
  # map_dbl(seq_len(nrow(x)), function(i) mean(sign(x[i, ]) == sign(y[i, ])))
  map_dbl(seq_len(nrow(x)), function(i) sum(x[i, ] * y[i, ])/(sqrt(sum(x[i, ]^2)) * sqrt(sum(y[i, ]^2))))
}

scores <- tibble(
  id = models$id
)

differences_future <- models$model[[which(models$method_id == "future")]]$differences
scores$cor <- map(models$differences, score, y = differences_future)

scores %>%
  unnest(cor) %>%
  ggplot(aes(id, cor)) + geom_boxplot()#ggbeeswarm::geom_quasirandom()











cell_id <- dataset$cell_ids[[5]]
models$difference %>%
  map(~.[cell_id, ]) %>%
  do.call(rbind, .) %>%
  melt(varnames = c("stage", "feature_id"), value.name = "difference") %>%
  ggplot(aes(feature_id, stage, fill = difference)) +
  geom_raster() +
  scale_fill_distiller(palette = "RdBu", limits = c(-5, 5))



















feature_id <- "B2_TF1"

difference_scvelo <- scale(expression_scvelo) - scale(expression_current)
difference_velocyto <- scale(expression_velocyto) - scale(expression_current)
difference_future <- scale(expression_future) - scale(expression_current)

# pseudotime <- dataset %>% add_root() %>% calculate_pseudotime() %>% enframe("cell_id", "pseudotime")
pseudotime <- model$experiment$cell_info %>% select(cell_id, pseudotime = sim_time)
expression_oi <- bind_rows(
  scale(expression_current)[, feature_id] %>% enframe("cell_id", "expression") %>% mutate(dataset = "current", type = "expression"),
  scale(expression_scvelo)[, feature_id] %>% enframe("cell_id", "expression") %>% mutate(dataset = "scvelo", type = "expression"),
  scale(expression_velocyto)[, feature_id] %>% enframe("cell_id", "expression") %>% mutate(dataset = "velocyto", type = "expression"),
  scale(expression_unspliced)[, feature_id] %>% enframe("cell_id", "expression") %>% mutate(dataset = "unspliced", type = "expression"),
  scale(expression_future)[, feature_id] %>% enframe("cell_id", "expression") %>% mutate(dataset = "future", type = "expression"),
  scale(difference_future)[, feature_id] %>% enframe("cell_id", "expression") %>% mutate(dataset = "difference_future", type = "difference"),
  scale(difference_velocyto)[, feature_id] %>% enframe("cell_id", "expression") %>% mutate(dataset = "difference_velocyto", type = "difference"),
  scale(difference_scvelo)[, feature_id] %>% enframe("cell_id", "expression") %>% mutate(dataset = "difference_scvelo", type = "difference"),
) %>%
  left_join(pseudotime)

expression_oi %>%
  ggplot(aes(pseudotime, expression, color = dataset)) +
  geom_smooth(se = F) +
  facet_grid(type~.)


expression_oi %>%
  ggplot(aes(pseudotime, expression, color = dataset)) +
  geom_point() +
  facet_grid(dataset~.)


dynplot::plot_heatmap(dataset)

difference_scvelo_scaled <- scale(difference_scvelo)
difference_scvelo_scaled[is.na(difference_scvelo_scaled)] <- 0
cor(as.numeric(difference_scvelo_scaled), as.numeric(scale(difference_future)))
cor(as.numeric(scale(difference_velocyto)), as.numeric(scale(difference_future)))


## Score based on similarity in signs (up or down)
changes_future <- sign(expression_future - expression_current)
changes_scvelo <- sign(dataset$expression_projected - dataset$expression)
changes_velocyto <- Matrix::t(sign(velocity$deltaE))[rownames(expression_current), colnames(expression_current)]

table(as.numeric(changes_future), as.numeric(changes_velocyto))
table(as.numeric(changes_future), as.numeric(changes_scvelo))

sum(as.logical(changes_future == changes_scvelo))
sum(as.logical(changes_future == changes_velocyto))
