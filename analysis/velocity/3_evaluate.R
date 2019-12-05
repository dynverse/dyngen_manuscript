library(tidyverse)

source("analysis/velocity/3_evaluate_functions.R")

list2env(design_velocity[1, ], .GlobalEnv)

scores <- design_velocity %>% pmap(function(method_id, params_id, dataset_id, ...) {
  model <- load_model(dataset_id)
  dataset <- load_dataset_velocity(dataset_id, method_id, params_id)
  groundtruth_velocity <- extract_groundtruth_velocity(model)

  velocity_differences <- dataset$expression_projected - dataset$expression
  velocity_differences[velocity_differences == 0] <- NA
  velocity_differences[is.na(velocity_differences)] <- runif(sum(is.na(velocity_differences)), -1e-10, 1e-10)

  scores <- pairwise_correlations(
    Matrix::t(velocity_differences),
    t(groundtruth_velocity)
  )
  scores <- tibble(
    feature_id = colnames(velocity_differences),
    score = scores
  )

  # scores <- pairwise_correlations(
  #   velocity_differences,
  #   groundtruth_velocity
  # )
  scores
}) %>%
  tibble(scores = .) %>%
  bind_cols(design_velocity)


scores %>%
  select(-params) %>%
  unnest(c(scores)) %>%
  # left_join(model$feature_info) %>%
  ggplot(aes(paste0(method_id, " ", params_id), score)) +
    ggbeeswarm::geom_beeswarm(aes(color = dataset_id)) +
    facet_wrap(~dataset_id)














scores_features <- scores %>%
  select(-params) %>%
  filter(dataset_id == "1") %>%
  unnest(scores)



feature_id <- scores_features %>%
  left_join(model$feature_info) %>%
  filter(method_id == "scvelo") %>%
  filter(is_tf) %>%
  pivot_wider(names_from = params_id, values_from = score) %>%
  mutate(difference = deterministic - dynamical) %>%
  arrange(desc(difference)) %>%
  pull(feature_id) %>%
  first()

feature_id <- find_updown_feature(model)[1]


velocity <- load_velocity(dataset_id, method_id, params_id)


dataset <- add_velocity(dataset, velocity = velocity)


dynplot_dimred(dataset) +
  geom_cell_point(aes(color = select_feature_velocity(feature_id, .data))) +
  scale_color_viridis_c()
dynplot_dimred(dataset) +
  geom_cell_point(aes(color = select_feature_expression(feature_id, .data))) +
  scale_color_viridis_c()


transform_groundtruth_velocity <- function(x) {
  log2(scales::squish(x, c(0.1, 10), only.finite = FALSE))
}

plot_velocity_gs <- function(
  design_velocity_oi,
  feature_id,
  dataset = load_dataset(design_velocity$dataset_id[[1]]),
  model = load_model(design_velocity$dataset_id[[1]])
) {
  dataset <- dataset %>% add_dimred(dyndimred::dimred_landmark_mds, pair_with_velocity = FALSE)
  groundtruth_velocity <- extract_groundtruth_velocity(model)

  gs_plot <- dynplot_dimred(dataset) +
    geom_cell_point(aes(color = transform_groundtruth_velocity(groundtruth_velocity[,feature_id]))) +
    scale_velocity_fillcolor(name = "")

  expression_plot <- dynplot_dimred(dataset) +
    geom_cell_point(aes(color = select_feature_expression(feature_id, .data))) +
    scale_expression_fillcolour(name = "")

  velocity_plots <- pmap(design_velocity_oi, function(dataset_id, method_id, params_id, ...) {
    velocity <- load_velocity(dataset_id, method_id, params_id)
    dataset <- dataset %>% add_velocity(velocity = velocity)
    dynplot_dimred(dataset) +
      geom_cell_point(aes(color = select_feature_velocity(feature_id, .data))) +
      scale_velocity_fillcolor(name = "") +
      ggtitle(paste0(method_id, " ", params_id))
  })

  patchwork::wrap_plots(
    c(
      list(gs_plot, expression_plot),
      velocity_plots
    )
  )
}

design_velocity_oi <- design_velocity %>% filter(dataset_id == "1")
plot_velocity_gs(design_velocity_oi, dataset$feature_ids[[6]])


compare_velocity_correlation <- function(model, dataset, velocity, ...) {
  model <- load_model(dataset_id)
  dataset <- load_dataset_velocity(dataset_id, method_id, params_id)
  groundtruth_velocity <- extract_groundtruth_velocity(model)

  velocity_differences <- dataset$expression_projected - dataset$expression
  velocity_differences[velocity_differences == 0] <- NA
  velocity_differences[is.na(velocity_differences)] <- runif(sum(is.na(velocity_differences)), -1e-10, 1e-10)

  plot(log2(groundtruth_velocity[, feature_id] %>% pmin(100) %>% pmax(0.01)), velocity_differences[, feature_id])
}

design_velocity %>% mutate(
  model = load_model(dataset_id),
  dataset = load_dataset(dataset_id),
  velocity = load_velocity(dataset_id, method_id, params_id)
)


design_velocity %>% filter(dataset_id == "1") %>% pmap(compare_velocity_correlation)









plot_velocity_grid <- function(
  design_velocity_oi,
  dataset = load_dataset(design_velocity$dataset_id[[1]])
) {
  dataset <- dataset %>% add_dimred(dyndimred::dimred_landmark_mds, pair_with_velocity = FALSE)

  velocity_plots <- pmap(design_velocity_oi, function(dataset_id, method_id, params_id, ...) {
    velocity <- load_velocity(dataset_id, method_id, params_id)
    dataset <- dataset %>% add_velocity(velocity = velocity)
    dataset <- dataset %>%
      add_dimred(dimred = dataset$dimred, dimred_projected = embed_velocity(velocity$scvelo, dataset$dimred))

    dynplot_dimred(dataset) +
      geom_cell_point(color = "grey") +
      geom_velocity_arrow(size = 1) +
      ggtitle(paste0(method_id, " ", params_id))
  })

  patchwork::wrap_plots(
    c(
      velocity_plots
    )
  )
}

plot_velocity_grid(design_velocity_oi)

































groundtruth_velocity <- extract_groundtruth_velocity(model)


plot(log2(groundtruth_velocity[, feature_id] %>% pmin(10) %>% pmax(0.1)), velocity_differences[, feature_id])


assertthat::assert_that(all(rownames(groundtruth_velocity) == rownames(velocity_differences)))
assertthat::assert_that(all(colnames(groundtruth_velocity) == colnames(velocity_differences)))

cors <- pairwise_correlations(
  Matrix::t(velocity_differences),
  t(groundtruth_velocity)
)


dataset <- dataset %>% add_dimred(dyndimred::dimred_landmark_mds, pair_with_velocity = FALSE)
dataset <- add_velocity(dataset, velocity = velocity)

qplot(cors)

feature_id <- find_updown_feature(model)[10]
feature_id <- dataset$feature_ids[10]
dynplot_dimred(dataset) +
  geom_cell_point(aes(color = pmin(log(groundtruth_velocity[,feature_id]), 5))) +
  scale_color_viridis_c()
dynplot_dimred(dataset) +
  geom_cell_point(aes(color = select_feature_expression(feature_id, .data))) +
  scale_color_viridis_c()
dynplot_dimred(dataset) +
  geom_cell_point(aes(color = select_feature_velocity(feature_id, .data))) +
  scale_color_viridis_c()
