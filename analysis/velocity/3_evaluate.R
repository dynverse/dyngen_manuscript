library(tidyverse)

source("analysis/velocity/3_evaluate_functions.R")

list2env(design_velocity[1, ], .GlobalEnv)

# Calculate scores -----
scores <- design_velocity %>% pmap(function(method_id, params_id, dataset_id, ...) {
  model <- load_model(dataset_id)
  dataset <- load_dataset_velocity(dataset_id, method_id, params_id)
  groundtruth_velocity <- extract_groundtruth_velocity(model)

  velocity_differences <- dataset$expression_projected - dataset$expression
  velocity_differences[velocity_differences == 0] <- NA
  velocity_differences[is.na(velocity_differences)] <- runif(sum(is.na(velocity_differences)), -1e-10, 1e-10)

  pairwise_correlations <- pairwise_correlations(
    Matrix::t(velocity_differences),
    t(groundtruth_velocity)
  )
  scores <- tibble(
    feature_id = colnames(velocity_differences),
    score = pairwise_correlations
  )

  # scores <- pairwise_correlations(
  #   velocity_differences,
  #   groundtruth_velocity
  # )
  scores
}) %>%
  tibble(scores = .) %>%
  bind_cols(design_velocity)


plot_scores <- scores %>%
  select(-params) %>%
  unnest(c(scores)) %>%
  group_by(dataset_id, method_id, params_id) %>%
  summarise(score = mean(score)) %>%
  left_join(dataset_design, c("dataset_id" = "id")) %>%
  ggplot(aes(paste0(method_id, "\n", params_id), score)) +
    ggbeeswarm::geom_quasirandom(aes(color = backbone))  +
    theme_bw() +
    scale_x_discrete("") +
    scale_y_continuous("Correlation (higher is better)") +
  theme_common()

plot_scores




dataset_id <- "bifurcating_2"

transform_groundtruth_velocity <- function(x) {
  log2(scales::squish(x, c(0.1, 10), only.finite = FALSE))
}

devtools::load_all('~/thesis/projects/dynverse/dynplot2/')
plot_velocity_gs <- function(
  design_velocity_oi,
  feature_id,
  dataset = NULL,
  model = NULL
) {
  if (is.null(dataset)) dataset = load_dataset(design_velocity_oi$dataset_id[[1]])
  if (is.null(model)) model = load_model(design_velocity_oi$dataset_id[[1]])

  if (is.null(dataset$dimred)) dataset <- dataset %>% add_dimred(dyndimred::dimred_landmark_mds)
  groundtruth_velocity <- extract_groundtruth_velocity(model)

  theme_common <- function(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust= 0.5),
    ...
    ) {
      ggplot2::theme(
        legend.position = legend.position,
        plot.title = plot.title,
        plot.subtitle = plot.subtitle,
        ...
    )
  }

  plot_trajectory <- dynplot_dimred(dataset) +
    geom_cell_point(aes(color = milestone_percentages)) +
    scale_milestones_fillcolour() +
    geom_trajectory_segments(size = 2, color = "#333333") +
    ggtitle(method_id) +
    theme_common(legend.position = "none") +
    ggtitle("Trajectory")
  plot_trajectory

  gs_plot <- dynplot_dimred(dataset) +
    geom_cell_point(aes(color = transform_groundtruth_velocity(groundtruth_velocity[,feature_id]))) +
    scale_velocity_color() +
    ggtitle("Ground truth velocity") +
    theme_common()
  gs_plot

  expression_plot <- dynplot_dimred(dataset) +
    geom_cell_point(aes(color = select_feature_expression(feature_id, .data))) +
    geom_trajectory_segments() +
    scale_expression_color(breaks = c(0, 1), labels = c("min", "max")) +
    theme_common() +
    ggtitle("Expression of a gene that goes up and down")
  expression_plot

  velocity_plots <- pmap(design_velocity_oi, function(dataset_id, method_id, params_id, ...) {
    velocity <- load_velocity(dataset_id, method_id, params_id)
    dataset <- dataset %>% add_velocity(velocity = velocity)
    dynplot_dimred(dataset) +
      geom_cell_point(aes(color = select_feature_velocity(feature_id, .data))) +
      scale_velocity_color(name = "", guide = "none") +
      ggtitle(method_id, subtitle = params_id) +
      theme_common()
  })


  plot_velocity_plots <- patchwork::wrap_plots(velocity_plots, nrow = 1)
  ((plot_trajectory | expression_plot | gs_plot) / plot_velocity_plots)
}

design_velocity_oi <- design_velocity %>% filter(dataset_id == !!dataset_id)
dataset = load_dataset(design_velocity_oi$dataset_id[[1]])
model <- load_model(design_velocity_oi$dataset_id[[1]])

# find feature that is most different in score between deterministic and dynamical
scores_features <- scores %>%
  select(-params) %>%
  filter(dataset_id == !!dataset_id) %>%
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

# find feature that goes up and down!
feature_id <- find_updown_feature(model)[1]

# or just select it yourself
feature_id <- "B1_TF1"

dataset <- dataset %>% add_dimred(dyndimred::dimred_landmark_mds)
plot_velocity_gs <- plot_velocity_gs(design_velocity_oi, feature_id, dataset = dataset)
plot_velocity_gs




plot_velocity_arrows <- function(
  design_velocity_oi,
  dataset = NULL,
  model = NULL
) {
  if (is.null(dataset)) dataset = load_dataset(design_velocity_oi$dataset_id[[1]])
  if (is.null(model)) model = load_model(design_velocity_oi$dataset_id[[1]])

  if (is.null(dataset$dimred)) dataset <- dataset %>% add_dimred(dyndimred::dimred_landmark_mds)

  velocity_plots <- pmap(design_velocity_oi, function(dataset_id, method_id, params_id, ...) {
    velocity <- load_velocity(dataset_id, method_id, params_id)
    dataset <- dataset %>% add_velocity(velocity = velocity)
    dataset <- dataset %>% add_dimred(
      dataset$dimred,
      dimred_projected = embed_velocity(dataset, expression_projected = velocity$expression_projected)
    )

    dynplot_dimred(dataset) +
      geom_cell_point(aes(color = milestone_percentages)) +
      scale_milestones_fillcolour() +
      geom_velocity_arrow(
        size = 1.2,
        color = "#333333",
        stat = stat_velocity_grid(grid_bandwidth = 1),
        arrow = arrow(length = unit(0.2, "cm"))
      ) +
      # geom_velocity_arrow(size = 1, color = "white") +
      ggtitle(method_id, subtitle = params_id) +
      theme_common() +
      theme(legend.position = "none")
  })


  theme_common <- function(...) ggplot2::theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust= 0.5),
    ...
  )

  patchwork::wrap_plots(velocity_plots, nrow = 1)
}


plot_velocity_arrows <- plot_velocity_arrows(design_velocity_oi, dataset = dataset)
plot_velocity_arrows


tag_first <- function(x, tag) {
  y <- x[[1]]

  if ("patchwork" %in% class(y)) {
    x[[1]] <- tag_first(x[[1]], tag = tag)
  } else {
    x[[1]] <- x[[1]] + labs(tag = tag)
  }
  x
}


pdf("results/velocity/usecase.pdf", height = 14, width = 10, useDingbats = F)
patchwork::wrap_plots(
  plot_scores,
  plot_velocity_gs,
  (plot_velocity_arrows & theme(plot.title = element_blank(), plot.subtitle = element_blank())),
  ncol = 1,
  heights = c(1, 2, 1)
) + patchwork::plot_annotation(tag_levels = c("A"))
graphics.off()

magick::image_read_pdf("results/velocity/usecase.pdf") %>%
  magick::image_convert("png") %>%
  magick::image_write("results/velocity/usecase.png") %>%
  invisible()










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
