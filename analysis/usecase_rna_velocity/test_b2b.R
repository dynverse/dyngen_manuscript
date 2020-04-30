library(tidyverse)
library(dyngen)
library(dyngen.manuscript)

exp <- start_analysis("usecase_rna_velocity")

# generate dataset
set.seed(1)
backbone <- backbone_cycle_simple()

kinetic_params <- kinetics_default()
kinetic_params$sampler_tfs <- function(feature_info, ...) {
  x <- kinetics_default()$sampler_tfs(feature_info, ...)
  x$transcription_rate <- x$transcription_rate * 100
  x
}

model <-
  initialise_model(
    id = "testb2b",
    num_tfs = nrow(backbone$module_info),
    num_targets = 0,
    num_hks = 0,
    backbone = backbone,
    num_cells = 1000,
    simulation_params = simulation_default(
      census_interval = 10,
      compute_rna_velocity = TRUE,
      store_reaction_propensities = TRUE,
      experiment_params = simulation_type_wild_type(
        num_simulations = 10
      )
    ),
    kinetics_params = kinetic_params,
    experiment_params = experiment_snapshot(
      sample_capture_rate = function(n) rep(10, n)
    ),
    num_cores = 7,
    download_cache_dir = "~/.cache/dyngen",
    verbose = TRUE
  )
out <- generate_dataset(model)
dataset <- out$dataset
model <- out$model

# GROUNDTRUTH -------------------------------------------------------------
reac_prop <- model$simulations$reaction_propensities[model$experiment$cell_info$step_ix, ]
rownames(reac_prop) <- model$experiment$cell_info$cell_id
transcription_prop <- reac_prop[, paste0("transcription_", dataset$feature_ids)]
degradation_prop <- reac_prop[, paste0("premrna_degradation_", dataset$feature_ids)] + reac_prop[, paste0("mrna_degradation_", dataset$feature_ids)]
colnames(transcription_prop) <- colnames(degradation_prop) <- dataset$feature_ids
groundtruth_velocity <- as(log2(transcription_prop + 1) - log2(degradation_prop + 1), "dgCMatrix")
# groundtruth_velocity <- transcription_prop - degradation_prop

# same:
# groundtruth_velocity <- dataset$rna_velocity


# PREDICTION --------------------------------------------------------------
velocity <- scvelo::get_velocity(
  spliced = dataset$expression,
  unspliced = dataset$expression_unspliced,
  mode = "dynamical", var_names = "all"
)
predicted_velocity <- velocity$velocity_vector


# COMPARE -----------------------------------------------------------------
plot_backbone_modulenet(model)
goi <- "M1_TF1"
df <-
  dataset$cell_info %>%
  left_join(dataset$progressions, by = "cell_id") %>%
  mutate(
    cycle_time = c("s1" = 0, "s2" = 1, "s3" = 2)[from] + percentage,
    expr_spliced = dataset$expression[,goi],
    expr_unspliced = dataset$expression_unspliced[,goi],
    transcription_prop = transcription_prop[,goi],
    degradation_prop = degradation_prop[,goi],
    groundtruth_velocity = groundtruth_velocity[,goi],
    predicted_velocity = predicted_velocity[,goi]
  )

plot(density(df$groundtruth_velocity))


gt_scales <- min(abs(range(df$groundtruth_velocity))) %>% c(-., .)
pr_scales <- min(abs(range(df$predicted_velocity))) %>% c(-., .)
patchwork::wrap_plots(
  ggplot(df) + geom_point(aes(expr_spliced, expr_unspliced)),
  ggplot(df) + geom_point(aes(expr_spliced, expr_unspliced, colour = scales::squish(groundtruth_velocity, gt_scales))) +
    scale_color_distiller(palette = "RdBu") + labs(col = "groundtruth velocity"),
  ggplot(df) + geom_point(aes(expr_spliced, expr_unspliced, colour = scales::squish(predicted_velocity, pr_scales))) +
    scale_color_distiller(palette = "RdBu") + labs(col = "predicted velocity"),
  ggplot(df) + geom_point(aes(degradation_prop, transcription_prop)),
  ggplot(df) + geom_point(aes(degradation_prop, transcription_prop, colour = scales::squish(groundtruth_velocity, gt_scales))) +
    scale_color_distiller(palette = "RdBu") + labs(col = "groundtruth velocity"),
  ggplot(df) + geom_point(aes(degradation_prop, transcription_prop, colour = scales::squish(predicted_velocity, pr_scales))) +
    scale_color_distiller(palette = "RdBu") + labs(col = "predicted velocity"),
  nrow = 2
)

ggplot(df) + geom_point(aes(groundtruth_velocity, predicted_velocity, col = cycle_time)) + viridis::scale_color_viridis()
ggplot(df %>% gather(var, val, expr_spliced:predicted_velocity) %>% mutate(var = forcats::fct_inorder(var))) +
  geom_point(aes(cycle_time, val, colour = var)) +
  facet_wrap(~var, ncol = 1, scales = "free_y") + labs(title = goi)

# simil metriek
psim <- paired_simil(
  predicted_velocity,
  groundtruth_velocity,
  method = "spearman",
  margin = 2
) %>%
  enframe("feature_id", "score")
mean(psim$score)
