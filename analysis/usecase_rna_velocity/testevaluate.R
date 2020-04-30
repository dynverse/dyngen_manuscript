library(tidyverse)
library(dyngen.manuscript)

exp <- start_analysis("usecase_rna_velocity")

design_datasets <- read_rds(exp$result("design_datasets.rds"))
design_velocity <- read_rds(exp$result("design_velocity.rds"))

method_id <- "scvelo"
params_id <- "dynamical"
dataset_id <- "cycle_1"


# GROUNDTRUTH -------------------------------------------------------------
dataset <- read_rds(exp$dataset_file(dataset_id))
model <- read_rds(exp$model_file(dataset_id))

reac_prop <- model$simulations$reaction_propensities[model$experiment$cell_info$step_ix, ]
rownames(reac_prop) <- model$experiment$cell_info$cell_id
transcription_prop <- reac_prop[, paste0("transcription_", dataset$feature_ids)]
degradation_prop <- reac_prop[, paste0("premrna_degradation_", dataset$feature_ids)] + reac_prop[, paste0("mrna_degradation_", dataset$feature_ids)]
colnames(transcription_prop) <- colnames(degradation_prop) <- dataset$feature_ids
groundtruth_velocity <- transcription_prop - degradation_prop

# same:
# groundtruth_velocity <- dataset$rna_velocity


# PREDICTION --------------------------------------------------------------
velocity <- scvelo::get_velocity(dataset$expression, dataset$expression_unspliced, mode = "dynamical", var_names = "all")

# same:
# velocity <- read_rds(exp$velocity_file(dataset_id, method_id, params_id))

predicted_velocity <- velocity$velocity_vector




# COMPARE -----------------------------------------------------------------
goi <- "A1_TF1"
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
gt_scales <- min(abs(range(df$groundtruth_velocity))) %>% c(-., .)
pr_scales <- min(abs(range(df$predicted_velocity))) %>% c(-., .)

patchwork::wrap_plots(
  ggplot(df) + geom_point(aes(expr_spliced, expr_unspliced)),
  ggplot(df) + geom_point(aes(expr_spliced, expr_unspliced, colour = scales::squish(groundtruth_velocity, gt_scales))) +
    scale_color_distiller(palette = "RdBu") + labs(col = "groundtruth velocity"),
  ggplot(df) + geom_point(aes(expr_spliced, expr_unspliced, colour = scales::squish(predicted_velocity, pr_scales))) +
    scale_color_distiller(palette = "RdBu") + labs(col = "predicted velocity"),
  ggplot(df) + geom_point(aes(transcription_prop, degradation_prop)),
  ggplot(df) + geom_point(aes(transcription_prop, degradation_prop, colour = scales::squish(groundtruth_velocity, gt_scales))) +
    scale_color_distiller(palette = "RdBu") + labs(col = "groundtruth velocity"),
  ggplot(df) + geom_point(aes(transcription_prop, degradation_prop, colour = scales::squish(predicted_velocity, pr_scales))) +
    scale_color_distiller(palette = "RdBu") + labs(col = "predicted velocity"),
  nrow = 2
)

paired_simil(
  predicted_velocity,
  groundtruth_velocity,
  method = "spearman",
  margin = 2
) %>%
  enframe("feature_id", "score") %>%
  mutate(method_id, params_id, dataset_id)




# poging tot een dimred / arrow metriek


dimred <- dataset$dimred
dataset$velocity_vector <- velocity$velocity_vector
dimred_future <- scvelo::embed_velocity(dataset, dimred)
dimred_diff <- dimred_future - dimred

progression_segments <- dataset$dimred_segment_progressions
dimred_segments <- dataset$dimred_segment_points

ranges <- apply(rbind(dimred, dimred_diff), 2, range)
diffrange <- apply(ranges, 2, diff)

grid_sd <- sqrt(mean((diffrange / 15)^2)) * sqrt(3) / 2


qplot(dimred[,1], dimred[,2]) + geom_point(x = segdr[[1]], y = segdr[[2]], colour = "red", size = 10)

segdr <- dimred_segments[1,]
segdr <- dimred_segments[2000,]
apply(dimred_segments, 1, function(segdr) {
  cd <- sqrt(rowSums(sweep(dimred, 2, segdr, "-")^2))
  cw <- cd < grid_sd

  qplot(dimred[,1], dimred[,2], col = ifelse(cw, "black", "gray")) + scale_colour_identity()  + geom_point(x = segdr[[1]], y = segdr[[2]], colour = "red", data = data.frame(x = 1))

  cws <- max(1, sum(cw))

  Matrix::colSums(dimred_diff[cw,]) / cws
})

# calculate for each gaussian the smoothed arrow using a gaussian kernel
garrows <- map_dfr(grid_x, function(x) {
  # cell distances and weights to each grid point
  cd <- sqrt(outer(cell_positions$y,-grid_y,'+')^2 + (x-cell_positions$x)^2)
  cw <- cd < grid_sd

  # calculate the actual arrow
  gw <- Matrix::colSums(cw)
  cws <- pmax(1,Matrix::colSums(cw))
  gxd <- Matrix::colSums(cw*cell_positions_difference$x)/cws
  gyd <- Matrix::colSums(cw*cell_positions_difference$y)/cws

  arrow_length <- sqrt(gxd^2+gyd^2)

  tibble(
    x = x,
    y = grid_y,
    x_difference = gxd,
    y_difference = gyd,
    length = arrow_length,
    angle = atan2(y_difference, x_difference),
    mass = gw
  )
})






dynwrap::calculate_geodesic_distances()
