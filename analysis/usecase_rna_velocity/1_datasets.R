library(tidyverse)
library(rlang)
library(dyngen)
library(dyngen.manuscript)

exp <- start_analysis("usecase_rna_velocity")

# file.remove(exp$result("design_datasets.rds"))

# setup dataset design
design_datasets <- exp$result("design_datasets.rds") %cache% {
  crossing(
    seed = 1:3,
    backbone_name = names(list_backbones()),
    difficulty = forcats::fct_inorder(c("easy", "medium", "hard"))
  ) %>%
    mutate(
      tr_rate_multiplier = c(easy = 25, medium = 5, hard = 1)[difficulty],
      id = paste0(backbone_name, "_", difficulty, "_seed", seed)
    ) %>%
    filter(
      difficulty == "hard" | !backbone_name %in% c("bifurcating_converging", "bifurcating_loop", "converging", "disconnected")
    )
}

#' @examples
#' design_datasets %>% dynutils::extract_row_to_list(1) %>% list2env(.GlobalEnv)

pwalk(design_datasets, function(id, seed, backbone_name, tr_rate_multiplier, num_cells, ...) {
  if (!file.exists(exp$dataset_file(id))) {

    cat("## Generating ", id, "\n", sep = "")
    set.seed(seed)

    kinetic_params <- kinetics_default()
    kinetic_params$sampler_tfs <- function(...) {
      x <- kinetics_default()$sampler_tfs(...)
      x$transcription_rate <- x$transcription_rate * tr_rate_multiplier
      x
    }

    backbone <- dyngen::list_backbones()[[backbone_name]]()
    model <-
      initialise_model(
        id = id,
        num_tfs = nrow(backbone$module_info),
        num_targets = 70,
        num_hks = 0,
        backbone = backbone,
        num_cells = num_cells,
        kinetics_params = kinetic_params,
        simulation_params = simulation_default(
          census_interval = 10,
          compute_rna_velocity = TRUE,
          store_reaction_propensities = TRUE,
          experiment_params = simulation_type_wild_type(
            num_simulations = 100
          )
        ),
        num_cores = 7,
        download_cache_dir = "~/.cache/dyngen",
        verbose = TRUE
      )
    generate_dataset(
      model,
      output_dir = exp$dataset_folder(id),
      make_plots = TRUE
    )

    gc()
  }
})


# pwalk(design_datasets, function(id, seed, backbone_name, tr_rate_multiplier) {
#   cat(id, "\n", sep = "")
#   if (!file.exists(exp$model_file(id))) return(NULL)
#   model <- read_rds(exp$model_file(id))
#   dataset <- read_rds(exp$dataset_file(id))
#
#   reac_prop <- model$simulations$reaction_propensities
#   transcription_prop <- reac_prop[, paste0("transcription_", dataset$feature_ids)]
#   degradation_prop <- reac_prop[, paste0("premrna_degradation_", dataset$feature_ids)] + reac_prop[, paste0("mrna_degradation_", dataset$feature_ids)]
#   colnames(transcription_prop) <- colnames(degradation_prop) <- dataset$feature_ids
#
#   groundtruth_velocity <- as(transcription_prop - degradation_prop, "dgCMatrix")
#   model$simulations$rna_velocity <- groundtruth_velocity
#   groundtruth_velocity_exp <- groundtruth_velocity[model$experiment$cell_info$step_ix, ]
#   rownames(groundtruth_velocity_exp) <- model$experiment$cell_info$cell_id
#   model$experiment$rna_velocity <- groundtruth_velocity_exp
#   dataset$rna_velocity <- groundtruth_velocity_exp
#
#   write_rds(dataset, exp$dataset_file(id), compress = "gz")
#   write_rds(model, exp$model_file(id), compress = "gz")
# })

library(dynplot2)
# id <- "bifurcating_cycle_medium_seed1_2500"
id <- "bifurcating_loop_hard_seed1_2500"
pwalk(design_datasets, function(id, ...) {
  cat(id, "\n", sep = "")
  if (!file.exists(exp$model_file(id))) return(NULL)
  dataset <- read_rds(exp$dataset_file(id))
  dataset <- dataset %>% dynwrap::add_waypoints(n_waypoints = 100, trafo = function(x) x * 0 + 1, recompute = TRUE)

  dimred <- dyndimred::dimred_mds(dataset$expression, ndim = 5)
  dataset <- dataset %>% dynwrap::add_dimred(dimred = dimred)
  out <- dynwrap::project_trajectory(dataset, dataset$dimred, trajectory_projection_sd = .025)
  dataset[names(out)] <- out
  # dynplot::plot_dimred(dataset)

  dynplot_dimred(dataset) +
    geom_cell_point(aes(color = milestone_percentages), size = 1) +
    scale_milestones_colour() +
    geom_trajectory_segments(size = 1, color = "#333333") +
    geom_milestone_label(aes(label = label), color = "black", fill = "#EEEEEE") +
    theme_common(legend.position = "none") +
    ggtitle("Trajectory")


  write_rds(dataset, exp$dataset_file(id), compress = "gz")
})

