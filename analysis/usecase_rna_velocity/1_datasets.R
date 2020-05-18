library(tidyverse)
library(rlang)
library(dyngen)
library(dyngen.manuscript)
library(dynplot2)

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

backbones <- dyngen::list_backbones()
backbones$disconnected <- function() dyngen::backbone_disconnected("linear", "cycle")

pwalk(design_datasets, function(id, seed, backbone_name, tr_rate_multiplier, ...) {
  dataset_file <- exp$dataset_file(id)
  plot_file <- gsub("\\.rds$", ".pdf", dataset_file)

  if (!file.exists(plot_file)) {

    cat("## Generating ", id, "\n", sep = "")
    set.seed(seed)

    kinetic_params <- kinetics_default()
    kinetic_params$sampler_tfs <- function(...) {
      x <- kinetics_default()$sampler_tfs(...)
      x$transcription_rate <- x$transcription_rate * tr_rate_multiplier
      x
    }

    backbone <- backbones[[backbone_name]]()
    model <-
      initialise_model(
        id = id,
        num_tfs = nrow(backbone$module_info),
        num_targets = 70,
        num_hks = 0,
        backbone = backbone,
        num_cells = 2500,
        kinetics_params = kinetic_params,
        simulation_params = simulation_default(
          burn_time = simtime_from_backbone(backbone, burn = TRUE) * 1.5,
          census_interval = ifelse(!backbone_name %in% c("linear_simple"), 10, 1),
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


    # add dimred
    dataset <- read_rds(dataset_file)

    dimred <- dyndimred::dimred_mds(dataset$expression, distance_method = "pearson", ndim = 3)

    dataset <- dataset %>%
      dynwrap::add_dimred(dimred = dimred) %>%
      dynwrap::add_waypoints()
    pl <- dynplot_dimred(dataset) +
      geom_cell_point(aes(color = milestone_percentages), size = 1) +
      scale_milestones_colour() +
      geom_trajectory_segments(size = 1, color = "#333333") +
      geom_milestone_label(aes(label = label), color = "black", fill = "#EEEEEE") +
      theme_common(legend.position = "none")
    ggsave(plot_file, pl, width = 6, height = 5)
    write_rds(dataset, dataset_file, compress = "gz")

    gc()
  }
})

# # recalculate dimreds
# pwalk(design_datasets, function(id, ...) {
#   cat(id, "\n", sep = "")
#   dataset_file <- exp$dataset_file(id)
#   plot_file <- gsub("\\.rds$", ".pdf", dataset_file)
#
#   if (!file.exists(dataset_file)) return(NULL)
#   dataset <- read_rds(dataset_file)
#
#   dimred <- dyndimred::dimred_mds(dataset$expression, distance_method = "pearson", ndim = 3)
#
#   dataset <- dataset %>%
#     dynwrap::add_dimred(dimred = dimred)
#   pl <- dynplot_dimred(dataset) +
#     geom_cell_point(aes(color = milestone_percentages), size = 1) +
#     scale_milestones_colour() +
#     geom_trajectory_segments(size = 1, color = "#333333") +
#     geom_milestone_label(aes(label = label), color = "black", fill = "#EEEEEE") +
#     theme_common(legend.position = "none")
#   ggsave(plot_file, pl, width = 6, height = 5)
#   write_rds(dataset, dataset_file, compress = "gz")
#
#   gc()
# })


# # write counts
# pwalk(design_datasets, function(id, ...) {
#   cat(id, "\n", sep = "")
#   dataset <- read_rds(exp$dataset_file(id))
#   dataset$counts %>% as.matrix %>% as.data.frame() %>% rownames_to_column("cell_id") %>% write_tsv(exp$temporary(paste0("counts/", id, "_spliced.tsv")))
#   dataset$counts_unspliced %>% as.matrix %>% as.data.frame() %>% rownames_to_column("cell_id") %>% write_tsv(exp$temporary(paste0("counts/", id, "_unspliced.tsv")))
#   gc()
# })

