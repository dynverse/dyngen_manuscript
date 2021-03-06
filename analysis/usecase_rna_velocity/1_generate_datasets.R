library(tidyverse)
library(rlang)
library(dyngen)
library(dyngen.manuscript)
library(dynplot2)

exp <- start_analysis("usecase_rna_velocity")

# file.remove(exp$result("design_datasets.rds"))

# setup dataset design
design_datasets <-
  crossing(
    seed = 1:3,
    backbone_name = names(list_backbones())) %>%
    mutate(
      id = paste0(backbone_name, "_seed", seed)
    )
write_rds(design_datasets, exp$result("design_datasets.rds"), compress = "gz")

#' @examples
#' design_datasets %>% dynutils::extract_row_to_list(1) %>% list2env(.GlobalEnv)

backbones <- dyngen::list_backbones()
backbones$disconnected <- function() dyngen::backbone_disconnected("linear", "cycle")

pwalk(design_datasets, function(id, seed, backbone_name, ...) {
  dataset_file <- exp$dataset_file(id)
  plot_file <- gsub("\\.rds$", ".pdf", dataset_file)

  if (!file.exists(plot_file)) {

    cat("## Generating ", id, "\n", sep = "")
    set.seed(seed)

    backbone <- backbones[[backbone_name]]()
    model <-
      initialise_model(
        id = id,
        num_tfs = nrow(backbone$module_info),
        num_targets = 70,
        num_hks = 0,
        backbone = backbone,
        num_cells = 2500,
        simulation_params = simulation_default(
          burn_time = simtime_from_backbone(backbone, burn = TRUE) * 1.5,
          census_interval = ifelse(!backbone_name %in% c("linear_simple"), 10, 1),
          compute_rna_velocity = TRUE,
          store_reaction_propensities = TRUE,
          experiment_params = simulation_type_wild_type(
            num_simulations = 100
          )
        ),
        verbose = FALSE
      )
    generate_dataset(
      model,
      format = "dyno",
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

