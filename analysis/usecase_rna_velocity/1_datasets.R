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
    backbone_name = names(list_backbones())
  ) %>%
    mutate(id = paste0(backbone_name, "_", seed))
}

#' @examples
#' design_datasets %>% dynutils::extract_row_to_list(1) %>% list2env(.GlobalEnv)

pwalk(design_datasets, function(id, seed, backbone_name) {
  if (!file.exists(exp$dataset_file(id))) {

    cat("## Generating ", id, "\n", sep = "")
    set.seed(seed)
    backbone <- dyngen::list_backbones()[[backbone_name]]()
    model <-
      initialise_model(
        id = id,
        num_tfs = nrow(backbone$module_info),
        num_targets = 70,
        num_hks = 30,
        backbone = backbone,
        num_cells = 1000,
        simulation_params = simulation_default(
          census_interval = 10,
          compute_rna_velocity = TRUE,
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
      make_plots = TRUE,
      store_rna_velocity = TRUE
    )

    gc()
  }
})
