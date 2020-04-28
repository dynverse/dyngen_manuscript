library(tidyverse)
library(dyngen)
library(dyngen.manuscript)
library(dynplot)

exp <- start_analysis("fig3_showcase_backbones")

backs <- list_backbones()

grid <- crossing(
  backbone_name = names(backs),
  seed = 1:3
) %>% mutate(
  id = paste0("bb_", backbone_name, "_", seed)
)

#' @examples
#' grid %>% dynutils::extract_row_to_list(1) %>% list2env(.GlobalEnv)

pwalk(grid, function(backbone_name, seed, id) {
  cat("=============== SIMULATING ", backbone_name, " seed ", seed, " ===============\n", sep = "")
  set.seed(seed)

  out_dir <- exp$dataset_folder(id)
  plot_file <- exp$result("backbone_plots/", id, ".pdf")

  if (!file.exists(plot_file)) {
    back <- backs[[backbone_name]]()
    model <-
      initialise_model(
        num_tfs = nrow(back$module_info),
        num_targets = 20,
        num_hks = 0,
        num_cells = 1000,
        backbone = back,
        verbose = TRUE,
        download_cache_dir = "~/.cache/dyngen",
        num_cores = 7,
        simulation_params = simulation_default(
          census_interval = 10,
          experiment_params = simulation_type_wild_type(
            num_simulations = 100
          )
        )
      )

    generate_dataset(
      model,
      make_plots = TRUE,
      output_dir = out_dir
    )

    dataset <- read_rds(exp$dataset_file(id))
    dimred <- dyndimred::dimred_landmark_mds(dataset$expression, distance_method = "pearson")
    g <- plot_dimred(dataset, dimred = dimred)
    ggsave(plot_file, g, width = 4, height = 4)

    gc()
  }
})
