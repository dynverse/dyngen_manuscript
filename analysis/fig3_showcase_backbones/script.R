library(tidyverse)
library(dyngen)
library(dyngen.manuscript)
library(dynplot)

RcppParallel::setThreadOptions(numThreads = 6)

exp <- start_analysis("fig3_showcase_backbones")

backs <- list_backbones()

grid <- crossing(
  backbone_name = names(backs) %>% setdiff("disconnected"),
  seed = seq_len(4)
) %>% mutate(
  id = paste0("bb_", backbone_name, "_", seed)
)

#' @examples
#' backbone_name <- "linear"
#' seed <- 1
#' id <- paste0("bb_", backbone_name, "_", seed)

pwalk(grid, function(backbone_name, seed, id) {
  cat("=============== SIMULATING ", backbone_name, " seed ", seed, " ===============\n", sep = "")
  set.seed(seed)

  out_dir <- exp$dataset_folder(id)

  if (!file.exists(paste0(out_dir, "traj_dimred.pdf"))) {
    back <- backs[[backbone_name]]()
    model <-
      initialise_model(
        num_tfs = nrow(back$module_info) * 1.5,
        num_targets = 20,
        num_hks = 0,
        num_cells = 1000,
        backbone = back,
        verbose = TRUE,
        download_cache_dir = "~/.cache/dyngen",
        num_cores = 6,
        simulation_params = simulation_default(
          experiment_params = simulation_type_wild_type(
            num_simulations = ifelse(backbone_name %in% c("binary_tree", "branching", "disconnected"), 40, 16)
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
    ggsave(exp$result("backbone_plots/", id, "/traj_dimred.pdf"), g, width = 4, height = 4)

    gc()
  }
})
