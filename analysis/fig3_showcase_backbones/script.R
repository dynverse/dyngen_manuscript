library(tidyverse)
library(dyngen)
library(dyngen.manuscript)
library(dynplot)
library(patchwork)

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
  plot_file <- exp$temporary("backbone_plots/", id, ".pdf")

  if (!file.exists(plot_file)) {
    if (backbone_name == "disconnected") {
      back <- backbone_disconnected(
        left_backbone = backbone_cycle(),
        right_backbone = backbone_linear()
      )
    } else {
      back <- backs[[backbone_name]]()
    }
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
    dimred <- dyndimred::dimred_mds(dataset$expression, distance_method = "pearson")
    g <- plot_dimred(dataset, dimred = dimred)
    ggsave(plot_file, g, width = 4, height = 4)

    gc()
  }
})


ids <- c(
  "bb_linear_1", "bb_bifurcating_1", "bb_converging_2",
  "bb_cycle_1", "bb_bifurcating_loop_1", "bb_bifurcating_converging_3",
  "bb_bifurcating_cycle_3", "bb_consecutive_bifurcating_3", "bb_disconnected_2"
)
all(ids %in% grid$id)
plots <- map(ids, function(id) {
  dataset <- read_rds(exp$dataset_file(id))
  dimred <- dyndimred::dimred_mds(dataset$expression, distance_method = "pearson")
  plot_dimred(dataset, dimred = dimred, size_milestones = 4, size_cells = 2)
})
g <- wrap_plots(plots, ncol = 3) + plot_annotation(tag_levels = "A")
ggsave(exp$result("overview.pdf"), g, width = 8, height = 6, useDingbats = FALSE)
ggsave(exp$result("overview.png"), g, width = 8, height = 6)
