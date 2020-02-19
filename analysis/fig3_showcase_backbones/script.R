library(tidyverse)
library(dyngen)
library(dynplot)

output_folder <- "results/fig3_showcase_backbones/"

backs <- list_backbones()

oks <- c("linear", "bifurcating", "converging", "cycle", "bifurcating_loop", "consecutive_bifurcating")
oks <- c(oks, "bifurcating_cycle", "bifurcating_converging", "binary_tree", "branching", "disconnected")
oks <- c("binary_tree", "branching", "disconnected")

for (seed in seq_len(10))
for (nb in oks) {
  cat("=============== SIMULATING ", nb, " seed ", seed, " ===============\n", sep = "")
  set.seed(seed)

  out_dir <- paste0(output_folder, "/backbone_plots/bb_", nb, "_", seed, "/")
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  if (!file.exists(paste0(out_dir, "traj_dimred.pdf"))) {
    back <- backs[[nb]]()
    model <-
      initialise_model(
        num_tfs = nrow(back$module_info) * 2,
        num_targets = 100,
        num_hks = 20,
        num_cells = 1000,
        backbone = back,
        verbose = TRUE,
        download_cache_dir = "~/.cache/dyngen",
        num_cores = 4,
        gold_standard_params = gold_standard_default(
          tau = .001,
          census_interval = .02
        ),
        simulation_params = simulation_default(
          ssa_algorithm = GillespieSSA2::ssa_etl(tau = .001),
          num_simulations = ifelse(nb %in% c("binary_tree", "branching", "disconnected"), 80, 16),
          census_interval = .05,
          burn_time = simtime_from_backbone(back),
          total_time = simtime_from_backbone(back)
        )
      )

    out <- generate_dataset(model, make_plots = TRUE)

    # write_rds(out, paste0(out_dir, "/out.rds"), compress = "gz")
    ggsave(paste0(out_dir, "plot.pdf"), out$plot, width = 20, height = 15)

    dataset <- out$dataset
    dimred <- dyndimred::dimred_landmark_mds(dataset$expression, distance_method = "pearson")
    g <- plot_dimred(dataset, dimred = dimred)
    ggsave(paste0(out_dir, "traj_dimred.pdf"), g, width = 4, height = 4)

    gc()
  }
}
