library(tidyverse)
library(dyngen)
library(dynplot)

backs <- list_backbones()

oks <- c("linear", "bifurcating", "converging", "cycle", "bifurcating_loop", "consecutive_bifurcating")
oks <- c(oks, "bifurcating_cycle", "bifurcating_converging", "binary_tree", "branching", "disconnected")
oks <- c("binary_tree", "branching", "disconnected")

for (seed in seq_len(10))
for (nb in oks) {
  cat("=============== SIMULATING ", nb, " seed ", seed, " ===============\n", sep = "")
  set.seed(seed)

  out_dir <- paste0("fig/example/bb_", nb, "_", seed, "/")
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


































# Old version -------------------------------------------------------------

# library(tidyverse)
# library(dyngen)
#
# backs <- list_backbones()
#
# for (nb in names(backs)) {
#   cat("=============== SIMULATING ", nb, " ===============\n", sep = "")
#   set.seed(1)
#
#   out_dir <- paste0("fig/example/bb_", nb, "/")
#   if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
#
#   back <- backs[[nb]]()
#   model <-
#     initialise_model(
#       num_tfs = nrow(back$module_info) * 2,
#       num_targets = 100,
#       num_hks = 100,
#       num_cells = 500,
#       backbone = back,
#       verbose = TRUE,
#       download_cache_dir = "~/.cache/dyngen",
#       num_cores = 8,
#       simulation_params = simulation_default(ssa_algorithm = GillespieSSA2::ssa_etl(), num_simulations = 50, census_interval = .01)
#     )
#
#   out <- generate_dataset(model, output_dir = out_dir, make_plots = TRUE)
#
#   # recalc dimred
#   out <- lst(
#     model = read_rds(paste0(out_dir, "model.rds")),
#     dataset = read_rds(paste0(out_dir, "dataset.rds"))
#   )
#   gs_counts <- out$model$gold_standard$counts %>% as.matrix
#   sim_counts <- out$model$simulations$counts[,colnames(gs_counts)] %>% as.matrix
#   counts <- rbind(gs_counts, sim_counts)
#   space <- dyndimred::dimred_landmark_mds(counts)
#   gs_ix <- seq_len(nrow(gs_counts))
#   out$model$gold_standard$dimred <- space[gs_ix,]
#   out$model$simulations$dimred <- space[-gs_ix,]
#
#   # plot_gold_mappings(out$model, do_facet = FALSE)
#
#   # save output
#   write_rds(out, paste0(out_dir, "/out.rds"), compress = "gz")
#
#   # plot model
#   model <- out$model
#
#   g <- plot_backbone_statenet(model)
#   ggsave(paste0(out_dir, "statenet.pdf"), g, width = 6, height = 6)
#   g <- plot_backbone_modulenet(model)
#   ggsave(paste0(out_dir, "/modulenet.pdf"), g, width = 8, height = 6)
#
#   g <- plot_feature_network(model, show_hks = TRUE)
#   ggsave(paste0(out_dir, "/featnet.pdf"), g, width = 8, height = 6)
#   g <- plot_feature_network(model, show_targets = FALSE)
#   ggsave(paste0(out_dir, "/featnet_onlytfs.pdf"), g, width = 8, height = 6)
#
#   g <- plot_gold_simulations(model) + scale_colour_brewer(palette = "Dark2")
#   ggsave(paste0(out_dir, "gold_simulations.pdf"), g, width = 8, height = 6)
#
#   g <- plot_gold_expression(model, what = "x") # mrna
#   ggsave(paste0(out_dir, "gold_mrna.pdf"), g, width = 8, height = 6)
#   g <- plot_gold_expression(model, label_changing = FALSE) # premrna, mrna, and protein
#   ggsave(paste0(out_dir, "gold_expression.pdf"), g, width = 8, height = 6)
#
#   g <- plot_simulations(model) + theme_classic() + theme(axis.ticks = element_blank(), axis.line = element_blank(), axis.text = element_blank(), axis.title = element_blank()) + labs(colour = "Simulation\ntime")
#   ggsave(paste0(out_dir, "simulations.pdf"), g, width = 8, height = 6)
#
#   g <- plot_gold_mappings(model, do_facet = FALSE) + scale_colour_brewer(palette = "Dark2")
#   ggsave(paste0(out_dir, "simulations_mapping.pdf"), g, width = 8, height = 6)
#
#   g <- plot_simulation_expression(model, 1:4, what = "x")
#   ggsave(paste0(out_dir, "simulation_expression.pdf"), g, width = 8, height = 6)
#
#   dataset <- out$dataset
#
#   library(dynplot)
#   g <- plot_dimred(dataset, dimred = dyndimred::dimred_landmark_mds)
#   ggsave(paste0(out_dir, "traj_dimred.pdf"), g, width = 4, height = 4)
#   g <- plot_graph(dataset)
#   ggsave(paste0(out_dir, "traj_graph.pdf"), g, width = 8, height = 6)
#   g <- plot_heatmap(dataset, features_oi = 40)
#   ggsave(paste0(out_dir, "traj_heatmap.pdf"), g, width = 8, height = 6)
#
#   pheatmap::pheatmap(
#     t(as.matrix(dataset$expression)),
#     filename = paste0(out_dir, "heatmap.pdf"),
#     width = 10,
#     height = 8,
#     border_color = NA,
#     show_colnames = FALSE,
#     show_rownames = FALSE,
#     treeheight_col = 0,
#     treeheight_row = 0
#   )
# }

