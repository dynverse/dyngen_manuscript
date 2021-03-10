library(tidyverse)
library(dyngen)
library(dyngen.manuscript)

exp <- start_analysis("usecase_trajectory_alignment")

# create several different backbone generators to choose from
backbone_funs <- list(
  linear_simple_short = function() {
    bblego(
      bblego_start("A", type = "simple", num_modules = 1),
      bblego_linear("A", "B", type = "simple", num_modules = 2),
      bblego_linear("B", "C", type = "simple", num_modules = 2),
      bblego_end("C", type = "simple", num_modules = 2)
    )
  },
  linear_simple = function() {
    bblego(
      bblego_start("A", type = "simple", num_modules = 2),
      bblego_linear("A", "B", type = "simple", num_modules = 3),
      bblego_linear("B", "C", type = "simple", num_modules = 3),
      bblego_end("C", type = "simple", num_modules = 2)
    )
  },
  linear_normal = function() {
    bblego(
      bblego_start("A", type = "simple"),
      bblego_linear("A", "B", type = "doublerep1"),
      bblego_linear("B", "C", type = "doublerep2"),
      bblego_end("C", type = "simple")
    )
  },
  linear_long = function() {
    bblego(
      bblego_start("A", type = "simple", num_modules = 5),
      bblego_linear("A", "B", num_modules = sample(6:10, 1), type = "doublerep2"),
      bblego_linear("B", "C", num_modules = sample(6:10, 1), type = "doublerep1"),
      bblego_end("C", type = "simple", num_modules = 6)
    )
  }
)

# setup dataset design
design_datasets <-
  crossing(
    seed = 1:10,
    backbone_name = names(backbone_funs)
  ) %>%
    mutate(id = paste0(backbone_name, "_", seed))
write_rds(design_datasets, exp$result("design_datasets.rds"))

#' @examples
#' design_datasets %>% dynutils::extract_row_to_list(1) %>% list2env(.GlobalEnv)
pwalk(design_datasets, function(id, seed, backbone_name) {
  if (!file.exists(exp$dataset_file(id))) {

    cat("## Generating ", id, "\n", sep = "")
    set.seed(seed)

    backbone <- backbone_funs[[backbone_name]]()

    # set number of cells
    num_cells <- 2000

    # divide up genes between tfs, targets and hks
    wanted_genes <- 500
    num_tfs <- nrow(backbone$module_info) * 2
    num_targets <- round((wanted_genes - num_tfs) / 2)
    num_hks <- wanted_genes - num_tfs - num_targets

    model <-
      initialise_model(
        id = id,
        num_tfs = num_tfs,
        num_targets = num_targets,
        num_hks = num_hks,
        backbone = backbone,
        num_cells = num_cells,

        # randomise kinetics more than usual
        kinetics_params = kinetics_random_distributions(),

        # set simulation params
        simulation_params = simulation_default(
          # shorten simulation a bit
          total_time = simtime_from_backbone(backbone) / 2,

          # perform 50 simulations per model
          experiment_params = bind_rows(
            simulation_type_wild_type(num_simulations = 50)
          ),

          # add extra noise to kinetics per simulation
          kinetics_noise_function = kinetics_noise_simple(sd = .1)
        ),
        verbose = TRUE
      ) %>%
      generate_tf_network() %>%
      generate_feature_network()

    # generate model once
    model1 <- model %>%
      generate_kinetics() %>%
      generate_gold_standard() %>%
      generate_cells()

    # generate another model with different kinetics
    model2 <- model %>%
      generate_kinetics() %>%
      generate_gold_standard() %>%
      generate_cells()

    # combine models and sample cells
    comb <- combine_models(
      list(left = model1, right = model2)
    ) %>%
      generate_experiment()

    # convert to dataset
    out <- as_dyno(comb, store_dimred = FALSE)

    # save dataset to file
    write_rds(out, exp$dataset_file(id), compress = "gz")

    # visual checks
    g1 <- plot_gold_simulations(comb)
    g2 <- plot_simulations(comb)
    g3 <- plot_experiment_dimred(comb)
    prm <- bind_rows(
      dyngen:::.kinetics_extract_parameters_as_df(model1$feature_info, model1$feature_network) %>% mutate(model = "model1"),
      dyngen:::.kinetics_extract_parameters_as_df(model2$feature_info, model2$feature_network) %>% mutate(model = "model2")
    ) %>%
      spread(model, value)
    g4 <- ggplot(prm) + geom_point(aes(model1, model2)) + facet_wrap(~param, scales = "free")
    g5 <- dynplot::plot_dimred(out)
    g6 <- dynplot::plot_heatmap(out, features_oi = 100)
    g <- patchwork::wrap_plots(g1, g2, g3, g4, g5, g6, byrow = TRUE, ncol = 3, widths = 1, heights = 1)
    ggsave(paste0(exp$dataset_folder(id), "/plot.pdf"), g, width = 3*10, height = 8*2)

    gc()
  }
})
