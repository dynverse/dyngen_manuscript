library(tidyverse)
library(dyngen)
library(dyngen.manuscript)

exp <- start_analysis("usecase_rna_velocity")

# setup dataset design
dataset_design <- tibble(
  seed = 1:3
) %>%
  crossing(backbone = c("linear", "linear_simple", "bifurcating", "cycle", "disconnected")) %>%
  mutate(id = paste0(backbone, "_", seed))

design_velocity <- tribble(
  ~method_id, ~params, ~params_id,
  "velocyto", list(assumption = "constant_velocity"), "constant_velocity",
  "velocyto", list(assumption = "constant_unspliced"), "constant_unspliced",
  "scvelo", list(mode = "deterministic"), "deterministic",
  "scvelo", list(mode = "dynamical"), "dynamical",
  "scvelo", list(mode = "stochastic"), "stochastic",
) %>%
  crossing(
    dataset_design %>% select(dataset_id = id)
  )

write_rds(lst(dataset_design, design_velocity), exp$result("design.rds"))

# generate datasets
pwalk(dataset_design, function(id, seed, backbone) {
  out_dir <- exp$temporary("datasets/", id, "/")
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  file_dataset <- paste0(out_dir, "dataset.rds")

  if (!file.exists(file_dataset)) {
    set.seed(seed)
    backbone <- dyngen::list_backbones()[[backbone]]()

    model <-
      initialise_model(
        id = id,
        num_tfs = 20,
        num_targets = 50,
        num_hks = 15,
        backbone = backbone,
        verbose = TRUE,
        num_cells = 1000,
        download_cache_dir = "~/.cache/dyngen",
        simulation_params = simulation_default(
          perform_dimred = FALSE,
          store_reaction_propensities = TRUE
        )
      )
    model <- generate_tf_network(model)
    model <- generate_feature_network(model)
    model <- generate_kinetics(model)
    model <- generate_gold_standard(model)
    model <- generate_cells(model)
    model <- generate_experiment(model)
    dataset <- wrap_dataset(model)

    write_rds(model, paste0(out_dir, "model.rds"))
    write_rds(dataset, file_dataset)
  }
})




#
# # characterise datasets
# # devtools::load_all('~/thesis/projects/dynverse/dynwrap/')
# # devtools::load_all('~/thesis/projects/dynverse/dynplot2/')
# library(dynwrap)
# library(dynplot2)
#
# dataset <- load_dataset("linear_simple_1")
# model <- load_model("linear_simple_1")
#
# dataset <- dataset %>% add_dimred(dyndimred::dimred_landmark_mds, pair_with_velocity = F)
#
# # dimred of dataset
# dynplot_dimred(dataset) +
#   geom_cell_point()
#
# # features
# feature_ids <- model$feature_info %>% group_by(module_id) %>% slice(1) %>% pull(feature_id)
#
# dynplot_dimred(dataset) +
#   geom_cell_point(aes(color = select_feature_expression(feature_id, .data))) +
#   facet_wrap_data(feature_id = feature_ids) +
#   scale_expression_fillcolour()
#
# # Plot the spliced vs unspliced changes
# # feature_id <- model$feature_info$feature_id[[2]]
# feature_id <- model$feature_info %>% filter(!burn) %>% pull(feature_id) %>% first()
# meta <- model$simulations$meta %>%
#   mutate(step_ix = row_number()) %>%
#   # filter(simulation_i == 1)
#   identity()
# counts <- model$simulations$counts[meta$step_ix, ]
# plotdata <- bind_rows(
#   tibble(
#     expression =  counts[,paste0("w_", feature_id)],
#     step_ix = meta$step_ix,
#     simulation_i = meta$simulation_i,
#     sim_time = meta$sim_time,
#     molecule = "unspliced"
#   ),
#   tibble(
#     expression =  counts[,paste0("x_", feature_id)],
#     step_ix = meta$step_ix,
#     simulation_i = meta$simulation_i,
#     sim_time = meta$sim_time,
#     molecule = "spliced"
#   )
# )
# plotdata %>%
#   pivot_wider(names_from = "molecule", values_from = "expression") %>%
#   ggplot(aes(spliced, unspliced)) + geom_path(aes(color = sim_time)) +
#   facet_wrap(~simulation_i)
#
#
#
# plotdata %>%
#   filter(simulation_i == 1) %>%
#   ggplot(aes(step_ix, expression)) + geom_line(aes(color = molecule))
