library(tidyverse)
library(dyngen)
library(dyngen.manuscript)

exp <- start_analysis("usecase_rna_velocity")

# setup dataset design
design_datasets <- exp$result("design_datasets.rds") %cache% {
  crossing(
    seed = 1:3,
    backbone = c("linear", "linear_simple", "bifurcating", "cycle", "disconnected")
  ) %>%
    mutate(id = paste0(backbone, "_", seed))
}

#' @examples
#' design_datasets %>% dynutils::extract_row_to_list(1) %>% list2env(.GlobalEnv)

pwalk(design_datasets, function(id, seed, backbone) {
  if (!file.exists(exp$dataset_file(id))) {
    cat("## Generating ", id, "\n", sep = "")
    set.seed(seed)
    backbone <- dyngen::list_backbones()[[backbone]]()
    model <-
      initialise_model(
        id = id,
        num_tfs = 20,
        num_targets = 50,
        num_hks = 15,
        backbone = backbone,
        num_cells = 1000,
        simulation_params = simulation_default(
          compute_propensity_ratios = TRUE
        ),
        num_cores = 4,
        download_cache_dir = "~/.cache/dyngen",
        verbose = TRUE
      )
    generate_dataset(
      model,
      output_dir = exp$dataset_folder(id),
      make_plots = TRUE,
      store_propensity_ratios = TRUE
    )
  }

  gc()
})




# make some plots
library(dynwrap)
library(dynplot2)

dataset <- read_rds(exp$dataset_file("linear_simple_1"))
model <- read_rds(exp$model_file("linear_simple_1"))

dataset <- dataset %>% add_dimred(dyndimred::dimred_landmark_mds, pair_with_velocity = FALSE)

# dimred of dataset
dynplot_dimred(dataset) +
  geom_cell_point()

# features
feature_ids <- model$feature_info %>% group_by(module_id) %>% slice(1) %>% pull(feature_id)

dynplot_dimred(dataset) +
  geom_cell_point(aes(color = select_feature_expression(feature_id, .data))) +
  facet_wrap_data(feature_id = feature_ids) +
  scale_expression_colour()

# Plot the spliced vs unspliced changes
feature_info <-
  model$feature_info %>%
  slice(80)
plotdata <-
  model$simulations$counts %>%
  as.matrix() %>%
  reshape2::melt(varnames = c("step_ix", "molecule"), value.name = "expression") %>%
  as_tibble() %>%
  inner_join(
    feature_info %>%
      select(feature_id, mol_mrna, mol_premrna) %>%
      pivot_longer(-feature_id, names_to = "mol_type", values_to = "molecule"),
    by = "molecule"
  ) %>%
  inner_join(
    model$simulations$meta %>%
      mutate(step_ix = row_number()),
    by = "step_ix"
  )
plotdata %>%
  select(-molecule) %>%
  pivot_wider(names_from = "mol_type", values_from = "expression") %>%
  ggplot(aes(mol_mrna, mol_premrna)) +
  geom_path(aes(color = sim_time, group = simulation_i)) +
  facet_wrap(~simulation_i) +
  viridis::scale_color_viridis()


plotdata %>%
  filter(simulation_i == 1) %>%
  ggplot(aes(sim_time, expression)) + geom_line(aes(color = molecule))
