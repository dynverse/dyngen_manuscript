library(tidyverse)
library(rlang)
library(dyngen)
library(dyngen.manuscript)

exp <- start_analysis("usecase_rna_velocity_legacydyngen")

# devtools::install_github("dynverse/dyngen@temp")

# setup dataset design
design_datasets <- exp$result("design_datasets.rds") %cache% {
  crossing(
    seed = 1:3,
    backbone_name = names(dyngentemp::list_backbones()) %>%
      setdiff(c("branching", "binary_tree", "consecutive_bifurcating", "disconnected", "trifurcating"))
  ) %>%
    mutate(id = paste0(backbone_name, "_", seed))
}

#' @examples
#' design_datasets %>% dynutils::extract_row_to_list(1) %>% list2env(.GlobalEnv)



sampler_tfs <-
  function(feature_info, feature_network, cache_dir = NULL, verbose = FALSE) {
    feature_info %>% mutate(
      transcription_rate = transcription_rate %|% rnorm_bounded(n(), 50, 10, min = 10),
      translation_rate = translation_rate %|% rnorm_bounded(n(), 5, 1, min = 1),
      mrna_halflife = mrna_halflife %|% rnorm_bounded(n(), .15, .03, min = .05),
      protein_halflife = protein_halflife %|% rnorm_bounded(n(), .25, .05, min = .1),
      independence = independence %|% runif(n(), 0, 1),
      splicing_rate = splicing_rate %|% rnorm_bounded(n(), 5, 1, min = 1)
    )
  }

sampler_interactions <-
  function(feature_info, feature_network, cache_dir = NULL, verbose = FALSE) {
    feature_network %>%
      mutate(
        effect = effect %|% sample(c(-1L, 1L), n(), replace = TRUE, prob = c(.25, .75)),
        strength = strength %|% 10 ^ runif(n(), log10(1), log10(100)),
        hill = hill %|% rnorm_bounded(n(), 2, 2, min = 1, max = 10)
      )
  }

pwalk(design_datasets, function(id, seed, backbone_name) {
  if (!file.exists(exp$dataset_file(id))) {

    cat("## Generating ", id, "\n", sep = "")
    set.seed(seed)
    backbone <- dyngentemp::list_backbones()[[backbone_name]]()
    model <-
      initialise_model(
        id = id,
        num_tfs = 50,
        num_targets = 70,
        num_hks = 30,
        backbone = backbone,
        num_cells = 1000,
        kinetics_params = lst(sampler_tfs, sampler_nontfs = sampler_tfs, sampler_interactions),
        gold_standard_params = gold_standard_default(
          tau = .001,
          census_interval = .05
        ),
        simulation_params = simulation_default(
          ssa_algorithm = GillespieSSA2::ssa_etl(tau = .001),
          census_interval = .05,
          compute_propensity_ratios = TRUE,
          experiment_params = simulation_type_wild_type(
            num_simulations = ifelse(backbone_name %in% c("binary_tree", "branching", "disconnected"), 40, 16)
          )
        ),
        num_cores = 6,
        download_cache_dir = "~/.cache/dyngen",
        verbose = TRUE
      )
    generate_dataset(
      model,
      output_dir = exp$dataset_folder(id),
      make_plots = TRUE,
      store_propensity_ratios = TRUE
    )

    gc()
  }
})




# make some plots
library(dynwrap)
library(dynplot2)

dataset <- read_rds(exp$dataset_file("cycle_1"))
model <- read_rds(exp$model_file("cycle_1"))

dataset <- dataset %>% add_dimred(dyndimred::dimred_landmark_mds, pair_with_velocity = FALSE)

# dimred of dataset
dynplot_dimred(dataset) +
  geom_cell_point()

# features
feature_ids <- model$feature_info %>% group_by(module_id) %>% sample_n(1) %>% pull(feature_id)

dynplot_dimred(dataset) +
  geom_cell_point(aes(color = select_feature_expression(feature_id, .data))) +
  facet_wrap_data(feature_id = feature_ids) +
  scale_expression_colour()

# Plot the spliced vs unspliced changes
feature_info <-
  model$feature_info %>%
  sample_n(1)
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
