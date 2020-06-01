library(tidyverse)
library(dyngen)
library(dyngen.manuscript)
library(dynplot)
library(patchwork)

set.seed(1)

exp <- start_analysis("summary")


# FIRST ROW: illustrate dyngen pipeline -----------------------------------
model1 <- exp$result("model_branching.rds") %cache% {
  backbone <- bblego(
    bblego_start("A", type = "simple", num_modules = 1),
    bblego_linear("A", "B", type = "simple", num_modules = 3),
    bblego_branching("B", c("C", "D"), num_steps = 1),
    bblego_end("C", type = "simple", num_modules = 2),
    bblego_end("D", type = "simple", num_modules = 2)
  )
  model1 <-
    initialise_model(
      num_tfs = 12,
      num_targets = 30,
      num_hks = 15,
      num_cells = 150,
      backbone = backbone,
      verbose = TRUE,
      download_cache_dir = "~/.cache/dyngen",
      num_cores = 8,
      distance_metric = "euclidean",
      simulation_params = simulation_default(
        kinetics_noise_function = kinetics_noise_none(),
        experiment_params = simulation_type_wild_type(num_simulations = 8),
        compute_cellwise_grn = TRUE,
        compute_rna_velocity = TRUE
      )
    ) %>%
    generate_tf_network() %>%
    generate_feature_network() %>%
    generate_kinetics() %>%
    generate_gold_standard() %>%
    generate_cells() %>%
    generate_experiment()

  model1
}

theme_dr <- function() {
  theme_bw() +
    theme(
      legend.position = "none",
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_blank()
    )
}

# feature network
g1a <-
  plot_feature_network(model1, show_targets = TRUE) +
  theme_dr()

# simulation expression
molecules <- model1$feature_info %>% filter(is_tf) %>% gather(mol, val, mol_mrna) %>% pull(val)
df <- bind_cols(
  model1$simulations$meta,
  as.data.frame(as.matrix(model1$simulations$counts)[,molecules])
) %>%
  filter(sim_time >= 0, simulation_i %in% 1:2) %>%
  gather(molecule, value, one_of(molecules)) %>%
  left_join(
    model1$feature_info %>%
      select(mol_mrna, module_id) %>%
      gather(type, molecule, mol_mrna),
    by = "molecule"
  ) %>%
  group_by(module_id, sim_time, simulation_i, type) %>%
  summarise(value = mean(value)) %>%
  ungroup() %>%
  mutate(module_group = gsub("[0-9]*$", "", module_id))

g1b <-
  ggplot(df, aes(sim_time, value)) +
  geom_step(aes(linetype = type, size = type, colour = module_id)) +
  scale_size_manual(values = c(mol_premrna = .5, mol_mrna = 1, mol_protein = .5)) +
  scale_colour_manual(values = model1$backbone$module_info %>% select(module_id, color) %>% deframe) +
  facet_wrap(~paste0("Cell ", simulation_i), ncol = 1, scales = "free_y") +
  theme_bw() +
  labs(x = "Simulation time", y = "Expression") +
  theme(legend.position = "none", axis.ticks = element_blank(), axis.text = element_blank(), strip.background = element_blank(), panel.grid = element_blank())

# plot simulations dimred
g1c <-
  plot_simulations(model1) +
  theme_classic() +
  theme(axis.ticks = element_blank(), axis.line = element_blank(), axis.text = element_blank(), axis.title = element_blank()) +
  labs(colour = "Simulation\ntime") +
  theme(legend.position = "none")

# plot trajectory dimred
dataset1 <- wrap_dataset(model1)
dimred <- dyndimred::dimred_mds(dataset1$expression, distance_method = "pearson")
g1d <- plot_dimred(dataset1, dimred = dimred, size_milestones = 4, size_cells = 2) +
  theme_dr()

# SECOND ROW: showcase backbones ------------------------------------------
exp2 <- start_analysis("fig3_showcase_backbones")

ids <- c(
  "bb_converging_2", "bb_cycle_1", "bb_bifurcating_converging_3", "bb_consecutive_bifurcating_3"
)
g2s <- map(ids, function(id) {
  dataset <- read_rds(exp2$dataset_file(id))
  dimred <- dyndimred::dimred_mds(dataset$expression, distance_method = "pearson")
  plot_dimred(dataset, dimred = dimred, size_milestones = 4, size_cells = 2) +
    theme_dr()
})


# THIRD ROW: illustrate dyngen pipeline -----------------------------------
exp3a <- start_analysis("usecase_trajectory_alignment")
g3a <- read_rds(exp3a$result("explanation_flat.rds"))[[4]]

exp3b <- start_analysis("usecase_rna_velocity")
g3b <- read_rds(exp3b$result("one_rna_velocity.rds"))

exp3c <- start_analysis("usecase_network_inference")
g3c <- read_rds(exp3c$result("cell1.rds")) + theme(legend.position = "none")

# Combine plots -----------------------------------------------------------
g2s[[1]] <- g2s[[1]] + labs(tag = "B")
summ <- wrap_plots(
  wrap_plots(g1a + labs(tag = "A"), g1b, g1c, g1d, nrow = 1, widths = rep(1, 4)),
  wrap_plots(g2s, nrow = 1, widths = rep(1, 4)),
  wrap_plots(g3a + labs(tag = "C"), g3b, g3c, plot_spacer(), nrow = 1, widths = rep(1, 4)),
  ncol = 1,
  heights = c(1,1,1)
)
overview_size <- 3.5
ggsave(exp$result("overview.pdf"), summ, width = overview_size * 4, height = overview_size * 3, useDingbats = FALSE)
ggsave(exp$result("overview.png"), summ, width = overview_size * 4, height = overview_size * 3)
