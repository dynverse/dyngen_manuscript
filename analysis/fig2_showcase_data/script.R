library(tidyverse)
library(dyngen)
library(dyngen.manuscript)

set.seed(1)

exp <- start_analysis("fig2_showcase_data")

backbone <- backbone(
  module_info = tribble(
    ~module_id, ~basal, ~burn, ~independence,
    "A", 1, TRUE, 1,
    "B", 1, TRUE, 1,
    "C", 1, TRUE, 1,
    "D", 0, TRUE, 1,
    "E", 0, TRUE, 1
  ),
  module_network = tribble(
    ~from, ~to, ~effect, ~strength, ~hill,
    "A", "B", -1L, 4, 2,
    "B", "C", -1L, 4, 2,
    "C", "D", 1L, 1, 2,
    "D", "E", 1L, 1, 2,
    "E", "A", -1L, 4, 2
  ),
  expression_patterns = tribble(
    ~from, ~to, ~module_progression, ~start, ~burn, ~time,
    "Burn0", "S1", "+A,+B,+C", TRUE, TRUE, 150,
    "S1", "S2", "+D,+E,-C", FALSE, FALSE, 120,
    "S2", "S3", "+C,-A,-B", FALSE, FALSE, 120,
    "S3", "S1", "+A,+B,-D,-E", FALSE, FALSE, 120
  )
)
census_interval <- 5
total_time <- 500
time_breaks <- c(0, 250, 500)
model <-
  initialise_model(
    num_tfs = nrow(backbone$module_info),
    num_targets = 0,
    num_hks = 0,
    num_cells = 300,
    backbone = backbone,
    verbose = TRUE,
    download_cache_dir = "~/.cache/dyngen",
    num_cores = 8,
    # kinetics_params = kinetics_default(
    #   sample_wpr = function(n) 50,
    #   sample_wsr = function(n) 5,
    #   sample_ypr = function(n) 2.5,
    #   sample_whl = function(n) .25,
    #   sample_xhl = function(n) .25,
    #   sample_yhl = function(n) .25
    # ),
    simulation_params = simulation_default(
      burn_time = 240,
      total_time = total_time,
      experiment_params = simulation_type_wild_type(num_simulations = 1),
      census_interval = census_interval,
      store_reaction_firings = TRUE,
      store_reaction_propensities = TRUE,
      compute_cellwise_grn = TRUE,
      compute_log_propensity_ratios = TRUE
    )
  )

model <- model %>%
  generate_tf_network() %>%
  generate_feature_network()

notf <- function(x) gsub("_TF1$", "", x)
model$feature_info <- model$feature_info %>% mutate_at(vars(feature_id), notf)
model$feature_network <- model$feature_network %>% mutate_at(vars(from, to), notf)

model <- model %>%
  generate_kinetics() %>%
  generate_gold_standard() %>%
  generate_cells()

plot_gold_expression(model)
plot_simulations(model)
plot_gold_mappings(model, do_facet = FALSE)

plot_backbone_modulenet(model)
# floor to nearest census
model$simulations$meta$sim_time <- floor(model$simulations$meta$sim_time / census_interval) * census_interval
time <- model$simulations$meta$sim_time

mol_map <- c(mol_premrna = "pre-mRNA", mol_mrna = "mRNA", mol_protein = "protein")
reac_map <- c(
  "transcription" = "Transcription", "splicing" = "Splicing", "translation" = "Translation",
  "premrna_degradation" = "pre-mRNA degradation", "mrna_degradation" = "mRNA degradation", "protein_degradation" = "Protein degradation"
)

expr_df <- data.frame(time = time, model$simulations$counts %>% as.matrix()) %>%
  gather(var, value, -time) %>%
  mutate(
    gene = paste0("Gene ", gsub(".*_", "", var)),
    molecule = factor(mol_map[gsub("_[^_]*$", "", var)], levels = mol_map)
  ) %>%
  filter(time > 0)
firings_df <- data.frame(time = time, model$simulations$reaction_firings %>% as.matrix(), check.names = FALSE) %>%
  gather(var, value, -time) %>%
  mutate(
    gene = paste0("Gene ", gsub(".*_", "", var)),
    reaction = factor(reac_map[sub("_[^_]*$", "", var)], levels = reac_map)
  ) %>%
  filter(time > 0)
reg_df <- data.frame(time = time, model$simulations$cellwise_grn %>% as.matrix(), check.names = FALSE) %>%
  gather(var, value, -time) %>%
  mutate(
    from = gsub("(.*)->.*", "\\1", var),
    to = gsub(".*->(.*)", "\\1", var)
  ) %>%
  left_join(model$feature_network, by = c("from", "to")) %>%
  mutate(name = paste0(from, ifelse(effect > 0, " → ", " ⊣ "), to)) %>%
  filter(time > 0)
state_df <- model$simulations$meta %>%
  mutate(cell_id = paste0("cell_", row_number())) %>%
  filter(sim_time > 0)
prop_df <- data.frame(time = time, model$simulations$reaction_propensities %>% as.matrix()) %>%
  gather(var, value, -time) %>%
  mutate(
    gene = paste0("Gene ", gsub(".*_", "", var)),
    reaction = factor(reac_map[sub("_[^_]*$", "", var)], levels = reac_map)
  ) %>%
  filter(time > 0)
velocity <- data.frame(time = time, model$simulations$net_propensity %>% as.matrix()) %>%
  gather(var, value, -time) %>%
  mutate(
    gene = paste0("Gene ", var)
  ) %>%
  filter(time > 0)


dummy <-
  dynwrap::wrap_data(
    cell_ids = state_df$cell_id
  ) %>%
  dynwrap::add_trajectory(
    milestone_network = model$backbone$expression_patterns %>% filter(!burn) %>% transmute(from, to, length = time, directed = TRUE),
    progressions = state_df %>% transmute(cell_id = cell_id, from, to, percentage = time)
  )
state_df <- state_df %>%
  left_join(dummy$milestone_percentages, by = "cell_id") %>%
  full_join(crossing(sim_time = state_df$sim_time, milestone_id = dummy$milestone_ids), by = c("milestone_id", "sim_time")) %>%
  mutate(percentage = ifelse(is.na(percentage), 0, percentage))

#######
# ABUNDANCE LEVELS
#######
g1 <- ggplot(expr_df, aes(time, value)) +
  geom_step(aes(colour = gene)) +
  facet_wrap(~molecule, ncol = 1, scales = "free_y") +
  theme_classic() +
  labs(x = "Simulation time", y = "Expression", colour = "Molecule") +
  scale_colour_brewer(palette = "Set1") +
  # scale_colour_manual(values = c("pre-mRNA" = "#4daf4a", "mRNA" = "#377eb8", "protein" = "#e41a1c", "gene" = "black")) +
  theme(
    text = element_text(family = "Helvetica"),
    strip.background = element_blank(),
    strip.text = element_blank(),
    legend.margin = margin()
  ) +
  scale_x_continuous(breaks = time_breaks) +
  geom_text(aes(x = time_breaks[[2]], y = max, label = molecule), expr_df %>% group_by(molecule) %>% summarise(max = max(value)), size = 3, hjust = 0.5, vjust = 1)

#######
# CELL STATE
#######
g2 <- ggplot(state_df) +
  geom_line(aes(sim_time, percentage, colour = milestone_id), size = 1.5) +
  theme_classic() +
  scale_colour_brewer(palette = "Accent") +
  labs(x = "Simulation time", y = "Percentage", colour = "Cell state") +
  theme(
    text = element_text(family = "Helvetica"),
    legend.margin = margin()
  ) +
  scale_x_continuous(breaks = time_breaks) +
  scale_y_continuous(breaks = c(0, .5, 1))


#######
# REACTION FIRINGS
#######
firings_df_A <- firings_df %>% filter(gene == "Gene A")
label_fir_df <- firings_df_A %>% group_by(reaction) %>% summarise(time = mean(time), value = max(value)*1.1) %>% mutate(reaction_text = paste0("Gene A ", reaction))

g3 <-
  ggplot(firings_df_A, aes(time, value)) +
  facet_wrap(~reaction, ncol = 1, scales = "free_y") +
  geom_area(aes(fill = forcats::fct_rev(reaction)), position = "stack") +
  theme_classic() +
  scale_fill_brewer(palette = "Set2") +
  labs(x = "Simulation time", y = "Number of reactions", fill = "Reaction type") +
  theme(
    text = element_text(family = "Helvetica"),
    strip.background = element_blank(),
    strip.text = element_blank(),
    legend.margin = margin()
  ) +
  scale_x_continuous(breaks = time_breaks) +
  scale_y_continuous(breaks = scales::breaks_extended(n = 3)) +
  geom_text(aes(label = reaction_text), label_fir_df, size = 3, hjust = 0.5, vjust = 1)



#######
# REACTION PROPENSITIES
#######
prop_df_A <- prop_df %>% filter(gene == "Gene A")
label_prop_df <- prop_df_A %>% group_by(reaction) %>% summarise(time = mean(time), value = max(value)*1.2) %>% mutate(reaction_text = paste0("Gene A ", reaction))

g4 <-
  ggplot(prop_df_A, aes(time, value)) +
  facet_wrap(~reaction, ncol = 1, scales = "free_y") +
  geom_area(aes(fill = forcats::fct_rev(reaction)), position = "stack") +
  theme_classic() +
  scale_fill_brewer(palette = "Set2") +
  labs(x = "Simulation time", y = "Propensity", fill = "Reaction type") +
  theme(
    text = element_text(family = "Helvetica"),
    strip.background = element_blank(),
    strip.text = element_blank(),
    legend.margin = margin()
  ) +
  scale_x_continuous(breaks = time_breaks) +
  scale_y_continuous(breaks = scales::breaks_extended(n = 3)) +
  geom_text(aes(label = reaction_text), label_prop_df, size = 3, hjust = 0.5, vjust = 1)


#######
# REGULATION
#######
g5 <- ggplot(reg_df) +
  geom_line(aes(time, abs(value), colour = name), size = 1) +
  theme_classic() +
  scale_colour_manual(values = RColorBrewer::brewer.pal(6, "Set3")[-2]) +
  labs(x = "Simulation time", y = "Regulation", colour = "Interaction") +
  theme(
    text = element_text(family = "Helvetica"),
    legend.margin = margin()
  ) +
  scale_x_continuous(breaks = time_breaks) +
  scale_y_continuous(breaks = c(0, .5, 1), limits = c(0, 1))


#######
# Prop diff
#######
g6 <- ggplot(velocity, aes(time, value)) +
  geom_step(aes(colour = gene)) +
  # facet_wrap(~gene, ncol = 1, scales = "free_y") +
  theme_classic() +
  labs(x = "Simulation time", y = "Velocity ground-truth", colour = "Molecule") +
  scale_colour_brewer(palette = "Set1") +
  # scale_colour_manual(values = c("pre-mRNA" = "#4daf4a", "mRNA" = "#377eb8", "protein" = "#e41a1c", "gene" = "black")) +
  theme(
    text = element_text(family = "Helvetica"),
    strip.background = element_blank(),
    strip.text = element_blank(),
    legend.margin = margin()
  ) +
  scale_x_continuous(breaks = time_breaks)



g <- patchwork::wrap_plots(
  g1,
  g2,
  g3,
  g5,
  heights = c(3, 1, 3, 1),
  ncol = 1
)

ggsave(exp$result("overview_unedited.pdf"), g, width = 8, height = 8, device = cairo_pdf)




g <- patchwork::wrap_plots(
  g1 + labs(title = "Expression"),
  g2 + labs(title = "State"),
  g3 + labs(title = "Reaction firings"),
  g4 + labs(title = "Propensity"),
  g5 + labs(title = "Regulation"),
  g6 + labs(title = "Velocity ground-truth"),
  heights = c(3, 1, 3),
  ncol = 2,
  byrow = FALSE
)

ggsave(exp$result("overview_alldata.pdf"), g, width = 16, height = 6, device = cairo_pdf)


