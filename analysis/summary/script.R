library(tidyverse)
library(dyngen)
library(dyngen.manuscript)
library(grImport)
library(patchwork)

set.seed(1)

exp <- start_analysis("summary")

# Read regulation svg from file
fil <- tempfile()
on.exit(file.remove(fil))
PostScriptTrace(exp$result("gene_regulation.eps"), outfilename = fil)
g0a <- pictureGrob(readPicture(fil))

# generate cyclic trajectory
census_interval <- 5
total_time <- 500
time_breaks <- c(0, 500)
model <- exp$temporary("model.rds") %cache% {
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
      simulation_params = simulation_default(
        burn_time = 240,
        total_time = total_time,
        experiment_params = simulation_type_wild_type(num_simulations = 1),
        census_interval = census_interval,
        store_reaction_firings = TRUE,
        store_reaction_propensities = TRUE,
        compute_cellwise_grn = TRUE,
        compute_rna_velocity = TRUE
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

  # floor to nearest census
  model$simulations$meta$sim_time <- floor(model$simulations$meta$sim_time / census_interval) * census_interval


  model
}

mol_map <- c(mol_premrna = "pre-mRNA", mol_mrna = "mRNA", mol_protein = "protein")
reac_map <- c(
  "transcription" = "Transcription", "splicing" = "Splicing", "translation" = "Translation",
  "premrna_degradation" = "pre-mRNA degradation", "mrna_degradation" = "mRNA degradation", "protein_degradation" = "Protein degradation"
)

time <- model$simulations$meta$sim_time
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
velocity <- data.frame(time = time, model$simulations$rna_velocity %>% as.matrix()) %>%
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
# BACKBONE PLOT
#######
library(ggraph)
set.seed(2)

nodes <- model$backbone$module_info %>% rename(name = module_id)
edges <- model$backbone$module_network %>% arrange(from == to)

gr <- tidygraph::tbl_graph(nodes = nodes, edges = edges)
layout <-
  gr %>%
  igraph::layout.graphopt(charge = .01, niter = 10000) %>%
  dynutils::scale_minmax() %>%
  magrittr::set_rownames(nodes$name) %>%
  magrittr::set_colnames(c("x", "y")) %>%
  as.data.frame()

r <- .1
cap <- ggraph::circle(4, "mm")
str <- .2
arrow_up <- grid::arrow(type = "closed", angle = 30, length = grid::unit(3, "mm"))
arrow_down <- grid::arrow(type = "closed", angle = 89, length = grid::unit(3, "mm"))

g0b <-
  ggraph(gr, layout = "manual", x = layout$x, y = layout$y) +
  geom_edge_fan(aes(filter = effect >= 0), arrow = arrow_up, start_cap = cap, end_cap = cap) +
  geom_edge_fan(aes(filter = effect < 0), arrow = arrow_down, start_cap = cap, end_cap = cap) +
  geom_node_point(aes(colour = name), size = 7) +
  geom_node_point(colour = "white", size = 6) +
  geom_node_text(aes(label = name)) +
  scale_edge_width_continuous(trans = "log10", range = c(.5, 3)) +
  scale_colour_brewer(palette = "Set1") +
  expand_limits(y = c(-.2, 1.2), x = c(-.05, 1.05)) +
  coord_cartesian() +
  theme_classic() +
  theme(
    legend.position = "none",
    text = element_text(family = "Helvetica"),
    axis.title = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank()
  )

#######
# ABUNDANCE LEVELS
#######
g1 <- ggplot(expr_df, aes(time, value)) +
  geom_step(aes(colour = gene)) +
  facet_wrap(~molecule, ncol = 1, scales = "free_y") +
  theme_classic() +
  labs(x = "Simulation time", y = "Molecule levels", colour = "Molecule") +
  scale_colour_brewer(palette = "Set1") +
  scale_y_continuous(n.breaks = 3, expand = c(0, 0)) +
  theme(
    text = element_text(family = "Helvetica"),
    strip.background = element_blank(),
    strip.text = element_blank(),
    legend.margin = margin(),
    legend.position = "none"
  ) +
  scale_x_continuous(breaks = time_breaks, expand = c(0.02, 0.02)) +
  geom_text(aes(x = mean(time_breaks), y = max, label = molecule), expr_df %>% group_by(molecule) %>% summarise(max = max(value)), size = 3, hjust = 0.5, vjust = 1)

#######
# CELL STATE
#######
labels_df <- state_df %>%
  group_by(sim_time) %>%
  arrange(milestone_id) %>%
  mutate(
    y2 = cumsum(percentage),
    y1 = lag(y2, 1, 0)
  ) %>%
  group_by(milestone_id) %>%
  slice(
    first(which(y2 - y1 > 0.8 & sim_time > 30))
  )

g2 <- ggplot(state_df) +
  geom_area(aes(sim_time, percentage, fill = milestone_id)) +
  geom_text(aes(sim_time, 1 - (y1 + (y2 - y1)/2), label = milestone_id), data = labels_df) +
  theme_classic() +
  scale_fill_brewer(palette = "Accent") +
  labs(x = "Simulation time", y = "Cell state", colour = "Cell state") +
  theme(
    text = element_text(family = "Helvetica"),
    legend.margin = margin(),
    legend.position = "none"
  ) +
  scale_x_continuous(breaks = time_breaks, expand = c(0.02, 0.02)) +
  scale_y_continuous(breaks = c(0, 1), labels = scales::percent, expand = c(0, 0))
g2

#######
# REACTION PROPENSITIES
#######
prop_df_A <- prop_df %>% filter(gene == "Gene A")
label_prop_df <- prop_df_A %>% group_by(reaction) %>% summarise(time = mean(time), value = max(value)*1.2) %>% mutate(reaction_text = paste0("Gene A ", reaction))

g3 <-
  ggplot(prop_df_A, aes(time, value)) +
  facet_wrap(~reaction, ncol = 1, scales = "free_y") +
  geom_area(aes(fill = forcats::fct_rev(reaction)), position = "stack") +
  theme_classic() +
  scale_fill_brewer(palette = "Set2") +
  labs(x = "Simulation time", y = "Reaction propensity\n(Gene A)", fill = "Reaction type") +
  theme(
    text = element_text(family = "Helvetica"),
    strip.background = element_blank(),
    strip.text = element_blank(),
    legend.margin = margin(),
    legend.position = "none"
  ) +
  scale_x_continuous(breaks = time_breaks, expand = c(0.02, 0.02)) +
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
  scale_x_continuous(breaks = time_breaks, expand = c(0.02, 0.02)) +
  scale_y_continuous(breaks = c(0, .5, 1), limits = c(0, 1.2)) +
  theme(legend.position = c(.5, .95)) +
  guides(colour = guide_legend(direction = "horizontal"))


#######
# Applications
#######


exp7d <- start_analysis("showcase_backbones")
dataset <- read_rds(exp7d$dataset_file("bb_consecutive_bifurcating_3"))
dimred <- dyndimred::dimred_mds(dataset$expression, distance_method = "pearson")
g7d <- dynplot::plot_dimred(dataset, dimred = dimred, size_milestones = 3, size_cells = 1.5) +
  coord_cartesian() +
  theme_classic() +
  theme(
    legend.position = "none",
    text = element_text(family = "Helvetica"),
    axis.title = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank()
  ) +
  labs(title = "Trajectory inference")

exp7a <- start_analysis("usecase_trajectory_alignment")
li7a <- read_rds(exp7a$result("explanation_plot_data.rds"))

g7a <- ggplot() +
  geom_segment(data = li7a$segm_traj_test, mapping = aes(x = color2, y = comp_1, xend = color22, yend = comp_3), alpha = 0.25) +
  geom_point(data = li7a$comb_traj_less %>% filter(!is.na(color)), mapping = aes(x = color2, y = comp_1, colour = factor(color))) +
  scale_colour_manual(values = c("#fd8d3c", "#6baed6")) +
  scale_x_continuous(breaks = c(0, 1), labels = c("Start", "End")) +
  scale_y_continuous(breaks = c(-.6, .6)) +
  theme_classic() +
  theme(
    legend.position = "none",
    text = element_text(family = "Helvetica"),
    axis.title.y = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank()
  ) +
  labs(x = "Simulation time", title = "Trajectory alignment") +
  coord_cartesian()
g7a

exp7b <- start_analysis("usecase_rna_velocity")
g7b <- read_rds(exp7b$result("one_rna_velocity.rds")) +
  # plot_spacer() +
  coord_cartesian() +
  theme_classic() +
  theme(
    legend.position = "none",
    text = element_text(family = "Helvetica"),
    axis.title = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank()
  ) +
  labs(title = "RNA velocity", subtitle = NULL)

exp7c <- start_analysis("usecase_network_inference")
g7c <- read_rds(exp7c$result("cell1.rds")) +
  coord_cartesian() +
  theme_classic() +
  theme(
    legend.position = "none",
    text = element_text(family = "Helvetica"),
    axis.title = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank()
  ) +
  labs(title = "Cell-specific network inference")

such_great_heights <- c(2.5, 3, 2, 4, 3)
g <- wrap_plots(
  wrap_plots(g0a, g0b, nrow = 1, widths = c(1, 1)),
  g1,
  g2,
  g3,
  g5,
  plot_spacer(),
  g7d,
  g7a,
  g7b,
  g7c,
  # wrap_plots(g0a, g0b, nrow = 1, widths = c(1, 1)), plot_spacer(),
  # g1, g7d,
  # g2, g7a,
  # g3, g7b,
  # g5, g7c,
  heights = such_great_heights,
  ncol = 2,
  widths = c(7, 3),
  byrow = FALSE
) +
  plot_annotation(tag_levels = c('A')) &
  theme(plot.tag.position = c(0, 1), plot.title = element_text(hjust = .5))
ggsave(exp$result("overview.pdf"), g, width = 10, height = 10, device = cairo_pdf)



