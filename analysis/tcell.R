library(tidyverse)
library(dyngen)
library(dynplot)

set.seed(1)

out_dir <- "fig/usecase/tcell/"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Read network ------------------------------------------------------------
feature_info <-
  read_tsv("fig/usecase/tcell/feature_info.tsv", col_types = c(id = "c", name = "c")) %>%
  transmute(
    feature_id = paste0("Gene", id),
    module_id = NA_character_,
    ba = 0,
    burn = TRUE,
    ind = 1,
    is_tf = TRUE,
    is_hk = FALSE,
    feature_name = name
  )
feature_network <-
  read_tsv("fig/usecase/tcell/feature_network.tsv", col_types = c(from = "c", to = "c", effect = "c")) %>%
  group_by(from, to) %>% # may have duplicates..
  slice(1) %>%
  ungroup() %>%
  transmute(
    from = paste0("Gene", from),
    to = paste0("Gene", to),
    effect = as.integer(ifelse(effect == "?", "1", effect)),
    strength = ifelse(effect == 1, sample.int(4, n(), replace = TRUE), sample.int(2, n(), replace = TRUE)),
    cooperativity = 2
  )
feature_network$from[!feature_network$from %in% feature_info$feature_id]
feature_network$to[!feature_network$to %in% feature_info$feature_id]
feature_info$feature_id[!feature_info$feature_id %in% c(feature_network$from, feature_network$to)]
feature_info %>% filter(duplicated(feature_info))


gr <- igraph::graph_from_data_frame(
  feature_network,
  vertices = feature_info,
  directed = TRUE
)
# igraph::cluster_edge_betweenness(gr)
cl <- igraph::cluster_infomap(gr, modularity = FALSE)

module_ids <- paste0("Mod", LETTERS[sort(unique(cl$membership))])
feature_info$module_id <- module_ids[cl$membership]
# igraph::cluster_optimal(gr)
# igraph::cluster_louvain(gr)
# plot(gr, vertex.label = igraph::V(gr)$feature_name)


# Create model ------------------------------------------------------------

back <- backbone(
  module_info = tibble(module_id = module_ids, ba = 0, burn = TRUE, ind = 1),
  module_network =
    left_join(
      feature_network %>% rename(fromf = from, tof = to),
      feature_info %>% select(fromf = feature_id, from = module_id),
      by = "fromf"
    ) %>%
    left_join(
      feature_info %>% select(tof = feature_id, to = module_id),
      by = "tof"
    ) %>%
    select(-fromf, -tof) %>%
    group_by(from, to) %>%
    summarise_all(mean) %>%
    ungroup(),
  expression_patterns = tibble(from = "s0", to = "s1", module_progression = paste0("+", module_ids, collapse = ","), start = TRUE, burn = TRUE, time = 10)
)
model <-
  initialise_model(
    num_tfs = nrow(feature_info),
    num_targets = 0,
    num_hks = 0,
    num_cells = 1000,
    backbone = back,
    verbose = TRUE,
    download_cache_dir = "~/.cache/dyngen",
    num_cores = 8,
    gold_standard_params = gold_standard_default(),
    simulation_params = simulation_default(
      ssa_algorithm = GillespieSSA2::ssa_etl(tau = .001),
      num_simulations = 100,
      census_interval = .5,
      burn_time = 10,
      total_time = 100
    )
  ) %>%
  generate_tf_network()

model$feature_info <- feature_info
model$feature_network <- feature_network


model <- model %>%
  generate_kinetics()

plot_feature_network(model, show_targets = FALSE)

initial_state <- model$simulation_system$initial_state
model$simulation_system$initial_state <- matrix(
  sample.int(100, length(model$simulation_system$initial_state) * model$simulation_params$num_simulations, replace = TRUE),
  ncol = length(initial_state),
  dimnames = list(NULL, names(initial_state)),
)
model <- model %>%
  generate_cells()

model <- dyngen:::calculate_dimred(model, dimred_gold = FALSE, dimred_premrna = FALSE)
model$simulations$meta$from <- "s0"
model$simulations$meta$to <- "s1"
model$simulations$meta$time <- dynutils::scale_minmax(model$simulations$meta$sim_time)
plot_simulations(model)
plot_simulation_expression(model)

lala <- model$feature_info %>% filter(grepl("Other|Genes|DN", feature_name))
df <- data.frame(
  model$simulations$meta,
  as.matrix(model$simulations$counts[,lala$x])
) %>%
  gather(x, value, one_of(lala$x)) %>%
  left_join(lala %>% select(feature_name, x), by = "x") %>%
  as_tibble()

ggplot(df) + geom_path(aes(sim_time, value, colour = simulation_i, group = simulation_i)) +
  facet_wrap(~feature_name) +
  viridis::scale_color_viridis()
