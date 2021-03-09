library(dyngen.manuscript)
library(dyngen)
library(tidyverse)

exp <- start_analysis("usecase_trajectory_alignment")

design_datasets <- exp$result("design_datasets.rds") %cache% {
  crossing(
    seed = 1:10,
    alpha = 10^seq(log10(0.001), log10(0.1), length.out = 10)
  ) %>%
    mutate(
      group = paste0("linear_seed", seed),
      id = sprintf("%s_alpha%.3f", group, alpha),
    )
}

#' design_datasets %>% dynutils::extract_row_to_list(1) %>% list2env(.GlobalEnv)
pwalk(design_datasets, function(seed, group, id, alpha) {
  if (!file.exists(exp$dataset_file(id))) {
    cat("## Generating ", id, "\n", sep = "")
    # generating common settings
    set.seed(seed)

    # get reference dataset
    data("realcounts", package = "dyngen")
    name_realcounts <- dyngen:::.sample_decent_realcount()
    url_realcounts <- realcounts %>% filter(name == name_realcounts) %>% pull(url)
    realcount <- dyngen:::.download_cacheable_file(url_realcounts, getOption("dyngen_download_cache_dir"), verbose = FALSE)

    # create backbone
    backbone <- bblego(
      bblego_start("A", type = "simple", num_modules = 4),
      bblego_linear("A", "B", type = "simple", num_modules = 6),
      bblego_linear("B", "C", type = "simple", num_modules = 6),
      bblego_end("C", type = "simple", num_modules = 4)
    )

    # generate model
    total_time <- 300

    model <-
      initialise_model(
        backbone = backbone,
        num_tfs = 50,
        num_targets = 300,
        num_hks = 200,
        num_cells = 1000,
        simulation_params = simulation_default(
          census_interval = 10,
          total_time = total_time,
          kinetics_noise_function = kinetics_noise_simple(sd = alpha),
          experiment_params = simulation_type_wild_type(
            num_simulations = 40
          )
        ),
        verbose = TRUE
      ) %>%
      generate_tf_network() %>%
      generate_feature_network()

    # generate kinetics for two bases
    model1 <- model %>%
      generate_kinetics() %>%
      generate_gold_standard() %>%
      generate_cells() %>%
      generate_experiment()

    # generate model2 with different kinetics
    model2 <- model %>%
      generate_kinetics() %>%
      generate_gold_standard() %>%
      generate_cells() %>%
      generate_experiment()

    # visual checks
    comb <- combine_models(list(model1 = model1, model2 = model2))
    plot_gold_simulations(comb)
    plot_experiment_dimred(comb)
    plot_simulations(comb)

    out <- list(
      left = as_dyno(model1),
      right = as_dyno(model2)
    )

    write_rds(out, exp$dataset_file(id), compress = "gz")
  }
})



