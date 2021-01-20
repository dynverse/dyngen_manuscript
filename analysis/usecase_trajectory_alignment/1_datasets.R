library(dyngen.manuscript)
library(dyngen)
library(readr)
library(tidyverse)
library(rlang)

exp <- start_analysis("usecase_trajectory_alignment")

selected_backbones <- list(
  list(
    name = "custom_linear",
    total_time = 300, # this backbone reaches endstate in about 300 time units
    backbone = bblego(
      bblego_start("A", type = "simple", num_modules = 4),
      bblego_linear("A", "B", type = "simple", num_modules = 6),
      bblego_linear("B", "C", type = "simple", num_modules = 6),
      bblego_end("C", type = "simple", num_modules = 4)
    )
  )
)
names(selected_backbones) <- map_chr(selected_backbones, "name")

alphas <- seq(0, 1, by = 0.1)

# setup dataset design
design_datasets <- exp$result("design_datasets.rds") %cache% {
  crossing(
    seed = 1:10,
    backbone_name = names(selected_backbones),
    alpha = alphas
  ) %>%
    mutate(
      group = paste0(backbone_name, "_seed", seed),
      base1 = paste0(group, "_base1"),
      base2 = paste0(group, "_base2"),
      id = sprintf("%s_alpha%.1f", group, alpha),
    )
}


# generate all datasets per base id
design_grouped <-
  design_datasets %>%
  group_by(group) %>%
  chop(c(id, alpha)) %>%
  ungroup()

# design_grouped %>% dynutils::extract_row_to_list(1) %>% list2env(.GlobalEnv)

pwalk(design_grouped, function(seed, backbone_name, total_time, group, base1, base2, id, alpha) {

  if (!file.exists(exp$dataset_file(group))) {

    # generating common settings
    set.seed(seed)

    backbone <- selected_backbones[[backbone_name]]$backbone
    total_time <- selected_backbones[[backbone_name]]$total_time

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
          experiment_params = simulation_type_wild_type(
            num_simulations = 100
          )
        ),
        verbose = TRUE
      ) %>%
      generate_tf_network() %>%
      generate_feature_network()

    # generate kinetics for two bases
    cat("## Generating ", base1, "\n", sep = "")
    model1 <- model %>%
      generate_kinetics() %>%
      generate_gold_standard() %>%
      generate_cells() %>%
      generate_experiment()
    model1$id <- base1

    # plot_simulations(model1)

    # generate model2 with different kinetics
    model2 <- model %>%
      generate_kinetics()
    model2$id <- base2

    # generate
    modelis <- map(seq_along(alpha), function(i) {
      alpha_ <- alpha[[i]]
      id_ <- id[[i]]

      cat("## Generating ", id_, "\n", sep = "")

      modeli <- model2
      modeli$id <- id_

      fii <- model2$feature_info
      fni <- model2$feature_network
      for (col in colnames(fii)) {
        if (is.numeric(model1$feature_info[[col]]) && !all(model1$feature_info[[col]] == model2$feature_info[[col]])) {
          fii[[col]] <- model1$feature_info[[col]] * (1-alpha_) + model2$feature_info[[col]] * alpha_
        }
      }
      for (col in colnames(fni)) {
        if (is.numeric(model1$feature_network[[col]]) && !all(model1$feature_network[[col]] == model2$feature_network[[col]])) {
          fni[[col]] <- model1$feature_network[[col]] * (1-alpha_) + model2$feature_network[[col]] * alpha_
        }
      }
      modeli$feature_info <- fii
      modeli$feature_network <- fni

      modeli <- modeli %>%
        generate_gold_standard() %>%
        generate_cells() %>%
        generate_experiment()

      # plot_gold_simulations(combine_models(list(model1 = model1, modeli = modeli)))
      # plot_simulations(combine_models(list(model1 = model1, modeli = modeli)))

      modeli
    })

    models <- lst(
      model1,
      model2,
      modelis,
    )
    write_rds(models, exp$dataset_file(group) %>% gsub("dataset.rds$", "model.rds", .), compress = "gz")

    datasets <- lst(
      dataset1 = as_dyno(model1),
      dataset2 = as_dyno(model2),
      dataseti = map(modelis, as_dyno)
    )
    write_rds(datasets, exp$dataset_file(group), compress = "gz")

    gc()
  }
})

