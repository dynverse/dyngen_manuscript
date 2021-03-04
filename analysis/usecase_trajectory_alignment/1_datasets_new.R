library(dyngen.manuscript)
library(dyngen)
library(readr)
library(tidyverse)
library(rlang)

exp <- start_analysis("usecase_trajectory_alignment")

# get reference dataset
data("realcounts", package = "dyngen")
name_realcounts <- "zenodo_1443566_real_silver_bone-marrow-mesenchyme-erythrocyte-differentiation_mca"
url_realcounts <- realcounts %>% filter(name == name_realcounts) %>% pull(url)
realcount <- dyngen:::.download_cacheable_file(url_realcounts, getOption("dyngen_download_cache_dir"), verbose = FALSE)

backbone <- bblego(
  bblego_start("A", type = "simple", num_modules = 4),
  bblego_linear("A", "B", type = "simple", num_modules = 6),
  bblego_linear("B", "C", type = "simple", num_modules = 6),
  bblego_end("C", type = "simple", num_modules = 4)
)

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
      experiment_params = simulation_type_wild_type(
        num_simulations = 20
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
  generate_cells()
model1$id <- base1

# generate model2 with different kinetics
model2 <- model %>%
  generate_kinetics() %>%
  generate_gold_standard() %>%
  generate_cells()
model2$id <- base2



alpha <- seq(0, 1, by = .05)

# generate
datasets <- map(seq_along(alpha), function(i) {
  alpha_ <- alpha[[i]]

  # create new reference
  total <- sum(realcount)
  newref <- realcount
  newref@x <- rbinom(length(realcount@x), realcount@x, alpha) %>% as.numeric
  newref <- newref %>% Matrix::drop0()

  model1i <- model1
  model2i <- model2

  model1i$experiment_params$realcount <- newref
  model2i$experiment_params$realcount <- newref

  model1i <- model1i %>%
    generate_experiment()
  model2i <- model2i %>%
    generate_experiment()

  list(
    alpha = alpha_,
    left = as_dyno(model1i),
    right = as_dyno(model2i)
  )
  # comb <- combine_models(list(model1 = model1i, model2 = model2i))
  # plot_gold_simulations(comb)
  #
  # plot_experiment_dimred(model1i)
  # plot_experiment_dimred(model2i)
  # plot_experiment_dimred(comb)
  #
  # plot_simulations(comb)
  #
  # modeli
})


write_rds(datasets, exp$dataset_file("test"), compress = "gz")

gc()
