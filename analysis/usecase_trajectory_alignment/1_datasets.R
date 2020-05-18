library(dyngen.manuscript)
library(dyngen)
library(readr)
library(tidyverse)
library(rlang)

exp <- start_analysis("usecase_trajectory_alignment")

backbone <- bblego(
  bblego_start("A", type = "simple", num_modules = 4),
  bblego_linear("A", "B", type = "simple", num_modules = 6),
  bblego_linear("B", "C", type = "simple", num_modules = 6),
  bblego_end("C", type = "simple", num_modules = 4)
)

noise_levels <- seq(from = 0.1, to = 1, by = 0.1)
only_noise <- TRUE

# setup dataset design
design_datasets <- exp$result("design_datasets.rds") %cache% {
  crossing(
    seed = 1:10,
    backbone_name = "linear",
    noise = noise_levels
  ) %>%
    mutate(
      base_id1 = paste0(backbone_name, seed, "_1"),
      base_id2 = paste0(backbone_name, seed, "_2"),
      id1 = paste0(base_id1, "_", noise),
      id2 = paste0(base_id2, "_", noise)
    )
}

# Generates only the base models. The noise is added onto these base models
pwalk(design_datasets, function(base_id1, base_id2, id1, id2, seed, backbone_name, noise) {

  if (!file.exists(exp$dataset_file(base_id1))) {

    cat("## Generating ", base_id1, "\n", sep = "")
    set.seed(seed)
    model <-
      initialise_model(
        id = base_id1,
        num_tfs = 50,
        num_targets = 70,
        num_hks = 30,
        backbone = backbone,
        num_cells = 1000,
        simulation_params = simulation_default(
          census_interval = 10,
          experiment_params = simulation_type_wild_type(
            num_simulations = 100
          )
        ),
        num_cores = 6,
        download_cache_dir = "~/.cache/dyngen",
        verbose = TRUE
      )
    generate_dataset(
      model,
      output_dir = exp$dataset_folder(base_id1),
      make_plots = TRUE
    )

    dataset1 <- read_rds(paste0(exp$dataset_folder(base_id1), "dataset.rds"))
    model1 <- read_rds(paste0(exp$dataset_folder(base_id1), "model.rds"))

    # Generate the second trajectory of the pair
    model2 <- model1 %>% generate_cells() %>% generate_experiment()
    dataset2 <- wrap_dataset(model2)

    write_rds(model2, paste0(exp$dataset_folder(base_id2), "model.rds"), compress = "gz")
    write_rds(dataset2, paste0(exp$dataset_folder(base_id2), "dataset.rds"), compress = "gz")

    max_dev1 <- dataset1$counts %>% as.vector() %>% Filter(f = function(x) x > 0) %>% quantile(0.75)
    max_dev2 <- dataset2$counts %>% as.vector() %>% Filter(f = function(x) x > 0) %>% quantile(0.75)

    # Add noise gradually to the datasets, according to the noiselevels
    walk(noise_levels, function(noise_perc, ds1=dataset1, ds2=dataset2){
      name1 <- paste0(base_id1, "_", noise_perc)
      name2 <- paste0(base_id2, "_", noise_perc)

      ds1$counts <- ds1$counts + rnorm(length(ds1$counts), mean = 0, sd = max_dev1 * noise_perc)
      ds2$counts <- ds2$counts + rnorm(length(ds2$counts), mean = 0, sd = max_dev2 * noise_perc)

      ds1$counts[ds1$counts<0] <- 0
      ds2$counts[ds2$counts<0] <- 0

      write_rds(ds1, paste0(exp$dataset_folder(name1), "dataset.rds"), compress = "gz")
      write_rds(ds2, paste0(exp$dataset_folder(name2), "dataset.rds"), compress = "gz")

    })

    gc()
  }
})

