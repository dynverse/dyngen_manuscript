library(tidyverse)
library(dynwrap)
library(dyneval)
library(dynutils)
library(dyngen.manuscript)

exp <- start_analysis("usecase_network_inference")

## Load all datasets as a tibble
dataset_files <- list.files(exp$dataset_folder(""), pattern = "dataset.rds", full.names = TRUE, recursive = TRUE)

datasets <- list_as_tibble(map(dataset_files, function(dataset_file) {
  dataset <- read_rds(dataset_file) %>% add_cell_waypoints()
  dataset$prior_information$regulators <- dataset$regulators
  dataset$prior_information$targets <- dataset$targets
  dataset
}))

## Load methods as a tibble
methods <- tribble(
  ~name, ~fun,
  "SSN*", cni_ssn(),
  # "pySCENIC GBM", cni_pyscenic_sgbm(subsample = 1, n_estimators = 500L),
  "pySCENIC", cni_pyscenic_sgbm(subsample = .9, n_estimators = 5000L),
  "LIONESS + Pearson", cni_lioness(method = "pearson")
) %>% mutate(
  id = name %>% tolower() %>% gsub("[^a-z]", "", .)
)

#' @examples
#' dataset <- datasets %>% extract_row_to_list(1)
#' priors <- dataset$prior_information
#' expression <- dataset$expression
#' parameters <- dynwrap::get_default_parameters(methods$fun[[2]])
#' parameters <- dynwrap::get_default_parameters(methods$fun[[3]])

## Run methods on all datasets

#' @examples
#' methods %>% mutate(i = row_number()) %>% dynutils::extract_row_to_list(2) %>% list2env(.GlobalEnv)
pwalk(
  methods %>% mutate(i = row_number()),
  function(name, fun, id, i) {
    cat("Evaluating ", i, "/", nrow(methods), ": ", name, "\n", sep = "")

    exp$temporary("eval/evaluation_", id, ".rds") %cache% {
      out <- evaluate_ti_method(
        dataset = datasets,
        method = fun,
        output_model = FALSE,
        metrics = list(cni_auc = cni_auc),
        parameters = NULL
      )
      out$summary %>% mutate(
        cni_method_id = id,
        cni_method_name = name
      )
    }
  }
)

## Summarise results
summaries <-
  map_df(list.files(exp$temporary("eval"), pattern = ".*\\.rds", full.names = TRUE), read_rds) %>%
  mutate(
    cni_method_name = factor(cni_method_name, levels = c("SSN*", "LIONESS + Pearson", "pySCENIC"))
  )

# check percentage errored
summaries %>%
  mutate(pct_errored = (!is.na(error)) + 0) %>%
  group_by(cni_method_name) %>%
  summarise_if(is.numeric, mean, na.rm = TRUE)

# compute metrics
aucs <-
  summaries %>%
  select(cni_method_id, cni_method_name, dataset_id, evals) %>%
  unnest(evals) %>%
  mutate(
    method_label = c("static_static" = "Regular Network Inference", "casewise_casewise" = "Casewise Network Inference")[method]
  ) %>%
  filter(method != "static_casewise")

# compute summaries
summaries %>% select(cni_method_name, dataset_id, starts_with("cc_"), starts_with("sc_"))

summ <- aucs %>%
  group_by(method, method_label, dataset_id, cni_method_id, cni_method_name) %>%
  summarise_if(is.numeric, mean) %>%
  ungroup()

# write results to file
write_rds(list(aucs = aucs, summ = summ), exp$result("scores.rds"), compress = "gz")

# preview results
ggstatsplot::grouped_ggwithinstats(
  summ %>% filter(method == "casewise_casewise") %>% gather(metric, score, auroc, aupr, F1) %>% mutate(metric = factor(metric, levels = c("auroc", "aupr", "F1"))),
  cni_method_name,
  score,
  grouping.var = metric,
  type = "np"
)
