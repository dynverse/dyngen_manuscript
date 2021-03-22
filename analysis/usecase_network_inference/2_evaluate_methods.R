library(tidyverse)
library(dynwrap)
library(dyneval)
library(dynutils)
library(dyngen.manuscript)

exp <- start_analysis("usecase_network_inference")

## might need to install an older version of miniconda because the latest version doesn't seem to want to work
# path = reticulate::miniconda_path()
# url <- "https://repo.anaconda.com/miniconda/Miniconda3-py38_4.8.3-Linux-x86_64.sh"
# installer <- reticulate:::miniconda_installer_download(url)
# reticulate:::miniconda_installer_run(installer, path)
# reticulate::py_install("pyscenic", pip = TRUE)
# reticulate::py_install("tbb", pip = TRUE)

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

aucs <-
  summaries %>%
  select(cni_method_id, cni_method_name, dataset_id, evals) %>%
  unnest(evals) %>%
  mutate(
    method_label = c("static_static" = "Regular Network Inference", "casewise_casewise" = "Casewise Network Inference")[method]
  ) %>%
  filter(method != "static_casewise")

summaries %>% select(cni_method_name, dataset_id, starts_with("cc_"), starts_with("sc_"))


summ <- aucs %>%
  group_by(method, method_label, dataset_id, cni_method_id, cni_method_name) %>%
  summarise_if(is.numeric, mean) %>%
  ungroup()
write_rds(list(aucs = aucs, summ = summ), exp$result("scores.rds"), compress = "gz")

ggstatsplot::grouped_ggwithinstats(
  summ %>% filter(method == "casewise_casewise") %>% gather(metric, score, auroc, aupr, F1) %>% mutate(metric = factor(metric, levels = c("auroc", "aupr", "F1"))),
  cni_method_name,
  score,
  grouping.var = metric,
  type = "np"
)


nam <- unique(summ$method)
summplots <- map(nam, function(meth) {
  s <- summ %>% filter(method == meth)
  meth_lab <- s$method_label[[1]]
  ggplot(s) +
    geom_point(aes(auroc, aupr), colour = "gray", function(df) df %>% select(-cni_method_id, -cni_method_name)) +
    geom_point(aes(auroc, aupr, colour = cni_method_name)) +
    theme_bw() +
    facet_wrap(~cni_method_name, nrow = 1) +
    labs(title = meth_lab, x = "mean AUROC", y = "mean AUPR", colour = "Method") +
    scale_color_brewer(palette = "Set1")
})
names(summplots) <- nam
ggsave(exp$result("score_summary.pdf"), patchwork::wrap_plots(summplots, ncol = 1), width = 10, height = 5)

write_rds(summplots, exp$result("score_summary.rds"), compress = "gz")

summplots$casewise_casewise
