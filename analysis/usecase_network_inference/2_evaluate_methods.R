library(tidyverse)
library(dynwrap)
library(dyneval)
library(dynutils)
library(dyngen.manuscript)

exp <- start_analysis("usecase_network_inference")

dataset_files <- list.files(exp$dataset_folder(""), pattern = "dataset.rds", full.names = TRUE, recursive = TRUE)
datasets <- list_as_tibble(map(dataset_files, function(dataset_file) {
  dataset <- read_rds(dataset_file) %>% add_cell_waypoints()
  dataset$prior_information$regulators <- dataset$regulators
  dataset$prior_information$targets <- dataset$targets
  dataset
}))

methods <- list(
  "SSN*" = cni_ssn(),
  "pySCENIC GBM" = cni_pyscenic_sgbm(subsample = 1, n_estimators = 500L),
  "pySCENIC SGBM" = cni_pyscenic_sgbm(subsample = .9, n_estimators = 5000L),
  "LIONESS" = cni_lioness()
)

out <- exp$result("evaluation.rds") %cache% function() {
  evals <- map(methods, function(method) {
    evaluate_ti_method(
      dataset = datasets,
      method = method,
      metrics = list(auc = cni_auc),
      parameters = NULL
    )
  })

  aucs <- map_df(seq_along(evals), function(i) {
    summ <- evals[[i]]$summary
    mapdf_dfr(summ, function(li) {
      li$evals %>% mutate(
        cni_method_id = names(methods)[[i]],
        dataset_id = li$dataset_id
      )
    })
  })
  int_ns <- mapdf_dfr(datasets, function(dat) dat$regulatory_network_sc %>% group_by(cell_id) %>% summarise(nints = n()) %>% mutate(dataset_id = dat$id))
  aucs <- aucs %>% left_join(int_ns, c("cell_id", "dataset_id"))

  lst(evals, aucs)
}

aucs <-
  out$aucs %>%
  filter(
    !cni_method_id %in% c("bred", "pySCENIC GBM")
  ) %>%
  mutate(
    method_label = c("static_static" = "Regular Network Inference", "casewise_casewise" = "Casewise Network Inference")[method],
    cni_method_id = ifelse(cni_method_id == "pySCENIC SGBM", "pySCENIC", cni_method_id)
  )

summ <- aucs %>%
  group_by(method, method_label, dataset_id, cni_method_id) %>%
  summarise_if(is.numeric, mean) %>%
  ungroup()

g <- ggplot(aucs %>% filter(method == "casewise_casewise")) +
  geom_point(aes(auroc, aupr), colour = "lightgray", function(df) select(df, -cni_method_id), size = 1) +
  geom_point(aes(auroc, aupr, colour = cni_method_id), size = 1) +
  facet_grid(dataset_id~cni_method_id) +
  theme_bw()
ggsave(exp$result("comparison_casewise_casewise.pdf"), g, width = 6, height = 15)

# g <- ggplot(aucs %>% filter(method == "casewise_casewise") %>% sample_n(n())) +
#   geom_point(aes(auroc, aupr), colour = "lightgray", function(df) select(df, -cni_method_id), size = 1) +
#   geom_point(aes(auroc, aupr, colour = nints), size = 1) +
#   facet_grid(dataset_id~cni_method_id) +
#   viridis::scale_color_viridis() +
#   theme_bw()
# ggsave(exp$result("comparison_casewise_casewise_nints.pdf"), g, width = 6, height = 15)

g <- ggplot(aucs %>% filter(method == "static_static")) +
  geom_point(aes(auroc, aupr), colour = "lightgray", function(df) select(df, -cni_method_id), size = 1) +
  geom_point(aes(auroc, aupr, colour = cni_method_id), size = 1) +
  facet_grid(~cni_method_id) +
  theme_bw()
ggsave(exp$result("comparison_static_static.pdf"), g, width = 10, height = 3)

ggplot(summ %>% filter(method == "casewise_casewise")) +
  geom_point(aes(auroc, aupr), colour = "gray", function(df) df %>% select(-cni_method_id)) +
  geom_point(aes(auroc, aupr, colour = cni_method_id)) +
  theme_bw() +
  facet_wrap(~cni_method_id)


nam <- unique(summ$method)
summplots <- map(nam, function(meth) {
  s <- summ %>% filter(method == meth)
  meth_lab <- s$method_label[[1]]
  ggplot(s) +
    geom_point(aes(auroc, aupr), colour = "gray", function(df) df %>% select(-cni_method_id)) +
    geom_point(aes(auroc, aupr, colour = cni_method_id)) +
    theme_bw() +
    facet_wrap(~cni_method_id, nrow = 1) +
    labs(title = meth_lab, x = "mean AUROC", y = "mean AUPR", colour = "Method") +
    scale_color_brewer(palette = "Set1")
})
names(summplots) <- nam
ggsave(exp$result("score_summary.pdf"), patchwork::wrap_plots(summplots, ncol = 1), width = 8, height = 5)

write_rds(summplots, exp$result("score_summary.rds"), compress = "gz")

