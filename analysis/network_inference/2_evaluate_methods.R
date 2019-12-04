library(tidyverse)
library(dynwrap)
library(dyneval)
library(dyngenanalysis)
library(dynutils)

dataset_files <- list.files("derived_files/network_inference/datasets/", pattern = "dataset.rds", full.names = TRUE, recursive = TRUE)
datasets <- list_as_tibble(map(dataset_files, function(dataset_file) {
  dataset <- read_rds(dataset_file) %>% add_cell_waypoints()
  dataset$prior_information$regulators <- dataset$regulators
  dataset$prior_information$targets <- dataset$targets
  dataset
}))

# expression <- dataset$expression
# regulators <- dataset$regulators
# targets <- dataset$targets
# parameters <- list(num_int_per_cell = 10000L)

methods <- list(
  "bred" = cni_bred(num_trees = 1000L),
  "deltacor" = cni_deltacor(),
  "pySCENIC GBM" = cni_pyscenic_sgbm(subsample = 1, n_estimators = 500L),
  "pySCENIC SGBM" = cni_pyscenic_sgbm(subsample = .9, n_estimators = 5000L)
)

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


ggplot(aucs) + geom_point(aes(auroc, aupr, colour = method)) + facet_wrap(~cni_method_id, nrow = 1)
ggplot(aucs) + geom_point(aes(auroc, aupr, colour = method)) + facet_grid(cni_method_id~dataset_id)
ggplot(aucs) + geom_point(aes(auroc, aupr, colour = cni_method_id)) + facet_wrap(~method, scales = "free")

g <- ggplot(aucs %>% filter(method == "casewise_casewise")) +
  geom_point(aes(auroc, aupr), colour = "lightgray", function(df) select(df, -cni_method_id), size = 1) +
  geom_point(aes(auroc, aupr, colour = cni_method_id), size = 1) +
  facet_grid(dataset_id~cni_method_id) +
  theme_bw()
ggsave(paste0("fig/network_inference/comparison_casewise_casewise.pdf"), g, width = 6, height = 15)


g <- ggplot(aucs %>% filter(method == "casewise_casewise") %>% sample_n(n())) +
  geom_point(aes(auroc, aupr), colour = "lightgray", function(df) select(df, -cni_method_id), size = 1) +
  geom_point(aes(auroc, aupr, colour = nints), size = 1) +
  facet_grid(dataset_id~cni_method_id) +
  viridis::scale_color_viridis() +
  theme_bw()
ggsave(paste0("fig/network_inference/comparison_casewise_casewise_nints.pdf"), g, width = 6, height = 15)

g <- ggplot(aucs %>% filter(method == "static_static")) +
  geom_point(aes(auroc, aupr), colour = "lightgray", function(df) select(df, -cni_method_id), size = 1) +
  geom_point(aes(auroc, aupr, colour = cni_method_id), size = 1) +
  facet_grid(~cni_method_id) +
  theme_bw()
ggsave(paste0("fig/network_inference/comparison_static_static.pdf"), g, width = 10, height = 3)




g <- ggplot(aucs %>% filter(method == "casewise_casewise")) +
  geom_density(aes(F1, colour = cni_method_id)) +
  facet_wrap(~dataset_id) +
  theme_bw()
g

summ <- aucs %>% filter(method == "casewise_casewise") %>% group_by(dataset_id, cni_method_id) %>% summarise_if(is.numeric, mean) %>% ungroup()
ggplot(summ) +
  geom_point(aes(auroc, aupr), colour = "gray", function(df) df %>% select(-cni_method_id)) +
  geom_point(aes(auroc, aupr, colour = cni_method_id)) +
  theme_bw() +
  facet_wrap(~cni_method_id)


summ <- aucs %>% group_by(dataset_id, cni_method_id, method) %>% summarise_if(is.numeric, mean) %>% ungroup()
g <- map(c("static_static", "casewise_casewise"), function(meth) {
  ggplot(summ %>% filter(method == meth)) +
    geom_point(aes(auroc, aupr), colour = "gray", function(df) df %>% select(-cni_method_id)) +
    geom_point(aes(auroc, aupr, colour = cni_method_id)) +
    theme_bw() +
    facet_wrap(~cni_method_id, nrow = 1) +
    labs(title = meth, x = "AUROC", y = "AUPR", colour = ifelse(meth == "static_static", "NI method", "CNI method")) +
    scale_color_brewer(palette = "Set1")
}) %>% patchwork::wrap_plots(ncol = 1)
ggsave("fig/network_inference/score_summary.pdf", g, width = 8, height = 5)



