library(tidyverse)
library(dyngen)
library(rlang)

derived_folder <- "derived_files/network_inference/"

dataset_files <- list.files(paste0(derived_folder, "datasets/"), pattern = "dataset.rds", recursive = TRUE, full.names = TRUE)

methods <- tribble(
  ~id, ~regressor_type, ~regressor_kwargs,
  "pyscenic_sgbm", "GBM", list(learning_rate = .01, n_estimators = 5000L, max_features = .1, subsample = .9),
  "pyscenic_gbm", "GBM", list(learning_rate = .01, n_estimators = 500L, max_features = .1),
  "pyscenic_et", "ET", list(n_jobs = 1L, n_estimators = 1000L, max_features = "sqrt"),
  "pyscenic_rf", "RF", list(n_jobs = 1L, n_estimators = 1000L, max_features = "sqrt")
)

cross <- crossing(
  dataset_i = seq_along(dataset_files),
  method_i = seq_len(nrow(methods))
)

desired_ints <- 10000

# python stuff
library(reticulate)
reticulate::use_python("/usr/bin/python3")
arboreto <- import("arboreto")
pyscenic <- import("pyscenic")
builtins <- import_builtins()

walk(seq_len(nrow(cross)), function(i) {
  cat("Run ", i, "/", nrow(cross), "\n", sep = "")
  method_params <- methods %>% dynutils::extract_row_to_list(cross$method_i[[i]])

  dataset_file <- dataset_files[[cross$dataset_i[[i]]]]
  dataset <- read_rds(dataset_file)

  out_file <- paste0(derived_folder, "/evals/eval_", gsub("_", "-", dataset$id), "_", gsub("_", "-", method_params$id), ".rds")
  if (!dir.exists(dirname(out_file))) dir.create(dirname(out_file), recursive = TRUE)

  if (!file.exists(out_file)) {

    expr <- as.data.frame(as.matrix(dataset$expression))
    regulators <- unique(dataset$feature_network$from)
    targets <- colnames(expr)
    samples <- rownames(expr)

    # run grnboost2
    adjacencies <- arboreto$algo$diy(
      expression_data = expr,
      regressor_type = method_params$regressor_type,
      regressor_kwargs = method_params$regressor_kwargs,
      tf_names = regulators,
      verbose = FALSE
    ) %>%
      as_tibble()

    pred_pyscenic_static <-
      adjacencies %>%
      transmute(from = factor(TF, levels = regulators), to = factor(target, levels = targets), importance) %>%
      arrange(desc(importance)) %>%
      head(desired_ints)

    # static eval
    ground_truth_static <-
      dataset$feature_network %>%
      transmute(from = factor(from, levels = regulators), to = factor(to, levels = targets), gold = 1, gold_value = strength, gold_effect = effect)

    eval_pyscenic_static <- with(
      pred_pyscenic_static %>%
        left_join(ground_truth_static, by = c("from", "to")) %>%
        mutate(gold = gold %|% 0, gold_value = gold_value %|% 0, gold_effect = gold_effect %|% 0) %>%
        as_tibble(),
      {
        GENIE4::evaluate_ranking_direct(
          values = importance,
          are_true = gold,
          num_positive_interactions = nrow(dataset$feature_network),
          num_possible_interactions = length(regulators) * length(targets)
        )
      }
    )

    # convert to regulons
    modules <- pyscenic$utils$modules_from_adjacencies(adjacencies, expr)

    # rename due to otherwise duplicate names
    for (i in seq_along(modules)) {
      modules[[i]] <- modules[[i]]$rename(paste0("regulon_", i))
    }

    # extract module information
    modules_info <- map_df(seq_along(modules), function(i) {
      mod <- modules[[i]]
      context <- builtins$list(mod$context)
      tibble(
        index = i,
        name = mod$name,
        from = factor(mod$transcription_factor, levels = regulators),
        to = list(match(unlist(mod$genes), targets)),
        weights = mod$weights %>% unlist %>% list,
        score = mod$score,
        context = paste0(context, collapse = ",")
      )
    }) %>%
      mutate(context = factor(context))

    # run aucell to get cell specific regulons
    auc_mtx <- pyscenic$aucell$aucell(expr, modules, num_workers = 8)

    # get 10000 edges per cell
    eval_pyscenic_sample <-
      pbapply::pblapply(
        samples,
        cl = 1,
        function(sample) {
          # ground truth
          ground_truth_sample <-
            with(
              dataset$feature_network_sc[sample,,drop=FALSE] %>% Matrix::summary(),
              tibble(
                from = factor(dataset$feature_network$from[j], levels = regulators),
                to = factor(dataset$feature_network$to[j], levels = targets),
                gold_effect = dataset$feature_network$effect[j],
                gold_value = x,
                gold = 1
              )
            )

          # pred
          scores <- auc_mtx[sample,] %>% as.matrix() %>% .[1,]
          pred_tmp <-
            modules_info %>%
            mutate(
              score = scores[name],
              n = map_int(to, length)
            ) %>%
            arrange(desc(score)) %>%
            mutate(ncs = cumsum(n))
          ix <- which(diff(pred_tmp$ncs < desired_ints) != 0)
          pred_pyscenic_sample <-
            pred_tmp %>%
            transmute(sample, from, to, importance = score) %>%
            unnest(to) %>%
            group_by(from, to) %>%
            summarise(importance = max(importance)) %>%
            ungroup() %>%
            mutate(to = factor(targets[to], levels = targets)) %>%
            head(desired_ints) %>%
            left_join(ground_truth_sample, by = c("from", "to")) %>%
            mutate(gold = gold %|% 0, gold_value = gold_value %|% 0, gold_effect = gold_effect %|% 0)


          eval_pyscenic_sample <- GENIE4::evaluate_ranking_direct(
            values = pred_pyscenic_sample$importance,
            are_true = pred_pyscenic_sample$gold,
            num_positive_interactions = nrow(ground_truth_sample),
            num_possible_interactions = length(regulators) * length(targets)
          )

          pred_pyscenic_static2 <-
            pred_pyscenic_static %>%
            left_join(ground_truth_sample, by = c("from", "to")) %>%
            mutate(gold = gold %|% 0, gold_value = gold_value %|% 0, gold_effect = gold_effect %|% 0) %>%
            as_tibble()

          eval_pyscenic_static2 <- GENIE4::evaluate_ranking_direct(
            values = pred_pyscenic_static2$importance,
            are_true = pred_pyscenic_static2$gold,
            num_positive_interactions = nrow(ground_truth_sample),
            num_possible_interactions = length(regulators) * length(targets)
          )

          bind_rows(
            eval_pyscenic_sample$area_under %>% mutate(type = "casewise"),
            eval_pyscenic_static2$area_under %>% mutate(type = "static")
          ) %>%
            mutate(
              sample = factor(sample, levels = samples),
              dataset = dataset$id,
              method = method_params$id
            )
        }
      ) %>%
      bind_rows()

    out <- lst(pred_pyscenic_static, modules_info, auc_mtx, eval_pyscenic_static = eval_pyscenic_static$area_under %>% mutate(dataset = dataset$id, method = method_params$id), eval_pyscenic_sample)
    write_rds(out, out_file, compress = "gz")
  }
})

outs <- map(seq_len(nrow(cross)), function(i) {
  dataset_id <- read_rds(dataset_files[[cross$dataset_i[[i]]]])$id
  method_id <- methods$id[[cross$method_i[[i]]]]
  out_file <- paste0(derived_folder, "/evals/eval_", gsub("_", "-", dataset_id), "_", gsub("_", "-", method_id), ".rds")
  cat(out_file, "\n", sep = "")
  if (file.exists(out_file)) {
    read_rds(out_file)
  } else {
    NULL
  }
})


comb <- bind_rows(
  map_df(outs, "eval_pyscenic_static") %>% mutate(type = "both_static", type_nice = "pySCENIC evaluated\nagainst static GRN"),
  map_df(outs, "eval_pyscenic_sample") %>% arrange(desc(type)) %>% mutate(type_nice = c(static = "pySCENIC evaluated\nagainst casewise GRN", casewise = "pySCENIC+AUCell evaluated\nagainst casewise GRN")[type])
) %>%
  mutate(type_nice = forcats::fct_inorder(type_nice))

g <-
  ggplot(comb %>% sample_n(n())) +
  geom_point(aes(auroc, aupr, colour = method, size = type)) +
  facet_grid(dataset ~ type_nice) +
  theme_bw() +
  scale_colour_brewer(palette = "Set3") +
  scale_size_manual(values = c("both_static" = 2, static = .75, casewise = .75))
ggsave("fig/network_inference/evaluation.pdf", g, width = 8, height = 6)

g <- ggplot(comb %>% sample_n(n())) +
  geom_point(aes(auroc, aupr, colour = type, size = type)) +
  facet_wrap(dataset ~ method, ncol = 4, scales = "free") +
  theme_bw() +
  scale_colour_brewer(palette = "Set1") +
  scale_size_manual(values = c("both_static" = 2, static = .75, casewise = .75))
ggsave("fig/network_inference/evaluation_2.pdf", g, width = 10, height = 10)


