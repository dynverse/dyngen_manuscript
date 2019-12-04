#' @export
cni_pyscenic_sgbm <- create_ti_method_r(
  package_required = c("reticulate"),
  package_loaded = c("dynutils", "dynwrap", "dplyr"),
  definition = definition(
    method = def_method(
      id = "pyscenic_sgbm",
      name = "pySCENIC SGBM"
    ),
    wrapper = def_wrapper(
      input_required = c("expression", "regulators", "targets"),
      input_optional = NULL
    ),

    # describe tuneable parameters
    parameters = parameter_set(
      numeric_parameter(
        id = "learning_rate",
        default = .01,
        distribution = expuniform_distribution(lower = .0001, upper = 1)
      ),
      integer_parameter(
        id = "n_estimators",
        default = 5000L,
        distribution = expuniform_distribution(lower = 100L, upper = 10000L)
      ),
      numeric_parameter(
        id = "max_features",
        default = .1,
        distribution = expuniform_distribution(lower = .0001, upper = 1)
      ),
      numeric_parameter(
        id = "subsample",
        default = .9,
        distribution = uniform_distribution(lower = 0, upper = 1)
      ),
      integer_parameter(
        id = "num_int_per_cell",
        default = 10000L,
        distribution = expuniform_distribution(lower = 100L, upper = 100000L),
        tuneable = FALSE
      )
    )
  ),

  run_fun = function(
    expression,
    parameters,
    priors,
    seed = NULL,
    verbose = TRUE
  ) {
    regulators <- priors$regulators
    targets <- priors$targets
    samples <- rownames(expression)

    exprdf <- as.data.frame(as.matrix(expression))
    regressor_type <- "GBM"
    regressor_kwargs <- parameters[c("learning_rate", "n_estimators", "max_features", "subsample")]

    reticulate::use_python("/usr/bin/python3")
    arboreto <- reticulate::import("arboreto")
    pyscenic <- reticulate::import("pyscenic")
    builtins <- reticulate::import_builtins()

    # TIMING: done with preproc
    checkpoints <- list(method_afterpreproc = as.numeric(Sys.time()))

    # run grnboost2
    adjacencies <- arboreto$algo$diy(
      expression_data = exprdf,
      regressor_type = regressor_type,
      regressor_kwargs = regressor_kwargs,
      tf_names = regulators,
      verbose = FALSE
    ) %>%
      as_tibble()
    regulatory_network <-
      adjacencies %>%
      transmute(regulator = TF, target, strength = importance)

    # convert to modules
    modules <- pyscenic$utils$modules_from_adjacencies(adjacencies, exprdf, top_n_regulators = list(10L), thresholds = list(), top_n_targets = list())

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
    auc_mtx <- pyscenic$aucell$aucell(exprdf, modules, num_workers = 1)

    # get 10000 edges per cell
    regulatory_network_sc <-
      map_df(
        samples,
        function(sample) {
          scores <- auc_mtx[sample,] %>% as.matrix() %>% .[1,]
          pred_tmp <-
            modules_info %>%
            mutate(
              score = scores[name],
              n = map_int(to, length)
            ) %>%
            arrange(desc(score)) %>%
            mutate(ncs = cumsum(n))
          pred_pyscenic_sample <-
            pred_tmp %>%
            transmute(cell_id = sample, regulator = from, target = to, strength = score) %>%
            unnest(target) %>%
            group_by(cell_id, regulator, target) %>%
            summarise(strength = max(strength)) %>%
            ungroup() %>%
            mutate(target = factor(targets[target], levels = targets)) %>%
            head(parameters$num_int_per_cell)
        }
      ) %>%
      arrange(desc(strength))


    # TIMING: done with method
    checkpoints$method_aftermethod <- as.numeric(Sys.time())

    # return output
    wrap_data(
      cell_ids = rownames(expression)
    ) %>%
      add_regulatory_network(
        regulatory_network = regulatory_network,
        regulatory_network_sc = regulatory_network_sc,
        regulators = regulators,
        targets = targets
      ) %>%
      add_timings(
        timings = checkpoints
      )
  },
  return_function = TRUE
)


