#' @export
cni_bred <- create_ti_method_r(
  package_required = c("bred"),
  package_loaded = c("dynutils", "dynwrap", "dplyr"),
  definition = definition(
    method = def_method(
      id = "bred",
      name = "bred"
    ),
    wrapper = def_wrapper(
      input_required = c("expression", "regulators", "targets"),
      input_optional = NULL
    ),

    # describe tuneable parameters
    parameters = parameter_set(
      integer_parameter(
        id = "num_trees",
        default = 1000L,
        distribution = expuniform_distribution(lower = 100L, upper = 100000L)
      ),
      integer_parameter(
        id = "min_node_size",
        default = 10L,
        distribution = expuniform_distribution(lower = 1L, upper = 100L)
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

    # TIMING: done with preproc
    checkpoints <- list(method_afterpreproc = as.numeric(Sys.time()))

    target_ix <- match(names(which(Matrix::colMeans(expression > 0) > 0)), targets)
    outs <-
      pbapply::pblapply(
        cl = 8,
      # map(
        target_ix,
        function(ix) {
          out <- bred::calculate_target_importance(
            target_ix = ix,
            expr_targets = expression,
            expr_regulators = expression,
            samples = samples,
            regulators = regulators,
            targets = targets,
            min_importance = -Inf,
            min_sc_importance = -Inf,
            num_trees = parameters$num_trees,
            min_node_size = parameters$min_node_size
          )
        }
      )

    importance <-
      map_df(outs, "importance") %>%
      arrange(desc(importance)) %>%
      head(parameters$num_int_per_cell) %>%
      rename(strength = importance)
    importance_sc <-
      map_df(outs, "importance_sc") %>%
      select(-importance_sc) %>%
      group_by(cell_id) %>%
      arrange(desc(importance)) %>%
      top_n(parameters$num_int_per_cell, importance) %>%
      ungroup() %>%
      rename(strength = importance)

    # TIMING: done with method
    checkpoints$method_aftermethod <- as.numeric(Sys.time())

    # return output
    # model <-
    wrap_data(
      cell_ids = rownames(expression)
    ) %>%
      add_regulatory_network(
        regulatory_network = importance,
        regulatory_network_sc = importance_sc,
        regulators = regulators,
        targets = targets
      ) %>%
      add_timings(
        timings = checkpoints
      )
  },
  return_function = TRUE
)


