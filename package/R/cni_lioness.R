#' @export
cni_lioness <- create_ti_method_r(
  package_required = c("lionessR"),
  package_loaded = c("dynutils", "dynwrap", "dplyr"),
  definition = definition(
    method = def_method(
      id = "lioness",
      name = "LIONESS"
    ),
    wrapper = def_wrapper(
      input_required = c("expression", "regulators", "targets"),
      input_optional = NULL
    ),

    # describe tuneable parameters
    parameters = parameter_set(
      character_parameter(
        id = "method",
        default = "pearson",
        values = c("pearson", "spearman", "cosine"),
        description = "Which similarity measure to use"
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
    # TIMING: done with preproc
    checkpoints <- list(method_afterpreproc = as.numeric(Sys.time()))

    # adapted from https://github.com/mararie/lionessR
    # should return the exact same results, but with a different implementation
    regulators <- priors$regulators
    targets <- priors$targets
    nr_samples <- nrow(expression)

    drop_same <- match(regulators, targets)
    drop_same[is.na(drop_same)] <- 0

    agg <- calculate_similarity(
      expression[, regulators, drop = FALSE],
      expression[, targets, drop = FALSE],
      method = parameters$method,
      margin = 2
    )

    regulatory_network <-
      reshape2::melt(agg, varnames = c("regulator", "target"), value.name = "strength") %>%
      filter(drop_same[regulator] != as.integer(target)) %>%
      arrange(desc(strength)) %>%
      head(parameters$num_int_per_cell) %>%
      as_tibble()

    regulatory_network_sc <-
      map_df(
        seq_len(nrow(expression)),
        function(i) {
          ss <- calculate_similarity(
            expression[-i, regulators, drop = FALSE],
            expression[-i, targets, drop = FALSE],
            method = parameters$method,
            margin = 2
          )

          (nr_samples * (agg - ss)  + ss) %>%
            reshape2::melt(varnames = c("regulator", "target"), value.name = "strength") %>%
            filter(drop_same[regulator] != as.integer(target)) %>%
            arrange(desc(strength)) %>%
            head(parameters$num_int_per_cell) %>%
            mutate(cell_id = factor(rownames(expression)[[i]], levels = rownames(expression))) %>%
            select(cell_id, regulator, target, strength) %>%
            as_tibble()
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


