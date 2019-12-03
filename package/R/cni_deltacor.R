#' @export
cni_deltacor <- create_ti_method_r(
  package_required = c("reshape2"),
  package_loaded = c("dynutils", "dynwrap"),
  definition = definition(
    method = def_method(
      id = "deltacor",
      name = "Delta Corr."
    ),
    wrapper = def_wrapper(
      input_required = c("expression", "regulators", "targets"),
      input_optional = NULL
    ),

    # describe tuneable parameters
    parameters = parameter_set(
      character_parameter(
        id = "method",
        default = "spearman",
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

    regulators <- priors$regulators
    targets <- priors$targets

    drop_same <- match(regulators, targets)
    drop_same[is.na(drop_same)] <- 0

    all_sim <- calculate_similarity(
      expression[, regulators, drop = FALSE],
      expression[, targets, drop = FALSE],
      method = parameters$method,
      margin = 2
    )

    regulatory_network <-
      reshape2::melt(all_sim, varnames = c("regulator", "target"), value.name = "strength") %>%
      filter(drop_same[regulator] != as.integer(target)) %>%
      arrange(desc(strength)) %>%
      head(parameters$num_int_per_cell) %>%
      as_tibble()

    regulatory_network_sc <-
      map_df(
        seq_len(nrow(expression)),
        function(i) {
          sub_sim <- calculate_similarity(
            expression[-i, regulators, drop = FALSE],
            expression[-i, targets, drop = FALSE],
            method = parameters$method,
            margin = 2
          )

          reshape2::melt(all_sim - sub_sim, varnames = c("regulator", "target"), value.name = "strength") %>%
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


