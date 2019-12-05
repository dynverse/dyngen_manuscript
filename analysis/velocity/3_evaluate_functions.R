pairwise_correlations <- function(x, y) {
  assertthat::assert_that(all(colnames(x) == colnames(y)))
  assertthat::assert_that(all(nrow(x) == nrow(y)))
  map_dbl(seq_len(nrow(x)), function(i) cor(x[i, ], y[i, ], method = "spearman"))
  # map_dbl(seq_len(nrow(x)), function(i) mean(sign(x[i, ]) == sign(y[i, ])))
  # map_dbl(seq_len(nrow(x)), function(i) sum(x[i, ] * y[i, ])/(sqrt(sum(x[i, ]^2)) * sqrt(sum(y[i, ]^2))))
}



run_metric <- function(
  dataset_id = "initial",
  method_id = "scvelo",
  params = list(),
  params_id = "default",
  metric_id = "pairwise_correlations"
) {
  dataset <- load_dataset_velocity(dataset_id, method_id, params_id)

  experiment_file <- dynamic_file(metric_file(dataset_id, method_id, params_id))

  reread(
    experiment_file("scores.rds"),
    function() {
      if(metric_id == "pairwise_correlations") {

      }

      write_rds(scores, experiment_file("scores.rds"))
    }
  )
}


metric_file <- function(dataset_id, method_id, params_id, metric_id) {
  dynamic_file(derived_file("metrics", dataset_id, method_id, params_id, metric_id))
}
