#' Filter cells
#'
#' @param dataset A dyno dataset
#' @param ... What to filter the dataset by. Can be one of the columns in `cell_info` or a column in `progressions`.
#' @param keep_dimred Whether to retain the dimred or not.
#'
#' @export
filter_cells <- function(dataset, ..., keep_dimred = FALSE) {
  df <-
    dataset$cell_info %>%
    left_join(dataset$progressions, by = "cell_id")

  df <- df %>% filter(...)

  cids <- df$cell_id

  dataset$cell_ids <- cids
  dataset$cell_info <- dataset$cell_info %>% slice(match(cids, .data$cell_id))
  dataset$counts_spliced <- dataset$counts_spliced[cids, , drop = FALSE]
  dataset$counts_protein <- dataset$counts_protein[cids, , drop = FALSE]
  dataset$counts_unspliced <- dataset$counts_unspliced[cids, , drop = FALSE]
  dataset$counts <- dataset$counts[cids, , drop = FALSE]
  dataset$expression <- dataset$expression[cids, , drop = FALSE]
  dataset$milestone_percentages <- dataset$milestone_percentages %>% filter(.data$cell_id %in% cids)
  dataset$progressions <- dataset$progressions %>% slice(match(cids, .data$cell_id))
  dataset$milestone_network <- dataset$milestone_network %>%
    inner_join(dataset$progressions %>% select(from, to) %>% unique, by = c("from", "to"))
  dataset$milestone_ids <- intersect(dataset$milestone_ids, c(dataset$progressions$from, dataset$progressions$to))

  if (keep_dimred && !is.null(dataset$dimred)) {
    dataset$dimred <- dataset$dimred[cids, , drop = FALSE]
    dataset$dimred_milestones <- dataset$dimred_milestones[dataset$milestone_ids, , drop = FALSE]
    progr <-
      dataset$dimred_segment_progressions %>%
      mutate(ix = row_number()) %>%
      semi_join(dataset$milestone_network, by = c("from", "to"))

    dataset$dimred_segment_progressions <- dataset$dimred_segment_progressions[progr$ix, , drop = FALSE]
    dataset$dimred_segment_points <- dataset$dimred_segment_points[progr$ix, , drop = FALSE]
  }

  dataset
}
