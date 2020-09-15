
#' Generate difference in abundance
#'
#' Multiplies the abundance of cells on particular edges by a given ratio
#'
#' @param model A dyngen model for which to alter the abundance in generated cells
#' @param froms The start milestone of each edge to alter
#' @param tos The end milestone of each edge to alter
#' @param ratios The ratio by which to alter the abundances
#'
#' @export
generate_diff_abundance <- function(model, froms, tos, ratios) {
  tib <- tibble(
    from = froms,
    to = tos,
    ratio = ratios
  )

  model$gold_standard$network <-
    model$gold_standard$network %>%
    full_join(tib, by = c("from", "to")) %>%
    mutate(
      ratio = ratio %|% 1,
      length = length * ratio
    ) %>%
    select(-ratio)


  model
}
