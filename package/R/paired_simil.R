#' Calculate paired similarity
#'
#' @inheritParams dynutils::calculate_distance
#'
#' @importFrom dynutils calculate_similarity
#'
#' @export
#' @examples
#' x <- matrix(runif(500), nrow = 10)
#' y <- x * rnorm(500, mean = 1, sd = .05)
#' margin <- 1; method <- "pearson"
paired_simil <- function(x, y, method = "spearman", margin = 1) {
  assertthat::assert_that(
    nrow(x) == nrow(y),
    ncol(x) == ncol(y)
  )
  out <- map_dbl(
    if (margin == 1) {seq_len(nrow(x))} else {seq_len(ncol(x))},
    function(i) {
      subx <- if (margin == 1) x[i,,drop = FALSE] else x[,i,drop = FALSE]
      suby <- if (margin == 1) y[i,,drop = FALSE] else y[,i,drop = FALSE]
      dynutils::calculate_similarity(subx, suby, method = method, margin = margin)[[1]]
    }
  )

  names(out) <- if (margin == 1) rownames(x) else colnames(x)

  out
}
