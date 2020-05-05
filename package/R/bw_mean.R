#' @export
bw_mean <- function(from, to, data, bandwidth = .05) {
  weights <- apply(to, 1, function(point) {
    weights <- LPCM::kernd(from, point, h = bandwidth)
    weights <- weights / sum(weights)
  })
  t(weights) %*% data
}
