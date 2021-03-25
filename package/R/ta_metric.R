#' ABWAP score
#'
#' @param pt1 Aligned pseudotime values for dataset 1
#' @param pt2 Aligned pseudotime values for dataset 2
#'
#' @importFrom pracma trapz
#' @export
ta_abwap <- function(pt1, pt2) {
  1 - pracma::trapz(pt1 + pt2, abs(pt1 - pt2))
}
