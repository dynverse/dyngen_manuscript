#' Obtain velocity with velocyto
#'
#' @importFrom reticulate import
#' @importFrom Matrix t sparseMatrix
#'
#' @importFrom reticulate py_dict np_array
#' @export
rnav_run_velocyto <- function(
  spliced,
  unspliced,
  assumption = "constant_velocity",
  velocyto = reticulate::import("velocyto")
) {
  # prepare loom file ------------------------------------------------------------------
  S <- Matrix::t(spliced)
  U <- Matrix::t(unspliced)
  A <- Matrix::sparseMatrix(integer(), integer(), dims = dim(S)) %>% as("dgCMatrix")

  cell_ids <- rownames(spliced)
  feature_ids <- colnames(spliced)


  loom_file <- tempfile()

  loompy <- import("loompy")
  loompy$create(
    loom_file,
    layers = py_dict(
      c("", "spliced", "unspliced", "ambiguous"),
      list(
        S,
        S,
        U,
        A
      )
    ),
    row_attrs = py_dict(
      c("Gene"),
      list(np_array(feature_ids))
    ),
    col_attrs = py_dict(
      c("CellID", "Cluster"),
      list(np_array(cell_ids), np_array(rep(1, length(cell_ids))))
    )
  )

  # run velocyto, based on http://velocyto.org/velocyto.py/tutorial/analysis.html
  vlm <- velocyto$VelocytoLoom(loom_file)
  vlm$`_normalize_S`(relative_size=rowSums(vlm$S), target_size=mean(rowSums(vlm$S)))
  vlm$`_normalize_U`(relative_size=rowSums(vlm$U), target_size=mean(rowSums(vlm$U)))
  vlm$perform_PCA("S")
  vlm$knn_imputation(n_pca_dims=20L, pca_space = TRUE, k=20L, balanced=TRUE, b_sight=50L, b_maxl=20L, n_jobs=16L, size_norm=FALSE)
  vlm$fit_gammas()
  vlm$predict_U()
  vlm$calculate_velocity()
  vlm$calculate_shift(assumption=assumption)
  vlm$extrapolate_cell_at_t(delta_t=1.)

  # extrapolated profile minus the smoothed spliced counts
  velocity_vector <- as(Matrix::t(vlm$Sx_sz_t - vlm$Sx), "dgCMatrix")
  dimnames(velocity_vector) <- dimnames(spliced)

  lst(
    velocity_vector,
    vlm
  )
}
