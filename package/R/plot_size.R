#' Helper function for reducing the file size of dynplot2 plots
#'
#' @param g A dynplot2 plot
#' @param verbose Whether to print a message
#' @export
reduce_size <- function(g, verbose = FALSE) {
  if (verbose) size1 <- pryr::object_size(g)

  attr(g$data, "data")$dataset <- NULL
  g$plot_env$p <- NULL
  attr(g$plot_env$cell_info, "data")$dataset <- NULL
  g$plot_env$data$dataset <- NULL
  g$plot_env$dataset <- NULL
  g$plot_env$trajectory <- NULL

  if (verbose) {
    size2 <- pryr::object_size(g)
    cat(
      "Reduced plot from ",
      utils:::format.object_size(size1, "auto"),
      " to ",
      utils:::format.object_size(size2, "auto"),
      "\n",
      sep = ""
    )
  }

  g
}
#' @export
recursive_check_size <- function(obj, path = "", thresh = 10 * 1024^2) {
  if (str_length(path) > 50) return()

  obj_size <- pryr::object_size(obj)
  if (obj_size <= thresh) return()

  attrs <- attributes(obj)
  attr_sizes <- map_dbl(attrs, pryr::object_size)

  # is subsettable (but not an environment)
  is_subsettable <- is.list(obj) || is(obj, "Layer") || is.environment(obj)
  if (is_subsettable) {
    li_sizes <- sapply(obj, pryr::object_size)
  } else {
    li_sizes <- NULL
  }

  if (any(c(li_sizes, attr_sizes) > thresh)) {
    walk(names(attrs), function(n) {
      path <- paste0("attr(", path, ", \"", n, "\")")
      recursive_check_size(attrs[[n]], path)
    })
    if (is_subsettable) {
      if (!is.null(names(obj))) {
        walk(names(obj), function(n) {
          path <- paste0(path, "$", n)
          recursive_check_size(obj[[n]], path)
        })
      } else {
        walk(seq_along(obj), function(i) {
          path <- paste0(path, "[[", i, "]]")
          recursive_check_size(obj[[i]], path)
        })
      }
    }
  } else {
    cat(path, " <- NULL # has size ", sep = "")
    print(obj_size)
  }
}
