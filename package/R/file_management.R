make_directory_function <- function(prefix, postfix = character(0)) {
  function(...) {
    file <- do.call(file.path, as.list(c(prefix, paste0(...), postfix)))
    folder <- gsub("[^/]*$", "", file)
    if (!file.exists(folder)) {
      dir.create(folder, recursive = TRUE)
    }
    file
  }
}

#' Helper function for knowing where to store files
#'
#' @param experiment_id The identifier of the experiment
#' @param ... File path
#'
#' @export
#'
#' @examples
#' \dontrun{
#' analysis <- start_analysis("usecase_network_inference")
#'
#' dataset_path <- analysis$temporary_file("dataset.rds")
#' temporary_file(experiment_id = "usecase_network_inference", "dataset.rds")
#'
#' figure_path <- analysis$result_file("figure.pdf")
#' result_file(experiment_id = "usecase_network_inference", "figure.pdf")
#'
#' }
start_analysis <- function(experiment_id) {
  list(
    temporary = make_directory_function(paste0("temporary_files/", experiment_id)),
    result = make_directory_function(paste0("result_files/", experiment_id)),
    dataset_folder = make_directory_function(paste0("temporary_files/", experiment_id, "/datasets"), postfix = ""),
    model_file = make_directory_function(paste0("temporary_files/", experiment_id, "/datasets"), postfix = "model.rds"),
    dataset_file = make_directory_function(paste0("temporary_files/", experiment_id, "/datasets"), postfix = "dataset.rds"),
    velocity_file = function(dataset_id, method_id, params_id) {
      make_directory_function(
        prefix = paste0("temporary_files/", experiment_id, "/velocity"),
        postfix = "velocity.rds"
      )(
        paste0(dataset_id, "-", method_id, "-", params_id)
      )
    }
  )
}

#' Obtain an object from cache, if it exists
#'
#' @param file The cache file
#' @param value A code block to obtain the object
#'
#' @export
#'
#' @examples
#' \dontrun{
#' obj <- "myfile.rds" %cache% {
#'   # long execution
#'   Sys.sleep(100)
#'   10
#' }
#'
#' }
`%cache%` <- function(file, value) {
  if (!file.exists(file)) {
    write_rds(value, file, compress = "gz")
    value
  } else {
    read_rds(file)
  }
}

