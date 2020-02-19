make_directory_function <- function(path) {
  function(...) {
    file <- file.path(path, paste0(...))
    folder <- fs::path_dir(file)
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
    result = make_directory_function(paste0("result_files/", experiment_id))
  )
}

#' @rdname start_analysis
temporary_file <- function(experiment_id, ...) {
  start_experiment(experiment_id = experiment_id)$temporary(...)
}

#' @rdname start_analysis
result_file <- function(experiment_id, ...) {
  start_experiment(experiment_id = experiment_id)$result(...)
}


#' Obtain an object from cache, if it exists
#'
#' @param file The cache file
#' @param func The function to obtain the object
#'
#' @examples
#' \dontrun{
#' long_execution_function <- function() {
#'   Sys.sleep(100)
#'   10
#' }
#'
#' obj <- "myfile.rds" %cache% long_execution_function
#' }
`%cache%` <- function(file, func) {
  if (!file.exists(file)) {
    x <- func()
    write_rds(x, file, compress = "gz")
    x
  } else {
    read_rds(file)
  }
}


#
#
#
#
#
# load_dataset <- function(id) {
#   dataset <- read_rds(dataset_file(id, "dataset.rds"))
#   dataset
# }
#
#
# load_model <- function(id) {
#   model <- read_rds(dataset_file(id, "model.rds"))
#   model
# }
#
#
#
#
#
# dataset_file <- dynamic_file(derived_file("datasets"))
#
# load_dataset <- function(id) {
#   dataset <- read_rds(dataset_file(id, "dataset.rds"))
#   dataset
# }
