
dynamic_file <- function(base) {
  function(...) {
    file <- file.path(base, ...)
    folder <- fs::path_dir(file)
    if(!file.exists(folder)) {dir.create(folder, recursive = TRUE)}

    file
  }
}

derived_file <- dynamic_file("derived/velocity/")

reread <- function(file, func, load = FALSE) {
  if(!file.exists(file)) {
    func()
  } else {
    if(load) {
      read_rds(file)
    } else {
      TRUE
    }
  }
}



