#' List directories in folder 
#'
#' @param p string path of directory.  Default = '.'
#' @param n integer value indicating the folder depth to search
#'
#' @return vector of strings containing the paths to the directories
#' @export
list.dirs.depth.n <- function(p, n) {
  res <- list.dirs(p, recursive = FALSE)
  if (n > 1) {
    add <- list.dirs.depth.n(res, n-1)
    c(res, add)
  } else {
    res
  }
}