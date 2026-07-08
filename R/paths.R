# Author: Alberto Quaini

#' Locate files in the AGCA paper replication compendium
#'
#' @param ... Path components below the replication-compendium root.
#' @return An absolute path.
#' @export
replication_path <- function(...) {
  root <- system.file(package = "replicateAGCApaper")
  if (basename(root) == "inst" && file.exists(file.path(dirname(root), "DESCRIPTION"))) {
    root <- dirname(root)
  }
  if (!nzchar(root)) {
    root <- normalizePath(file.path(getwd()), mustWork = FALSE)
  }
  file.path(root, ...)
}

#' Check that a replication-compendium file exists
#'
#' @param ... Path components below the replication-compendium root.
#' @return An absolute path.
#' @export
replication_file <- function(...) {
  path <- replication_path(...)
  if (!file.exists(path)) {
    stop("File not found in replication compendium: ", path, call. = FALSE)
  }
  path
}
