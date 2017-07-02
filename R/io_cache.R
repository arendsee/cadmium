#' Store and retrieve data using groups and labels
#'
#' Currently I am caching by saving the files in as RData objects in a local
#' folder. However, I may replace this at some point with a database. This API
#' should be general enough to accomadate such changes. Hence I refer to the
#' data 'group', rather than directory and 'label' rather than basename.
#'
#' @param x Data to be cached
#' @param group a caching group (e.g. a species name)
#' @param label the base label under which to store the data
#' @name fagin_cache
NULL

#' @rdname fagin_cache
#' @export
to_cache <- function(x, group, label) {

  # TODO: save x

}

#' @rdname fagin_cache
#' @export
from_cache <- function(group, label) {

  # TODO: load x

}

