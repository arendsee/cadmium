#' Merge named vectors into a data.frame
#'
#' If the input list of vectors has names, these names will be used as the
#' column names. Otherwise automatically generated names will be used (e.g.
#' 'value.x').
#'
#' @param ... named vectors
#' @param .fill value to replace NAs with
#' @export
merge_named_vectors <- function(.fill=NA, ...){
  xs <- list(...)
  ds <- lapply(xs, function(x){
    data.frame(
      name = names(x),
      value = unname(x)
    )
  })
  x <- Reduce(function(...) { merge(..., by='name', all=TRUE) }, ds[-1], ds[[1]]) 
  x[is.na(x)] <- .fill
  if(! is.null(names(xs)))
    names(x) <- c('name', names(xs))
  x
}
