#' Store and retrieve data using groups and labels
#'
#' Currently I am caching by saving the files in as RData objects in a local
#' folder. However, I may replace this at some point with a database. This API
#' should be general enough to accomadate such changes. Hence I refer to the
#' data 'group', rather than directory and 'label' rather than basename.
#'
#' @param x Data to be cached
#' @param group A caching group (e.g. a species name)
#' @param label The base label under which to store the data
#' @param cache_dir The directory in which to cache the results
#' @name fagin_cache
NULL

.get_cached_filename <- function(cache_dir, group, label, ext){
  if(is.null(group)){
    file.path(cache_dir, paste(label, ext, sep="."))
  } else {
    file.path(cache_dir, group, paste(label, ext, sep="."))
  }
}

.setup_cache_folders <- function(group, cache_dir){
  if(! dir.exists(cache_dir)){
    dir.create(cache_dir)
  }
  if(! is.null(group)){
    grpdir <- file.path(cache_dir, group)
    if(! dir.exists(grpdir)){
      dir.create(grpdir)
    }
  }
}

#' @rdname fagin_cache
#' @export
to_cache <- function(x, label, group=NULL, cache_dir=".fagin_cache") {

  .setup_cache_folders(group, cache_dir)

  if(class(x) == 'TxDb'){

    cached_filename <- .get_cached_filename(cache_dir, group, label, ext="sqlite")
    AnnotationDbi::saveDb(x, file=cached_filename)

  } else {

    cached_filename <- .get_cached_filename(cache_dir, group, label, ext="RData")
    save(x, file=cached_filename)

  }

  cached_filename

}

#' @rdname fagin_cache
#' @export
from_cache <- function(file, type='Rdata') {

  if(file.exists(file)){
    if(type == 'sqlite'){
      AnnotationDbi::loadDb(file)
    } else {
      # Return x or NULL
      load(file)
      if(exists("x")){
        base::get("x")
      } else {
        NULL
      }
    }
  } else {
    NULL
  }

}
