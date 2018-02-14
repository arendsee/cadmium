#' Build the cacher that will be used inside rmonad
#'
#' @param archive_dir The main archive directory
#' @param cache_dir The rmonad cache directory (relative to the archive directory)
make_fagin_cacher <- function(archive_dir, cache_dir){
  rmonad::make_recacher(rmonad::make_local_cacher(
    path = file.path(archive_dir, cache_dir),
    save = function(x, filename) {
      if(class(x) == 'TxDb') {
        AnnotationDbi::saveDb(x, filename)
      } else {
        saveRDS(x, filename)
      }
    },
    get = function(filename) {
      ext <- sub(".*\\.", "", filename)
      if(ext == 'sqlite') {
        AnnotationDbi::loadDb(filename)
      } else if(ext == 'Rdata') {
        readRDS(filename)
      } else {
        stop("Illegal cache extension: ", filename)
      }
    },
    ext = function(cls) {
      if(cls == 'TxDb') {
        '.sqlite'
      } else {
        '.Rdata'
      }
    }
  ))
}
