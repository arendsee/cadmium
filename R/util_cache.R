#' Build the cacher that will be used inside rmonad
#'
#' @return rmonad::Cacher object
make_fagin_cacher <- function(){
  rmonad::make_cacher(
    f_save = function(x, filename) {
      if(class(x) == 'TxDb') {
        AnnotationDbi::saveDb(x, filename)
      } else {
        saveRDS(x, filename)
      }
    },
    f_get = function(filename) {
      ext <- sub(".*\\.", "", filename)
      if(ext == 'sqlite') {
        AnnotationDbi::loadDb(filename)
      } else if(ext == 'Rdata') {
        readRDS(filename)
      } else {
        stop("Illegal cache extension: ", filename)
      }
    },
    f_ext = function(cls) {
      if(cls == 'TxDb') {
        '.sqlite'
      } else {
        '.Rdata'
      }
    }
  )
}
