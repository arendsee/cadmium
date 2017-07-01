#' Naming stuff
#'
#' Spaces in species names will be automatically converted to underscores.
#'
#' @param species_name The base name for a species (spaces allowed)
#' @param dir The base directory
#' @param focal_name The name of the focal species
#' @param target_name The name of the target species
#' @name fagin_naming
NULL

make_filename <- function(x, dir=NULL, ext=NULL, spaceToUnderscore=TRUE){
  if(spaceToUnderscore){
    x <- gsub(" ", "_", x)
  }
  if(!is.null(dir)){
    x <- file.path(dir, x)
  }
  if(!is.null(ext)){
    x <- paste0(x, ".", ext)
  }
  x
}

#' @rdname fagin_naming
#' @export
get_genome_filename <- function(species_name, dir) {
  make_filename(species_name, dir=dir, ext="fna")
}

#' @rdname fagin_naming
#' @export
get_gff_filename <- function(species_name, dir) {
  make_filename(species_name, dir=dir, ext="gff3")
}

#' @rdname fagin_naming
#' @export
get_synmap_filename <- function(focal_name, target_name, dir) {
  make_filename(paste0(focal_name, ".vs.", target_name), dir=dir, ext="syn")
}
