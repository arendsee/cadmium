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

.check_name <- function(x){
  if(is.na(x) || is.null(x) || x == ""){
    stop("Expected a species name, found nothing")
  }
}

#' @rdname fagin_naming
#' @export
get_genome_filename <- function(species_name, dir) {
  .check_name(species_name)
  make_filename(species_name, dir=dir, ext="fna")
}

#' @rdname fagin_naming
#' @export
get_gff_filename <- function(species_name, dir) {
  .check_name(species_name)
  # TODO: make this hack less repulsive
  gff3_filename <- make_filename(species_name, dir=dir, ext="gff3")
  gff_filename  <- make_filename(species_name, dir=dir, ext="gff")
  if(file.exists(gff_filename)){
    gff_filename
  } else {
    gff_filename3
  }
}

#' @rdname fagin_naming
#' @export
get_synmap_filename <- function(focal_name, target_name, dir) {
  .check_name(focal_name)
  .check_name(target_name)
  make_filename(paste0(focal_name, ".vs.", target_name), dir=dir, ext="syn")
}
