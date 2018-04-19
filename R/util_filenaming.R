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

.label <- function(label, ...){
  label <- if(is.null(label)){
    ""
  } else {
    glue::glue("In {label}: ")
  }
}

.check_name <- function(x){

  "Assert that this is a reasonable string"

  stopifnot(is.character(x))
  stopifnot(length(x) == 1)
  stopifnot(!is.na(x))
  stopifnot(!is.null(x))
  stopifnot(nchar(x) > 0)

  x

}

get_readable_filename <- function(x, dir=NULL, ext=NULL){

  "Make a filename from a base string, a directory path, and list of legal
  extensions. Spaces in the base string are replaced with underscores. Die if
  the resulting filename is not readable. If multiple extensions are given, the
  first one to produces an existing filename is used. If no file is found,
  raise an error"

  x <- gsub(" ", "_", x)

  .check_name(x)
  .check_name(dir)

  if(!is.null(dir)){
    x <- file.path(dir, x)
  }
  if(!is.null(ext)){
    x <- paste0(x, ".", ext)
  }

  readable_files <- x[file.exists(x)]

  n_valid_files <- length(readable_files)

  if(n_valid_files == 0){
    stop(sprintf(
      "No readable file found, checked [%s]",
      paste(x, collapse=", ")
    ))
  }

  if(n_valid_files > 1){
    stop(sprintf(
      "Multiple valid files found: [%s]",
      paste(readable_files, collapse=", "),
      readable_files[1]
    ))
  }

  readable_files[1]
}

#' @rdname fagin_naming
#' @export
get_gff_filename <- function(species_name, dir) {

  "Get a GFF filename from a species name. Filename should have the format:
  'Genus_species.gff3'. 'gff' extension is also allowed."

  species_name %v>%
  .check_name %v>%
  get_readable_filename(dir=dir, ext=c("gff3","gff"))

}

#' @rdname fagin_naming
#' @export
get_synmap_filename <- function(focal_name, target_name, dir) {

  "Get a synteny map filename from the focal and target species names. The
  filename should have the format:
  <focal-genus>_<focal-species>.vs.<target-genus>_<target-species>.syn"

  .check_name(focal_name)
  .check_name(target_name)
  get_readable_filename(
    paste0(focal_name, ".vs.", target_name),
    dir = dir,
    ext = "syn"
  )
}
