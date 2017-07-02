#' Input and output
#'
#' \code{load_dna} strips the metadata from the FASTA headers, leaving only the
#' sequence identifier. The headers are assumed to have the format '>SEQID DESC'.
#'
#' @param filename filename
#' @name fagin_io
NULL

#' @rdname fagin_io
#' @export
load_dna <- function(filename) {
  d <- Biostrings::readDNAStringSet(filename)
  # strip metadata from FASTA headers
  names(d) <- gsub(" .*", "", names(d))
  d
}

#' @rdname fagin_io
#' @export
load_queries <- function(filename){
  readr::read_table(
    filename,
    col_names=FALSE,
    col_types="c",

    comment="#"
  )[[1]]
}

#' @rdname fagin_io
#' @export
load_gff <- function(filename) {
  # NOTE: this can fail
  # synder's gff loader will pick up basic errors in the gff format
  gff <- synder::read_gff(filename)

  # NOTE: this can fail
  # fail on missing parent/child relations or malformed attribute columns
  make_GI_with_parent_child_relations(gff)
}

#' @rdname fagin_io
#' @export
load_synmap <- function(filename) {
  synder::read_synmap(filename)
}

#' @rdname fagin_io
#' @export
load_tree <- function(filename){
  ape::read.tree(filename)
}
