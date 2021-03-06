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

  "Open a connection to the indexed fasta file, this will not slurp the whole
  thing into memory. Build an indexed fasta file if one does not exist"

  Rsamtools::indexFa(filename)
  Rsamtools::FaFile(filename)
}

#' @rdname fagin_io
#' @export
load_gene_list <- function(filename){

  "Read a single column list of gene IDs"

  readr::read_table(
    filename,
    col_names=FALSE,
    col_types="c",
    comment="#"
  )[[1]]
}

#' @rdname fagin_io
#' @export
load_tree <- function(filename){

  "Read a phylogenetic tree as a `phylo` object"

  ape::read.tree(filename)
}
