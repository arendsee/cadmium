#' Load a species GFF
#'
#' @param species species name 
#' @param gfffilename gene feature format file
#' @param genome_lengths tab delimited file with columns [ seqid, length ]
#' @return GI bioconductor object
#' @export
make_load_gff <- function(species, genome_lengths, gfffilename) {
  function(genome_lengths) {
    xxx(species, gfffilename)
  }
}

#' Load a species DNA genome
#'
#' @param species species name
#' @param genome a genomic fasta file name
#' @export
make_load_genome <- function(species, genome_filename) {
  function(genome_filename){
    x <- xxx(genome_filename) 
    xxx(x) <- name
    x
  }
}

#' Load a species predicted proteins
#'
#' @param name the species scientific name
#' @param proteome a proteomic fasta file name
#' @export
make_load_proteome <- function(name, faa_filename) {
  function(faa_filename){
    x <- xxx(faa_filename) 
    xxx(x) <- name
    x
  }
}

#' Load a synteny map
#'
#' @param qname name of query species
#' @param tname name of target species
#' @param synmap synteny map file name
#' @export
load_synmap <- function(qname, tname, qlen, tlen, synmap) {
  function(synmap, qlen, tlen) {
    x <- xxx(synmap, qlen, tlen, qname, tname) 
    x
  }
}

#' Load the species tree
#' 
#' @param treefile a Newick format tree file
#' @return tree species tree
#' @export
load_species_tree <- function(treefile) {
  function(treefile) {

  }
}
