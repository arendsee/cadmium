#' Calculate the lengths of each scaffold in a given genome
#' 
#' @param x Genome 
#' @return ([Conid], [Integer]) data frame
#' @export
get_genome_lengths <- function(x){
  NULL
}

#' Get all ORFs for input genome
#' 
#' @param x Genome
#' @return GFF
#' @export
get_genome_orfs <- function(x){
  NULL
}

#' Get all ORFs from transcripts
#' 
#' @param x Genome 
#' @return FFN
#' @export
get_transcript_orfs <- function(x){
  NULL
}

#' Get all ORFs from annotated mRNAs
#' 
#' @param x Genome
#' @return FFN
#' @export
get_annotated_mRNA_orfs <- function(x){
  NULL
}

#' Get proteins from all predicted genomic models
#' 
#' @param x ([Genome], [GFF])
#' @return [Protein]
#' @export
get_proteins <- function(x){
  NULL
}

#' Parse parent/child relationships from GFF file
#' 
#' @param x GFF
#' @return [(Character, Character)]
#' @export
parse_feature_graph <- function(x){
  NULL
}

#' Predict search intervals from synteny map
#' 
#' @param x (Synmap, Seqids)
#' @return SearchIntervals
#' @export
get_search_intervals <- function(x){
  NULL
}

