#' Derive N-interval table from a genome
#'
#' @param dna DNAStringSet genome
#' @return GenomicRanges object
#' @export
derive_nstring <- function(dna) { }

#' Derive mono-exonic ORF intervals from a genome
#'
#' Note, here I am doing the stupidest thing possible: just finding all
#' '<START>(...)+<STOP>' intervals. I should replace this with a real ORF
#' finder.
#'
#' @param dna DNAStringSet genome
#' @return GenomicRanges object
derive_orfgff <- function(dna) { }

#' Derive ORF protein sequences from genome and ORF ranges
#'
#' @param dna DNAStringSet genome
#' @param orfgff GenomicRanges object
#' @return AAStringSet object
derive_orffaa <- function(dna, orfgff) { }

#' Derive protein sequence from genome and gene model GFF
#'
#' @param dna DNAStringSet genome
#' @param gff GenomicRanges object annotated with \code{type} and \code{parent}
#' @return AAStringSet object
derive_aa <- function(gff, dna) { }

#' Derive mRNA sequence from genome and gene model GFF
#'
#' @param dna DNAStringSet genome
#' @param gff GenomicRanges object annotated with \code{type} and \code{parent}
#' @return DNAStringSet object
derive_trans <- function(gff, dna) { }

#' Derive ORF intervals from a transcriptome
#'
#' @param trans DNAStringSet genome
#' @return GenomicRanges object
derive_transorf <- function(trans) { }
