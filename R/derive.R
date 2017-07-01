#' Derive N-interval table from a genome
#'
#' @param dna DNAStringSet genome
#' @return IRanges object
#' @export
derive_nstring <- function(dna) {

  # TODO: add tests for this 

  # gregexpr search for a pattern, returning the start and width on each given string
  lapply(dna, function(s) gregexpr("N+", as.character(s))) %>%
  {
    data.frame(
      names = lapply(., '[[', 1) %>% lapply(length) %>% rep.int(names(.), times=.),
      start = unlist(.) %>% as.vector,
      width = lapply(., '[[', 1) %>% lapply(attr, "match.length") %>% unlist %>% as.vector
    )
  } %>%
  # gregexpr stores absence of a match as -1 in both start and width
  dplyr::filter(.data$start != -1) %>%
  {
    IRanges::IRanges(
      start=.$start,
      width=.$width,
      names=.$names
    )
  }

} 

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
